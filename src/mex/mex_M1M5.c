/*------------------------------------------------------------------------
 * mex_M1M5.c
 *------------------------------------------------------------------------
 *Fast Multipole code for particles in the plane designed for use in 
 *Matlab. At each particle z(k) (corresponding to a discretization point on
 *a parameterized boundary in the complex plane) the sums
 *    M1 = dens(m)*zp(m)/(z(m)-z(k)) 
 *    M5 = -dens(m)*real(zp(m)/(z(m)-z(k))) +
 *         i*conj(dens(m))*imag(zp(m)*conj(z(m)-z(k)))/(conj(z(m)-z(k)))^2
 *are computed simultaneously for all m != k, where dens is the density and
 *zp is the derivative of the boundary with respect to parameter.
 *
 *The Matlab syntax is
 *
 *[M1,M5,maxlevel] = mex_M1M5(z,zp,dens,tol,levmax,numthreads)
 * 
 *where
 *  z           - Positions of the boundary particles, complex.
 *  zp          - Tangents of the boundary, complex.
 *  dens        - Density, real/complex.
 *  tol         - The requested accuracy(a bit pessimistic, tol = 1e-13
 *                is sufficient for machine epsilon accuracy)
 *  levmax      - The max number of box subdivisions before we resort to
 *                slow direct evaluation. This is to be able to handle
 *                (albeit slowly) highly non-uniform distributions of 
 *                particles.
 *  numthreads  - The number of computational threads to use. 
 *  M1,M5       - The output sums
 *  maxlevel    - The final depth of the refinement recursion.
 *
 *Some details:
 *
 *SSE3       - This code uses SSE3-instructions in some inner loops, so
 *             a Pentium 4 processor or later is required. For this to 
 *             be efficient, some arrays need to be aligned on a 16 byte
 *             boundary. To accomplish this, gcc's aligned malloc was
 *             rewritten to use Matlab's mxMalloc (the code is in 
 *             mm_mxmalloc.h). 
 *Threads    - The code is multithreaded. To be sure no race conditions 
 *             occur, mutexes are used. In order not to lock up the entire 
 *             target arrays, a number of mutexes are used. Based on some 
 *             quick testing on only two machines, 64 mutexes are used for 
 *             2 threads and 2048 for 4 threads or above. This can most 
 *             likely be tuned.
 *Memory     - It is hard to exactly calculate exactly how much memory
 *             the code uses but an estimate is:
 *                  ints    : 2*n_particles+5*n_boxes+326
 *                  doubles : 7*n_particles+n_terms*(n_terms+4*ntboxes)+
 *                            10*n_threads*n_terms
 *             Where n_particles is the number of particles, n_boxes 
 *             is the total number of boxes at the finest grid, n_terms is 
 *             the number of terms in the multipole expansion, ntboxes is 
 *             the number of boxes containing > n_terms/2 particles at the 
 *             finest grid and n_threads is the number of threads used.
 *Efficiency - This code runs in O(nlogn) time, for n particles, so 
 *             technically it is not an implementation of the fast
 *             multipole method, which is usually O(n). A simple 
 *             adaptive scheme is used, but there will be problems with
 *             efficiency for highly non-uniform distributions of 
 *             particles.
 *
 *By Rikard Ojala, 2012-04-04.
 */
#include "mex.h"
#include "mex_FMM.h"

/*Error messages unique for this function, the rest of the messages are in
 *mex_FMM.h.*/
#define STD_ERROR "Error in mex_M1M5: "
#define STD_WARN "Warning in mex_M1M5: "
#define SYNTAX_ERROR "incorrect syntax. \n\n usage : [M1,M5,maxlevel] = mex_M1M5(z,zp,dens,tol,levmax,numthreads)"
#define ZP_CLASS_ERROR "zp must be real or complex."
#define DENS_CLASS_ERROR "dens must be real or complex."
#define Z_ZP_SIZE_ERROR "length(z) is not equal to length(zp)."
#define Z_DENS_SIZE_ERROR "length(z) is not equal to length(dens)."
#define ZP_DIM_ERROR "zp must be a column-vector."
#define DENS_DIM_ERROR "dens must be a column-vector."

/*Most of the working data is bunched up in larger arrays depending on the
 *type of data and where it is used. These are the offsets of the
 *sub-arrays in the main arrays.*/
/*----------------*/
/*realloc_data(ints):*/
#define BOX_OFFSETS_OFFSET         (0)
#define NPARTICLES_IN_BOX_OFFSET   (nside*nside+1)
#define TAYLOR_BOXES_OFFSET        (2*nside*nside+1)
#define TAYLOR_BOX_NBR_OFFSET      (3*nside*nside+1)
#define ASSIGN_TMP_OFFSET          (4*nside*nside+1)
/*int_data:*/
#define IN_BOX_OFFSET              (0)
#define PARTICLE_OFFSETS_OFFSET    (n_particles)
#define ILIST_X_OFFSET             (2*n_particles)
#define ILIST_Y_OFFSET             (2*n_particles+108)
#define INTERACTION_LIST_OFFSET    (2*n_particles+216)
/*double_data:*/
#define M1_C_OFFSET                (0)
#define M2_C_OFFSET                (2*n_particles)
#define QA_C_OFFSET                (4*n_particles)
#define QB_OFFSET                  (6*n_particles)
#define CC_OFFSET                  (7*n_particles + (n_particles&1))
#define DBL_THDATA_OFFSET          (7*n_particles+(n_particles&1)+n_terms*n_terms+(16-((n_terms*n_terms)&0xf)))
#define MPOLE1_C_OFFSET            (0)
#define MPOLE2_C_OFFSET            (2*n_terms)
#define TAYLOREXP1_C_OFFSET        (4*n_terms)
#define TAYLOREXP2_C_OFFSET        (6*n_terms)
#define ZSH_POWERS_C_OFFSET        (8*n_terms)
#define DBL_THDATA_SIZE            (10*n_terms)

/*The structure that is passed as the input parameter to MpolesWorker
 *and MpolesWorkerSum.*/
typedef struct {
    double *z_re, *z_im;
    double *thread_data;
    double *double_data, *localexp;
    int *realloc_data, *int_data, *cursquare;
    int n_terms, n_particles, ntboxes, nside, taylor_threshold;
    int maxparticles_in_box;
} MpolesWorkerStruct;

/*The structure that is passed as the input parameter to DirectWorker.*/
typedef struct {
    double *z_re, *z_im;
    double *double_data;
    int *int_data, *realloc_data, *ilist_x, *ilist_y, *cursquare;
    int nside, n_particles;
} DirectWorkerStruct;

/*------------------------------------------------------------------------
 *The main multipole driver function. The function spawns threads that: 
 *
 *1) Computes the multipole expansion of the particles in each box
 *2) Traverses the interaction list of these boxes and, depending on
 *   the number of target particles, either evaluates the interactions 
 *   directly at these targets, evaluates the multipole expansion at
 *   these targets, or converts the multipole expansion into a taylor
 *   expansion centered around the center of the target box, and sums these
 *   together.
 *3) After 1 and 2 are finished for all boxes, a new set of threads 
 *   evaluates the local taylor series in each box with enough particles.
 *------------------------------------------------------------------------
 */
void Mpoles(double *z_re, double *z_im,
            double* double_data, int* int_data, int* realloc_data, 
            int n_terms, int nside, int taylor_threshold, int n_particles, 
            int num_threads);

/*------------------------------------------------------------------------
 *When we have reached the level where all boxes contain less than n_terms
 *particles, this function computes direct interaction between nearest
 *neighbors as well as box self-interaction.
 *------------------------------------------------------------------------
 */
void Direct(double *z_re, double *z_im, 
            double* double_data, int* int_data, int* realloc_data, 
            int nside, int n_particles, int num_threads);

/*------------------------------------------------------------------------
 *This function assigns particles to boxes on the current grid.
 *------------------------------------------------------------------------
 */
void Assign(double *z_re, double *z_im, int n_particles, int nside,
            int *maxparticles_in_box, int *int_data, int *realloc_data);

/*-----------------------------------------
 *The main entry function
 *-----------------------------------------*/
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
 
    double *double_data;
    double *M1_re,*M1_im,*M2_re,*M2_im,*z_re,*z_im,*qa_c,*qb;
    double *M1_c,*M2_c,*C,*CC,*zp_re,*zp_im,*dens_re,*dens_im;
    double tol;
    int *int_data,*realloc_data;
    int *ilist_x,*ilist_y;
    int lev_max,n_particles,n_terms,taylor_threshold,csize,i,j,k;
    int current_level,nside,ptr,numthreads,maxparticles_in_box;

    /*Quit with error if the number of input parameters is not equal
     *to six(all inputs are required) and the number of outputs are two
     *or three(current_level is optional).*/
    if (nrhs != 6 || nlhs < 2 || nlhs > 3) 
        mexErrMsgTxt(STD_ERROR SYNTAX_ERROR);

    /*Quit with error if z is not of class double*/
    if(!mxIsClass(prhs[0],"double"))
        mexErrMsgTxt(STD_ERROR Z_CLASS_ERROR);

    /*Quit with error if zp is not of class double*/
    if(!mxIsClass(prhs[1],"double"))
        mexErrMsgTxt(STD_ERROR ZP_CLASS_ERROR);
    
    /*Quit with error if dens is not of class double*/
    if(!mxIsClass(prhs[2],"double"))
        mexErrMsgTxt(STD_ERROR DENS_CLASS_ERROR);
    
    /*Quit with error if tol is not a double scalar*/
    if(mxGetM(prhs[3]) != 1 || mxGetN(prhs[3]) != 1 || !mxIsClass(prhs[3],"double"))
        mexErrMsgTxt(STD_ERROR TOL_CLASS_ERROR);

    /*Quit with error if levmax is not a double scalar*/
    if(mxGetM(prhs[4]) != 1 || mxGetN(prhs[4]) != 1 || !mxIsClass(prhs[4],"double"))
        mexErrMsgTxt(STD_ERROR LEVMAX_CLASS_ERROR);

    /*Quit with error if numthreads is not a double scalar*/
    if(mxGetM(prhs[5]) != 1 || mxGetN(prhs[5]) != 1 || !mxIsClass(prhs[5],"double"))
        mexErrMsgTxt(STD_ERROR NUMTHREADS_CLASS_ERROR);

    /*Get some constants*/
    tol = mxGetScalar(prhs[3]);
    lev_max = (int)mxGetScalar(prhs[4]);
    numthreads = (int)mxGetScalar(prhs[5]);
    n_particles = mxGetM(prhs[0]);

    /*Quit with error if the number of threads is less than 1.*/
    if(numthreads < 1)
        mexErrMsgTxt(STD_ERROR NUM_THREADS_TOO_LOW_ERROR);
    
    /*Quit with error if z is not a column vector*/
    if((n_particles == 1 && n_particles < (int)mxGetN(prhs[0])) || 
       (n_particles > 1 && mxGetN(prhs[0])> 1))
        mexErrMsgTxt(STD_ERROR Z_DIM_ERROR);

    /*Quit with error if zp is not a column vector*/
    if((n_particles == 1 && n_particles < (int)mxGetN(prhs[1])) || 
       (n_particles > 1 && mxGetN(prhs[1])> 1))
        mexErrMsgTxt(STD_ERROR ZP_DIM_ERROR);
    
    /*Quit with error if dens is not a column vector*/
    if((mxGetM(prhs[2]) == 1 && mxGetM(prhs[2]) < mxGetN(prhs[2])) || 
       (mxGetM(prhs[2]) > 1 && mxGetN(prhs[2])>1))
        mexErrMsgTxt(STD_ERROR DENS_DIM_ERROR);
    
    /*Quit with error if z is not the same size as zp*/
    if(n_particles != mxGetM(prhs[1])) {
        mexErrMsgTxt(STD_ERROR Z_ZP_SIZE_ERROR);
    }
    
    /*Quit with error if z is not the same size as dens*/
    if(n_particles != mxGetM(prhs[2])) {
        mexErrMsgTxt(STD_ERROR Z_DENS_SIZE_ERROR);
    }

    /*If z is empty then return two empty matrices and current_level = 0 if
     *current_level was requested.*/
    if(mxGetM(prhs[0]) == 0 || mxGetN(prhs[0]) == 0) {
        plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]), mxGetN(prhs[0]), mxREAL);
        plhs[1] = mxCreateDoubleMatrix(mxGetM(prhs[0]), mxGetN(prhs[0]), mxREAL);
        if(nlhs==3)
            plhs[2] = mxCreateDoubleScalar(0);
        return;
    }
        
    /*Warn user if lev_max is larger than 16.*/
    if(lev_max > 16) {
        lev_max = 16;
        mexWarnMsgTxt(STD_WARN LEVMAX_WARNING);
    }
    /*We need at least 3 levels.*/
    if(lev_max < 3) {
        lev_max = 3;
    }
    
    /*Get the positions of the particles.*/
    z_re = mxGetPr(prhs[0]);
    z_im = mxGetPi(prhs[0]);
    
    /*The number of terms required in the multipole expansion to
     *get an accuracy of tol. Actually, this formula is a bit pessimistic,
     *to get a relative error around machine epsilon it suffices to 
     *specify tol = 1e-13.*/
    n_terms = -(int)ceil(log(tol)/LOG2);
    /*The breakpoint in terms of number of particles in a box deciding if
     *to use summed up local taylor expansions or direct evaluated 
     *multipole expansions.*/
    taylor_threshold = n_terms/2;
 
    /*The main double data vector. It is aligned on a 16 byte boundary
     *to allow for SSE instructions. Contains:
     * M1_c[2*n_particles]     Output field. Complex. Read/Write. Mutexed
     * M2_c[2*n_particles]     Output field. Complex. Read/Write. Mutexed
     * QA_c[2*n_particles]     Temp density. Complex. Read-only
     * QB[n_particles]         Temp density. Real. Read-only
     * CC[n_terms*n_terms]     Binomial matrix. Real. Read-only
     * --- One for each thread --- 
     * mpole1_c[2*n_terms]     Temp data. Separated from other threads
     * mpole2_c[2*n_terms]     Temp data. Separated from other threads
     * taylorexp1_c[2*n_terms] Temp data. Separated from other threads
     * taylorexp2_c[2*n_terms] Temp data. Separated from other threads
     * zsh_powers_c[2*n_terms] Temp data. Separated from other threads
     */
    double_data = _mm_mxMalloc((7*n_particles+(n_particles&1)+n_terms*n_terms+(16-((n_terms*n_terms)&0xf))+numthreads*10*n_terms)*sizeof(double),16);
    M1_c = &double_data[M1_C_OFFSET];
    M2_c = &double_data[M2_C_OFFSET];
    qa_c = &double_data[QA_C_OFFSET];
    qb = &double_data[QB_OFFSET];
    CC = &double_data[CC_OFFSET];
    
    /*Clear the output*/
    for(j=0;j<2*n_particles;j++) {
        M1_c[j] = 0;
        M2_c[j] = 0;
    }
    
    /*Temporary pointers to the derivatives.*/
    zp_re = mxGetPr(prhs[1]);
    zp_im = mxGetPi(prhs[1]);
  
    /*Temporary pointers to the density.*/
    dens_re = mxGetPr(prhs[2]);
    dens_im = mxGetPi(prhs[2]);
    
    /*The code assumes complex derivatives.*/
    if(zp_re == NULL)
        zp_re = mxCalloc(n_particles,sizeof(double));    
    if(zp_im == NULL)
        zp_im = mxCalloc(n_particles,sizeof(double));
    
    /*The code assumes complex density.*/
    if(dens_re == NULL)
        dens_re = mxCalloc(n_particles,sizeof(double));    
    if(dens_im == NULL)
        dens_im = mxCalloc(n_particles,sizeof(double));
    
    for(j = 0;j<n_particles;j++) {
        qa_c[2*j] = -dens_re[j]*zp_re[j] + dens_im[j]*zp_im[j];
        qa_c[2*j+1] = -dens_re[j]*zp_im[j] - dens_im[j]*zp_re[j];
        /*This vector is actually purely imaginary, but that will only
          become apparent later when multiplying with it.*/
        qb[j] = -2*(dens_re[j]*zp_im[j] - dens_im[j]*zp_re[j]);
    }
    
    /*In the off chance that the input particles are purely real or 
     *imaginary, we will run into trouble. Fix this.*/
    if(z_re==NULL) {
       z_re = mxCalloc(n_particles,sizeof(double));
    }
    if(z_im==NULL) {
       z_im = mxCalloc(n_particles,sizeof(double));
    }

    /*Check that all particles are inside the computational box, that is 
     *the unit box [-0.5 0.5 -0.5 0.5]. Actually the code allows for points
     *that are slightly outside this box, to allow for numerical artefacts
     *resulting from scaling and translation.*/
    for(j=0;j<n_particles;j++) {
        int re = (int)floor(0.99999999*z_re[j] + 0.5);
        int im = (int)floor(0.99999999*z_im[j] + 0.5);
        /*Quit with error if particles are outside the computational box*/
        if(re*re + im*im > 0)
            mexErrMsgTxt(STD_ERROR OUTSIDE_ERROR);
    }

    /*Compute a matrix containing binomial numbers. It is used when
     *converting from multipole to taylor expansions.*/
    csize = 2*n_terms-1;
    C = mxCalloc(csize*csize,sizeof(double));

    for(i = 0;i<csize;i++)
        C[i] = 1;

    for(i = 2;i<=csize;i++)
        for(j=2;j<=i;j++)
            C[(j-1)*csize +i-1] = C[(j-2)*csize +i-2]+C[(j-1)*csize +i-2];
    
    k = 1;
    for(i = 1;i<=n_terms;i++) {
        for (j = 1;j<=n_terms;j++)
            CC[(j-1)*n_terms+i-1]=k*C[(i-1)*csize+j+i-2];
        k*=-1;
    }
    mxFree(C);

    /*We start with a level-2 grid, that is 4x4 boxes. On coarser grids
     *the interaction list is empty, so a bit of work would be wasted
     *setting up for these cases.*/
    current_level = 2;
    nside = 4;

    /*int_data. Contains data which is read-only in the threads. Contains:
     *
     *in_box[n_particles]           The box point j is assigned to
     *particle_offsets[n_particles] The addresses of the particles in each 
     *                              box.
     *ilist_x[108]                  x-coordinates of the boxes in the 
     *                              interaction list.
     *ilist_y[108]                  y-coordinates of the boxes in the 
     *                              interaction list.
     *interactionlist[108]          interaction list in terms of relative
     *                              box numbers.
     */
    int_data = mxCalloc((2*n_particles+3*108),sizeof(int));
    ilist_x = &int_data[ILIST_X_OFFSET];
    ilist_y = &int_data[ILIST_Y_OFFSET];
    /*Initialize the interaction list. See mex_FMM.h.*/
    Intlist(ilist_x,ilist_y);
    
    /*Set up the mutexes. output_mutex[] and localexp_mutex[] are mutexes 
     *protecting the output and localexp_c arrays.*/
    if(numthreads <= 2)
        num_mutexes = 64;
    else
        num_mutexes = 4*512;

    output_mutex = mxMalloc(num_mutexes*sizeof(MUTEX_TYPE));
    localexp_mutex = mxMalloc(num_mutexes*sizeof(MUTEX_TYPE));
    for(i=0;i<num_mutexes;i++) {
        INIT_MUTEX(&output_mutex[i]);
        INIT_MUTEX(&localexp_mutex[i]);
    }
    INIT_MUTEX(&mpoleSqTblmutex);
    INIT_MUTEX(&directSqTblmutex);
    
    /*The realloc_data array contains read-only data shared by the 
     *threads. The array is realloc'd each level in the hierarchy since
     *it contains members that depend on the number of sides of the current
     *grid. Calloc-ing this once and for all would be bad since for 
     *uneven distributions of particles many layers may be needed to 
     *resolve tight clusters of particles. A safe size would thus be 
     *prohibitively large. Instead we use a more moderate initial size and
     *start reallocing when this memory is no longer sufficient. 
     *
     *The array contains: 
     *box_offsets[nside*nside+1]          The addresses of the boxes in 
     *                                    particle_offsets
     *nparticles_in_box[nside*nside]      The number of particles in each 
     *                                    box
     *particlenum_in_box[nside*nside+1]   Temporary array in Assign
     *taylor_boxes[nside*nside]           Contains the boxes with more than 
     *                                    taylor_threshold particles. 
     *taylor_box_nbr[nside*nside]         contains the "local numbering" of 
     *                                    the taylor boxes*/
    realloc_data = mxCalloc((5*PREALLOCNSIDE*PREALLOCNSIDE+2),sizeof(int));
    
    maxparticles_in_box = n_terms+1;
    /*This is the main loop. Work our way down finer grids until the 
     *maximum number of particles in a box is less than a threshold, or  
     *until the maximum number of refinements is reached, in 
     *which case we break and compute the rest of the interactions via  
     *direct evaluation.*/
    while(current_level < lev_max && maxparticles_in_box > n_terms) {
 
        if(nside > PREALLOCNSIDE)
            realloc_data = mxRealloc(realloc_data,(5*nside*nside+2)*sizeof(int));

        Assign(z_re,z_im,n_particles,nside,
                &maxparticles_in_box,int_data,realloc_data);
        Mpoles(z_re,z_im,double_data,int_data,realloc_data,
                n_terms,nside,taylor_threshold,n_particles,numthreads);

        current_level++;
        nside *= 2;
        
    }
    nside /=2;
    /*Compute the last interactions via direct evaluation.*/
    Direct(z_re,z_im,double_data,int_data,realloc_data,nside,n_particles,numthreads);
    
    /*Clean up mutexes*/
    for(i=0;i<num_mutexes;i++) {
        DESTROY_MUTEX(&output_mutex[i]);
        DESTROY_MUTEX(&localexp_mutex[i]);
    }
    DESTROY_MUTEX(&mpoleSqTblmutex);
    DESTROY_MUTEX(&directSqTblmutex);
    
    /*Free some of the allocated memory*/
    mxFree(output_mutex);
    mxFree(localexp_mutex);
    mxFree(int_data);
    mxFree(realloc_data);
    
    /*Create the output vectors.*/
    plhs[0] = mxCreateDoubleMatrix(n_particles, 1, mxCOMPLEX);
    M1_re = mxGetPr(plhs[0]);
    M1_im = mxGetPi(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(n_particles, 1, mxCOMPLEX);
    M2_re = mxGetPr(plhs[1]);
    M2_im = mxGetPi(plhs[1]);
    
    /*Convert vector back to matlab format*/
    ptr = 0;
    for(j = 0;j<n_particles;j++, ptr+=2) {
        M1_re[j] = -M1_c[ptr];
        M1_im[j] = -M1_c[ptr+1];
        M2_re[j] = (M1_c[ptr]+M2_c[ptr])/2;
        M2_im[j] = (M1_c[ptr+1]+M2_c[ptr+1])/2;
    }
    /*If requested, the last level of the hierarchy is returned.*/
    if(nlhs==3)
        plhs[2] = mxCreateDoubleScalar(current_level-1);
    
    /*Free the last of the allocated memory. (Aligned free)*/
    _mm_mxFree(double_data);
}
/*------------------------------------------------------------------------
 *The main multipole driver function.
 *The function spawns threads that:
 *1) Computes the multipole expansion of the particles in each box
 *2) Traverses the interaction list of these boxes and, depending on
 *   the number of target particles, either evaluates the interactions 
 *   directly at these targets, evaluates the multipole expansion at
 *   these targets, or converts the multipole expansion into a taylor
 *   expansion centered around the center of the target box, and sums these
 *   together.
 *3) After 1 and 2 are finished for all boxes, a new set of threads 
 *   evaluates the local taylor series in each box with enough particles.
 *------------------------------------------------------------------------
 */
void Mpoles(double *z_re, double *z_im,
            double* double_data, int* int_data, int* realloc_data, 
            int n_terms, int nside, int taylor_threshold, int n_particles, 
            int num_threads) {

    MpolesWorkerStruct* arguments;
    /*See below for explanations of these pointers.*/
    int *interaction_list = &int_data[INTERACTION_LIST_OFFSET];
    int *taylor_box_nbr = &realloc_data[TAYLOR_BOX_NBR_OFFSET];
    int *taylor_boxes = &realloc_data[TAYLOR_BOXES_OFFSET];
    double *localexp_c;
    double *double_offset;
    /*The interaction list.*/
    int *ilist_x = &int_data[ILIST_X_OFFSET];
    int *ilist_y = &int_data[ILIST_Y_OFFSET];
    /*The number of particles in each box.*/
    int *nparticles_in_box = &realloc_data[NPARTICLES_IN_BOX_OFFSET];
    /*ntboxes is the number of Taylor boxes and cursquare is the next
     *square to be treated by the threads.*/
    int ntboxes,cursquare=0,i;

    /*Get boxes that contain enough particles to qualify for taylor series
     *conversion. Store the numbering of these in taylor_boxes. ntboxes 
     *is the number of such boxes.*/
    ntboxes=0;
    for(i=0;i<nside*nside;i++)
        if(nparticles_in_box[i] > taylor_threshold)
            taylor_boxes[ntboxes++] = i;

    /*taylor_box_nbr contains the "local numbering" in terms of 
     *taylor_boxes.*/
    for(i=0;i<ntboxes;i++)
        taylor_box_nbr[taylor_boxes[i]] = i;           

    /*localexp_c. Contains the local taylor expansions which are common
     *to all threads. Actually contains localexp1_c and localexp2_c*/
    if(ntboxes > 0)
        localexp_c = _mm_mxMalloc(4*n_terms*ntboxes*sizeof(double),16);
    else
        localexp_c = NULL;

    /*Clear the local Taylor expansions.*/
    for(i=0;i<4*n_terms*ntboxes;i++)
        localexp_c[i] = 0;
    
    /*Set up the interaction list.*/
    for(i=0;i<108;i++)  {
        interaction_list[i] = ilist_y[i]*nside+ilist_x[i];
    }
    
    /*Alloc some arrays. mpoleWorkerThd holds the thread instances and 
     *arguments holds the call arguments to the threads.*/
    mpoleWorkerThd = mxMalloc(num_threads*sizeof(THREAD_TYPE));
    arguments = mxMalloc(num_threads*sizeof(MpolesWorkerStruct));
    
    /*Keep track of the offsets of the thread's private data.*/
    double_offset = &double_data[DBL_THDATA_OFFSET];
    
    /*Fill the arguments structs and spawn the threads that are 
     *responsible for step 1 and 2 in the description above.*/
    for(i=0;i<num_threads;i++) {
        
        arguments[i].z_re = z_re;
        arguments[i].z_im = z_im;
        arguments[i].thread_data = double_offset;
        arguments[i].double_data = double_data;
        arguments[i].realloc_data = realloc_data;
        arguments[i].localexp = localexp_c;
        arguments[i].int_data = int_data;
        arguments[i].n_terms = n_terms;
        arguments[i].n_particles = n_particles;
        arguments[i].ntboxes = ntboxes;
        arguments[i].nside = nside;
        arguments[i].cursquare = &cursquare;
        arguments[i].taylor_threshold = taylor_threshold;
        
        /*Spawn thread*/
        THREAD_CREATE(mpoleWorkerThd[i],MpolesWorker,(void*) &arguments[i]);
        
        /*Step forward in the private array*/
        double_offset = &double_offset[DBL_THDATA_SIZE];
    }

    /*Wait for all threads to complete */
	for(i = 0;i<num_threads;i++)
        THREAD_JOIN(mpoleWorkerThd[i]);

    /*Start again from the first box*/
    cursquare = 0;
    /*Spawn the threads that do step 3 above.*/
    for(i=0;i<num_threads;i++)
        THREAD_CREATE(mpoleWorkerThd[i],MpolesWorkerSum,(void*) &arguments[i]);
    /*Wait for all threads to complete */
	for(i = 0;i<num_threads;i++)
        THREAD_JOIN(mpoleWorkerThd[i]);

    /*Clean up and free the memory used.*/
    _mm_mxFree(localexp_c);
    mxFree(mpoleWorkerThd);
    mxFree(arguments);
  
}
/*------------------------------------------------------------------------
 *Threaded worker function that computes multipole expansions for the
 *particles in each box and then traverses the interaction list doing one
 *of three things: either the interactions are computed directly, or the 
 *multipole expansion is evaluated directly in the target box, or the 
 *multipole expansions are converted into local taylor expansions around 
 *the centers of the target boxes and summed up.
 *------------------------------------------------------------------------
 */
THREAD_FUNC_TYPE MpolesWorker(void *argument) {
    /*Get the variables from the argument structure.*/
    MpolesWorkerStruct* arg = (MpolesWorkerStruct*) argument;
    /*The number of terms in the expansions.*/
    int n_terms = arg->n_terms;
    /*The total number of particles.*/
    int n_particles = arg->n_particles;
    /*The number of boxes on a side. That is we have nside*nside boxes in
     *total.*/
    int nside = arg->nside;
    /*The number of boxes where we transform into Taylor series.*/
    int ntboxes = arg->ntboxes;
    /*The threshold number of particles when we switch to Taylor series.*/
    int taylor_threshold = arg->taylor_threshold;

    /*The "global numbering" of the Taylor boxes.*/
    int *taylor_box_nbr = &(arg->realloc_data[TAYLOR_BOX_NBR_OFFSET]);
    /*The relative positions of the interaction list.*/
    int *interaction_list = &(arg->int_data[INTERACTION_LIST_OFFSET]);
    /*The current box to be treated. Common to all threads.*/
    int *cursquare = arg->cursquare;
    
    /*Pointers to particle position and charges.*/
    double *z_re = arg->z_re;
    double *z_im = arg->z_im;
    /*Pointers to the temporary densities.*/
    double *qa_c = &(arg->double_data[QA_C_OFFSET]);
    double *qb = &(arg->double_data[QB_OFFSET]);
    /*Pointer to the output sums.*/
    double *M1_c = &(arg->double_data[M1_C_OFFSET]);
    double *M2_c = &(arg->double_data[M2_C_OFFSET]);
    /*A binomial matrix. Used when going from multipole to Taylor.*/
    double *CC = &(arg->double_data[CC_OFFSET]);
    /*Matrices containing the Taylor series for each Taylor box. Common to
     *all threads.*/
    double *localexp1_c = arg->localexp;
    double *localexp2_c = &(arg->localexp[2*ntboxes*n_terms]);
    /*--- The following arrays are unique to each thread. ---*/
    /*The multipole expansions for the boxes.*/
    double *mpole1_c = &(arg->thread_data[MPOLE1_C_OFFSET]);
    double *mpole2_c = &(arg->thread_data[MPOLE2_C_OFFSET]);
    /*Temporary taylor arrays.*/
    double *taylorexp1_c = &(arg->thread_data[TAYLOREXP1_C_OFFSET]);
    double *taylorexp2_c = &(arg->thread_data[TAYLOREXP2_C_OFFSET]);
    /*Temporary array containing powers of 1/z.*/
    double *zsh_powers_c = &(arg->thread_data[ZSH_POWERS_C_OFFSET]);
    /*Buffer arrays for multipole evaluation. Points to the same memory
     *as taylorexp_c and zsh_powers_c, but while they're used only for
     *Taylor boxes, temp_c is used only for multipole boxes.*/
    double *temp1_c = &(arg->thread_data[TAYLOREXP1_C_OFFSET]);
    double *temp2_c = &(arg->thread_data[TAYLOREXP2_C_OFFSET]);
    
    /*The offsets of the particles in the arrays. Sorted by box.*/
    int *particle_offsets = &(arg->int_data[PARTICLE_OFFSETS_OFFSET]);
    /*The starting offsets in particle_offsets for each box.*/
    int *box_offsets = &(arg->realloc_data[BOX_OFFSETS_OFFSET]);
    /*The number of particles in each box.*/
    int *nparticles_in_box = &(arg->realloc_data[NPARTICLES_IN_BOX_OFFSET]);
    /*The interaction list.*/
    int *ilist_x = &(arg->int_data[ILIST_X_OFFSET]);
    int *ilist_y = &(arg->int_data[ILIST_Y_OFFSET]);

    double box_center_re,box_center_im;
    
    /*Loop through the boxes.*/
    for(;;) {
        int j,current_box,i;

        /*Find the next untreated box. This needs to be mutexed so we
         *don't run the risk of employing several threads to the same box.
         */
        LOCK_MUTEX(&mpoleSqTblmutex);
        i = *cursquare;
        
        while(i < nside*nside && nparticles_in_box[i] == 0)
            i++;
        
        /*If we are done we exit immediately.*/
        if(i==nside*nside) {
            UNLOCK_MUTEX(&mpoleSqTblmutex);
            break;
        }

        *cursquare = i+1;
        UNLOCK_MUTEX(&mpoleSqTblmutex);
        
        current_box = i;
        
        /*The center of the current box.*/
        box_center_re = (double)(current_box%nside-(double)(nside-1)/2)/nside;
        box_center_im = (double)(current_box/nside-(double)(nside-1)/2)/nside;

        /*We need to reset the multipole expansions for this box.*/
        for(j=0;j<2*n_terms;j++){
            mpole1_c[j] = 0;
            mpole2_c[j] = 0;
            zsh_powers_c[j] = 0;
        }
        /*Compute the multipole expansions around the center of this box
         *This loop shows up as a hotspot in Shark, so SIMD instructions
         *are used to speed it up.*/
         for(j=0;j<nparticles_in_box[current_box];j++) {
            int k;
            __m128d tz_re = _mm_set1_pd(z_re[particle_offsets[box_offsets[current_box]+j]]);
            __m128d tz_im = _mm_set1_pd(z_im[particle_offsets[box_offsets[current_box]+j]]);
            __m128d zshift_re = _mm_sub_pd(tz_re,_mm_set1_pd(box_center_re));
            __m128d zshift_im = _mm_sub_pd(tz_im,_mm_set1_pd(box_center_im));
            
            __m128d zpow1 = _mm_load_pd(&qa_c[2*particle_offsets[box_offsets[current_box]+j]]);
            __m128d zpow2 = _mm_setr_pd(0,qb[particle_offsets[box_offsets[current_box]+j]]);
            
            _mm_store_pd(mpole1_c,_mm_sub_pd(_mm_load_pd(mpole1_c),zpow1));
            _mm_store_pd(mpole2_c,_mm_add_pd(_mm_load_pd(mpole2_c),zpow2));

            for(k=2;k<2*n_terms;k+=2) {
                __m128d tmp;
                /*Multiply conj(z) by zpow1.*/
                tmp = _mm_mul_pd(zpow1, tz_re);
                tmp = _mm_addsub_pd(_mm_shuffle_pd(tmp,tmp,1),_mm_mul_pd(zpow1, tz_im));
                tmp = _mm_shuffle_pd(tmp,tmp,1);
                /*The result is in tmp, add it to zsh_powers_c*/
                _mm_store_pd(&zsh_powers_c[k],_mm_add_pd(_mm_load_pd(&zsh_powers_c[k]),tmp));
                
                /*We multiply zpow1 by zshift here.*/
                tmp = _mm_mul_pd(zpow1, zshift_re);
                zpow1 = _mm_shuffle_pd(zpow1, zpow1, 1);
                zpow1 = _mm_addsub_pd(tmp,_mm_mul_pd(zshift_im, zpow1));
                /*Subtract zpow1 from mpole1.*/
                _mm_store_pd(&mpole1_c[k],_mm_sub_pd(_mm_load_pd(&mpole1_c[k]),zpow1));
                
                /*We multiply zpow2 by zshift here.*/
                tmp = _mm_mul_pd(zpow2, zshift_re);
                zpow2 = _mm_shuffle_pd(zpow2, zpow2, 1);
                zpow2 = _mm_addsub_pd(tmp,_mm_mul_pd(zshift_im, zpow2));
                /*Add zpow2 to mpole2.*/
                _mm_store_pd(&mpole2_c[k],_mm_add_pd(_mm_load_pd(&mpole2_c[k]),zpow2));
            }
        }
        
        for(j=1;j<n_terms;j++) {
            __m128d tmp = _mm_mul_pd(_mm_load_pd(&zsh_powers_c[2*j]),_mm_set1_pd(j));
            _mm_store_pd(&mpole2_c[2*j],_mm_sub_pd(_mm_load_pd(&mpole2_c[2*j]),tmp));
        }
        
        /*Loop though the interaction list.*/
        for(j=0;j<27;j++) {
            /*Make sure that the target box is inside the computational
             *grid. The ci_flag variable reflects the fact that the inter-
             *action list of adjacent boxes may be similar if they are both
             *in a 2x2 superbox properly aligned. See Beatson and 
             *Greengard page 19. current_box%nside and current_box/nside 
             *are the box coordinates.*/
            int ci_flag = 2*((current_box%nside)&1)+((current_box/nside)&1);
         
            if((current_box%nside + ilist_x[j+ci_flag*27])<nside &&
               (current_box/nside + ilist_y[j+ci_flag*27])<nside &&
               (current_box%nside + ilist_x[j+ci_flag*27])>=0 &&
               (current_box/nside + ilist_y[j+ci_flag*27])>=0) {
                
                /*The number of the target box.*/
                int target_box = current_box+interaction_list[j+ci_flag*27];

                /*Check the number of sources of the target box.*/
                if(nparticles_in_box[target_box] > taylor_threshold) {
                    /*If larger than taylor_threshold we convert the
                     *multipole expansions into local taylor expansions 
                     *centered around the target box center.*/
                    int k;
                    double* tptr2 = &localexp1_c[2*taylor_box_nbr[target_box]*n_terms];
                    double* tptr3 = &localexp2_c[2*taylor_box_nbr[target_box]*n_terms];
                    /*zsh_re + i*zsh_im is the reciprocal of the coordinate 
                     *of the target box center.*/ 
                    double t1 = ilist_x[j+ci_flag*27]*ilist_x[j+ci_flag*27]+
                            ilist_y[j+ci_flag*27]*ilist_y[j+ci_flag*27];
                    double zsh_re = (ilist_x[j+ci_flag*27]*nside)/t1;
                    double zsh_im = -(ilist_y[j+ci_flag*27]*nside)/t1;
                    
                    /*Compute the local taylor expansions.*/
                    for(k=0;k<n_terms;k++){
                        taylorexp1_c[2*k] = mpole1_c[0]*CC[k];
                        taylorexp1_c[2*k+1] = mpole1_c[1]*CC[k];
                        taylorexp2_c[2*k] = mpole2_c[0]*CC[k];
                        taylorexp2_c[2*k+1] = mpole2_c[1]*CC[k];
                    }
                    /*zsh_powers_c is an array of increasing powers of the
                     *reciprocal zsh_re+i*zsh_im*/
                    zsh_powers_c[0] = zsh_re;
                    zsh_powers_c[1] = zsh_im;

                    for(k=1;k<n_terms;k++){
                        unsigned int l;
                        double* tptr = &CC[k*n_terms];
                        __m128d m1,m2;

                        m1 = _mm_setr_pd(zsh_powers_c[2*k-2]*mpole1_c[2*k]-zsh_powers_c[2*k-1]*mpole1_c[2*k+1],
                                         zsh_powers_c[2*k-2]*mpole1_c[2*k+1]+zsh_powers_c[2*k-1]*mpole1_c[2*k]);
                        m2 = _mm_setr_pd(zsh_powers_c[2*k-2]*mpole2_c[2*k]-zsh_powers_c[2*k-1]*mpole2_c[2*k+1],
                                         zsh_powers_c[2*k-2]*mpole2_c[2*k+1]+zsh_powers_c[2*k-1]*mpole2_c[2*k]);
                        for(l=0;l<n_terms;l++) {
                            __m128d cc = _mm_set1_pd(tptr[l]);
                            _mm_store_pd(&taylorexp1_c[2*l],_mm_add_pd(_mm_load_pd(&taylorexp1_c[2*l]),_mm_mul_pd(cc,m1)));
                            _mm_store_pd(&taylorexp2_c[2*l],_mm_add_pd(_mm_load_pd(&taylorexp2_c[2*l]),_mm_mul_pd(cc,m2)));
                        }
                        
                        zsh_powers_c[2*k] = zsh_powers_c[2*k-2]*zsh_re - zsh_powers_c[2*k-1]*zsh_im;
                        zsh_powers_c[2*k+1] = zsh_powers_c[2*k-2]*zsh_im + zsh_powers_c[2*k-1]*zsh_re;
                    }

                    /*Sum the local taylor expansions; these are common to
                     *all threads, hence the mutex. Actually, there are
                     *several mutexes(tunable), one for each group of
                     *boxes, so that we lock up as little of the localexp-
                     *matrix as possible.*/
                    LOCK_MUTEX(&localexp_mutex[target_box&(num_mutexes-1)]);    
                    for(k=0;k<2*n_terms;k+=2) {
                        __m128d zsh_re = _mm_set1_pd(zsh_powers_c[k]);
                        __m128d zsh_im = _mm_set1_pd(zsh_powers_c[k+1]);
                        /*Multiply zsh_powers_c[k] by taylorexp1_c[k]*/
                        __m128d tmp = _mm_load_pd(&taylorexp1_c[k]);
                        __m128d tmp2 = _mm_mul_pd(tmp, zsh_re);
                        tmp = _mm_shuffle_pd(tmp, tmp, 1);
                        tmp = _mm_addsub_pd(tmp2, _mm_mul_pd(zsh_im, tmp));
                        /*Store the series coefficient.*/
                        _mm_store_pd(&tptr2[k], _mm_add_pd(_mm_load_pd(&tptr2[k]), tmp));
                        /*Multiply zsh_powers_c[k] by taylorexp2_c[k]*/
                        tmp = _mm_load_pd(&taylorexp2_c[k]);
                        tmp2 = _mm_mul_pd(tmp, zsh_re);
                        tmp = _mm_shuffle_pd(tmp, tmp, 1);
                        tmp = _mm_addsub_pd(tmp2, _mm_mul_pd(zsh_im, tmp));
                        /*Store the series coefficient.*/
                        _mm_store_pd(&tptr3[k], _mm_add_pd(_mm_load_pd(&tptr3[k]), tmp));
                    }
                    UNLOCK_MUTEX(&localexp_mutex[target_box&(num_mutexes-1)]);
                    
                }else if(nparticles_in_box[target_box] > 0) {
                    
                    if(nparticles_in_box[current_box] < 12) {
                        /*There are relatively few particles in the target
                         *box, and if there are few particles in the 
                         *current box, we may as well evaluate interactions 
                         *directly.*/
                        unsigned int k;
                        int* tptr = &particle_offsets[box_offsets[current_box]]; 
                        int* tptr2 = &particle_offsets[box_offsets[target_box]]; 
                        LOCK_MUTEX(&output_mutex[target_box&(num_mutexes-1)]);
                        
                        /*This loop computes the direct interaction between
                         *all the particles in current_box and the
                         *particles in target_box.*/
                        for(k=0;k<(unsigned int)nparticles_in_box[target_box];k++) {
                            unsigned int l;
                            for(l=0;l<(unsigned int)nparticles_in_box[current_box];l++) {
                                double tz_re = z_re[tptr[l]]-z_re[tptr2[k]];
                                double tz_im = z_im[tptr[l]]-z_im[tptr2[k]];
                                double tmp = 1.0/(tz_re*tz_re+tz_im*tz_im);
                                double t1 = (qa_c[2*tptr[l]]*tz_re
                                        + qa_c[2*tptr[l]+1]*tz_im)*tmp;
                                double t2 = (qa_c[2*tptr[l]+1]*tz_re
                                        - qa_c[2*tptr[l]]*tz_im)*tmp;
                                double t3;
                                M1_c[2*tptr2[k]] += t1;
                                M1_c[2*tptr2[k]+1] += t2;
                                
                                t3 = t1*tz_re+t2*tz_im;
                                t2 = t1*tz_im-t2*tz_re-qb[tptr[l]];
                                
                                M2_c[2*tptr2[k]] += (t3*tz_re-t2*tz_im)*tmp;
                                M2_c[2*tptr2[k]+1] += (t3*tz_im+t2*tz_re)*tmp;
                                
                            }
                            
                        }
                        UNLOCK_MUTEX(&output_mutex[target_box&(num_mutexes-1)]);
                    }else{
                        /*There are too many particles in the current box,
                         *and direct evaluation would be expensive. Instead
                         *evaluate the multipole expansion at each particle
                         *in the target box.*/
                        unsigned int k;
                        int* tptr = &particle_offsets[box_offsets[target_box]];
                     
                        /*Evaluate the multipole expansions using Horner's rule.
                         *Shark hotspot, so we go SIMD.*/
                        for(k=0;k<(unsigned int)nparticles_in_box[target_box];k++) {
                            int m;
                            double t1 = (z_re[tptr[k]]-box_center_re);
                            double t2 = (z_im[tptr[k]]-box_center_im);
                            double tmpd = 1.0/(t1*t1+t2*t2);
                            double zsh_re = t1*tmpd;
                            double zsh_im = -t2*tmpd;
                            __m128d zshift_re = _mm_set1_pd(zsh_re);
                            __m128d zshift_im = _mm_set1_pd(zsh_im);
                            __m128d tmp, tmp2;
                            __m128d zpow1 = _mm_load_pd(&mpole1_c[2*n_terms-2]);
                            __m128d zpow2 = _mm_load_pd(&mpole2_c[2*n_terms-2]);
                            __m128d zpow3 = zpow1;
                            zpow3 = _mm_mul_pd(zpow3, _mm_set1_pd(n_terms));
                            
                            t1 = zsh_re*zsh_re-zsh_im*zsh_im;
                            t2 = 2*zsh_re*zsh_im;
                            /*This is the main hotspot in shark.*/
                            for(m=n_terms-2;m>=0;m--) {
                                /*We multiply zpow1 by zshift here. Do it via SSE3*/
                                __m128d tmp = _mm_mul_pd(zpow1, zshift_re);
                                zpow1 = _mm_shuffle_pd(zpow1, zpow1, 1);
                                zpow1 = _mm_addsub_pd(tmp, _mm_mul_pd(zshift_im, zpow1));
                                /*Add the multipole coefficient.*/
                                zpow1 = _mm_add_pd(_mm_load_pd(&mpole1_c[2*m]), zpow1);
                                
                                /*We multiply zpow2 by zshift here. Do it via SSE3*/
                                tmp = _mm_mul_pd(zpow2, zshift_re);
                                zpow2 = _mm_shuffle_pd(zpow2, zpow2, 1);
                                zpow2 = _mm_addsub_pd(tmp, _mm_mul_pd(zshift_im, zpow2));
                                /*Add the multipole coefficient.*/
                                zpow2 = _mm_add_pd(_mm_load_pd(&mpole2_c[2*m]), zpow2);
                                
                                /*We multiply zpow3 by zshift here. Do it via SSE3*/
                                tmp = _mm_mul_pd(zpow3, zshift_re);
                                zpow3 = _mm_shuffle_pd(zpow3, zpow3, 1);
                                zpow3 = _mm_addsub_pd(tmp, _mm_mul_pd(zshift_im, zpow3));
                                /*Add the multipole coefficient.*/
                                zpow3 = _mm_add_pd(_mm_mul_pd(_mm_load_pd(&mpole1_c[2*m]), _mm_set1_pd(m+1)), zpow3);
                                
                            }
                            /*We multiply zpow1 by zshift here. Do it via SSE3*/
                            tmp2 = _mm_mul_pd(zpow1, zshift_re);
                            zpow1 = _mm_shuffle_pd(zpow1, zpow1, 1);
                            zpow1 = _mm_addsub_pd(tmp2, _mm_mul_pd(zshift_im, zpow1));
                            
                            _mm_store_pd(&temp1_c[2*k], zpow1);
                            /*We multiply zpow2 by zshift here. Do it via SSE3*/
                            tmp2 = _mm_mul_pd(zpow2, zshift_re);
                            zpow2 = _mm_shuffle_pd(zpow2, zpow2, 1);
                            zpow2 = _mm_addsub_pd(tmp2, _mm_mul_pd(zshift_im, zpow2));

                            tmp = _mm_mul_pd(zpow3, _mm_set1_pd(t1));
                            zpow3 = _mm_shuffle_pd(zpow3, zpow3, 1);
                            zpow3 = _mm_addsub_pd(tmp, _mm_mul_pd(_mm_set1_pd(t2), zpow3));
                            
                            /*Multiply conj(zpow3) by z.*/
                            tmp = _mm_mul_pd(zpow3, _mm_set1_pd(z_re[tptr[k]]));
                            tmp = _mm_addsub_pd(_mm_mul_pd(zpow3, _mm_set1_pd(z_im[tptr[k]])), _mm_shuffle_pd(tmp, tmp, 1));
                            tmp = _mm_shuffle_pd(tmp, tmp, 1);
                            /*Add/Sub zpow2.*/
                            tmp = _mm_addsub_pd(tmp, zpow2);
                            _mm_store_pd(&temp2_c[2*k], tmp);
                        }
                        
                        /*Add the contributions, this is safe since we
                         *know that there is less than taylor_threshold
                         *particles in the target box. taylor_threshold is 
                         *n_terms/2 by default, and a value of n_terms
                         *would be very slow(and now also crash the program).*/
                        LOCK_MUTEX(&output_mutex[target_box&(num_mutexes-1)]);
                        for(k=0;k<(unsigned int)nparticles_in_box[target_box];k++) {
                            __m128d tmp = _mm_load_pd(&temp1_c[2*k]);
                            _mm_store_pd(&M1_c[2*tptr[k]], _mm_add_pd(tmp, _mm_load_pd(&M1_c[2*tptr[k]])));
                            tmp = _mm_load_pd(&temp2_c[2*k]);
                            _mm_store_pd(&M2_c[2*tptr[k]], _mm_add_pd(tmp, _mm_load_pd(&M2_c[2*tptr[k]])));
                        }
                        UNLOCK_MUTEX(&output_mutex[target_box&(num_mutexes-1)]);
                    }/*Direct or multipole in target box.*/
                } /*Taylor or Direct/Multipole in target box.*/
            }/*Box is inside the computational domain.*/
        }/*Interaction list loop.*/
    }/*Box loop.*/
    THREAD_EXIT();
}
/*------------------------------------------------------------------------
 *Threaded worker function that evaluates and sums the local taylor 
 *expansions that were computed in MpolesWorker.
 *------------------------------------------------------------------------
 */
THREAD_FUNC_TYPE MpolesWorkerSum(void* argument) {
    /*Get the variables from the argument structure.*/
    MpolesWorkerStruct* arg = (MpolesWorkerStruct*) argument;
    /*Number of terms in the expansions.*/
    int n_terms = arg->n_terms;
    /*The total number of particles.*/
    int n_particles = arg->n_particles;
    /*The number of boxes on a side. That is we have nside*nside boxes in
     *total.*/
    int nside = arg->nside;
    /*Number of Taylor boxes.*/
    int ntboxes = arg->ntboxes;
    
    /*Pointers to particle positions.*/
    double *z_re = arg->z_re,*z_im = arg->z_im;
    /*Pointer to matrix of Taylor series for boxes with enough particles.*/
    double *localexp1_c = arg->localexp;
    double *localexp2_c = &(arg->localexp[2*ntboxes*n_terms]);
    /*Pointers to the output data.*/
    double *M1_c = &(arg->double_data[M1_C_OFFSET]);
    double *M2_c = &(arg->double_data[M2_C_OFFSET]);
    
    /*The next box to be treated. Is common to all threads. Has mutex guard.*/
    int *cursquare = arg->cursquare;
    /*The "global numbering" of the boxes with Taylor series.*/
    int *taylor_boxes = &(arg->realloc_data[TAYLOR_BOXES_OFFSET]);
    /*The offsets of the particles in the arrays. Sorted by box.*/
    int *particle_offsets = &(arg->int_data[PARTICLE_OFFSETS_OFFSET]);
    /*The starting offsets in particle_offsets for each box.*/
    int *box_offsets = &(arg->realloc_data[BOX_OFFSETS_OFFSET]);
    /*The number of particles in each box.*/
    int *nparticles_in_box = &(arg->realloc_data[NPARTICLES_IN_BOX_OFFSET]);
    
    double box_center_re,box_center_im;
    
    /*Loop through the boxes with more than taylor_threshold points in them.*/
    for(;;) {
        int j,current_box,i;
        int* tptr;

        /*Find the next untreated box. This needs to be mutexed so we
         *don't run the risk of employing several threads to the same box.
         */
        LOCK_MUTEX(&mpoleSqTblmutex);
        i = *cursquare;
        /*If we are done we exit immediately.*/
        if(i==ntboxes) {
            UNLOCK_MUTEX(&mpoleSqTblmutex);
            break;
        }

        *cursquare = i+1;
        UNLOCK_MUTEX(&mpoleSqTblmutex);

        /*The number of the current box.*/
        current_box = taylor_boxes[i];
        
        /*The center of the current box.*/
        box_center_re = (double)(current_box%nside-(double)(nside-1)/2)/nside;
        box_center_im = (double)(current_box/nside-(double)(nside-1)/2)/nside;
        
        /*A temporary pointer to the pre-offsets of this box's particles.*/
        tptr = &particle_offsets[box_offsets[current_box]];

        /*Evaluate the taylor series using Horner's rule.
         *(Minor)Shark hotspot, so we go SIMD.*/
        for(j=0;j<nparticles_in_box[current_box];j++) {
            int m;
            __m128d tmp;
            __m128d zshift_re = _mm_set1_pd(z_re[tptr[j]]-box_center_re);
            __m128d zshift_im = _mm_set1_pd(z_im[tptr[j]]-box_center_im);
            
            __m128d zpow1 = _mm_load_pd(&localexp1_c[2*(n_terms-1+i*n_terms)]);
            __m128d zpow2 = _mm_load_pd(&localexp2_c[2*(n_terms-1+i*n_terms)]);
            __m128d zpow3 = zpow1;
            zpow3 = _mm_mul_pd(zpow3,_mm_set1_pd(n_terms-1));
            
            for(m=n_terms-2;m>=1;m--) {
               int tadr = 2*(m+i*n_terms);
                /*We multiply zpow1 by zshift here. Do it via SSE3*/
                __m128d tmp = _mm_mul_pd(zpow1, zshift_re);
                zpow1 = _mm_shuffle_pd(zpow1, zpow1, 1);
                zpow1 = _mm_addsub_pd(tmp, _mm_mul_pd(zshift_im, zpow1));
                /*Add the expansion coefficient.*/
                zpow1 = _mm_add_pd(_mm_load_pd(&localexp1_c[tadr]), zpow1);
 
                /*We multiply zpow2 by zshift here. Do it via SSE3*/
                tmp = _mm_mul_pd(zpow2, zshift_re);
                zpow2 = _mm_shuffle_pd(zpow2, zpow2, 1);
                zpow2 = _mm_addsub_pd(tmp, _mm_mul_pd(zshift_im, zpow2));
                /*Add the expansion coefficient.*/
                zpow2 = _mm_add_pd(_mm_load_pd(&localexp2_c[tadr]), zpow2);

                /*We multiply zpow3 by zshift here. Do it via SSE3*/
                tmp = _mm_mul_pd(zpow3, zshift_re);
                zpow3 = _mm_shuffle_pd(zpow3, zpow3, 1);
                zpow3 = _mm_addsub_pd(tmp, _mm_mul_pd(zshift_im, zpow3));
                /*Add the expansion coefficient.*/
                zpow3 = _mm_add_pd(_mm_mul_pd(_mm_set1_pd(m),_mm_load_pd(&localexp1_c[tadr])), zpow3);
            }
            tmp = _mm_mul_pd(zpow1, zshift_re);
            zpow1 = _mm_shuffle_pd(zpow1, zpow1, 1);
            zpow1 = _mm_addsub_pd(tmp, _mm_mul_pd(zshift_im, zpow1));
            /*Add the expansion coefficient.*/
            zpow1 = _mm_add_pd(_mm_load_pd(&localexp1_c[2*i*n_terms]), zpow1);
            _mm_store_pd(&M1_c[2*tptr[j]], _mm_add_pd(_mm_load_pd(&M1_c[2*tptr[j]]), zpow1));

            tmp = _mm_mul_pd(zpow2, zshift_re);
            zpow2 = _mm_shuffle_pd(zpow2, zpow2, 1);
            zpow2 = _mm_addsub_pd(tmp, _mm_mul_pd(zshift_im, zpow2));
            /*Add the expansion coefficient.*/
            zpow2 = _mm_add_pd(_mm_load_pd(&localexp2_c[2*i*n_terms]), zpow2);

            /*Multiply conj(zpow3) by z.*/
            tmp = _mm_mul_pd(zpow3, _mm_set1_pd(-z_re[tptr[j]]));
            tmp = _mm_addsub_pd(_mm_mul_pd(zpow3, _mm_set1_pd(-z_im[tptr[j]])), _mm_shuffle_pd(tmp, tmp, 1));
            tmp = _mm_shuffle_pd(tmp, tmp, 1);
            /*Add/Sub zpow2.*/
            tmp = _mm_addsub_pd(tmp, zpow2);
            _mm_store_pd(&M2_c[2*tptr[j]], _mm_add_pd(_mm_load_pd(&M2_c[2*tptr[j]]), tmp));
        }        

    }
    THREAD_EXIT();
    
}

/*------------------------------------------------------------------------
 *When we have reached the level where all boxes contain less than n_terms
 *particles, this function computes direct interaction between nearest
 *neighbors as well as box self-interaction.
 *------------------------------------------------------------------------
 */
void Direct(double *z_re, double *z_im,
            double* double_data, int* int_data, int* realloc_data, 
            int nside, int n_particles, int num_threads) {
    
    /*cursquare is the next box to be treated. Common to all threads.*/
    int i,cursquare=0;
    DirectWorkerStruct *arguments;

    /*The nearest neighbor interaction list.*/
    int ilist_x[8] = {-1,-1,-1,0,0,1,1,1};
    int ilist_y[8] = {-1,0,1,-1,1,-1,0,1};

    /*Allocate memory for the threads and thread
     *argument structs.*/
    directWorkerThd = mxMalloc(num_threads*sizeof(THREAD_TYPE));
    arguments = mxMalloc(num_threads*sizeof(DirectWorkerStruct));
    
    /*Spawn the threads*/
    for(i=0;i<num_threads;i++) {
        arguments[i].ilist_x = ilist_x;
        arguments[i].ilist_y = ilist_y;
        arguments[i].z_re = z_re;
        arguments[i].z_im = z_im;
        arguments[i].double_data = double_data;
        arguments[i].realloc_data = realloc_data;
        arguments[i].int_data = int_data;
        arguments[i].n_particles = n_particles;
        arguments[i].nside = nside;
        arguments[i].cursquare = &cursquare;

        THREAD_CREATE(directWorkerThd[i],DirectWorker,(void *) &arguments[i]);
    }

    /*Wait for all threads to complete */
	for(i = 0;i<num_threads;i++)
        THREAD_JOIN(directWorkerThd[i]);

    /*Clean up*/
    mxFree(directWorkerThd);
    mxFree(arguments);

}
/*------------------------------------------------------------------------
 *Threaded worker routine for computing the remaining interactions 
 *directly. Accesses to the output vector is completely parallel, so no
 *mutexes are necessary.
 *------------------------------------------------------------------------
 */
THREAD_FUNC_TYPE DirectWorker(void* in_struct) {
    /*Get the structure arrays and variables.*/
    DirectWorkerStruct *arg = (DirectWorkerStruct*)in_struct;
    /*The number of particles.*/
    int n_particles = arg->n_particles;
    /*The number of boxes on a side. That is we have nside*nside boxes in
     *total.*/
    int nside = arg->nside;
    
    /*Pointers to the output sums.*/
    double *M1_c = &(arg->double_data[M1_C_OFFSET]);
    double *M2_c = &(arg->double_data[M2_C_OFFSET]);
    /*Pointers to the temporary densities.*/
    double *qa_c = &(arg->double_data[QA_C_OFFSET]);
    double *qb = &(arg->double_data[QB_OFFSET]);
    
    /*Pointers to the positions and charges of the particles.*/
    double *z_re = arg->z_re;
    double *z_im = arg->z_im;
    
    /*The relative offsets of the nearest neighbors.*/
    int *ilist_x = arg->ilist_x;
    int *ilist_y = arg->ilist_y;
    
    /*The offsets of the particles in the arrays. Sorted by box.*/
    int *particle_offsets = &(arg->int_data[PARTICLE_OFFSETS_OFFSET]);
    /*The starting offsets in particle_offsets for each box.*/
    int *box_offsets = &(arg->realloc_data[BOX_OFFSETS_OFFSET]);
    /*The number of particles in each box.*/
    int *nparticles_in_box = &(arg->realloc_data[NPARTICLES_IN_BOX_OFFSET]);
    /*The next box to be treated. Common to all threads.*/
    int *cursquare = arg->cursquare;
    
    /*Loop through the boxes.*/
    for(;;) {
        int i,j,current_box;
        int *tptr;
        
        /*Find the next untreated box. This needs to be mutexed so we
         *don't run the risk of employing several threads to the same box.
         */
        LOCK_MUTEX(&directSqTblmutex);
        i = *cursquare;
        while(i < nside*nside && nparticles_in_box[i] == 0)
            i++;
        if(i==nside*nside) {
            UNLOCK_MUTEX(&directSqTblmutex);
            break;
        }
        current_box = i;
        *cursquare = i+1;
        UNLOCK_MUTEX(&directSqTblmutex);        
        
        /*Temporary pointer to the particles of current box.*/
        tptr = &particle_offsets[box_offsets[current_box]];

        /*Compute the box self-interactions.*/
        for(j=0;j<nparticles_in_box[current_box];j++) {
            int k;
            for(k=0;k<nparticles_in_box[current_box];k++) {
                if(k!=j) {
                    double tz_re = z_re[tptr[k]]-z_re[tptr[j]];
                    double tz_im = z_im[tptr[k]]-z_im[tptr[j]];
                    double tmp = 1.0/(tz_re*tz_re+tz_im*tz_im);
                    double t1 = (qa_c[2*tptr[k]]*tz_re 
                         + qa_c[2*tptr[k]+1]*tz_im)*tmp;
                    double t2 = (qa_c[2*tptr[k]+1]*tz_re 
                         - qa_c[2*tptr[k]]*tz_im)*tmp;
                    double t3;
                    M1_c[2*tptr[j]] += t1;
                    M1_c[2*tptr[j]+1] += t2;
                    t3 = t1*tz_re+t2*tz_im;
                    t2 = t1*tz_im-t2*tz_re-qb[tptr[k]];
                    M2_c[2*tptr[j]] += (t3*tz_re-t2*tz_im)*tmp;
                    M2_c[2*tptr[j]+1] += (t3*tz_im+t2*tz_re)*tmp;

                }
            }
        }
        
        /*Compute interactions from the nearest neighbors.*/
        for(j=0;j<8;j++) {
            /*Make sure that the box is inside the computational grid.*/
            if((current_box%nside+ilist_x[j]) < nside && 
               (current_box/nside+ilist_y[j]) < nside &&
                (current_box%nside+ilist_x[j])>= 0 &&
                (current_box/nside+ilist_y[j])>= 0) {

                int k;
                /*The number of the source nearest neighbor box.*/
                int source_box = current_box+ilist_y[j]*nside+ilist_x[j];
                
                if(nparticles_in_box[source_box] > 0) {
                    /*Pointer to the particles of the source box.*/
                    int* tptr2 = &particle_offsets[box_offsets[source_box]];
                    for(k=0;k<nparticles_in_box[current_box];k++) {
                        int l;
                        
                        /*This loop computes the direct interaction between
                         *all the particles in current_box and the source  
                         *particles.*/
                        for(l=0;l<nparticles_in_box[source_box];l++) {
                            double tz_re = z_re[tptr2[l]]-z_re[tptr[k]];
                            double tz_im = z_im[tptr2[l]]-z_im[tptr[k]];
                            double tmp = 1.0/(tz_re*tz_re+tz_im*tz_im);
                            double t1 = (qa_c[2*tptr2[l]]*tz_re
                                    + qa_c[2*tptr2[l]+1]*tz_im)*tmp;
                            double t2 = (qa_c[2*tptr2[l]+1]*tz_re
                                    - qa_c[2*tptr2[l]]*tz_im)*tmp;
                            double t3;
                            M1_c[2*tptr[k]] += t1;
                            M1_c[2*tptr[k]+1] += t2;
                            
                            t3 = t1*tz_re+t2*tz_im;
                            t2 = t1*tz_im-t2*tz_re-qb[tptr2[l]];
                            
                            M2_c[2*tptr[k]] += (t3*tz_re-t2*tz_im)*tmp;
                            M2_c[2*tptr[k]+1] += (t3*tz_im+t2*tz_re)*tmp;
                            
                        }
                        
                    }
                }
            }
        }
        
    }
    THREAD_EXIT();
}

/*------------------------------------------------------------------------
 *This function assigns particles to boxes on the current grid.
 *------------------------------------------------------------------------
 */
void Assign(double *z_re, double *z_im, int n_particles, int nside,
            int *maxparticles_in_box, int *int_data, int *realloc_data) {
            
    /*The offsets of the particles in the arrays. Sorted by box.*/
    int* particle_offsets = &int_data[PARTICLE_OFFSETS_OFFSET];
    /*The starting offsets in particle_offsets for each box.*/
    int* box_offsets = &realloc_data[BOX_OFFSETS_OFFSET];
    /*The number of particles in each box.*/
    int* nparticles_in_box = &realloc_data[NPARTICLES_IN_BOX_OFFSET];
    /*Temporary arrays.*/
    int* in_box = &int_data[IN_BOX_OFFSET];
    int* particlenum_in_box = &realloc_data[ASSIGN_TMP_OFFSET];
    /*The total number of boxes.*/
    int number_of_boxes = nside*nside;
    int i,j;
    
    /*Clear the nparticles_in_box array. It is realloc:d for each grid.*/
    for(j = 0;j<number_of_boxes;j++){
        nparticles_in_box[j] = 0;
    }

    /*Assign the particles to boxes.*/
    for(j = 0;j<n_particles;j++) {
        int box_x = (int)floor(nside*(z_re[j]+0.5));
        int box_y = (int)floor(nside*(z_im[j]+0.5));
        if(box_x < 0) box_x = 0;
        if(box_x >= nside) box_x = nside-1;
        if(box_y < 0) box_y = 0;
        if(box_y >= nside) box_y = nside-1;
        in_box[j] =  box_y*nside + box_x;
        nparticles_in_box[in_box[j]]++;
    }

    /*Keep track of the maximum number of particles in a box. When this
     *quantity falls below n_terms, we do one last sweep of multipole
     *calculations and then compute the rest of the interactions directly.
     */
    maxparticles_in_box[0] = -1;
    for(j=0;j<number_of_boxes;j++) {
        if(nparticles_in_box[j] > *maxparticles_in_box)
            *maxparticles_in_box = nparticles_in_box[j];
    }
    
    /*box_offsets is the offsets of the boxes in particle_offsets to be 
     *computed below. particlenum_in_box is a temporary array.*/
    box_offsets[0] = particlenum_in_box[0] = 0;
    i = 0;
    for(j=1;j<number_of_boxes+1;j++) {
        i += nparticles_in_box[j-1];
        particlenum_in_box[j] = box_offsets[j] = i;
    }
    
    /*The offsets of the particles in each box in the z and q arrays.*/
    for(j=0;j<n_particles;j++) 
        particle_offsets[particlenum_in_box[in_box[j]]++] = j;
}


