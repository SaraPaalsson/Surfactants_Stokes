#include "mex.h"
#include "Complex.h"
#include "pthread.h"

/*Define thread related macros and types. This is a crude wrapper to make
 *the code somewhat portable.*/
#ifdef WIN32
typedef CRITICAL_SECTION MUTEX_TYPE; 
typedef HANDLE THREAD_TYPE; 
/*Ugly!*/
#define THREAD_FUNC_TYPE DWORD WINAPI 
#define INIT_MUTEX(MUTEX) InitializeCriticalSection(MUTEX)
#define DESTROY_MUTEX(MUTEX) DeleteCriticalSection(MUTEX)
#define LOCK_MUTEX(MUTEX) EnterCriticalSection(MUTEX)
#define UNLOCK_MUTEX(MUTEX) LeaveCriticalSection(MUTEX)
#define THREAD_CREATE(THREAD,FUNC,ARGS) THREAD=CreateThread(NULL, 0, FUNC, ARGS, 0, NULL)
#define THREAD_JOIN(THREAD) WaitForSingleObject(THREAD,INFINITE)
#define THREAD_EXIT() return 0
#else
typedef pthread_mutex_t MUTEX_TYPE;
typedef pthread_t THREAD_TYPE; 
/*Ugly!*/
#define THREAD_FUNC_TYPE void*
#define INIT_MUTEX(MUTEX) pthread_mutex_init(MUTEX,NULL)
#define DESTROY_MUTEX(MUTEX) pthread_mutex_destroy(MUTEX)
#define LOCK_MUTEX(MUTEX) pthread_mutex_lock(MUTEX)
#define UNLOCK_MUTEX(MUTEX) pthread_mutex_unlock(MUTEX)
#define THREAD_CREATE(THREAD,FUNC,ARGS) pthread_create(&THREAD, NULL, FUNC, ARGS)
#define THREAD_JOIN(THREAD) pthread_join(THREAD,NULL)
#define THREAD_EXIT() pthread_exit((void*) 0)
#endif

const double pi = 3.1415926535897932385;
const int startsize = 32;
void vandernewton(Complex *T, Complex *b, int N);
void vandernewtonT(double *T, double *b, int N);

typedef struct grid_struct grid_panel;

struct grid_struct {
	int panel_nbr;
	grid_panel* next;
};

typedef struct {
    Complex *corrs,*icorrs;
    int mynumber,ncores,Ncorrs,XBoxes,YBoxes,Nbubbles,N;
    int* cidx;
    double xmin,ymin,meanlen;
    double *z_re,*z_im,*zp_re,*zp_im,*W,*pe_re,*pe_im,*bubble,*idx,*out_real,*beta;
    grid_panel** grid;
    
} WorkerStruct;

THREAD_FUNC_TYPE Worker(void* in_struct);

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
    
    double *z_re,*z_im,*zp_re,*zp_im,*W,*pe_re,*pe_im,*bubble,*idx,*beta;
    double *out_re,*out_im,*out_real,*corridx;
//     int* cidx = new int[startsize*2];
    int Nbubbles,N;
//     Complex* corrs =new Complex[startsize*32];
//     Complex* icorrs =new Complex[startsize*32];
    int num_threads;
//     int* mynumber;
    
    if(nrhs != 8)
        mexErrMsgTxt("mex_dospecquad3: incorrect number of input arguments.\n");
    z_re = mxGetPr(prhs[0]);
    z_im = mxGetPi(prhs[0]);
    N = mxGetM(prhs[0]);
    zp_re = mxGetPr(prhs[1]);
    zp_im = mxGetPi(prhs[1]);
    W = mxGetPr(prhs[2]);
    pe_re = mxGetPr(prhs[3]);
    pe_im = mxGetPi(prhs[3]);
    bubble = mxGetPr(prhs[4]);
    idx = mxGetPr(prhs[5]);
    Nbubbles = mxGetM(prhs[5]);
    beta = mxGetPr(prhs[6]);
//     omega_re = mxGetPr(prhs[6]);
//     omega_im = mxGetPi(prhs[6]);   
    num_threads = static_cast<int>(mxGetScalar(prhs[7]));
    
    plhs[2] = mxCreateDoubleMatrix(N,1,mxREAL);
    out_real = mxGetPr(plhs[2]);
    
    
    double xmax=-1e100,xmin=1e100,ymax=-1e100,ymin=1e100;
    
    for(int j = 0;j<N;j++) {
        if(z_re[j] > xmax)
            xmax = z_re[j];
        if(z_re[j] < xmin)
            xmin = z_re[j];
        if(z_im[j] > ymax)
            ymax = z_im[j];
        if(z_im[j] < ymin)
            ymin = z_im[j];
    }
    
    //Compute the mean panel length
    int ptr=0;
    double meanlen = 0;
    for(int j = 0;j<Nbubbles;j++) {
        int np = ((int)idx[j+Nbubbles]-(int)idx[j]+1)/8;
        
        for(int l = 0;l<np;l++,ptr++) 
            meanlen += abs(Complex(pe_re[ptr+j+1]-pe_re[ptr+j],pe_im[ptr+j+1]-pe_im[ptr+j]));
    }
    
    meanlen/=N/16;

    int XBoxes = (int)ceil((xmax-xmin)/meanlen);
    int YBoxes = (int)ceil((ymax-ymin)/meanlen);
    
    grid_panel** grid = new grid_panel*[XBoxes*YBoxes];
    for(int j = 0;j < XBoxes*YBoxes;j++)
        grid[j] = NULL;
    
    ptr = 0;
    for(int j = 0;j<Nbubbles;j++) {
        int np = ((int)idx[j+Nbubbles]-(int)idx[j]+1)/8;
        for(int l = 0;l<np;l++,ptr++) {
            
            Complex z1 = Complex(pe_re[ptr+j],pe_im[ptr+j]);
            Complex z2 = Complex(pe_re[ptr+j+1],pe_im[ptr+j+1]);
            Complex mid = 0.5*(z1+z2);
            Complex zrel = (mid - Complex(xmin,ymin))/meanlen;
            int midx = floor(real(zrel));
            int midy = floor(imag(zrel));
            /*Dangerous for very curvy panels.*/
            int radius = ceil(abs(z1-z2)/meanlen);
            for(int x = midx-radius;x <= midx+radius;x++)
                for(int y = midy-radius;y <= midy+radius;y++)
                    if(x >= 0 && x < XBoxes && y >= 0 && y < YBoxes) {
                        grid_panel* tmp = new grid_panel;
                        tmp->panel_nbr = ptr;
                
                        if(grid[x*YBoxes+y] == NULL) {
                            tmp->next = NULL;
                        }else{
                            tmp->next = grid[x*YBoxes+y];
                        }
                        grid[x*YBoxes+y] = tmp;
                    }
        }
    }
    
   
    
    WorkerStruct *arguments = new WorkerStruct[num_threads];
    

    THREAD_TYPE* workerThd = new THREAD_TYPE[num_threads];
    for(int j = 0;j<num_threads;j++) {
        arguments[j].z_re = z_re;
        arguments[j].z_im = z_im;
        arguments[j].zp_re = zp_re;
        arguments[j].zp_im = zp_im;
        arguments[j].W = W;
        arguments[j].pe_re = pe_re;
        arguments[j].pe_im = pe_im;
        arguments[j].bubble = bubble;
        arguments[j].beta = beta;
        arguments[j].idx = idx;
        arguments[j].grid = grid;
        arguments[j].xmin = xmin;
        arguments[j].ymin = ymin;
        arguments[j].meanlen = meanlen;
        arguments[j].mynumber = j;
        arguments[j].ncores = num_threads;
        arguments[j].XBoxes = XBoxes;
        arguments[j].YBoxes = YBoxes;
        arguments[j].Nbubbles = Nbubbles;
        arguments[j].N = N;
        arguments[j].out_real = out_real;
        arguments[j].corrs = new Complex[startsize*32];
        arguments[j].icorrs = new Complex[startsize*32];
        arguments[j].cidx = new int[startsize*32];
        arguments[j].Ncorrs = 0;
        
        /*Spawn thread*/
        THREAD_CREATE(workerThd[j],Worker,(void*) &arguments[j]);
        
    }
    
    /*Wait for all threads to complete */
	for(int j = 0;j<num_threads;j++)
        THREAD_JOIN(workerThd[j]);
   
    int totsize = 0;
    for(int j = 0;j<num_threads;j++)
        totsize += arguments[j].Ncorrs;
    
//     mexPrintf("Totsize : %d\n",totsize);
    
    plhs[0] = mxCreateDoubleMatrix(totsize,1,mxCOMPLEX);
    out_re = mxGetPr(plhs[0]);
    out_im = mxGetPi(plhs[0]);
    
    ptr = 0;
    for(int k = 0;k<num_threads;k++)
        for(int j = 0;j < arguments[k].Ncorrs;j++,ptr++) {
            out_re[ptr] = real(arguments[k].corrs[j]);
            out_im[ptr] = imag(arguments[k].corrs[j]);
        }
    

    plhs[1] = mxCreateDoubleMatrix(totsize,1,mxCOMPLEX);
    out_re = mxGetPr(plhs[1]);
    out_im = mxGetPi(plhs[1]);
    
    ptr = 0;
    for(int k = 0;k<num_threads;k++)
        for(int j = 0;j < arguments[k].Ncorrs;j++,ptr++) {
            out_re[ptr] = real(arguments[k].icorrs[j]);
            out_im[ptr] = imag(arguments[k].icorrs[j]);
        }
    
    
    plhs[3] = mxCreateDoubleMatrix(totsize/32,2,mxREAL);
    corridx = mxGetPr(plhs[3]);
    ptr = 0;
    for(int k = 0;k<num_threads;k++)
        for(int j = 0;j < arguments[k].Ncorrs/32;j++,ptr++) {
                corridx[ptr] = arguments[k].cidx[2*j];
                corridx[ptr+totsize/32] = arguments[k].cidx[2*j+1];
            }
    
    //Clean up
    for(int j = 0;j < XBoxes*YBoxes;j++) {
        grid_panel* tmpgrid = grid[j];
        while(tmpgrid != NULL) {
            grid_panel* t = tmpgrid->next;
            delete tmpgrid;
            tmpgrid = t;
        }
    }
    delete grid;
    for(int j = 0;j<num_threads;j++) {
        
        delete arguments[j].corrs;
        delete arguments[j].icorrs;
        delete arguments[j].cidx;
    }
}
    
THREAD_FUNC_TYPE Worker(void* in_struct_in) {
    WorkerStruct* in_struct = static_cast<WorkerStruct*>(in_struct_in);
    int N = in_struct->N;
    int Nbubbles = in_struct->Nbubbles;
    Complex* corrs = in_struct->corrs;
    Complex* icorrs = in_struct->icorrs;
    int* cidx = in_struct->cidx;
    double* z_re = in_struct->z_re;
    double* z_im = in_struct->z_im;
    double* zp_re = in_struct->zp_re;
    double* zp_im = in_struct->zp_im;
    double* W = in_struct->W;
    double* pe_re = in_struct->pe_re;
    double* pe_im = in_struct->pe_im;
    double* bubble = in_struct->bubble;
    double* beta = in_struct->beta;
    double* idx = in_struct->idx;
    grid_panel** grid = in_struct->grid;
    double xmin = in_struct->xmin;
    double ymin = in_struct->ymin;
    double meanlen = in_struct->meanlen;
    int mynumber = in_struct->mynumber;
    int ncores = in_struct->ncores;
//     int XBoxes = in_struct->XBoxes;
    int YBoxes = in_struct->YBoxes;
    int cptr = 0;
    int size = startsize*32;
    double* out_real = in_struct->out_real;
    
    
    for(int k = mynumber;k<N;k+=ncores) {
        Complex pq[32],old[16],nzpan[16],tz[16],tzp[16];
        double tmpT[16],tmpb[16];
        int x = floor((z_re[k]-xmin)/meanlen);
        int y = floor((z_im[k]-ymin)/meanlen);
        
        grid_panel* tmpgrid = grid[x*YBoxes+y];
        while(tmpgrid != NULL) {
            int pp = tmpgrid->panel_nbr;
            int b1 = (int)bubble[pp*16];
            Complex mid = Complex((pe_re[pp+1+b1]+pe_re[pp+b1])/2,(pe_im[pp+1+b1]+pe_im[pp+b1])/2);
            Complex len = Complex(pe_re[pp+1+b1]-pe_re[pp+b1],pe_im[pp+1+b1]-pe_im[pp+b1]);
            Complex zk = Complex(z_re[k],z_im[k]);
            //Close enough
            if(abs(zk-mid) < abs(len) ) {
                //Not close in parameter
                if(bubble[k] != bubble[pp*16] || (fabs(pp-(int)(k/16)) >1 && fabs(pp-(int)(k/16)) != ((int)idx[b1+Nbubbles]-(int)idx[b1]+1)/8-1)) {
                    Complex nz = 2*(zk-mid)/len;
                    Complex oldsum = 0;
                    pq[0] = log(1-nz)-log(-1-nz);
                    
                    for(int j = 0;j<16;j++) {
                        tz[j] = Complex(z_re[pp*16+j],z_im[pp*16+j]);
                        nzpan[j] = 2*(tz[j]-mid)/len;
                    }
                    //Is the point between the panel and the real axis?
                    if(real(nz) > -1 && real(nz) < 1) {
                        if(imag(nz) > 0) {
                            //Above the real axis, check if nz is enclosed
                            //by the panel and the real axis.
                            int furthercheck = 0;
                            for(int j = 0;j<16;j++)
                                if(imag(nzpan[j])>imag(nz)) {
                                furthercheck = 1;
                                break;
                                }
                            if(furthercheck) {
                                for(int j = 0;j<16;j++) {
                                    tmpT[j] = real(nzpan[j]);
                                    tmpb[j] = imag(nzpan[j]);
                                }
                                vandernewtonT(tmpT,tmpb,16);
                                double test = tmpb[15];
                                for(int j = 14;j>=0;j--)
                                    test = test*real(nz) + tmpb[j];
                                
                                if(test > imag(nz)) {
                                    //Yes, it is, assuming a reasonably well refined mesh
                                    //Correct the value of the integral.
                                    pq[0] += 2*pi*_i;
//                                     mexPrintf("Huhu1\n");
                                }
                            }
                        }
                        if(imag(nz) < 0) {
                            //Below the real axis, check if nz is enclosed
                            //by the panel and the real axis.
                            int furthercheck = 0;
                            for(int j = 0;j<16;j++)
                                if(imag(nzpan[j])<imag(nz)) {
                                furthercheck = 1;
                                break;
                                }
                            if(furthercheck) {
                                for(int j = 0;j<16;j++) {
                                    tmpT[j] = real(nzpan[j]);
                                    tmpb[j] = imag(nzpan[j]);
                                }
                                vandernewtonT(tmpT,tmpb,16);
                                double test = tmpb[15];
                                for(int j = 14;j>=0;j--)
                                    test = test*real(nz) + tmpb[j];
                                
                                if(test <imag(nz)) {
                                    //Yes, it is, assuming a reasonably well refined mesh
                                    //Correct the value of the integral.
                                    pq[0] -= 2*pi*_i;
//                                      mexPrintf("Huhu2\n");
                                }
                            }
                        }
                    }


                    for(int j = 0;j<16;j++) {
                        tzp[j] = Complex(zp_re[pp*16+j],zp_im[pp*16+j]);
//                         omega[j] = Complex(omega_re[pp*16+j],omega_im[pp*16+j]);
                        old[j] = W[pp*16+j]*tzp[j]/(tz[j]-zk);
                        oldsum += old[j];
                    }
                    //Does standard quadrature suffice?
                    if(abs(pq[0]-oldsum) > 1e-14) {
                        pq[16] = -1/(1-nz)-1/(1+nz);
                        
                        for(int kk = 1;kk<14;kk+=2) {
                            pq[kk] = nz*pq[kk-1] + 2.0/kk;
                            pq[kk+16] = nz*pq[kk+15] + pq[kk-1];
                            pq[kk+1] = nz*pq[kk];
                            pq[kk+17] = nz*pq[kk+16] + pq[kk];
                        }
                        pq[15] = nz*pq[14] + 2.0/15;
                        pq[31] = nz*pq[30] + pq[14];
                        out_real[k] += real(pq[0]-oldsum)/pi;
//                         Complex out = Complex(omega_re[k],omega_im[k])*real(pq[0]-oldsum);
//                         Complex out;
                        vandernewton(nzpan,pq,16);
                        
//                         for(int j = 0;j<16;j++) {
//                             out -= real(pq[j]-old[j])*omega[j];
//                             Complex tmp = conj((pq[j]-old[j])/tzp[j])*tzp[j]-
//                                     conj(pq[j+16]*2/len-old[j]/(tz[j]-zk))*(tz[j]-zk);
//                             out += tmp*conj(omega[j])/2;
//                         }
//                         allout += out;
                        cidx[cptr/16] = k;
                        cidx[cptr/16+1] = pp*16;
                        for(int j = 0;j<16;j++) {
                            corrs[cptr+j] = -real(pq[j]-old[j]);
                            Complex tmp = conj((pq[j]-old[j])/tzp[j])*tzp[j]-
                                    conj(pq[j+16]*2/len-old[j]/(tz[j]-zk))*(tz[j]-zk);
                            corrs[cptr+j+16] = tmp/2;                            
                            icorrs[cptr+j] = beta[k]*imag(pq[j]-old[j]);
			    tmp = -_i*(conj((pq[j]-old[j])/tzp[j])*tzp[j]-
                                    conj(pq[j+16]*2/len-old[j]/(tz[j]-zk))*(tz[j]-zk));
                            
                            //tmp = conj((pq[j]-old[j])/tzp[j])*tzp[j]+
                            //        conj(pq[j+16]*2/len-old[j]/(tz[j]-zk))*(tz[j]-zk);
                            icorrs[cptr+j+16] = beta[k]*tmp/2;                            
                        }
                        
//                         nmodifs[0]++;
                        cptr += 32;
                        if(cptr >= size) {
                            //Need to increase size of vectors.
                            size *= 2;
                            if(size > 10000000)
                                break;
                            Complex* newcorr = new Complex[size];
                            Complex* newicorr = new Complex[size];
                            int* newcidx = new int[size/16];
                            for(int j=0;j<cptr;j++) 
                                newcorr[j] = corrs[j];
                            for(int j=0;j<cptr;j++) 
                                newicorr[j] = icorrs[j];
                            for(int j=0;j<cptr/16;j++)
                                newcidx[j] = cidx[j];
                            delete corrs;
                            delete icorrs;
                            delete cidx;
                            corrs = newcorr;
                            icorrs = newicorr;
                            cidx = newcidx;
                        }
                    }
                }
            }
            tmpgrid = tmpgrid->next;
        }
        if(size > 10000000)
            break;
//         out_re[k] = real(allout)/pi;
//         out_im[k] = imag(allout)/pi;
    }
    
    in_struct->corrs = corrs;
    in_struct->icorrs = icorrs;
    in_struct->cidx =cidx;
    in_struct->Ncorrs = cptr;
    
    THREAD_EXIT();
}


void vandernewton(Complex *T, Complex *b, int N) {
	for(int k = 1;k < N;k++) 
		for(int j = N-1;j>=k;j--) {
			b[j] -= T[k-1]*b[j-1];
            b[j+16] -= T[k-1]*b[j+15];
        }
	
	for(int k = N-1;k >= 1;k--) 
		for(int j=k;j<N;j++) {
			b[j] /= T[j]-T[j-k];
            b[j+16] /= T[j]-T[j-k];
			b[j-1] -= b[j];
            b[j+15] -= b[j+16];
		}
}

void vandernewtonT(double *T, double *b, int N) {
	for(int k = 0;k < N-1;k++) 
		for(int j = N-1;j>=k+1;j--)
			b[j] = (b[j]-b[j-1])/(T[j]-T[j-k-1]);
	
	for(int k = N-1;k >= 0;k--) 
		for(int j=k;j<N-1;j++) 
			b[j] -= T[k]*b[j+1];
}
