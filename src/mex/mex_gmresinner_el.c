#include "mex.h"
#include <math.h>

double mynorm(double* in_re,double* in_im,int N)
{
    int i;
    double t=0;
    for(i=0;i<N;i++)
        t += in_re[i]*in_re[i]+in_im[i]*in_im[i];
    return sqrt(t);
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    double *rhs_re,*rhs_im;
    
    register int i,j,k;
    int err_code;
    
    double *tptr_re,*tptr_im,*output_re,*output_im;
    double a,b,c,d,error,bnrm2,tol;
    double *cs,*sn,*H,*s,*w_re,*w_im,*Vin_re,*Vin_im;
    double **V_re,**V_im;
    double *numiter;
    double *counter,*relerr;
    mxArray *mxVin,*mxCounter,*mxRelerr,*fhandle,*cellptr;
    
   
    mxArray *mxrhs[3],*mxlhs[1],*mxPrintlhs[3];
    
    int M,N,HM;
/*
    //Check for incorrect number of in and out arguments.
*/
    if (nrhs != 7) {
        mexErrMsgTxt("\nmex_gmresinner_el : Incorrect number of input arguments.");
    } else if (nlhs != 2) {
        mexErrMsgTxt("\nmex_gmresinner_el : Incorrect number of output arguments.");
    }

    M = floor(mxGetScalar(prhs[0]));
    N = mxGetM(prhs[1]);

    rhs_re = mxGetPr(prhs[1]);
    rhs_im = mxGetPi(prhs[1]);

    if(rhs_re == NULL)
        rhs_re = mxCalloc(N,sizeof(double));
    if(rhs_im == NULL)
        rhs_im = mxCalloc(N,sizeof(double));

    
    tol = mxGetScalar(prhs[2]);
    fhandle = (mxArray*)prhs[5];
    cellptr = (mxArray*)prhs[6];

    mxrhs[0] = fhandle;
    mxrhs[2] = cellptr;
    
/*
    //Matlab-variables for printout.
*/
    mxCounter = mxCreateDoubleMatrix(1,1, mxREAL);
    mxRelerr = mxCreateDoubleMatrix(1,1, mxREAL);
    
    mxPrintlhs[0] = (mxArray*)prhs[4];
    mxPrintlhs[1] = mxCounter;
    mxPrintlhs[2] = mxRelerr;

/*
    //mxPrintlhs[0] = mxCounter;
    //mxPrintlhs[1] = mxRelerr;
*/

    counter = mxGetPr(mxPrintlhs[1]);
    relerr = mxGetPr(mxPrintlhs[2]);

    H = mxCalloc((M+1)*M,sizeof(double));
    
    mxVin = mxCreateDoubleMatrix(N,1, mxCOMPLEX);
    Vin_re = mxGetPr(mxVin);
    Vin_im = mxGetPi(mxVin);
    
    cs = mxCalloc(M,sizeof(double));
    sn = mxCalloc(M,sizeof(double));
    s = mxCalloc((M+1),sizeof(double));
    
    V_re = mxCalloc(M+1,sizeof(double*));
    V_im = mxCalloc(M+1,sizeof(double*));
    
    V_re[0] = mxCalloc(N,sizeof(double));
    V_im[0] = mxCalloc(N,sizeof(double));
    
/*
    //Calculate original rhs-norm in the driver method to
    //fix problems involving initial guesses.
    //bnrm2 = mynorm(rhs,N);
*/
    bnrm2 = mxGetScalar(prhs[3]);
    
    tptr_re = V_re[0];
    tptr_im = V_im[0];
    for(i=0;i<N;i++) {
        tptr_re[i] = rhs_re[i]/bnrm2;
        tptr_im[i] = rhs_im[i]/bnrm2;
    }    
    s[0] = bnrm2;

    HM = M+1;

    for(i=0;i<M;i++) {
        
        tptr_re = V_re[i];
        tptr_im = V_im[i];
        for(k = 0;k<N;k++) {
            Vin_re[k] = tptr_re[k];
            Vin_im[k] = tptr_im[k];
        }
        mxrhs[1] = mxVin;

/*
        //Perform call to multiplication function.
*/
        err_code = mexCallMATLAB(1,mxlhs,3,mxrhs, "feval");

        w_re = mxGetPr(mxlhs[0]);
        w_im = mxGetPi(mxlhs[0]);
        
        if(w_re == NULL)
            w_re = mxCalloc(N,sizeof(double));
        if(w_im == NULL)
            w_im = mxCalloc(N,sizeof(double));
        
        for(k=0;k<=i;k++) {
            tptr_re = V_re[k];
            tptr_im = V_im[k];
            a = 0;
            for(j=0;j<N;j++) 
                a += w_re[j]*tptr_re[j] + w_im[j]*tptr_im[j];
            H[i*(HM) + k] = a;
            for(j=0;j<N;j++) {
                w_re[j] -= a*tptr_re[j];
                w_im[j] -= a*tptr_im[j];
            }
        }
 
        H[i*(HM+1)]+=(double)1;

        c = mynorm(w_re,w_im,N);
        H[i*(HM+1)+ 1] = c;
            
        V_re[i+1] = mxMalloc(N*sizeof(double));
        V_im[i+1] = mxMalloc(N*sizeof(double));

        tptr_re = V_re[i+1];
        tptr_im = V_im[i+1];
        a = (double)1/c;
        for(k = 0;k<N;k++) {
            tptr_re[k] = w_re[k]*a;
            tptr_im[k] = w_im[k]*a;
        }

        mxDestroyArray(mxlhs[0]);

        for(k = 0;k<=i-1;k++) {
            d = cs[k]*H[i*(HM)+k] + sn[k]*H[i*(HM)+k+1];
            H[i*(HM)+k+1] = -sn[k]*H[i*(HM)+k] + cs[k]*H[i*(HM)+k+1];
            H[i*(HM)+k] = d;
        }

        a = H[i*(HM+1)];
        b = H[i*(HM)+i+1];

        if(b==0) {
            cs[i] = 1;
            sn[i] = 0;
            
        }else
            if(fabs(b) > fabs(a)) {
                c = a/b;
                sn[i] = ((double)1)/sqrt(1+c*c);
                cs[i] = c*sn[i];
            }else {
                c = b/a;
                cs[i] = ((double)1)/sqrt(1+c*c);
                sn[i] = c*cs[i];
        }

        H[i*(HM+1)] = cs[i]*H[i*(HM+1)] + sn[i]*H[i*(HM)+i+1];
        
        s[i+1] = -(sn[i]*s[i]);
        s[i] = cs[i]*s[i];
  
        error = fabs(s[i+1])/bnrm2;
/*
//        error = fabs(s[i+1]);
*/

        relerr[0] = error;
        counter[0] = i;
/*
//        err_code = mexCallMATLAB(0,NULL,2,mxPrintlhs, "do_printout");
*/
        err_code = mexCallMATLAB(0,NULL,3,mxPrintlhs, "feval");
        if(error<=tol) {
                        
            plhs[0] = mxCreateDoubleMatrix(N,1, mxCOMPLEX);
            output_re = mxGetPr(plhs[0]);
            output_im = mxGetPi(plhs[0]);

            plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
            numiter = mxGetPr(plhs[1]);

            numiter[0] = i+1;

            for(j=i;j>=0;j--) {
                s[j] /= H[j*(HM+1)];
                a = s[j];
                for(k=0;k<j;k++)
                    s[k] -= H[j*HM+k]*a;
            }

            for(j=0;j<=i;j++) {
                tptr_re = V_re[j];
                tptr_im = V_im[j];
                for(k=0;k<N;k++) {
                    output_re[k] += tptr_re[k]*s[j];
                    output_im[k] += tptr_im[k]*s[j];
                }
                    
            }
            return;
        }
        
    }
    i--;
    plhs[0] = mxCreateDoubleMatrix(N,1, mxCOMPLEX);
    output_re = mxGetPr(plhs[0]);
    output_im = mxGetPi(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
    numiter = mxGetPr(plhs[1]);
    
    numiter[0] = i+1;
    
    for(j=i;j>=0;j--) {
        s[j] /= H[j*(HM+1)];
        a = s[j];
        for(k=0;k<j;k++)
            s[k] -= H[j*HM+k]*a;
    }
    
    for(j=0;j<=i;j++) {
        tptr_re = V_re[j];
        tptr_im = V_im[j];
        for(k=0;k<N;k++) {
            output_re[k] += tptr_re[k]*s[j];
            output_im[k] += tptr_im[k]*s[j];
        }
        
    }
    mexWarnMsgTxt("\nmex_gmresinner_el : did not converge.");
}
