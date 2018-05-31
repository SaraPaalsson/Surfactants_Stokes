#include "mex.h"

const double pi = 3.1415926535897932385;

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
    
    double* tar_re = mxGetPr(prhs[0]);
    double* tar_im = mxGetPi(prhs[0]);
    double* src_re = mxGetPr(prhs[1]);
    double* src_im = mxGetPi(prhs[1]);
    double* mod_re = mxGetPr(prhs[2]);
    double* mod_im = mxGetPi(prhs[2]);
    double* cidx = mxGetPr(prhs[3]);
    int N = mxGetM(prhs[3]);
    int cptr = 0;
    for(int j = 0;j<N;j++) {
        int kidx = static_cast<int>(cidx[j]);
        int jidx = static_cast<int>(cidx[j+N]);
        double corr_re=0;
        double corr_im=0;
        for(int k = 0;k<16;k++) {
            double c = src_re[jidx+k];
            double d = src_im[jidx+k];
            double a1 = mod_re[cptr+k];
            double b1 = mod_im[cptr+k];
            double a2 = mod_re[cptr+k+16];
            double b2 = mod_im[cptr+k+16];
            corr_re += a1*c-b1*d + a2*c+b2*d;
            corr_im += b1*c+a1*d + b2*c-a2*d;
        }
        tar_re[kidx] += corr_re/pi;
        tar_im[kidx] += corr_im/pi;
        cptr += 32;
    }
    
}