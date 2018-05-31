#include "mex.h"
#include <math.h>

#define pi 3.1415926535897932385

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    if(nrhs != 2)
        mexErrMsgTxt("mex_evalspline : Incorrect number of input parameters");

    double *c_re = mxGetPr(prhs[0]);
    double *c_im = mxGetPi(prhs[0]);
    double *ts = mxGetPr(prhs[1]);
    
    int cs = mxGetM(prhs[0]);
    int N = mxGetN(prhs[0]);
    
    int M = mxGetM(prhs[1]);
    
    plhs[0] = mxCreateDoubleMatrix(M, 1, mxCOMPLEX);
    plhs[1] = mxCreateDoubleMatrix(M, 1, mxCOMPLEX);
    plhs[2] = mxCreateDoubleMatrix(M, 1, mxCOMPLEX);
    
    double* z_re = mxGetPr(plhs[0]);    
    double* z_im = mxGetPi(plhs[0]);    
    
    double* zp_re = mxGetPr(plhs[1]);    
    double* zp_im = mxGetPi(plhs[1]);    
    
    double* zpp_re = mxGetPr(plhs[2]);    
    double* zpp_im = mxGetPi(plhs[2]);    
    double ti;
    for(int j = 0;j<M;j++) {
        double t = N*(ts[j]+pi)/pi;
        
        t = modf(t/2/N,&ti)*2*N;
       
        int k = floor(0.5*t);
        t = t-2*floor(0.5*t)-1;
        z_re[j] = c_re[k*cs];
        z_im[j] = c_im[k*cs];
        for(int l = 1;l<cs;l++) {
            z_re[j] = t*z_re[j]+c_re[l+k*cs];
            z_im[j] = t*z_im[j]+c_im[l+k*cs];            
        }
        
        zp_re[j] = (cs-1)*c_re[k*cs];
        zp_im[j] = (cs-1)*c_im[k*cs];
        for(int l = 1;l<cs-1;l++) {
            zp_re[j] = t*zp_re[j]+(cs-1-l)*c_re[l+k*cs];
            zp_im[j] = t*zp_im[j]+(cs-1-l)*c_im[l+k*cs];            
        }
        zp_re[j] *= N/pi;
        zp_im[j] *= N/pi;
        
        zpp_re[j] = (cs-1)*(cs-2)*c_re[k*cs];
        zpp_im[j] = (cs-1)*(cs-2)*c_im[k*cs];
        for(int l = 1;l<cs-2;l++) {
            zpp_re[j] = t*zpp_re[j]+(cs-1-l)*(cs-2-l)*c_re[l+k*cs];
            zpp_im[j] = t*zpp_im[j]+(cs-1-l)*(cs-2-l)*c_im[l+k*cs];            
        }
        zpp_re[j] *=N*N/pi/pi;
        zpp_im[j] *=N*N/pi/pi;
    }
    
}