#include "mex.h"
#include <math.h>
#define pi 3.1415926535897932385

#define spr 12
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
 
    int M,Mr,j,Mtarg;
    double tau;
    double *fmtau_re,*fmtaup_re,*fmtaupp_re;
    /*
     double *fmtaup_im,*fmtau_im,*fmtaupp_im;
     double *z_im,*zp_im,*zpp_im;
    */
    double *z_re,*zp_re,*zpp_re;
    double *T,*E3;
    
    Mr = mxGetM(prhs[0]);
    Mtarg = mxGetM(prhs[3]);
    
    tau = 40.0/Mr/Mr;
    fmtau_re = mxGetPr(prhs[0]);
    fmtaup_re = mxGetPr(prhs[1]);
    fmtaupp_re = mxGetPr(prhs[2]);
    /*
    fmtau_im = mxGetPi(prhs[0]);  
    fmtaup_im = mxGetPi(prhs[1]);
    fmtaupp_im = mxGetPi(prhs[2]);
    */
    T = mxGetPr(prhs[3]);
    
    
    plhs[0] = mxCreateDoubleMatrix(Mtarg, 1, mxCOMPLEX);
    z_re = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(Mtarg, 1, mxCOMPLEX);
    zp_re = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(Mtarg, 1, mxCOMPLEX);
    zpp_re = mxGetPr(plhs[2]);
    /*
    z_im = mxGetPi(plhs[0]); 
    zp_im = mxGetPi(plhs[1]); 
    zpp_im = mxGetPi(plhs[2]); 
    */
    /*
    if (fmtau_im != NULL) 
    {
    } else { 
        fmtau_im = mxCalloc(Mr,sizeof(double));
    }

    if (fmtaup_im != NULL) 
    {
    } else {
        fmtaup_im = mxCalloc(Mr,sizeof(double));
    }
    
    if (fmtaupp_im != NULL) 
    {
     mexPrintf("fmtaupp_im imaginary\n");   
    } else {
        fmtaupp_im = mxCalloc(Mr,sizeof(double));
    }
     */
    
    
     E3 = mxCalloc(2*spr+2,sizeof(double)); 
/* RICKARD COMMENTED
//    E3 = exp(-(pi*(-spr-1:spr)'/Mr).^2/tau);
*/
    
    
    for(j = 0;j<2*spr+2;j++) 
        E3[j] = exp(-(pi*(-spr-1+j)/Mr)*(pi*(-spr-1+j)/Mr)/tau);
     

    for(j = 0;j<Mtarg;j++) {
       
        int idx = (int)ceil(Mr*T[j]/2.0/pi);
        int fidx = (Mr+idx-spr-1)%Mr;
        int k;
        double pos = 2.0*pi*idx/Mr;
        double E1 = exp(-(T[j]-pos)*(T[j]-pos)/4.0/tau)/Mr;
        double E2 = exp((T[j]-pos)*pi/Mr/tau);
        double E4 = exp((-spr-1)*(T[j]-pos)*pi/Mr/tau);

        for(k=0;k<2*spr+2;k++) {
            double E5 = E4*E3[k];
            z_re[j] += fmtau_re[fidx]*E5;
            zp_re[j] += fmtaup_re[fidx]*E5;
            zpp_re[j] += fmtaupp_re[fidx]*E5;
            /*
            z_im[j] += fmtau_im[fidx]*E5; 
            zp_im[j] += fmtaup_im[fidx]*E5;
            zpp_im[j] += fmtaupp_im[fidx]*E5; 
            */
            fidx = (fidx+1)%Mr;
            E4 *= E2;
            
        }

        z_re[j] *= E1;
        zp_re[j] *= E1;
        zpp_re[j] *= E1;
        /*
        z_im[j] *= E1; 
        zp_im[j] *= E1;
        zpp_im[j] *= E1;
         */ 
    }
     
/* RIKARD COMMENTED  
//      idx = floor(Mr*T/2/pi);
//     zp = zeros(size(z));
// zpp = zeros(size(z));
// for j = 1:M
//     fidx = mod(Mr+(idx(j)-spr-1:idx(j)+spr)',Mr)+1;
//     E4 = E2_1(j).*E2_2(j).^(0:2*spr+1)'.*E3;
//     z(j) = E1(j)*sum(fmtau(fidx).*E4)/Mr;
//     zp(j) = E1(j)*sum(fmtaup(fidx).*E4)/Mr;
//     zpp(j) = E1(j)*sum(fmtaupp(fidx).*E4)/Mr;
// end
*/


    
}
