% TEST STRESSLET

% x = ones(size(x));

% x = ones(size(x));
% x = 1i*x;
% %
x = ones(size(x));
x = x + 1i*x;

%  FMM computations. ON INTERFACE

disp('-----------------------');
disp('POINT ON INTERFACE');
disp('-----------------------');

[~,y_onint] = mex_M1M4(z_sc,zp_sc,W.*x/pi,1e-13,8,4);
y_onint = -1i*y_onint;
y_onint = y_onint + W.*x.*imag(zpp./zp)/2/pi + W.*conj(x).*imag(zpp.*conj(zp))./conj(zp).^2/2/pi;
% Computes the integrals for all interface points

% ON INTEFACE, MATLAB COMP
i_comp = 750;
y_onintMi = 0;
for j=1:length(x)
    if j ~= i_comp
       y_onintMi = y_onintMi + x(j).*W(j).*imag(zp(j)./(z(j)-z(i_comp))) + conj(x(j)).*W(j).*imag(zp(j).*conj(z(j)-z(i_comp)))./(conj(z(j)-z(i_comp))^2); 
    else
        y_onintMi = y_onintMi + x(j).*W(j).*imag(zpp(j)/(2*zp(j))) + conj(x(j)).*W(j).*imag(zpp(j)*conj(zp(j)))/(2*conj(zp(j))^2);
    end
end
y_onintMi = y_onintMi/pi;

% Compare on interface:
disp(['By FMM we have T = ' num2str(y_onint(i_comp))])
disp(['Our computations give T = ' num2str(y_onintMi)])
disp(['They differ by: ' num2str(abs(y_onint(i_comp)-y_onintMi))])
disp(' ')

% INSIDE, FMM computations
zin = 0.1;
disp('-----------------------');
disp('POINT INSIDE');
disp('-----------------------');
%
z2a = [z; zin];
W2 = [W; 0];
zp2a = [zp; 0];
zpp2a = [zpp; 0];
[z2,zp2,zpp2] = scaletrans(z2a,zp2a,zpp2a);
x2 = [x; 0];

[~,y_in] = mex_M1M4(z2,zp2,W2.*x2/pi,1e-13,8,4);
y_in = -1i*y_in;
y_in = y_in(end);
disp(['By FMM we have T = ' num2str(y_in)])

% OUR COMPUTATIONS
y_inM = 0;
for j=1:length(x)
       y_inM = y_inM + x(j).*W(j).*imag(zp(j)./(z(j)-zin)) + conj(x(j)).*W(j).*imag(zp(j).*conj(z(j)-zin))./(conj(z(j)-zin)^2);
end
y_inM = y_inM/pi;
disp(['Our computations give T = ' num2str(y_inM)])


% Compare on interface:
disp(['They differ by: ' num2str(abs(y_in-y_inM))])
disp(' ')

% OUTSIDE, FMM computations
zout = 5;

disp('-----------------------');
disp('POINT OUTSIDE');
disp('-----------------------');
%
z2a = [z; zout];
W2 = [W; 0];
zp2a = [zp; 0];
zpp2a = [zpp; 0];
[z2,zp2,zpp2] = scaletrans(z2a,zp2a,zpp2a);
x2 = [x; 0];

[~,y_out] = mex_M1M4(z2,zp2,W2.*x2/pi,1e-13,8,4);
y_out = -1i*y_out;
y_out = y_out(end);
disp(['By FMM we have T = ' num2str(y_out)])

% OUR COMPUTATIONS
y_outM = 0;
for j=1:length(x)
       y_outM = y_outM + x(j).*W(j).*imag(zp(j)./(z(j)-zout)) + conj(x(j)).*W(j).*imag(zp(j).*conj(z(j)-zout))./(conj(z(j)-zout)^2);
end
y_outM = y_outM/pi;
disp(['Our computations give T = ' num2str(y_outM)])


% Compare on interface:
disp(['They differ by: ' num2str(abs(y_onint(i_comp)-y_onintMi))])
disp(' ')
