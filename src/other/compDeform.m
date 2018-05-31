function D = compDeform(z)

x = real(z);
y = imag(z);
alpha = linspace(0,2*pi,length(x))';

[~,Ix] = max(x);
if Ix == 1 || Ix == length(x)
    xvec(:,1) = [x(end); x(1:2)];
    xvec(:,2) = [alpha(end-1)-2*pi; alpha(1:2)];
else
    xvec(:,1) = x(Ix-1:Ix+1);
    xvec(:,2) = alpha(Ix-1:Ix+1);
end
px = polyfit(xvec(:,2),xvec(:,1),2);
nux = linspace(min(xvec(:,2)),max(xvec(:,2)))';
xval = polyval(px,nux);
xmax = max(xval);
Rmax_interp = xmax;
[~,Iy] = max(y);
if Iy == 1 || Iy == length(y)
    yvec(:,1) = [y(end); y(1:2)];
    yvec(:,2) = [alpha(end-1)-2*pi; alpha(1:2)];
else
    yvec(:,1) = y(Iy-1:Iy+1);
    yvec(:,2) = alpha(Iy-1:Iy+1);
end
py = polyfit(yvec(:,2),yvec(:,1),2);
nuy = linspace(min(yvec(:,2)),max(yvec(:,2)))';
yval = polyval(py,nuy);
ymax = max(yval);
Rmin_interp = ymax;

D = (Rmax_interp-Rmin_interp)/(Rmax_interp+Rmin_interp);


end