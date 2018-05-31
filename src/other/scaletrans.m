function [z,zp,zpp] = scaletrans(z,zp,zpp)
xmin = min(real(z));
ymin = min(imag(z));
xmax = max(real(z));
ymax = max(imag(z));
if xmax-xmin > ymax-ymin
    z = (z-(xmax+xmin)/2-1i*(ymax+ymin)/2)/(xmax-xmin);
    zp = zp/(xmax-xmin);
    zpp = zpp/(xmax-xmin);
else
    z = (z-(xmax+xmin)/2-1i*(ymax+ymin)/2)/(ymax-ymin);
    zp = zp/(ymax-ymin);
    zpp = zpp/(ymax-ymin);
end