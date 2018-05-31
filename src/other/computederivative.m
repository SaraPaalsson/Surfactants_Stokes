function zp = computederivative(z,D,idx,pointsPerPanel)
%Gauss-Legendre derivative with respect to parameter.
%Assumes equisized panels from 0 to 2*pi.
zp = zeros(length(z),1);
if pointsPerPanel == 16
    for j = 1:size(idx,1)
        i1 = idx(j,1);
        i2 = idx(j,2);
        N = i2-i1+1;
        pe = 2*pi/(N/16);
        for k = 1:N/16
            zp((k-1)*16+(1:16)+i1-1) = D*z((k-1)*16+(1:16)+i1-1)*2/pe;
        end
    end
else
    for j = 1:size(idx,1)
        i1 = idx(j,1);
        i2 = idx(j,2);
        N = i2-i1+1;
        pe = 2*pi/(N/32);
        for k = 1:N/32
            zp((k-1)*32+(1:32)+i1-1) = D32*z((k-1)*32+(1:32)+i1-1)*2/pe;
        end
    end  
end