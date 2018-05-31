function [y] = gmres_KX(x,A,wt)
% Computes (I+A)x for
% wt = 1 : MATLAB's gmres solver
% else Ax for Rickard's gmres

    if wt == 1
        y = x + A*x;
    else 
        y = A*x;
    end
        
        
end

