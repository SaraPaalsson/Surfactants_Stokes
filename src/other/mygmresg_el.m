function [x,trueres,i,iterhist] = mygmresg_el(rhs,initial,tol,m,p_style,f_handle,varargin)
%[x,trueres,i,iterhist] = mygmresg_el(rhs,initial,tol,m,p_style,f_handle,varargin)
%
%REMARK! Special version for elasticity. Will work for purely real systems
%but may be slower.
%Solves (I+K)x=rhs using GMRES. The matrix-vector product K*x
%is computed in the function f_handle, which can have any number of
%input parameters - all of which are passed at the end of the parameter
%list.
%
%initial - an initial guess for the solution.
%tol - the desired error tolerance
%m - the maximum number of iterations. If a solution has not been found
%    inside m iterations, an error is raised.
%p_style - printout style. Valid values are :
%          0 - no printout.
%          1 - printout on one updated line only.
%          2 - printout on one line per iteration.
%f_handle - the function returning K*x.

%In order to easily change printout style, a global variable is used.
global printstyle;
global tmpiterhist;

printstyle = p_style;
tmpiterhist = [];

%Compute the norm of the right hand side.
bnrm2=norm(rhs);

    if bnrm2==0
%         fprintf('rhs equals zero');
        x = rhs;
        trueres = 0;
        i = 0;
        iterhist = 0;
        return;
    end


%Print start of iteration message if one line printout is requested.

tic
if norm(initial) == 0
    if printstyle == 1
        fprintf('GMRES pseudo-residual at iteration ');
    end

    %Perform MEX-call in the no-initial-guess-case
    [x,i] = mex_gmresinner_el(m,rhs,...
        tol,bnrm2,@do_printout,f_handle,varargin);
else
%     fprintf('Initial residual norm : %e\n',norm(rhs - feval(f_handle,initial,varargin) - initial)/bnrm2)

    if printstyle == 1
        fprintf('GMRES pseudo-residual at iteration ');
    end

    %Perform MEX-call with initial guess
 
    %Add initial guess to get proper solution of the system.
    x = x + initial;
end
%Compute the relative residual.
trueres=norm(rhs-x-feval(f_handle,x,varargin))/bnrm2;

%Print stats.
if printstyle == 1
    fprintf('\n');
end
if printstyle ~= 0
    fprintf('System solved in %f sec. Relres : %e, Iter : %d\n',toc,trueres,i);
end
iterhist = tmpiterhist;




function do_printout(i,relerr)

global printstyle
global tmpiterhist

tmpiterhist(i+1) = relerr;

if printstyle == 0
    return
end

if printstyle == 1
    if i == 0
        fprintf('%d : %e',i+1,relerr);
    else
        iterl = floor(log10(i));
        for j = 1:(16+iterl)
            fprintf('\b');
        end
        fprintf('%d : %e',i+1,relerr);
    end
    return
end

if printstyle == 2
    fprintf('GMRES pseudo-residual at iteration %d : %e\n',i+1,relerr);
    return
end
