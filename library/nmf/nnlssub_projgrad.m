function [X,grad,iter,subIter] = nnlssub_projgrad(A,B,tol,init,MAX_ITER,option)
% Nonnegative Leastsquare for Multiple Righthandside : minimize |Ax-B|_F
% Jingu Kim (jingu@cc.gatech.edu)
%
% modified from:
%   Chih-Jen Lin,
%   Projected Gradient Methods for Nonnegative Matrix Factorization
%   Neural Computation, MIT Press, 2007, 19, 2756-2779
%
% This function must be only used as a subroutine for NMF. Look for 'tol' to see why.
%
% Updated 2008.06.18
% Updated 2011.03.31: iter

[m,n]=size(A);, [m,k]=size(B);
if nargin < 6, option=0;, end
if nargin < 5, MAX_ITER = 100;, end
if nargin < 4, X=rand(n,k);, else, X = init;, end

MAX_LINE_SEARCH = 20;
alpha = 1;, beta = 0.1;, sigma=0.01;

if option
    AtA=A;, AtB=B;
else
    AtB = A'*B;
    AtA = A'*A;
end

subIter = 0;
for iter=0:MAX_ITER,  
    % convergence test
    grad = AtA*X - AtB;
    projgrad = norm(grad(grad < 0 | X >0));
    % fprintf(1,'projgrad: norm(%3.2e) , tol:%3.2e\n',projgrad,tol);
    if projgrad < tol, break, end
    
    % search step size 
    for inner_iter=1:MAX_LINE_SEARCH,
        X_next = max(X - alpha*grad, 0); d = X_next - X;
        gradd=sum(sum(grad.*d)); dQd = sum(sum((AtA*d).*d));
        suff_decr = (1-sigma)*gradd + 0.5*dQd < 0;
        if inner_iter==1,
            decr_alpha = ~suff_decr; X_prev = X;
        end
        if decr_alpha,
            if suff_decr,
                X = X_next; break;
            else
                alpha = alpha * beta;
            end
        else
            if ~suff_decr | X_prev == X_next,
                X = X_prev; break;
            else
                alpha = alpha/beta; X_prev = X_next;
            end
        end
    end
    subIter = subIter + inner_iter;
end

%if iter==MAX_ITER,
%  fprintf('Max iter in nnlsm_projgrad\n');
%end
