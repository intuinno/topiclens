function [X,grad,iter] = nnlssub_projnewton_mod(A, B, tol, init, MAX_ITER,option) 
% Nonnegative Leastsquare for Multiple Righthandside : minimize |Ax-B|_F
% Jingu Kim (jingu@cc.gatech.edu)
% last updated 2008.06.27
% modified from:
%   Dongmin Kim and Suvrit Sra and Inderjit S. Dhillon
%   Fast Newton-type Methods for the Least Squares Nonnegative Matrix
%   Approximation Problem
%   Proceedings of the 2007 SIAM International Conference on Data Mining,
%   2007
%
% This function must be only used as a subroutine for NMF. Look for 'tol' to see why.
    [m,n]=size(A);
    [m,k]=size(B);
	if nargin < 6, option = 1;, end
	if nargin < 5, MAX_ITER = 100;, end
	if nargin < 4, X=rand(n,k);, else, X = init;, end

    grad = zeros(n,k);
    iter = 0;
    
	if option
    	AtA = A; AtB = B;
	else
    	AtA = A'*A; AtB = A'*B;
	end
    
    subtol = tol/sqrt(k);
    for i=1:k
        [X(:,i),grad(:,i),column_iter]=nnls1_sub_mod(AtA, AtB(:,i), subtol,X(:,i),MAX_ITER);
        iter = iter + column_iter;
    end
end

function [x,grad,iter] = nnls1_sub_mod(AtA, Atb, tol, init, MAX_ITER) 
    [temp,n] = size(AtA);
    if nargin < 4, x=rand(n,1);, else, x = init;, end

    AtAG = AtA;
    AtbG = Atb;
    grad = AtAG * x - AtbG;  % gradient
    invHess = eye(n);
    
    % convergence test
    projgrad = norm(grad(grad < 0 | x >0));
    if projgrad < tol, 
        iter = 0;
        return
    end
    
    for iter=1:MAX_ITER
        xOld = x;
        gp = find((grad > 0) & ~x);
        grad(gp) = 0;                  % compute projected gradient
        srchDir = -invHess * grad;
        ns = find((srchDir < 0) & ~x);
        srchDir(gp) = 0;               % do not update fixed variables
        srchDir(ns) = 0;               % do not go below zero

        % limited minimization rule
        % A = Aglb;
        % A(:,gp) = 0;
        % Ad = A * srchDir;
        % alphaLnSrch = -Ad' * ((A * x) - b) / (Ad' * Ad);
        AtA = AtAG;, AtA(gp,:) = 0; AtA(:,gp) = 0;
        Atb = AtbG;, Atb(gp) = 0;
        alphaLnSrch = (- srchDir' * AtA * x + srchDir' * Atb) / (srchDir' * AtA * srchDir);
        x = x + alphaLnSrch * srchDir;
        x(x < 0) = 0;

        % convergence test
        grad = AtAG * x - AtbG;
        projgrad = norm(grad(grad < 0 | x >0));
        if projgrad < tol, break, end
        
        % update the quasi-Newton approximate inverse Hessian matrix
        % Au = A * deltaX;
        % nAu = norm(Au)^2;
        % AtAu = A' * Au;
        % DAtAu = invHess * AtAu;
        % invHess = invHess + ((1 + AtAu' * DAtAu / nAu) * ...
        %                    deltaX * deltaX' - ...
        %                    (DAtAu * deltaX' + deltaX * DAtAu')) / nAu;
        dx = x - xOld;
        nAu = dx' * AtA * dx;
        AtAu = AtA * dx;
        DAtAu = invHess * AtAu;
        Ddx = DAtAu * dx';
        invHess = invHess + ((1 + AtAu' * DAtAu / nAu) * dx * dx' - (Ddx + Ddx')) / nAu;
    end
end
