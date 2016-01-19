function [W, H, iter] = orthogonalnmf(A, Winit, Hinit, params)
%%Orthogonal NMF with Frobenius norm:
%http://users.cis.fiu.edu/~taoli/tenure/p126-DLPH-KDD05.pdf
% Input parameters
% A: data matrix (size m*n)
% Winit: initialization of W (size m*k)
% Hinit: initialization of H (size k*n)
% (k is the number of clusters)
% params (optional)
% params.vec_norm (default=2): indicates which norm to use for the normalization of W or H,
%                              e.g. vec_norm=2 means Euclidean norm; vec_norm=0 means no normalization.
% params.normW (default=true): true if normalizing columns of W; false if normalizing rows of H.
% params.tol (default=1e-4): a positive number (<1) for checking the stopping criterion;
%                            smaller tol value means longer running time.
% params.maxiter (default=10000): maximum number of iteration times, a positive integer
%
% Output:
% W, H: result of NMF
% iter: actual number of iterations
%
% Typical usage:
% [m, n] = size(A);
% Winit = rand(m, k);
% Hinit = rand(k, n);
% [W, H, iter] = nmfsh_comb_kl(A, Winit, Hinit);

if ~exist('params', 'var')
    vec_norm = 2.0;
    normW = true;
    tol = 1e-4;
    maxiter = 200;
else
    if isfield(params, 'vec_norm')
        vec_norm = params.vec_norm;
    else
        vec_norm = 2.0;
    end
    if isfield(params, 'normW')
        normW = params.normW;
    else
        normW = true;
    end
    if isfield(params, 'tol')
        tol = params.tol;
    else
        tol = 1e-4;
    end
    if isfield(params, 'maxiter')
        maxiter = params.maxiter;
    else
        maxiter = 200;
    end
end

[m, n] = size(A);
W = Winit;
H = Hinit;
EPSILON=1e-6;
k=size(W,2);
for iter = 1 : maxiter
    %update H
    AtW = A'*W;
    HtWtW = H'*(W'*W);
    HtWtW(HtWtW <=0) = EPSILON;    
    H = H.*(AtW./(HtWtW))';    
    %update W
    AHt = A*H';
    WWtAH = W*(W'*AHt);
    WWtAH(WWtAH <=0) = EPSILON;
    W = W.*((AHt./WWtAH).^(1/2));
    %vec_norm == 0 as we normalizing to sum to 1;
    W = bsxfun(@rdivide,W,sum(W));
    if vec_norm ~= 0
        if normW
            norms = sum(W.^vec_norm) .^ (1/vec_norm);
            W = bsxfun(@rdivide, W, norms);
            H = bsxfun(@times, H, norms');
        else
            norms = sum(H.^vec_norm, 2) .^ (1/vec_norm);
            W = bsxfun(@times, W, norms');
            H = bsxfun(@rdivide, H, norms);
        end
    end    
end
end