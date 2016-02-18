function [W, H, w1, h1, iter] = noveltynmf(A, Winit, Hinit, l, vec_norm, normW, tol, maxiter)
%%Affine NMF with Frobenius norm
% min ||A-[w1,W][e^T;H]||_F^2 subject to w1>=0,W>=0,H>=0
%
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
% [W, H, iter] = nmfsh_comb_mixture(A, Winit, Hinit);

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
%w1 = full(mean(A,2));
%w1 = w1/norm(w1);
%h1 = rand(1, n);
if (exist('l','var'))
    [w1,h1]=nmfsh_comb_fro(A,rand(m,l),rand(l,n));
else
    l = 1;
    w1 = full(mean(A,2));
    h1 = rand(1, n);    
end
Wlarge = [w1, W];
Hlarge = [h1; H];
left = Wlarge' * Wlarge;
right = A' * Wlarge;
for iter = 1 : maxiter
    Hlarge = anls_entry_blockpivot_precompute(left, right, Hlarge')';
    h1 = Hlarge(1:l,:);
    H = Hlarge(l+1:end,:);	
	left = H * H';
    %H(2:end, :) * A'+ gamma * Z_hat' - H(2:end, :) * H(1,:)' * avg_vector?;
	right = A*H' - w1 * (h1 * H');    
	W = anls_entry_blockpivot_precompute(left, right, W);
	Wlarge = [w1 W];
    gradW = W*left-right;
    left = Wlarge' * Wlarge;
    right = A' * Wlarge;
	gradHlarge = left*Hlarge-right';
	if iter == 1
		initgrad = sqrt(norm(gradW(gradW<=0|W>0))^2 + norm(gradHlarge(gradHlarge<=0|Hlarge>0))^2);
		continue;
	else
		projnorm = sqrt(norm(gradW(gradW<=0|W>0))^2 + norm(gradHlarge(gradHlarge<=0|Hlarge>0))^2);
	end
	if projnorm < tol * initgrad
		break;
    end
    [iter projnorm]
end

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