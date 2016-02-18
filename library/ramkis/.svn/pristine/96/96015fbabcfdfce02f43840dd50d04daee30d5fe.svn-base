function [W, H, w1, h1, iter] = shiftednmf(A, Winit, Hinit, vec_norm, normW, tol, maxiter)
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
% [W, H, iter] = nmfsh_comb_kl(A, Winit, Hinit);

if ~exist('params', 'var')
	vec_norm = 2.0;
	normW = true;
	tol = 1e-4;
	maxiter = 10000;
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
		maxiter = 10000;
	end
end

[m, n] = size(A);
W = Winit;
H = Hinit;
w1 = rand(m, 1);
h1 = rand(1, n);
Wlarge = [w1, W];
Hlarge = [h1; H];
left = Hlarge * Hlarge';
right = A * Hlarge';
for iter = 1 : maxiter
	Wlarge = anls_entry_blockpivot_precompute(left, right, Wlarge);
	w1 = Wlarge(:, 1);
	W = Wlarge(:, 2:end);
	%temp = W'*w1;
	left = Wlarge' * Wlarge;
	right = A' * Wlarge;
    Hlarge = anls_entry_blockpivot_precompute(left, right, Hlarge')';	
	%Hlarge = [h1; H];
    h1 = Hlarge(1,:);
    H = Hlarge(2:end,:);
	gradH = left*Hlarge-right';%+temp(:,ones(1,n));
	left = Hlarge * Hlarge';
	right = A * Hlarge';    
	gradWlarge = Wlarge*left-right;
	if iter == 1
		initgrad = sqrt(norm(gradWlarge(gradWlarge<=0|Wlarge>0))^2 + norm(gradH(gradH<=0|Hlarge>0))^2);
		continue;
	else
		projnorm = sqrt(norm(gradWlarge(gradWlarge<=0|Wlarge>0))^2 + norm(gradH(gradH<=0|Hlarge>0))^2);
	end
	if projnorm < tol * initgrad
		break;
	end
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
