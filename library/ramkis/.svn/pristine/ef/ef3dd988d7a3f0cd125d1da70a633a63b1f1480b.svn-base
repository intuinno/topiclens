function [W, H, iter] = nmfsh_comb_kl(A, Winit, Hinit, params)
%%NMF with KL-divergence
% min (summing over i=1~m, j=1~n) [(WH)_ij - A_ij*log(WH)_ij] subject to W>=0,H>=0
% (Here, the algorithm of choice is multiplicative update rule
% in order to produce consistent experiment results with Brunet et al. PNAS 2004).
%
% Input parameters
% A: data matrix (size m*n)
% Winit: initialization of W (size m*k)
% Hinit: initialization of H (size k*n)
% (k is the number of clusters)
% params (optional)
% params.vec_norm (default=0): indicates which norm to use for the normalization of W or H,
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
	vec_norm = 0;
	normW = true;
	tol = 1e-4;
	maxiter = 10000;
else
	if isfield(params, 'vec_norm')
		vec_norm = params.vec_norm;
	else
		vec_norm = 0;
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

Htrace = zeros(1, maxiter);
Hprev = zeros(size(H));
for iter = 1 : maxiter
	W = anls_entry_kl(H', A', W);
	H = anls_entry_kl(W, A, H')';
	Htrace(iter) = norm(H - Hprev, 'fro') / norm(H, 'fro');
	if (Htrace(iter) < tol)
		break;
	end
	Hprev = H;
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
