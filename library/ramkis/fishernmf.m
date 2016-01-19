function [W, H, iter] = fishernmf(A, Winit, Hinit, params)
%%Standard NMF with Frobenius norm:
% min ||A-WH||_F^2 + alpha*S_w - alpha*S_b subject to W>=0,H>=0
% http://www.cs.ucsb.edu/~mturk/pubs/ACCVa2004.pdf
% Implementation follows https://github.com/marinkaz/nimfa/blob/master/nimfa/methods/factorization/lfnmf.py
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
    maxiter = 50;
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
        maxiter = 50;
    end
end

[m, n] = size(A);

W = Winit;
H = Hinit;
c = size(Winit,2);
threshold=0.2;
EPSILON=1e-6;
for iter = 1 : maxiter
    %update H
    tic
    H_temp = H;
    H_temp(H_temp<threshold)=0;
    [n_i,avgs,classlabels]=encoding(H_temp);
    for k=1:c
        for l=1:n
            currentClassLabel=classlabels(l);
            currentn_i=n_i(k);
            muki = avgs(k,currentClassLabel);
            hkl = H(k,l);
            b = 4/(currentn_i*c*(c-1))*sum(avgs(:,k)-(muki-hkl/currentn_i))-(2/(currentn_i*c))*muki+1;
            Whl = W*H(:,l);
            Whl(Whl==0)=EPSILON;
            firsttermsqrt=A(:,l)'*(W(:,k).*(hkl./Whl));
            %for i=1:n
            %    firsttermsqrt = firsttermsqrt + A(i,l)*(W(i,k)*hkl/Whl(i));
            %end            
            secondtermsqrt = 2/(currentn_i*c)-4/(currentn_i.^2*(c-1));
            H(k,l)=-b + sqrt(b.^2+4*firsttermsqrt*secondtermsqrt);
        end
    end
    for k=1:m
        for l=1:c                        
            %for j = 1:n
                %w_1 = sum(self.H[k, j] * self.V[i, j] / (dot(self.W[i, :], self.H[:, j])[0, 0] + 1e-5)
            %    nr = nr + A(k,j)*(H(l,j)/sum(W(:,l)*H(l,j)));
            %end
            hnorm = zeros(1,n);
            wl = W(:,l);
            for j=1:n
                hnorm(j)=sum(wl*H(l,j));
            end
            hnorm(hnorm==0)=EPSILON;
            nr = A(k,:)*(H(l,:)./hnorm)';
            W(k,l)=(W(k,l)*nr)/sum(H(l,:));
        end
    end
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
    [iter toc]
end
end

%H is of size kxm
%n_i is a of size k
%avgs is of size kxm
function [n_i, avgs,classlabels]=encoding(H)
H_tmp=false(size(H));
H_tmp(find(H~=0))=true;
n_i = sum(H_tmp,2);
avgs = zeros(size(H));
for i=1:size(H,1)
    h_tmp = H(i,:);
    idxs = find(h_tmp);
    avgs(i,idxs)=sum(H(:,idxs))/numel(idxs);
    [~,classlabels]=max(H);
end
end
