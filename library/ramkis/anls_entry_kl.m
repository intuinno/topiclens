function H = anls_entry_kl(W, A, H)

% W: m*k
% A: m*n
% H: n*k
% Returning H of size n*k also

[n, k] = size(H);
sumW = sum(W);
H = H .* ((A./(W*H'))'*W) ./ (sumW(ones(1, n), :) + 1e-9);
