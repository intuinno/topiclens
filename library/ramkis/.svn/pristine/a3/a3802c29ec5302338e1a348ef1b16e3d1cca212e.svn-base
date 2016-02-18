function H = anls_entry_blockpivot_precompute(WtW, AtW, H)

% WtW: k*k
% AtW: n*k
% H: n*k
% Returning H of size n*k also

H = nnlsm_blockpivot(WtW, AtW', 1, H')';
