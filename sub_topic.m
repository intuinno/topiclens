%% running standard nmf of subset
% function [mappedX_sub, cl_idx_sub, Wtopk_idx_sub] = sub_topic(idx, cl_idx)

addpath('./library/nmf');     
addpath('./library/ramkis');
addpath('./library/peripheral');
addpath('./library/discnmf');
addpath('./library/tSNE_matlab');

%% importing term-doc matrix
% in_txt_dir = 'text_files';
% out_fname = 'out.txt';
% 
% is_stemming = 1;
% 
% % merge multiple text files in a particular directory to a single file
% %[raw_txt,fnames] = merge_to_single_text(in_txt_dir,out_fname);
% 
% if ~exist('tmp1','dir')
%     mkdir('tmp1');
% end
% 
% delete(['tmp1' filesep '*.*']);
% 
% % run python function to generate a term-document matrix
% load('vis_paper_ver201508.mat');
% fid = fopen('out.txt', 'w');
% for i=1:size(PaperTitle,1)
%     fprintf(fid, '%s %s\n', PaperTitle{i}, Abstract{i});  
% end
% fclose(fid);
% dos(sprintf('python txt2mtx_fast.py %s',out_fname));
% 
% % importing dictionary and term-document matrix
% dict = import_dictionary(['tmp1' filesep 'vocabulary.txt']);
% tdm = import_tdm(['tmp1' filesep 'tmp.mtx']);
% A = sparse(tdm(:,1),tdm(:,2),tdm(:,3),max(tdm(:,1)),max(tdm(:,2)));
% clear tdm;
% 
% if is_stemming
%     tmp = cell(size(dict));
%         tmp{i}=cell(size(dict));
%     for i=1:length(dict)
%          tmp{i} = porterStemmer(dict{i});
%     end
%     [dict_stemmed,ia,ic] = unique(tmp);
%     A_stemmed = sparse([],[],[],length(dict_stemmed),size(A,2));
%     A_stemmed = A_stemmed';
%     A_transpose = A';
%     for i=1:length(dict_stemmed)
%         A_stemmed(:,i) = sum(A_transpose(:,ic==i),2);
%     end
%     A_stemmed = A_stemmed';
% 
%     A = A_stemmed;
%     dict = dict_stemmed;
%     clear A_transpose A_stemmed dict_stemmed;
% end
% 
% % [A w] = tfidf2(A);
% 
% %% additional stopwords
% addl_stopwords = {'visualization','visual','information','analysis','analysi','data','approach','process','based','technique','techniques','paper'};
% 
% idxs =[];
% 
% for i=1:length(addl_stopwords)
%     idxs = [idxs find(strcmp(dict,addl_stopwords{i}))];
% end
% 
% idxs = setdiff(1:size(A,1),idxs);
% A = A(idxs,:);
% dict = dict(idxs);

%% sub_A set
load tdm;
if iscell(idx)
    idx = cell2mat(idx);
end

A_sub = A(:,idx);
cluster_idx=unique(cl_idx(idx));
disp(cluster_idx);
sub_k=length(cluster_idx);
% no of topics
k_sub = max(min([floor(length(idx)/20) 5]), sub_k);
if (k_sub<2) 
    k_sub = 2;
end
% no of top keywords
topk_sub = 5;

% normalization
A_norm_sub = bsxfun(@rdivide,A_sub,sqrt(sum(A_sub.^2)));  

% choosing one among different preprocessings
target_A_sub = A_norm_sub;     % replaced by below code (9/10) <- original
% target_A = A_norm;

%% hier8_neat1 (initial run)
cnt = 0;
if(cnt == 0)
    if ~exist('params', 'var')
        trial_allowance = 3;
        unbalanced = 0.1;
        vec_norm = 2.0;
        normW = true;
        anls_alg = @anls_entry_rank2_precompute;
        tol = 1e-4;
        maxiter = 10000;
    else
        if isfield(params, 'trial_allowance')
            trial_allowance = params.trial_allowance;
        else
            trial_allowance = 3;
        end
        if isfield(params, 'unbalanced')
            unbalanced = params.unbalanced;
        else
            unbalanced = 0.1;
        end
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
        if isfield(params, 'anls_alg')
            anls_alg = params.anls_alg;
        else
            anls_alg = @anls_entry_rank2_precompute;
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

    params = [];
    params.vec_norm = vec_norm;
    params.normW = normW;
    params.anls_alg = anls_alg;
    params.tol = tol;
    params.maxiter = maxiter;

    [m, n] = size(target_A_sub);

    cluster_idx=unique(cl_idx(idx));                    %
    sub_k=length(cluster_idx);                          % number of cluster selected
    ctrary = zeros(size(target_A_sub,1),sub_k);         %
    timings_sub = zeros(1, k_sub-sub_k);                
    clusters_sub = cell(1, 2*(k_sub-sub_k));
    Ws_sub = cell(1, 2*(k_sub-sub_k));
    W_buffer = cell(1, 2*(k_sub-sub_k));
    H_buffer = cell(1, 2*(k_sub-sub_k));
    priorities_sub = zeros(1, 2*(k_sub-sub_k));
    is_leaf_sub = -1 * ones(1, 2*(k_sub-sub_k));
    tree_sub = zeros(2, 2*(k_sub-sub_k));
    splits_sub = -1 * ones(1, k_sub-sub_k);
    min_priority = 1e308;


    for i=1:sub_k
        subidx = find(cl_idx(idx)==cluster_idx(i));
        ctrary(:,i) = mean(target_A_sub(:, subidx),2);
        Ws_sub{i} = ctrary(:,i);
        is_leaf_sub(i) = 1;
        [~, W_buffer_one, H_buffer_one, priority_one]...
            = trial_split(trial_allowance, unbalanced, min_priority, target_A_sub, subidx, Ws_sub{i}, params);
        clusters_sub{i} = subidx;
        W_buffer{i} = W_buffer_one;
        H_buffer{i} = H_buffer_one;
        priorities_sub(i) = priority_one;
    end
    cnt = 1;
    result_used = sub_k;
end
%% hier8_neat2 (run interatively)
if k_sub == sub_k
    W_sub = cell2mat(Ws_sub(find(is_leaf_sub)));
    [H_sub,temp,suc_H,numChol_H,numEq_H] = nnlsm_blockpivot(W_sub'*W_sub,W_sub'*target_A_sub,1,rand(k_sub,size(target_A_sub,2)));

else k_sub > sub_k
    for i = 1:k_sub-sub_k
        leaves = find(is_leaf_sub == 1);
        temp_priority = priorities_sub(leaves);
        min_priority = min(temp_priority(temp_priority > 0));
        [max_priority, split_node] = max(temp_priority);
        
        if max_priority < 0
            fprintf('Cannot generate all %d leaf clusters_sub\n', k);
            return;
        end
        split_node = leaves(split_node);
        is_leaf_sub(split_node) = 0;
        W_sub = W_buffer{split_node};
        H_sub = H_buffer{split_node};
        split_subset = clusters_sub{split_node};
        new_nodes = [result_used+1 result_used+2];
        tree_sub(1, split_node) = new_nodes(1);
        tree_sub(2, split_node) = new_nodes(2);
        
        result_used = result_used + 2;
        [max_val, cluster_subset] = max(H_sub);
        clusters_sub{new_nodes(1)} = split_subset(cluster_subset == 1);
        clusters_sub{new_nodes(2)} = split_subset(cluster_subset == 2);
        Ws_sub{new_nodes(1)} = W_sub(:, 1);
        Ws_sub{new_nodes(2)} = W_sub(:, 2);
        splits_sub(i) = split_node;
        is_leaf_sub(new_nodes) = 1;
        
        subset = clusters_sub{new_nodes(1)};
        [subset, W_buffer_one, H_buffer_one, priority_one] = trial_split(trial_allowance, unbalanced, min_priority, target_A_sub, subset, W_sub(:, 1), params);
        clusters_sub{new_nodes(1)} = subset;
        W_buffer{new_nodes(1)} = W_buffer_one;
        H_buffer{new_nodes(1)} = H_buffer_one;
        priorities_sub(new_nodes(1)) = priority_one;
        
        subset = clusters_sub{new_nodes(2)};
        [subset, W_buffer_one, H_buffer_one, priority_one] = trial_split(trial_allowance, unbalanced, min_priority, target_A_sub, subset, W_sub(:, 2), params);
        clusters_sub{new_nodes(2)} = subset;
        W_buffer{new_nodes(2)} = W_buffer_one;
        H_buffer{new_nodes(2)} = H_buffer_one;
        priorities_sub(new_nodes(2)) = priority_one;
        
        W_sub = cell2mat(Ws_sub(find(is_leaf_sub)));
        [H_sub,temp,suc_H,numChol_H,numEq_H] = nnlsm_blockpivot(W_sub'*W_sub,W_sub'*target_A_sub,1,rand(i+sub_k,size(target_A_sub,2)));
    end
end

%%
% displaying top keywords for each topic
[Wtopk_sub,Htopk_sub,DocTopk_sub,Wtopk_idx_sub] = parsenmf(W_sub,H_sub,dict,topk_sub);

Wtopk_idx_sub = Wtopk_idx_sub';
[~,cl_idx_sub] = max(H_sub);

mappedX_sub = 'undefined';

%% t-sne visualization
% Run t-SNE
% no_dims=2;
% initial_dims_sub=min([30, size(A_sub,2)]);
% shrink_factor = .3;
% perplexity = 30;
% mappedX_sub = tsne_sup(target_A_sub', cl_idx_sub, shrink_factor, no_dims, initial_dims_sub, perplexity); 
% Run t?SNE
% mappedX = tsne(target_A', cl_idx, no_dims, initial_dims, perplexity);
%mappedX_sub = tsne(target_A_sub', cl_idx_sub, no_dims, initial_dims, perplexity);
