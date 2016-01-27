%% running standard nmf of subset
function [mappedX_sub, cl_idx_sub, Wtopk_idx_sub] = sub_topic_tsne(idx)

addpath('./library/nmf');     
addpath('./library/ramkis');
addpath('./library/peripheral');
addpath('./library/discnmf');
addpath('./library/tSNE_matlab');

%% importing term-doc matrix
%in_txt_dir = 'text_files';
out_fname = 'out.txt';

is_stemming = 1;

% merge multiple text files in a particular directory to a single file
[raw_txt,fnames] = merge_to_single_text(in_txt_dir,out_fname);

if ~exist('tmp1','dir')
    mkdir('tmp1');
end

delete(['tmp1' filesep '*.*']);

% run python function to generate a term-document matrix
load('vis_paper_ver201508.mat');
fid = fopen('out.txt', 'w');
for i=1:size(PaperTitle,1)
    fprintf(fid, '%s %s\n', PaperTitle{i}, Abstract{i});  
end
fclose(fid);
dos(sprintf('python txt2mtx_fast.py %s',out_fname));

% importing dictionary and term-document matrix
dict = import_dictionary(['tmp1' filesep 'vocabulary.txt']);
tdm = import_tdm(['tmp1' filesep 'tmp.mtx']);
A = sparse(tdm(:,1),tdm(:,2),tdm(:,3),max(tdm(:,1)),max(tdm(:,2)));
clear tdm;

if is_stemming
    tmp = cell(size(dict));
    for i=1:length(dict)
        tmp{i} = porterStemmer(dict{i});
    end
    [dict_stemmed,ia,ic] = unique(tmp);
    A_stemmed = sparse([],[],[],length(dict_stemmed),size(A,2));
    A_stemmed = A_stemmed';
    A_transpose = A';
    for i=1:length(dict_stemmed)
        A_stemmed(:,i) = sum(A_transpose(:,ic==i),2);
    end
    A_stemmed = A_stemmed';

    A = A_stemmed;
    dict = dict_stemmed;
    clear A_transpose A_stemmed dict_stemmed;
end

% [A w] = tfidf2(A);

%% additional stopwords
addl_stopwords = {'visualization','visual','information','analysis','analysi','data','approach','process','based','technique','techniques','paper'};

idxs =[];

for i=1:length(addl_stopwords)
    idxs = [idxs find(strcmp(dict,addl_stopwords{i}))];
end

idxs = setdiff(1:size(A,1),idxs);

A = A(idxs,:);
dict = dict(idxs);


idx;
idx = cell2mat(idx);
A_sub = A(:,idx);
% no of topics
k_sub = min([floor(length(idx)/10) 10]);
% no of top keywords
topk_sub = 5;

% normalization
A_norm_sub = bsxfun(@rdivide,A_sub,sqrt(sum(A_sub.^2)));  

% choosing one among different preprocessings
target_A_sub = A_norm_sub;     % replaced by below code (9/10) <- original
% target_A = A_norm;

%%
% target_A = A_idf;
size(target_A_sub)
[tree_sub, splits_sub, is_leaf_sub, clusters_sub, timings_sub, Ws_sub, priorities_sub, W_sub, H_sub] = hier8_neat(target_A_sub, k_sub);

%%
% displaying top keywords for each topic
[Wtopk_sub,Htopk_sub,DocTopk_sub,Wtopk_idx_sub] = parsenmf(W_sub,H_sub,dict,topk_sub);

Wtopk_idx_sub = Wtopk_idx_sub';
[~,cl_idx_sub] = max(H_sub);
Wlen=size(Wtopk_sub,2);

%% t-sne visualization
% Run t-SNE
no_dims=2;
initial_dims_sub=min([30, size(A_sub,2)]);
shrink_factor = .3;
perplexity = 30;
mappedX_sub = tsne_sup(target_A_sub', cl_idx_sub, shrink_factor, no_dims, initial_dims_sub, perplexity); 
% Run t?SNE
% mappedX = tsne(target_A', cl_idx, no_dims, initial_dims, perplexity);
%mappedX_sub = tsne(target_A_sub', cl_idx_sub, no_dims, initial_dims, perplexity);

end
