function [mappedX, cl_idx, Wtopk_idx, dict] = main_topic_tsne()

addpath('./library/nmf');     
addpath('./library/ramkis');
addpath('./library/peripheral');
addpath('./library/discnmf');
addpath('./library/tSNE_matlab');

%% importing term-doc matrix
in_txt_dir = 'text_files';
out_fname = 'out.txt';
is_stemming = 1;

% merge multiple text files in a particular directory to a single file
%[raw_txt,fnames] = merge_to_single_text(in_txt_dir,out_fname);

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


%% running standard nmf
% no of topics
k = 10 ;
% no of top keywords
topk =  5;

% normalization
A_norm = bsxfun(@rdivide,A,sqrt(sum(A.^2)));  

% choosing one among different preprocessings
target_A = A;     % replaced by below code (9/10) <- original
% target_A = A_norm;

%%
% target_A = A_idf;
[tree, splits, is_leaf, clusters, timings, Ws, priorities, W, H] = hier8_neat(target_A, k);

%%
% displaying top keywords for each topic
[Wtopk,Htopk,DocTopk,Wtopk_idx] = parsenmf(W,H,dict,topk);

Wtopk_idx = Wtopk_idx';
[~,cl_idx] = max(H);


%% t-sne visualization
no_dims = 2;
initial_dims = 50;
perplexity = 30;
% Run t?SNE
% mappedX = tsne(target_A', cl_idx, no_dims, initial_dims, perplexity);

mappedX = tsne_sup(target_A', cl_idx, .7, no_dims, initial_dims, perplexity);
% Run t?SNE
%mappedX = tsne(target_A', cl_idx, no_dims, initial_dims, perplexity);

end

