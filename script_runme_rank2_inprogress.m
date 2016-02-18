clear all;
close all;

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
[raw_txt,fnames] = merge_to_single_text(in_txt_dir,out_fname);

if ~exist('tmp1','dir')
    mkdir('tmp1');
end

delete('tmp1\\*.*');

% run python function to generate a term-document matrix
dos(sprintf('python txt2mtx_fast.py %s',out_fname));

% importing dictionary and term-document matrix
dict = import_dictionary('tmp1\\vocabulary.txt');
tdm = import_tdm('tmp1\\tmp.mtx');
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
tic
[tree, splits, is_leaf, clusters, timings, Ws, priorities, W, H] = hier8_neat(target_A, k);
toc

tic
[W_nmf,H_nmf]=nmf(target_A, k); 
% nmf() is matrix decomposition on A to get W,H (i.e. A=W*H); num of topic = k ; =
toc


%%
% displaying top keywords for each topic
[Wtopk,Htopk,DocTopk,Wtopk_idx] = parsenmf(W,H,dict,topk);
Wtopk

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

%%
% Plot results
clf;
gscatter(mappedX(:,1), mappedX(:,2), cl_idx);
title ('tnse')
fig1=gcf;
for i=1:k
    tmp = mean(mappedX(cl_idx==i,:));
    tmp_str =[];
    for j=1:topk
        tmp_str = [tmp_str sprintf(' %s',Wtopk{j,i})];
        if mod(j,5)==0
            tmp_str = [tmp_str sprintf('\n')];
        end
    end
%     tmp_str

    text(tmp(1),tmp(2),sprintf('%02d. %s',i,tmp_str),'HorizontalAlignment','center');
    
end

%% making the subset
%find the neighborhood of picked point
[x y] = ginput(1);              
X= mappedX;                     %all the coordinates
Y = [x y];                      %point picked
r = 5;                          %distance
idx = rangesearch(X, Y, r);     %neighborhood

%%
% make term-doc of neighborhood   
idx = idx{1};
A_sub = A(:,idx);
%% running standard nmf of subset
% no of topics
k_sub = min([floor(length(idx)/10) 10]) ;
% no of top keywords
topk_sub = 5;

% normalization
A_norm_sub = bsxfun(@rdivide,A_sub,sqrt(sum(A_sub.^2)));  

% choosing one among different preprocessings
target_A_sub = A_norm_sub;     % replaced by below code (9/10) <- original
% target_A = A_norm;

%%
% target_A = A_idf;
tic
[tree_sub, splits_sub, is_leaf_sub, clusters_sub, timings_sub, Ws_sub, priorities_sub, W_sub, H_sub] = hier8_neat(target_A_sub, k_sub);
toc
%%
% displaying top keywords for each topic
[Wtopk_sub,Htopk_sub,DocTopk_sub,Wtopk_idx_Sub] = parsenmf(W_sub,H_sub,dict,topk);
Wtopk_sub

[~,cl_idx_sub] = max(H_sub);
Wlen=size(Wtopk_sub,2);


%% t-sne visualization
% Run t-SNE
figure;
no_dims=2;
initial_dims_sub=min([30, size(A_sub,2)]);
shrink_factor = .3;

mappedX_sub = tsne_sup(target_A_sub', cl_idx_sub, shrink_factor, no_dims, initial_dims_sub, perplexity); 
% Run t?SNE
% mappedX = tsne(target_A', cl_idx, no_dims, initial_dims, perplexity);
%mappedX_sub = tsne(target_A_sub', cl_idx_sub, no_dims, initial_dims, perplexity);
%%
gscatter(mappedX_sub(:,1), mappedX_sub(:,2), cl_idx_sub);
% title ('tnse subset')
fig2=gcf;
for i=1:Wlen
    tmp_sub = mean(mappedX_sub(cl_idx_sub==i,:));
    tmp_str_sub =[];
    for j=1:topk_sub
        tmp_str_sub = [tmp_str_sub sprintf(' %s',Wtopk_sub{j,i})];
        if mod(j,5)==0
            tmp_str_sub = [tmp_str_sub sprintf('\n')];
        end
    end
%     tmp_str
    text(tmp_sub(1),tmp_sub(2),sprintf('%02d. %s',i,tmp_str_sub),'HorizontalAlignment','center');   
end
%%
[h_m, h_i]=inset(fig1,fig2, 0.3, x, y);
%set(h_i,'position', [x y 0.5 0.5])

