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


%% running standard nmf
% no of topics
k = 10 ;
% no of top keywords
topk = 15;

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
mappedX = tsne(target_A', cl_idx, no_dims, initial_dims, perplexity);

%%
% Plot results
clf;
fig = gscatter(mappedX(:,1), mappedX(:,2), cl_idx);
for i=1:k
    tmp = mean(mappedX(cl_idx==i,:));
    tmp_str =[];
    for j=1:3
        tmp_str = [tmp_str sprintf(' %s',Wtopk{j,i})];
        if mod(j,5)==0
            tmp_str = [tmp_str sprintf('\n')];
        end
    end
%     tmp_str
    text(tmp(1),tmp(2),sprintf('%02d. %s',i,tmp_str),'HorizontalAlignment','center');
    
end

%%
% find the index of neighborhood of picked point
[x y] = ginput(1);              
X= mappedX;                     %all the coordinates
Y = [x y];                      %point picked
r = 5;                          %distance
idx = rangesearch(X, Y, r);     %neighborhood
idx=cat(1,idx{:});

% make term-doc of neighborhood   
A_sub=[];
for i=1:size(idx,2)  
    A_sub = [A_sub, A(:,idx)];
end



    







