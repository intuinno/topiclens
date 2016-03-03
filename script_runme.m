clear all;
close all;

addpath('./library/nmf');     
addpath('./library/ramkis');
addpath('./library/peripheral');
addpath('./library/discnmf');


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
[W,H]=nmf(target_A, k); % nmf() is matrix decomposition on A to get W,H (i.e. A=W*H); num of topic = k ; =
toc

%%
% displaying top keywords for each topic
[Wtopk,Htopk,DocTopk,Wtopk_idx] = parsenmf(W,H,dict,topk);
Wtopk
