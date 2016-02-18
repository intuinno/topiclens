function compare(AAA,BBB,par,data_name_str,i)
    disp(data_name_str);
    
    k=10; kd=5;
    par.k = k;
% for kd = 1:k    
    [m,n1] = size(AAA);
	[~,n2] = size(BBB);
    init.W=rand(m,k); init.U=rand(m,k); init.H=rand(k,n1); init.V=rand(k,n2);
    num = par.num; num2 = par.num2;
    if size(num,2)>size(num2,2)
        n = size(num,2)-size(num2,2);
        num2 = [repmat(' ',size(num,1),n) num2];
    end
    if size(num,2)<size(num2,2)
        n = size(num2,2)-size(num,2);
        num = [repmat(' ',size(num,1),n) num];
    end
    par.num = num; par.num2 = num2;
    
    % measures
    recon = zeros(3,3);

    % 1. baseline : separate NMF
    [W,H,~,iter,~]=nmf(AAA,k,'method','hals','verbose',0,'init',init);
    init2.W = init.U; init2.H = init.V;
    [U,V,~,iter,~]=nmf(BBB,k,'method','hals','verbose',0,'init',init2);
    recon(1,1:2) = [norm(AAA-W*H,'fro')/n1, norm(BBB-U*V,'fro')/n2];
    sep.W=W;sep.H=H;sep.U=U;sep.V=V;
    
    
    % 2. baseline : batch NMF
    [W,H,U,V,~,iter,~]=nmf2(AAA,BBB,k,kd,1000,100,'method','comp','verbose',0,'init',init);
    recon(2,1:2) = [norm(AAA-W*H,'fro')/n1, norm(BBB-U*V,'fro')/n2];
    batch.W=W;batch.H=H;batch.U=U;batch.V=V;
    
    
    % 3. disc. NMF
    [W,H,U,V,~,iter,~]=nmf4(AAA,BBB,k,kd,1,1000,100,100,'method','disc','verbose',0,'init',init);
    recon(3,1:2) = [norm(AAA-W*H,'fro')/n1, norm(BBB-U*V,'fro')/n2];
    disc.W=W;disc.H=H;disc.U=U;disc.V=V;
    
    % reorder
    [assignment,cost] = munkres(-sep.U'*sep.W);
    sep.W = sep.W(:,assignment);
    sep.H = sep.H(assignment,:);
    [val ix] = sort(diag(sep.U'*sep.W));
    sep.W = sep.W(:,ix);
    sep.H = sep.H(ix,:);
    sep.U = sep.U(:,ix);
    sep.V = sep.V(ix,:);
    [assignment,cost] = munkres(-disc.W(:,1:kd)'*sep.W(:,1:kd));
    sep.W(:,1:kd) = sep.W(:,assignment);
    sep.H(1:kd,:) = sep.H(assignment,:);
    [assignment,cost] = munkres(-disc.U(:,1:kd)'*sep.U(:,1:kd));
    sep.U(:,1:kd) = sep.U(:,assignment);
    sep.V(1:kd,:) = sep.V(assignment,:);
    discW = [disc.W(:,kd+1:k);disc.U(:,kd+1:k)]; sepW = [sep.W(:,kd+1:k);sep.U(:,kd+1:k)];
    [assignment,cost] = munkres(-discW'*sepW);
    sep.W(:,kd+1:k) = sep.W(:,assignment+kd);
    sep.H(kd+1:k,:) = sep.H(assignment+kd,:);
    sep.U(:,kd+1:k) = sep.U(:,assignment+kd);
    sep.V(kd+1:k,:) = sep.V(assignment+kd,:);    
    [assignment,cost] = munkres(-disc.W(:,1:kd)'*batch.W(:,1:kd));
    batch.W(:,1:kd) = batch.W(:,assignment);
    batch.H(1:kd,:) = batch.H(assignment,:);
    [assignment,cost] = munkres(-disc.U(:,1:kd)'*batch.U(:,1:kd));
    batch.U(:,1:kd) = batch.U(:,assignment);
    batch.V(1:kd,:) = batch.V(assignment,:);
    discW = [disc.W(:,kd+1:k);disc.U(:,kd+1:k)]; batchW = [batch.W(:,kd+1:k);batch.U(:,kd+1:k)];
    [assignment,cost] = munkres(-discW'*batchW);
    batch.W(:,kd+1:k) = batch.W(:,assignment+kd);
    batch.H(kd+1:k,:) = batch.H(assignment+kd,:);
    batch.U(:,kd+1:k) = batch.U(:,assignment+kd);
    batch.V(kd+1:k,:) = batch.V(assignment+kd,:);
    
    % print
    disp('separate NMF');
    sep.topics = print(sep.W,sep.U,par);
    disp('batch NMF');
    batch.topics = print(batch.W,batch.U,par);
    disp('disc NMF');
    disc.topics = print(disc.W,disc.U,par);

    recon(:,3) = sum(recon(:,1:2),2);
    disp(recon);
    
    fname = sprintf('.\\%s-k%i_%i', data_name_str, kd, i);
    save(fname, 'AAA','BBB','par','init','recon','sep','batch','disc');
% end    
%     if strcmp(data_name_str,'infovast_norm') == 1
%         fig = 5;
%     elseif strcmp(data_name_str,'infovast_tfidf') == 1
%         fig = 6;
%     elseif strcmp(data_name_str,'infovast_tfidf_norm') == 1
%         fig = 7;
%     else
%         fig = 4;
%     end
%     figure(fig);subplot(4,3,1);imagesc(sep.W);
%     figure(fig);subplot(4,3,2);imagesc(sep.U);
%     figure(fig);subplot(4,3,3);imagesc(sep.H);
%     figure(fig);subplot(4,3,4);imagesc(sep.V);
%     figure(fig);subplot(4,3,5);imagesc(batch.W);
%     figure(fig);subplot(4,3,6);imagesc(batch.U);
%     figure(fig);subplot(4,3,7);imagesc(batch.H);
%     figure(fig);subplot(4,3,8);imagesc(batch.V);
%     figure(fig);subplot(4,3,9);imagesc(disc.W);
%     figure(fig);subplot(4,3,10);imagesc(disc.U);
%     figure(fig);subplot(4,3,11);imagesc(disc.H);
%     figure(fig);subplot(4,3,12);imagesc(disc.V); 
end
function topics = print(W,U,par)
    num = par.num; num2 = par.num2; dic = par.dic;
    [val,IX] = sort([W U],'descend'); topw = IX(1:10,:);
    len = (size(dic,2)+1+size(num,2));
    topics = [];
    k=par.k;
    for i = 1:k
        topics = [topics [num(topw(:,i),:),repmat(' ',10,1),dic(topw(:,i),:); repmat('-',1,len); num2(topw(:,i+k),:),repmat(' ',10,1),dic(topw(:,i+k),:)]];
    end
    topics
%     topics = [num(topw(:,1),:),repmat(' ',10,1),dic(topw(:,1),:),num(topw(:,2),:),repmat(' ',10,1),dic(topw(:,2),:),num(topw(:,3),:),repmat(' ',10,1),dic(topw(:,3),:),num(topw(:,4),:),repmat(' ',10,1),dic(topw(:,4),:),num(topw(:,5),:),repmat(' ',10,1),dic(topw(:,5),:),num(topw(:,6),:),repmat(' ',10,1),dic(topw(:,6),:),num(topw(:,7),:),repmat(' ',10,1),dic(topw(:,7),:),num(topw(:,8),:),repmat(' ',10,1),dic(topw(:,8),:),num(topw(:,9),:),repmat(' ',10,1),dic(topw(:,9),:),num(topw(:,10),:),repmat(' ',10,1),dic(topw(:,10),:);repmat('-',1,len);num2(topw(:,11),:),repmat(' ',10,1),dic(topw(:,11),:),num2(topw(:,12),:),repmat(' ',10,1),dic(topw(:,12),:),num2(topw(:,13),:),repmat(' ',10,1),dic(topw(:,13),:),num2(topw(:,14),:),repmat(' ',10,1),dic(topw(:,14),:),num2(topw(:,15),:),repmat(' ',10,1),dic(topw(:,15),:),num2(topw(:,16),:),repmat(' ',10,1),dic(topw(:,16),:),num2(topw(:,17),:),repmat(' ',10,1),dic(topw(:,17),:),num2(topw(:,18),:),repmat(' ',10,1),dic(topw(:,18),:),num2(topw(:,19),:),repmat(' ',10,1),dic(topw(:,19),:),num2(topw(:,20),:),repmat(' ',10,1),dic(topw(:,20),:)]
    disp(W'*U)
end