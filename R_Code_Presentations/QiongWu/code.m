clear
%% Load correlation matrix between exposome and metabolic variables
load('meta_merge.mat')
load('expos_merge.mat')
[cor P] = corrcoef([meta' expos']);
cor_bi = cor(1:221,222:390);
figure;imagesc(cor_bi);colormap jet;colorbar;snapnow

%% Detect a dense bicluster in this n*m matrix
W=abs(cor_bi);
lambda_vec = 0.5:0.02:1.4;
c_vec = [1./(1:1:20) 1:1:20];
W_vec = W(:);
W_vec = W_vec(W_vec>0);
r = median(W_vec);

[s_rev,t_rev, c_rev, lambda_rev,max_lik]=greedy_lik_fun(W,c_vec,lambda_vec,r,30);
[N M]= size(W);
s_in = setdiff(1:N,s_rev);
t_in = setdiff(1:M,t_rev);
figure;imagesc(W(s_in,t_in));colormap jet;colorbar;
figure;imagesc(W([s_in sort(s_rev)],[t_in sort(t_rev)]));colormap jet;colorbar;


%% Refine the pattern through similarity matrices
% similarity among rows
W_in = cor_bi(s_in,t_in);
S1 = [];
norm_r = sqrt(sum(abs(W_in).^2,2));
for i = 1:size(W_in,1)
    for j = i:size(W_in,1)
        S1(i,j) = dot(W_in(i,:), W_in(j,:)) / (norm_r(i) * norm_r(j));
        S1(j,i) = S1(i,j);
    end
        S1(i,i)=0;
end

figure;imagesc(S1);colormap jet;colorbar;

[Cindx1,CID1,Clist1,K_final1]=NICE(squareform(S1), 0.5,0,5);
S1_sort = S1(Clist1,Clist1);
figure;imagesc(S1_sort);colormap jet;colorbar;snapnow

% similarity among columns
S2 = [];
norm_r = sqrt(sum(abs(W_in).^2,1));
for i = 1:size(W_in,2)
    for j = i:size(W_in,2)
        S2(i,j) = dot(W_in(:,i), W_in(:,j)) / (norm_r(i) * norm_r(j));
        S2(j,i) = S2(i,j);
    end
        S2(i,i)=0;
end

figure;imagesc(S2);colormap jet;colorbar;

[Cindx2,CID2,Clist2,K_final2]=NICE(squareform(S2), 0.3,0,5);
S2_sort = S2(Clist2,Clist2);
figure;imagesc(S2_sort);colormap jet;colorbar;snapnow

figure;imagesc(W_in(Clist1,Clist2));colormap jet;colorbar;snapnow

%% Save results
save('expos_meta_res.mat','s_in','t_in','Cindx1','CID1','Clist1','K_final1','Cindx2','CID2','Clist2','K_final2')
