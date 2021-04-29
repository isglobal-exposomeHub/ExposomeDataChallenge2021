clear
load('meta_merge.mat')
load('expos_merge.mat')
%% Draw original matrix
[cor P] = corrcoef([meta' expos']);
cor_bi = cor(1:221,222:390);
figure;imagesc(cor_bi);colormap jet;colorbar;snapnow
set(gca,'FontSize',12)
%print('cor_bi_raw','-dpng','-r300');

%% Draw matrix with detected structure
load('expos_meta_res.mat')
s_out = setdiff(1:221,s_in);
t_out = setdiff(1:169,t_in);

list_row = [s_in(Clist1) sort(s_out)];
list_col = [t_in(Clist2) sort(t_out)];
figure;imagesc(cor_bi(list_row,list_col));colormap jet;colorbar;
set(gca,'FontSize',12)
%daspect([1 1 1])
%print('cor_bi_sorted','-dpng','-r300');
figure;imagesc(cor_bi(s_in(Clist1),t_in(Clist2)));colormap jet;colorbar;
pbaspect([size(t_in,2) size(s_in,2) 1])
%print('W_in_sorted','-dpng','-r300');

%% Add variable names 
expos_names = replace(expos_names, '_', '.');
row_names = meta_names(s_in(Clist1));
col_names = expos_names(t_in(Clist2));
meta_group1 = meta_names(s_in(find(Cindx1==CID1(1))))
meta_group2 = meta_names(s_in(find(Cindx1==CID1(2))))
expos_group1 = expos_names(t_in(find(Cindx2==CID2(1))))
expos_group2 = expos_names(t_in(find(Cindx2==CID2(2))))
figure;imagesc(cor_bi(s_in(Clist1)',t_in(Clist2)));colormap jet;colorbar;
set(gca, 'XTick',1:size(t_in,2),'XTickLabel', col_names,'YTick',1:size(s_in,2),'YTickLabel', row_names,'FontSize',7)
xtickangle(45)
pbaspect([size(t_in,2) size(s_in,2) 1])
%print('expos2meta_names.png','-dpng','-r300');


%% Sort the correlation matrics by multi-to-multi associations
cor_meta = cor(1:221,1:221);
cor_expos = cor(222:390,222:390);

figure;imagesc(cor_meta(s_in(Clist1),s_in(Clist1)));colormap jet;colorbar;snapnow;
set(gca,'FontSize',12)
%print('meta_cor_rotate.png','-dpng','-r300');
figure;imagesc(cor_meta(list_row ,list_row));colormap jet;colorbar;snapnow;
set(gca,'FontSize',12)
%print('meta_cor_all_rotate.png','-dpng','-r300');

figure;imagesc(cor_expos(t_in(Clist2),t_in(Clist2)));colormap jet;colorbar;snapnow;
set(gca,'FontSize',12)
%print('expos_cor_rotate.png','-dpng','-r300');
figure;imagesc(cor_expos(list_col,list_col));colormap jet;colorbar;snapnow;
set(gca,'FontSize',12)
%print('expos_cor_all_rotate.png','-dpng','-r300');