function [s_out t_out c_max]=greedy_bipar(W,c_vec,lambda,k)
W1 = W;
clear remove_s_vec remove_t_vec dST_vec deleteS_list deleteT_list
for l=1:size(c_vec,2)
    %l
    W=W1;
    c=c_vec(l);
    deleteS_list={};
    deleteT_list={};
    density_list = zeros((size(W,1)+size(W,2))-1,8);
    for i=1:(size(W,1)+size(W,2))
        if W==0
            break
        end
        
        %if i==1
        %    col=size(W,2);
        %    row=size(W,1);
        %elseif density_list(i-1,3)==1
        %        row=row-size(deleteS_list{i-1},2);
        %    else
        %        col=col-size(deleteT_list{i-1},2);
        %end
       
        row = size(W,1)-sum(sum(W,2)==0);
        col = size(W,2)-sum(sum(W,1)==0);
        
        k1 = ceil(row/k);
        k2 = ceil(col/k);
  
        % mean of rows
        R=sum(W,2);
        % mean of columns
        C=sum(W,1);
        % min k1 sum in row 
        [~,index]=sort(R);
        IndS=sort(index(sum(R==0)+1:sum(R==0)+k1));
        resultR=R(IndS);
        dS=mean(resultR);
        %z=ismember(R,resultR);
        % min k2 sum in column 
        [~,index]=sort(C);
        IndT=sort(index(sum(C==0)+1:sum(C==0)+k2));
        resultC=C(IndT);
        dT=mean(resultC);
        %z=ismember(C,resultC);
        
        if sqrt(c)*dS <= dT/sqrt(c)
            W(IndS,:)=zeros(k1,size(W,2));
            density_list(i,1)=sum(R~=0);
            density_list(i,2)=sum(C~=0);
            density_list(i,3)=1;
        else
            W(:,IndT)=zeros(size(W,1),k2);
            density_list(i,1)=sum(R~=0);
            density_list(i,2)=sum(C~=0);
            density_list(i,3)=2;
        end
        deleteS_list{i}=IndS;
        deleteT_list{i}=IndT;
        %density_list(i,5)=dS;
        %density_list(i,7)=dT;
        density_list(i,8)=(sum(R))/(sqrt(density_list(i,1)*density_list(i,2)))^lambda;
    end
    
    %ouput W
    [dST, indST]=max(density_list(:,8));
    remove_s=[];
    remove_t=[];
    for j=1:indST
        if density_list(j,3)==1
            remove_s = union(remove_s, deleteS_list{j});
            %remove_s =[remove_s,deleteS_list{j}];
        else
            remove_t = union(remove_t, deleteT_list{j});
            %remove_t =[remove_t,deleteT_list{j}];
        end
    end
remove_s_vec{l}=remove_s;
remove_t_vec{l}=remove_t;
dST_vec(l)=dST;
end
c_max = find(dST_vec == max(dST_vec));
s_out0 = remove_s_vec{c_max(1)}';
t_out0 = remove_t_vec{c_max(1)}';

W = W1;
W(s_out0,:) = zeros(size(s_out0,2),size(W,2));
W(:,t_out0) = zeros(size(W,1),size(t_out0,2));
t_out = find(sum(W,1)==0);
s_out = find(sum(W,2)==0)';
end
