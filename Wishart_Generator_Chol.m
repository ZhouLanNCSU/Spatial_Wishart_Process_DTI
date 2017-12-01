function [z11 z22 z33 z21 z31 z32]=Wishart_Generator_Chol(X,n_S,n_R,dof,Scale_Status)

if Scale_Status== true
    Data=arrayfun(@(index) X(:,(1+n_S*(index-1)):n_S*index)'*X(:,(1+n_S*(index-1)):n_S*index)/dof,[1:n_R],'UniformOutput', false);
else
    Data=arrayfun(@(index) X(:,(1+n_S*(index-1)):n_S*index)'*X(:,(1+n_S*(index-1)):n_S*index),[1:n_R],'UniformOutput', false);
end


result = cellfun(@chol, Data, 'UniformOutput', false);
result=cell2mat(result);

z11=(result(1,1:n_S:n_S*n_R).^2)';
z22=(result(2,2:n_S:n_S*n_R).^2)';
z33=(result(3,3:n_S:n_S*n_R).^2)';
z21=(result(1,2:n_S:n_S*n_R))';
z31=(result(1,3:n_S:n_S*n_R))';
z32=(result(2,3:n_S:n_S*n_R))';



end

