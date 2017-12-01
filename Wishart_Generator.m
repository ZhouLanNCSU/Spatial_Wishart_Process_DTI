function Data=Wishart_Generator(X,n_S,n_R,dof,Scale_Status)

if Scale_Status== true
    Data=arrayfun(@(index) X(:,(1+n_S*(index-1)):n_S*index)'*X(:,(1+n_S*(index-1)):n_S*index),[1:n_R],'UniformOutput', false);
else
    Data=arrayfun(@(index) X(:,(1+n_S*(index-1)):n_S*index)'*X(:,(1+n_S*(index-1)):n_S*index)/dof,[1:n_R],'UniformOutput', false);
end



end

