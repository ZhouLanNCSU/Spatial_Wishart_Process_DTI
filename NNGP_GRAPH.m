function Graph = NNGP_GRAPH(m,ori_grid,ndim )

Graph=sparse(repmat(0,[ndim ndim]));
for i=1:ndim-1
    index=knnsearch(ori_grid(i:end,:),ori_grid(i,:),'k',m+1);
    Graph(i,setdiff(index,1)+(i-1))=1;
    
end


end

