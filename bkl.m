function t_for_blk=bkl(t_all,graph,ndim,nsub)

t_for_blk=cell(nsub,1);

for s=1:nsub
    t_all_par=sparse(t_all(:,s));
    

%     a=arrayfun(@(d) t_all_par(graph(d,:)~=0), [1:ndim-1]','UniformOutput',false);
    
    
    
    a = cell(ndim-1,1);
    for d = 1:ndim-1
        a{d}= t_all_par(graph(d,:)~=0);
    end
    t_for_blk{s}=a;
end




end
