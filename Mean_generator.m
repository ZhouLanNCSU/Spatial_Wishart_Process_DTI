function t_central=Mean_generator(s,B_all,t_all,a,graph,ndim) 



%t_central=t_all(:,s)-[arrayfun(@(d) B_all{d}*t_all((graph(d,:)~=0),s) ,[1:ndim-1]'); 0];
MYB=(blkdiag(B_all{1:ndim-1}));
% 
% a = cell(ndim-1,1);
% for d = 1:ndim-1
%     a{d}= sparse(t_all((graph(d,:)~=0),s));
% end


%t_all_par=sparse(t_all(:,s));
%a=arrayfun(@(d) (t_all_par(graph(d,:)~=0)),[1:ndim-1], 'UniformOutput',false);
% a = cell(ndim-1,1);
% for d = 1:ndim-1
%     a{d}= t_all_par(graph(d,:)~=0);
% end
MYT=(blkdiag(a{:}));




mean=diag(MYB*MYT);


%mean=arrayfun(@(d) B_all{d}*t_all((graph(d,:)~=0),s) ,[1:ndim-1]');
mean=[mean; 0];

t_central=t_all(:,s)-mean;

end