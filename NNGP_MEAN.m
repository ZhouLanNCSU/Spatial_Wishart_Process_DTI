function t_central_all=NNGP_MEAN(B_all,t_all,t_all_for_blk,graph,nsub,ndim)


t_central_all=cell2mat(arrayfun(@(s) Mean_generator(s,B_all,t_all,t_all_for_blk{s},graph,ndim) ,[1:nsub],'UniformOutput',0));





end