function [B_all F_all]=NNGP_COV(graph,R,ndim)


B_all=[arrayfun(@(s) B_generator(s,graph,R),[1:ndim-1],'UniformOutput',0) 1];
F_all=[arrayfun(@(s) F_generator(s,graph,R,B_all),[1:ndim-1]) 1];








end