function [B]=B_generator(s,graph,R)

B=sparse(R(s,find(graph(s,:)~=0))*inv(R(find(graph(s,:)~=0),find(graph(s,:)~=0))));



end