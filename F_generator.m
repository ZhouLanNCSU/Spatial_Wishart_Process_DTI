function [F]=F_generator(s,graph,R,B_all)


F=R(s,s)-B_all{s}*R(find(graph(s,:)~=0),s);


end