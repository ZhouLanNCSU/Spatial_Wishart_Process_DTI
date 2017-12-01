function [A_inv_partial] = off_cross_diag( A_inv,dim_index ,index,ndim,nsub)

A_inv_partial=A_inv(index);

%A_inv(dim_index,:)=[];

%A_inv_cross=repmat(A_inv_partial, [ndim-1 1]).*A_inv;



% A_inv_partial=A_inv(index);
% 
% AA=A_inv([1:dim_index-1 dim_index+1:ndim],:);
% 
% A_inv_cross=repmat(A_inv_partial, [ndim-1 1]).*AA;




% A_inv_partial=A_inv(index);
% 
% 
% A = zeros(size(A_inv) - [numel(index) 0]);
% r = true(size(A_inv,1),1);
% c = true(size(A_inv,2),1);
% r(index) = false;
% c([]) = false;
% A = A_inv(r,c);
% 
% A_inv_cross=repmat(A_inv_partial, [ndim-1 1]).*A;


end

