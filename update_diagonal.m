% out: updated matrix, n: size of the input matrix, shift: diagonal selection
function out = update_diagonal(n,shift)
A=magic(n)
result=diag(A,shift);
% Create matrix of ones along specified diagonal
iden = diag(ones(n-abs(shift),1),shift);
% Modify the elements at location of ones in 'iden'. Here, increment all elements
% along the specified diagonal by 2
 A(iden(:,:)~=0) = A(iden(:,:)~=0)+2;
 out = A;