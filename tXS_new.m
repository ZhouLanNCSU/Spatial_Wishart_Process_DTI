function v = tXS_new(z_X_sub,A_inv_sub,RR_w_inv,dof,ndim,npred)

index=1:ndim;
A_inv_sub_diag=sparse(index,index,A_inv_sub);
Siginv=mmtimes(A_inv_sub_diag,RR_w_inv,A_inv_sub_diag);

z_X_sub=z_X_sub';
z_X_sub_vec=repmat(z_X_sub(:),  [1 ndim]);

v=z_X_sub_vec.*kron(Siginv,ones(npred,1))*dof;


end

