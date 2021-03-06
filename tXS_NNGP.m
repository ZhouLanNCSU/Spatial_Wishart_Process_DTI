function v = tXS(z_BigX_sub,A_inv_sub,RR_w_F_all,dof,sigma_w,ndim)

index=1:ndim;
A_inv_sub_diag=sparse(index,index,A_inv_sub);
Siginv=mmtimes(A_inv_sub_diag,sparse(diag(1./RR_w_F_all)),A_inv_sub_diag);

v=z_BigX_sub'*Siginv*dof/sigma_w;


end

