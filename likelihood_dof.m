function value=likelihood_dof( s,t_central, R_Xi_Choleksy_Lower_inv,R_Xi_det,sigma_Xi,ndim)


t_central_sub=t_central(:,s);
Rchol_inv_t_central=R_Xi_Choleksy_Lower_inv*t_central_sub;

value=-0.5*Rchol_inv_t_central'*Rchol_inv_t_central/sigma_Xi-0.5*ndim*log(sigma_Xi);



end