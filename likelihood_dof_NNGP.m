function value=likelihood_dof( s,t_central, R_Xi_F_all,sigma_Xi,ndim)


t_central_sub=t_central(:,s);


value=-0.5*sum(t_central_sub.^2./R_Xi_F_all'/sigma_Xi)-0.5*ndim*log(sigma_Xi);



end