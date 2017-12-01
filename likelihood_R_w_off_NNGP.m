function value=likelihood_R_w_off_NNGP( s,t_central,A_inv,RR_w_F_all ,dof,sigma_w)

t_central_A_inv=t_central(:,s).*A_inv(:,s);

value=-0.5*sum(t_central_A_inv.^2./RR_w_F_all')*dof/sigma_w-0.5*log(prod(RR_w_F_all));



end