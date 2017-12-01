function value = likelihood_dof_off_NNGP( s,t_central,A_inv, RR_w_F_all)

t_central_A_inv=t_central(:,s).*A_inv(:,s);

value=-0.5*sum(t_central_A_inv.^2./RR_w_F_all');

end