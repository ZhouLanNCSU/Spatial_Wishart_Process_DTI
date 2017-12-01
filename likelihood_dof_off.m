function value = likelihood_dof_off( s,t_central,A_inv, RR_w_Choleksy_Lower_inv)

t_central_A_inv=t_central(:,s).*A_inv(:,s);
RRchol_inv_t_central_A_inv=RR_w_Choleksy_Lower_inv*t_central_A_inv;

value=-0.5*RRchol_inv_t_central_A_inv'*RRchol_inv_t_central_A_inv;

end

