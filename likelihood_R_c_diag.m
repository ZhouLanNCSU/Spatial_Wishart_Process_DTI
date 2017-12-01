function value= likelihood_R_c_diag(s,u,R_c_det,R_c_Choleksy_Lower_inv)

u_ind=u(:,s);
Rchol_inv_u=R_c_Choleksy_Lower_inv*u_ind;

value=-0.5*log(R_c_det)-0.5*(Rchol_inv_u'*Rchol_inv_u);

end

