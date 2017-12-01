function value = likelihood_dof_diag(s,u,R_c_Choleksy_Lower_inv,Jocobian )

u_ind=u(:,s);
uu=u_ind'*u_ind;
Rchol_inv_u=R_c_Choleksy_Lower_inv*u_ind;


value=-0.5*(Rchol_inv_u'*Rchol_inv_u-uu)+sum(Jocobian(:,s));

end

