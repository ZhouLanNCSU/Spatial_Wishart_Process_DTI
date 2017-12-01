function value = likelihood_dof_diag_NNGP(s,u,u_NNGP,R_c_F_all,Jocobian ,ndim)

u_ind=u(1:ndim-1,s);
uu=u_ind'*u_ind;


u_ind_NNGP=u_NNGP(:,s);

value=-0.5*(sum(u_ind_NNGP.^2./R_c_F_all')-uu)+sum(Jocobian(:,s));

end

