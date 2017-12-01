function value = likihood_R_Xi_diag2(s,t_central,R_Xi_Choleksy_Lower_inv,R_Xi_det,sigma_Xi)
t_central=t_central(:,s);
Rchol_inv_t_central=R_Xi_Choleksy_Lower_inv*t_central;

value=-0.5*log(R_Xi_det)-0.5*Rchol_inv_t_central'*Rchol_inv_t_central/sigma_Xi;

end