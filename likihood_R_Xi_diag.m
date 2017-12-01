function value = likihood_R_Xi_diag(s,t_central,R_Xi_Choleksy_Lower_inv)
t_central=t_central(:,s);
Rchol_inv_t_central=R_Xi_Choleksy_Lower_inv*t_central;

value=-0.5*Rchol_inv_t_central'*Rchol_inv_t_central;

end

