function value = likelihood_R_b( beta_all,R_b_det,R_b_Choleksy_Lower_inv)

[n p]=size(beta_all);

half= arrayfun(@(pp)  R_b_Choleksy_Lower_inv*beta_all(:,pp),[1:p],'UniformOutput',false);
half_prod_sum=sum(cellfun(@inner,half));

value=-0.5*log(R_b_det)*p-0.5*half_prod_sum;


end

