function [Bookkeeping_Current_mapObj]=MCMC_R_c_ModifiedModel(it,Bookkeeping_Input_mapObj,Bookkeeping_Data_mapObj, Bookkeeping_Priors_mapObj, Bookkeeping_Current_mapObj);

     %% Extracting Variables
     %Input
     grid=Bookkeeping_Input_mapObj('grid');
     ndim=Bookkeeping_Input_mapObj('ndim');
     nsub=Bookkeeping_Input_mapObj('nsub');
     npred=Bookkeeping_Input_mapObj('npred');
     X=Bookkeeping_Input_mapObj('X');
     BigX=Bookkeeping_Input_mapObj('BigX');
     %Data
     t11=Bookkeeping_Data_mapObj('t11');
     t22=Bookkeeping_Data_mapObj('t22');
     t33=Bookkeeping_Data_mapObj('t33');
     t21=Bookkeeping_Data_mapObj('t21');
     t31=Bookkeeping_Data_mapObj('t31');
     t32=Bookkeeping_Data_mapObj('t32');
     %Prior
     mean_nu=Bookkeeping_Priors_mapObj('mean_nu');
     sd_nu=Bookkeeping_Priors_mapObj('sd_nu');
     mean_range=Bookkeeping_Priors_mapObj('mean_range');
     sd_range=Bookkeeping_Priors_mapObj('sd_range');
     mean_r=Bookkeeping_Priors_mapObj('mean_r');
     sd_r=Bookkeeping_Priors_mapObj('sd_r');
     a_var=Bookkeeping_Priors_mapObj('a_var');
     b_var=Bookkeeping_Priors_mapObj('b_var');
 
     
     
     %Current
     Current_dof_mapObj=Bookkeeping_Current_mapObj('Current_dof_mapObj');
     Current_theta_w_mapObj=Bookkeeping_Current_mapObj('Current_theta_w_mapObj');
     Current_theta_c_mapObj=Bookkeeping_Current_mapObj('Current_theta_c_mapObj');
     Current_beta_mapObj=Bookkeeping_Current_mapObj('Current_beta_mapObj');
     Current_scale_mapObj=Bookkeeping_Current_mapObj('Current_scale_mapObj');
     Current_quantile_mapObj=Bookkeeping_Current_mapObj('Current_quantile_mapObj');
     Current_norminv_mapObj=Bookkeeping_Current_mapObj('Current_norminv_mapObj');
     Current_Jocobian_mapObj=Bookkeeping_Current_mapObj('Current_Jocobian_mapObj');
     Current_R_c_mapObj=Bookkeeping_Current_mapObj('Current_R_c_mapObj');
     Current_RR_w_mapObj=Bookkeeping_Current_mapObj('Current_RR_w_mapObj');
     Current_others_mapObj=Bookkeeping_Current_mapObj('Current_others_mapObj');
     %dof
     dof=Current_dof_mapObj('dof');
     %theta_w
     rho_w=Current_theta_w_mapObj('rho_w');
     nu_w=Current_theta_w_mapObj('nu_w');
     sigma_w=Current_theta_w_mapObj('sigma_w');
     nugget=Current_theta_w_mapObj('nugget');
     %theta_c
     rho_c=Current_theta_c_mapObj('rho_c');
     nu_c=Current_theta_c_mapObj('nu_c');
     %beta
     beta11=Current_beta_mapObj('beta11');
     beta22=Current_beta_mapObj('beta22');
     beta33=Current_beta_mapObj('beta33');
     beta21=Current_beta_mapObj('beta21');
     beta31=Current_beta_mapObj('beta31');
     beta32=Current_beta_mapObj('beta32');
     %scale
     scale11=Current_scale_mapObj('scale11');
     scale22=Current_scale_mapObj('scale22');
     scale33=Current_scale_mapObj('scale33');
     %quantile
     q_t11=Current_quantile_mapObj('q_t11');
     q_t22=Current_quantile_mapObj('q_t22');
     q_t33=Current_quantile_mapObj('q_t33');
     %norminv
     u_t11=Current_norminv_mapObj('u_t11');
     u_t22=Current_norminv_mapObj('u_t22');
     u_t33=Current_norminv_mapObj('u_t33');
     %Jocobian 
     Jocobian_t11=Current_Jocobian_mapObj('Jocobian_t11');
     Jocobian_t22=Current_Jocobian_mapObj('Jocobian_t22');
     Jocobian_t33=Current_Jocobian_mapObj('Jocobian_t33');
     %R_c
     R_c_Choleksy_Lower=Current_R_c_mapObj('R_c_Choleksy_Lower');
     R_c_Choleksy_Lower_inv=Current_R_c_mapObj('R_c_Choleksy_Lower_inv');
     R_c_det=Current_R_c_mapObj('R_c_det');
     %RR_w
     RR_w_Choleksy_Lower=Current_RR_w_mapObj('RR_w_Choleksy_Lower');
     RR_w_Choleksy_Lower_inv=Current_RR_w_mapObj('RR_w_Choleksy_Lower_inv');
     RR_w_inv=Current_RR_w_mapObj('RR_w_inv');
     RR_w_det=Current_RR_w_mapObj('RR_w_det');
     %others
     A_inv22=Current_others_mapObj('A_inv22');
     A_inv33=Current_others_mapObj('A_inv33');
     z11=Current_others_mapObj('z11');
     z22=Current_others_mapObj('z22');
     t21_central=Current_others_mapObj('t21_central');
     t31_central=Current_others_mapObj('t31_central');
     t32_central=Current_others_mapObj('t32_central');
     
     %% Update rho_c nu_c
     rho_c_can=exp(normrnd(log(rho_c),0.1));
     nu_c_can=0.5;
     if rho_c_can>=5
         rho_c_can=5;
     end
     %nu_c_can=exp(normrnd(log(nu_c),0.1));
     R_c_can=correlation_matern(grid, nu_c_can, rho_c_can);
     R_c_Choleksy_Lower_can=sparse(chol(R_c_can,'lower'));
     R_c_Choleksy_Lower_inv_can=inv(R_c_Choleksy_Lower_can);
     R_c_det_can=prod(diag(R_c_Choleksy_Lower_can))^2;
     
     can11=sum(arrayfun(@(s) likelihood_R_c_diag(s,u_t11,R_c_det_can,R_c_Choleksy_Lower_inv_can),[1:nsub]));
     can22=sum(arrayfun(@(s) likelihood_R_c_diag(s,u_t22,R_c_det_can,R_c_Choleksy_Lower_inv_can),[1:nsub]));
     can33=sum(arrayfun(@(s) likelihood_R_c_diag(s,u_t33,R_c_det_can,R_c_Choleksy_Lower_inv_can),[1:nsub]));
     
     old11=sum(arrayfun(@(s) likelihood_R_c_diag(s,u_t11,R_c_det,R_c_Choleksy_Lower_inv),[1:nsub]));
     old22=sum(arrayfun(@(s) likelihood_R_c_diag(s,u_t22,R_c_det,R_c_Choleksy_Lower_inv),[1:nsub]));
     old33=sum(arrayfun(@(s) likelihood_R_c_diag(s,u_t33,R_c_det,R_c_Choleksy_Lower_inv),[1:nsub]));

     % Compute prob
     can=can11+can22+can33+log(normpdf(log(rho_c_can),mean_range,sd_range));
     old=old11+old22+old33+log(normpdf(log(rho_c),mean_range,sd_range));
     prob=min(1,exp(can-old));
     Accept=randsample([0,1],1,true,[1-prob,prob]);
     
     if Accept==1
         rho_c=rho_c_can;
         nu_c=nu_c_can;
         R_c=R_c_can;
         R_c_Choleksy_Lower=R_c_Choleksy_Lower_can;
         R_c_Choleksy_Lower_inv=R_c_Choleksy_Lower_inv_can;
         R_c_det=R_c_det_can;
     end
     
     
%      %% Update nu_c
%      nu_c_can=exp(normrnd(log(nu_c),0.1));
%      R_c_can=correlation_matern(grid, nu_c_can, rho_c);
%      R_c_Choleksy_Lower_can=sparse(chol(R_c_can,'lower'));
%      R_c_Choleksy_Lower_inv_can=inv(R_c_Choleksy_Lower_can);
%      R_c_det_can=prod(diag(R_c_Choleksy_Lower_can))^2;
%      
%      can11=sum(arrayfun(@(s) likelihood_R_c_diag(s,u_t11,R_c_det_can,R_c_Choleksy_Lower_inv_can),[1:nsub]));
%      can22=sum(arrayfun(@(s) likelihood_R_c_diag(s,u_t22,R_c_det_can,R_c_Choleksy_Lower_inv_can),[1:nsub]));
%      can33=sum(arrayfun(@(s) likelihood_R_c_diag(s,u_t33,R_c_det_can,R_c_Choleksy_Lower_inv_can),[1:nsub]));
%      
%      old11=sum(arrayfun(@(s) likelihood_R_c_diag(s,u_t11,R_c_det,R_c_Choleksy_Lower_inv),[1:nsub]));
%      old22=sum(arrayfun(@(s) likelihood_R_c_diag(s,u_t22,R_c_det,R_c_Choleksy_Lower_inv),[1:nsub]));
%      old33=sum(arrayfun(@(s) likelihood_R_c_diag(s,u_t33,R_c_det,R_c_Choleksy_Lower_inv),[1:nsub]));
% 
%      %Compute prob
%      can=can11+can22+can33+log(normpdf(log(nu_c_can),mean_nu,sd_nu));
%      old=old11+old22+old33+log(normpdf(log(nu_c),mean_nu,sd_nu));
%      prob=min(1,exp(can-old));
%      Accept=randsample([0,1],1,true,[1-prob,prob]);
%      
%      if Accept==1
%          nu_c=nu_c_can;
%          R_c=R_c_can;
%          R_c_Choleksy_Lower=R_c_Choleksy_Lower_can;
%          R_c_Choleksy_Lower_inv=R_c_Choleksy_Lower_inv_can;
%          R_c_det=R_c_det_can;
%      end
%      
     
     %% Finalizing
     Current_theta_c_mapObj('rho_c')=rho_c;
    
     Current_theta_c_mapObj('nu_c')=nu_c;
     
     if mod(it,20)==0
      ['rho_c']
     rho_c
     ['nu_c']
     nu_c
     end
     
     Current_R_c_mapObj('R_c_Choleksy_Lower')=R_c_Choleksy_Lower;
     Current_R_c_mapObj('R_c_Choleksy_Lower_inv')=R_c_Choleksy_Lower_inv;
     Current_R_c_mapObj('R_c_det')=R_c_det;
     Current_R_c_mapObj('R_c_inv')=full(R_c_Choleksy_Lower_inv'*R_c_Choleksy_Lower_inv);
     
     Bookkeeping_Current_mapObj('Current_theta_c_mapObj')=Current_theta_c_mapObj;
     Bookkeeping_Current_mapObj('Current_R_c_mapObj')=Current_R_c_mapObj;
     
     

end

