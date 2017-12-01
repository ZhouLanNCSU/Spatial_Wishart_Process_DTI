function [Bookkeeping_Current_mapObj]=Log_Gamma_V1_MCMC_R_Xi_ModifiedModel_NNGP(it,Bookkeeping_Input_mapObj,Bookkeeping_Data_mapObj, Bookkeeping_Priors_mapObj, Bookkeeping_Current_mapObj)
     %% Extracting Variables
     %Input
     grid=Bookkeeping_Input_mapObj('grid');
     ndim=Bookkeeping_Input_mapObj('ndim');
     nsub=Bookkeeping_Input_mapObj('nsub');
     npred=Bookkeeping_Input_mapObj('npred');
     X=Bookkeeping_Input_mapObj('X');
     BigX=Bookkeeping_Input_mapObj('BigX');
     graph=Bookkeeping_Input_mapObj('graph');
     %Data
     t11=Bookkeeping_Data_mapObj('t11');
     t22=Bookkeeping_Data_mapObj('t22');
     t33=Bookkeeping_Data_mapObj('t33');
     log_t11=Bookkeeping_Data_mapObj('log_t11');
     log_t22=Bookkeeping_Data_mapObj('log_t22');
     log_t33=Bookkeeping_Data_mapObj('log_t33');
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
     sd_beta=Bookkeeping_Priors_mapObj('sd_beta');
     
     
     
     %Current
     Current_dof_mapObj=Bookkeeping_Current_mapObj('Current_dof_mapObj');
     Current_theta_w_mapObj=Bookkeeping_Current_mapObj('Current_theta_w_mapObj');
     Current_theta_Xi_mapObj=Bookkeeping_Current_mapObj('Current_theta_Xi_mapObj');
     Current_beta_mapObj=Bookkeeping_Current_mapObj('Current_beta_mapObj');
     Current_scale_mapObj=Bookkeeping_Current_mapObj('Current_scale_mapObj');
     Current_R_Xi_mapObj=Bookkeeping_Current_mapObj('Current_R_Xi_mapObj');
     Current_RR_w_mapObj=Bookkeeping_Current_mapObj('Current_RR_w_mapObj');
     Current_others_mapObj=Bookkeeping_Current_mapObj('Current_others_mapObj');
     %dof
     dof=Current_dof_mapObj('dof');
     %theta_w
     rho_w=Current_theta_w_mapObj('rho_w');
     nu_w=Current_theta_w_mapObj('nu_w');
     sigma_w=Current_theta_w_mapObj('sigma_w');
     nugget=Current_theta_w_mapObj('nugget');
     %theta_Xi
     rho_Xi=Current_theta_Xi_mapObj('rho_Xi');
     nu_Xi=Current_theta_Xi_mapObj('nu_Xi');
     sigma_Xi11=Current_theta_Xi_mapObj('sigma_Xi11');
     sigma_Xi22=Current_theta_Xi_mapObj('sigma_Xi22');
     sigma_Xi33=Current_theta_Xi_mapObj('sigma_Xi33');
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
     %R_Xi
     R_Xi_B_all=Current_R_Xi_mapObj('R_Xi_B_all');
     R_Xi_F_all=Current_R_Xi_mapObj('R_Xi_F_all');
     
     log_t11_central=Current_R_Xi_mapObj('log_t11_central');
     log_t22_central=Current_R_Xi_mapObj('log_t22_central');
     log_t33_central=Current_R_Xi_mapObj('log_t33_central');
     
     log_t11_central_for_blk=Current_R_Xi_mapObj('log_t11_central_for_blk');
     log_t22_central_for_blk=Current_R_Xi_mapObj('log_t22_central_for_blk');
     log_t33_central_for_blk=Current_R_Xi_mapObj('log_t33_central_for_blk');
     
     log_t11_central_NNGP=Current_R_Xi_mapObj('log_t11_central_NNGP');
     log_t22_central_NNGP=Current_R_Xi_mapObj('log_t22_central_NNGP');
     log_t33_central_NNGP=Current_R_Xi_mapObj('log_t33_central_NNGP');

     %RR_w
     RR_w_B_all=Current_RR_w_mapObj('RR_w_B_all');
     RR_w_F_all=Current_RR_w_mapObj('RR_w_F_all');
     t21_central_NNGP=Current_RR_w_mapObj('t21_central_NNGP');
     t31_central_NNGP=Current_RR_w_mapObj('t31_central_NNGP');
     t32_central_NNGP=Current_RR_w_mapObj('t32_central_NNGP');
     
     %others
     A_inv22=Current_others_mapObj('A_inv22');
     A_inv33=Current_others_mapObj('A_inv33');
     z11=Current_others_mapObj('z11');
     z22=Current_others_mapObj('z22');
     t21_central=Current_others_mapObj('t21_central');
     t31_central=Current_others_mapObj('t31_central');
     t32_central=Current_others_mapObj('t32_central');
     
     
     
     %% Update rho_Xi nu_Xi
     rho_Xi_can=exp(normrnd(log(rho_Xi),0.1));
     nu_Xi_can=exp(normrnd(log(nu_Xi),0.1));
     
     if rho_Xi_can>=5
         rho_Xi_can=5;
     end
     
     if nu_Xi_can>=2
         nu_Xi_can=2;
     end
     
     R_Xi_can=correlation_matern(grid, nu_Xi_can, rho_Xi_can);
     
     [R_Xi_B_all_can R_Xi_F_all_can]=NNGP_COV(graph,R_Xi_can,ndim);
     
     log_t11_central_can_NNGP=NNGP_MEAN(R_Xi_B_all_can,log_t11_central,log_t11_central_for_blk,graph,nsub,ndim);
     log_t22_central_can_NNGP=NNGP_MEAN(R_Xi_B_all_can,log_t22_central,log_t22_central_for_blk,graph,nsub,ndim);
     log_t33_central_can_NNGP=NNGP_MEAN(R_Xi_B_all_can,log_t33_central,log_t33_central_for_blk,graph,nsub,ndim);

     

     % log central t11 t22 t33
     t_central_all=cell(3,1);
     t_central_all(1)={log_t11_central_NNGP};
     t_central_all(2)={log_t22_central_NNGP};
     t_central_all(3)={log_t33_central_NNGP};
     
     t_central_all_can=cell(3,1);
     t_central_all_can(1)={log_t11_central_can_NNGP};
     t_central_all_can(2)={log_t22_central_can_NNGP};
     t_central_all_can(3)={log_t33_central_can_NNGP};
     
     sigma_Xi_all=[sigma_Xi11 sigma_Xi22 sigma_Xi33];
     
     cankk=[0 0 0];
     oldkk=[0 0 0];
     
     for dindex=1:3
        t_central=cell2mat(t_central_all(dindex));
        t_central_can=cell2mat(t_central_all_can(dindex));
        cankk(dindex)=sum(arrayfun(@(s) likelihood_dof_NNGP( s,t_central_can, R_Xi_F_all_can,sigma_Xi_all(dindex),ndim),[1:nsub]));
        oldkk(dindex)=sum(arrayfun(@(s) likelihood_dof_NNGP( s,t_central, R_Xi_F_all,sigma_Xi_all(dindex),ndim),[1:nsub]));         
     end

     
     % Compute prob
     can=sum(cankk)+log(normpdf(log(rho_Xi_can),mean_range,sd_range));
     old=sum(oldkk)+log(normpdf(log(rho_Xi),mean_range,sd_range));
     prob=min(1,exp(can-old));
     Accept=randsample([0,1],1,true,[1-prob,prob]);
     if Accept==1
         rho_Xi=rho_Xi_can;
         nu_Xi=nu_Xi_can;
         R_Xi=R_Xi_can;
         R_Xi_B_all=R_Xi_B_all_can;
         R_Xi_F_all=R_Xi_F_all_can;
         log_t11_central_NNGP=log_t11_central_can_NNGP;
         log_t22_central_NNGP=log_t22_central_can_NNGP;
         log_t33_central_NNGP=log_t33_central_can_NNGP;
     end
     
     
%      %% Update nu_Xi
%      nu_Xi_can=exp(normrnd(log(nu_Xi),0.1));
%      R_Xi_can=correlation_matern(grid, nu_Xi_can, rho_Xi);
%      R_Xi_Choleksy_Lower_can=sparse(chol(R_Xi_can,'lower'));
%      R_Xi_Choleksy_Lower_inv_can=inv(R_Xi_Choleksy_Lower_can);
%      R_Xi_det_can=prod(diag(R_Xi_Choleksy_Lower_can))^2;
%      
%      % log central t11 t22 t33
%      t_central_all=cell(3,1);
%      t_central_all(1)={log_t11_central};
%      t_central_all(2)={log_t22_central};
%      t_central_all(3)={log_t33_central};
%      
%      sigma_Xi_all=[sigma_Xi11 sigma_Xi22 sigma_Xi33];
%      
%      can=[0 0 0];
%      old=[0 0 0];
%      
%      for dindex=1:3
%         t_central=cell2mat(t_central_all(dindex));
%         can(dindex)=sum(arrayfun(@(s) likelihood_R_w_off( s,t_central,repmat(1 ,[ndim nsub]), R_Xi_Choleksy_Lower_inv_can,R_Xi_det_can,1,sigma_Xi_all(dindex)),[1:nsub]));
%         old(dindex)=sum(arrayfun(@(s) likelihood_R_w_off( s,t_central,repmat(1 ,[ndim nsub]), R_Xi_Choleksy_Lower_inv,R_Xi_det,1,sigma_Xi_all(dindex)),[1:nsub]));         
%      end
%      
%      % Compute prob
%      can=sum(can)+log(normpdf(log(nu_Xi_can),mean_range,sd_range));
%      old=sum(old)+log(normpdf(log(nu_Xi),mean_range,sd_range));
%      prob=min(1,exp(can-old));
%      Accept=randsample([0,1],1,true,[1-prob,prob]);
%      if Accept==1
%          nu_Xi=nu_Xi_can;
%          R_Xi=R_Xi_can;
%          R_Xi_Choleksy_Lower=R_Xi_Choleksy_Lower_can;
%          R_Xi_Choleksy_Lower_inv=R_Xi_Choleksy_Lower_inv_can;
%          R_Xi_det=R_Xi_det_can;
%      end
%      
%      
%      
     
     %% Finalizing
     Current_theta_Xi_mapObj('rho_Xi')=rho_Xi;    
     Current_theta_Xi_mapObj('nu_Xi')=nu_Xi;

     
     if mod(it,20)==0
     ['rho_Xi']
     rho_Xi
     ['nu_Xi']
     nu_Xi
     end
     
      Current_R_Xi_mapObj('R_Xi_B_all')=R_Xi_B_all;
      Current_R_Xi_mapObj('R_Xi_F_all')=R_Xi_F_all;
      Current_R_Xi_mapObj('log_t11_central_NNGP')=log_t11_central_NNGP;
      Current_R_Xi_mapObj('log_t22_central_NNGP')=log_t22_central_NNGP;
      Current_R_Xi_mapObj('log_t33_central_NNGP')=log_t33_central_NNGP;
      
      
     Bookkeeping_Current_mapObj('Current_theta_Xi_mapObj')=Current_theta_Xi_mapObj;
     Bookkeeping_Current_mapObj('Current_R_Xi_mapObj')=Current_R_Xi_mapObj;


 



end