function [Bookkeeping_Current_mapObj]=Log_Gamma_V1_MCMC_R_dof_ModifiedModel(it,Bookkeeping_Input_mapObj,Bookkeeping_Data_mapObj, Bookkeeping_Priors_mapObj, Bookkeeping_Current_mapObj)
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
     log_t22_central_NNGP=Current_R_Xi_mapObj('log_t22_central_NNGP');
     log_t33_central_NNGP=Current_R_Xi_mapObj('log_t33_central_NNGP');
     
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

     
     %% Generating candidate dof
     dof_can=normrnd(dof,0.1);
     if dof_can<=9.5
       dof_can=9.5;
     elseif dof_can>=10.5
        dof_can=10.5;
     end
     %% t21 t31 t32

     t_central_all=cell(3,1);
     t_central_all(1)={t21_central_NNGP};
     t_central_all(2)={t31_central_NNGP};
     t_central_all(3)={t32_central_NNGP};

     beta_all=cell(3,1);
     beta_all(1)={beta21};
     beta_all(2)={beta31};
     beta_all(3)={beta32};
     
     A_all=cell(3,1);
     A_all(1)={A_inv22};
     A_all(2)={A_inv33};
     A_all(3)={A_inv33};
     
     middle_part=[0 0 0];
     
     for dindex=1:3
        t_central=cell2mat(t_central_all(dindex));
        A_inv=cell2mat(A_all(dindex));
        middle_part(dindex)=sum(arrayfun(@(s) likelihood_dof_off_NNGP(s,t_central,A_inv, RR_w_F_all),[1:nsub]));         
     end
     
     log_dof_part_can=-0.5*nsub*ndim*log(1/dof_can);
     can21=log_dof_part_can+middle_part(1)*dof_can/sigma_w;
     can31=log_dof_part_can+middle_part(2)*dof_can/sigma_w;
     can32=log_dof_part_can+middle_part(3)*dof_can/sigma_w;
     
     log_dof_part_old=-0.5*nsub*ndim*log(1/dof);
     old21=log_dof_part_old+middle_part(1)*dof/sigma_w;
     old31=log_dof_part_old+middle_part(2)*dof/sigma_w;
     old32=log_dof_part_old+middle_part(3)*dof/sigma_w;
     
     %% log t11 t22 t33
     sigma_Xi11_can=psi(1,(dof_can-0)/2);
     sigma_Xi22_can=psi(1,(dof_can-1)/2);
     sigma_Xi33_can=psi(1,(dof_can-2)/2);
     
     sigma_Xi_all_can=[sigma_Xi11_can sigma_Xi22_can sigma_Xi33_can];
     sigma_Xi_all=[sigma_Xi11 sigma_Xi22 sigma_Xi33];
     
     log_t11_central_can=log_t11-psi(0,(dof_can-0)/2)-log(2)+log(dof_can)-log(reshape(scale11,[ndim nsub]));
     log_t22_central_can=log_t22-psi(0,(dof_can-1)/2)-log(2)+log(dof_can)-log(reshape(scale22,[ndim nsub]));
     log_t33_central_can=log_t33-psi(0,(dof_can-2)/2)-log(2)+log(dof_can)-log(reshape(scale33,[ndim nsub]));
     
     log_t11_central_can_for_blk=bkl(log_t11_central_can,graph,ndim,nsub);
     log_t22_central_can_for_blk=bkl(log_t22_central_can,graph,ndim,nsub);
     log_t33_central_can_for_blk=bkl(log_t33_central_can,graph,ndim,nsub);
     
     log_t11_central_can_NNGP=NNGP_MEAN(R_Xi_B_all,log_t11_central_can,log_t11_central_can_for_blk,graph,nsub,ndim);
     log_t22_central_can_NNGP=NNGP_MEAN(R_Xi_B_all,log_t22_central_can,log_t22_central_can_for_blk,graph,nsub,ndim);
     log_t33_central_can_NNGP=NNGP_MEAN(R_Xi_B_all,log_t33_central_can,log_t33_central_can_for_blk,graph,nsub,ndim);

     
     % log central t11 t22 t33
     t_central_all=cell(3,1);
     t_central_all(1)={log_t11_central_NNGP};
     t_central_all(2)={log_t22_central_NNGP};
     t_central_all(3)={log_t33_central_NNGP};
     
     t_central_all_can=cell(3,1);
     t_central_all_can(1)={log_t11_central_can_NNGP};
     t_central_all_can(2)={log_t22_central_can_NNGP};
     t_central_all_can(3)={log_t33_central_can_NNGP};
     
     cankk=[0 0 0];
     oldkk=[0 0 0];
     
     for dindex=1:3
        t_central=cell2mat(t_central_all(dindex));
        t_central_can=cell2mat(t_central_all_can(dindex));
        cankk(dindex)=sum(arrayfun(@(s) likelihood_dof_NNGP( s,t_central_can, R_Xi_F_all,sigma_Xi_all_can(dindex),ndim),[1:nsub]));
        oldkk(dindex)=sum(arrayfun(@(s) likelihood_dof_NNGP( s,t_central, R_Xi_F_all,sigma_Xi_all(dindex),ndim),[1:nsub]));         
     end

     %% Compute prob
     can=sum(cankk)+can21+can31+can32;
     old=sum(oldkk)+old21+old31+old32;
     prob=min(1,exp(can-old));
     Accept=randsample([0,1],1,true,[1-prob,prob]);
     
     if Accept==1
         Current_dof_mapObj('dof')=dof_can;
         Current_theta_Xi_mapObj('sigma_Xi11')=sigma_Xi11_can;
         Current_theta_Xi_mapObj('sigma_Xi22')=sigma_Xi22_can;
         Current_theta_Xi_mapObj('sigma_Xi33')=sigma_Xi33_can;
         
         Current_R_Xi_mapObj('log_t11_central_for_blk')=log_t11_central_can_for_blk;
         Current_R_Xi_mapObj('log_t22_central_for_blk')=log_t22_central_can_for_blk;
         Current_R_Xi_mapObj('log_t33_central_for_blk')=log_t33_central_can_for_blk;
         
         Current_R_Xi_mapObj('log_t11_central_NNGP')=log_t11_central_can_NNGP;
         Current_R_Xi_mapObj('log_t22_central_NNGP')=log_t22_central_can_NNGP;
         Current_R_Xi_mapObj('log_t33_central_NNGP')=log_t33_central_can_NNGP;
         
         
  
         
         Bookkeeping_Current_mapObj('Current_dof_mapObj')=Current_dof_mapObj;
         Bookkeeping_Current_mapObj('Current_R_Xi_mapObj')=Current_R_Xi_mapObj;

     end
     
     if mod(it,20)==0
         ['dof']
         Current_dof_mapObj('dof')
     end

     



end