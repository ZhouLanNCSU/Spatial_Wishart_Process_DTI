function [Bookkeeping_Current_mapObj]=MCMC_dof_ModifiedModel_NNGP(it,Bookkeeping_Input_mapObj,Bookkeeping_Data_mapObj, Bookkeeping_Priors_mapObj, Bookkeeping_Current_mapObj)

     %% Extracting Variables
     %Input
     grid=Bookkeeping_Input_mapObj('grid');
     ndim=Bookkeeping_Input_mapObj('ndim');
     nsub=Bookkeeping_Input_mapObj('nsub');
     npred=Bookkeeping_Input_mapObj('npred');
     X=Bookkeeping_Input_mapObj('npred');
     graph=Bookkeeping_Input_mapObj('graph');
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
     R_c_B_all=Current_R_c_mapObj('R_c_B_all');
     R_c_F_all=Current_R_c_mapObj('R_c_F_all');
     u_t11_central_NNGP=Current_R_c_mapObj('u_t11_central_NNGP');
     u_t22_central_NNGP=Current_R_c_mapObj('u_t22_central_NNGP');
     u_t33_central_NNGP=Current_R_c_mapObj('u_t33_central_NNGP');
     
     u_t11_for_blk=Current_R_c_mapObj('u_t11_for_blk');
     u_t22_for_blk=Current_R_c_mapObj('u_t22_for_blk');
     u_t33_for_blk=Current_R_c_mapObj('u_t33_for_blk');
     
 
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
     
     
     
          %% Generating candidate dof
     dof_can=normrnd(dof,0.1);
     if dof_can<=3
       dof_can=3;
     end
     
     %% t11 t22 t33
     q_t11_can=gamcdf(t11,repmat(dof_can-0,[ndim,nsub])/2,1/dof_can*2*reshape(scale11,[ndim,nsub]));
     q_t22_can=gamcdf(t22,repmat(dof_can-1,[ndim,nsub])/2,1/dof_can*2*reshape(scale22,[ndim,nsub]));
     q_t33_can=gamcdf(t33,repmat(dof_can-2,[ndim,nsub])/2,1/dof_can*2*reshape(scale33,[ndim,nsub]));     
     q_t11_can(find(q_t11_can==1))=0.99999999999999994;
     q_t22_can(find(q_t22_can==1))=0.99999999999999994;
     q_t33_can(find(q_t33_can==1))=0.99999999999999994;
     u_t11_can=norminv(q_t11_can);
     u_t22_can=norminv(q_t22_can);
     u_t33_can=norminv(q_t33_can);
     Jocobian_t11_can=log(gampdf(t11,repmat(dof_can-0,[ndim,nsub])/2,1/dof_can*2*reshape(scale11,[ndim,nsub])));
     Jocobian_t22_can=log(gampdf(t22,repmat(dof_can-1,[ndim,nsub])/2,1/dof_can*2*reshape(scale22,[ndim,nsub])));
     Jocobian_t33_can=log(gampdf(t33,repmat(dof_can-2,[ndim,nsub])/2,1/dof_can*2*reshape(scale33,[ndim,nsub])));

     t=cell(3,1);
     t(1)={u_t11_can};
     t(2)={u_t22_can};
     t(3)={u_t33_can};
     myblk=cell(3,1);
     parfor i=1:3
          myblk(i)={bkl(t{i},graph,ndim,nsub)};
     end
     
u_t11_for_blk_can=myblk{1};
u_t22_for_blk_can=myblk{2};
u_t33_for_blk_can=myblk{3};

%      u_t11_for_blk_can=bkl(u_t11_can,graph,ndim,nsub);
%      u_t22_for_blk_can=bkl(u_t22_can,graph,ndim,nsub);
%      u_t33_for_blk_can=bkl(u_t33_can,graph,ndim,nsub);

     
     u_t11_central_NNGP_can=NNGP_MEAN(R_c_B_all,u_t11,u_t11_for_blk_can,graph,nsub,ndim);
     u_t22_central_NNGP_can=NNGP_MEAN(R_c_B_all,u_t22,u_t22_for_blk_can,graph,nsub,ndim);
     u_t33_central_NNGP_can=NNGP_MEAN(R_c_B_all,u_t33,u_t33_for_blk_can,graph,nsub,ndim);
     
    
     %candidate
     can11=sum(arrayfun(@(s) likelihood_dof_diag_NNGP(s,u_t11_can,u_t11_central_NNGP_can,R_c_F_all,Jocobian_t11_can,ndim),[1:nsub]));
     can22=sum(arrayfun(@(s) likelihood_dof_diag_NNGP(s,u_t22_can,u_t22_central_NNGP_can,R_c_F_all,Jocobian_t22_can,ndim),[1:nsub]));
     can33=sum(arrayfun(@(s) likelihood_dof_diag_NNGP(s,u_t33_can,u_t33_central_NNGP_can,R_c_F_all,Jocobian_t33_can,ndim),[1:nsub]));
     %old
     old11=sum(arrayfun(@(s) likelihood_dof_diag_NNGP(s,u_t11,u_t11_central_NNGP,R_c_F_all,Jocobian_t11,ndim),[1:nsub]));
     old22=sum(arrayfun(@(s) likelihood_dof_diag_NNGP(s,u_t22,u_t22_central_NNGP,R_c_F_all,Jocobian_t22,ndim),[1:nsub]));
     old33=sum(arrayfun(@(s) likelihood_dof_diag_NNGP(s,u_t33,u_t33_central_NNGP,R_c_F_all,Jocobian_t33,ndim),[1:nsub]));
     
%      can11=0;
%      can22=0;
%      can33=0;
%      old11=0;
%      old22=0;
%      old33=0;
%      
     
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
     
     
     
     %% Compute prob
     can=can11+can22+can33+can21+can31+can32;
     old=old11+old22+old33+old21+old31+old32;
     prob=min(1,exp(can-old));
     Accept=randsample([0,1],1,true,[1-prob,prob]);
     
     if Accept==1
         Current_dof_mapObj('dof')=dof_can;
         
         Current_quantile_mapObj('q_t11')=q_t11_can;
         Current_quantile_mapObj('q_t22')=q_t22_can;
         Current_quantile_mapObj('q_t33')=q_t33_can;
         
         Current_norminv_mapObj('u_t11')=u_t11_can;
         Current_norminv_mapObj('u_t22')=u_t22_can;
         Current_norminv_mapObj('u_t33')=u_t33_can;
         
         Current_R_c_mapObj('u_t11_central_NNGP')=u_t11_central_NNGP_can;
         Current_R_c_mapObj('u_t22_central_NNGP')=u_t22_central_NNGP_can;
         Current_R_c_mapObj('u_t33_central_NNGP')=u_t33_central_NNGP_can;
         
         Current_R_c_mapObj('u_t11_for_blk')=u_t11_for_blk_can;
         Current_R_c_mapObj('u_t22_for_blk')=u_t22_for_blk_can;
         Current_R_c_mapObj('u_t33_for_blk')=u_t33_for_blk_can;
         
    
         
         
         Current_Jocobian_mapObj('Jocobian_t11')=Jocobian_t11_can;
         Current_Jocobian_mapObj('Jocobian_t22')=Jocobian_t22_can;
         Current_Jocobian_mapObj('Jocobian_t33')=Jocobian_t33_can;
         
         Bookkeeping_Current_mapObj('Current_dof_mapObj')=Current_dof_mapObj;
         Bookkeeping_Current_mapObj('Current_quantile_mapObj')=Current_quantile_mapObj;
         Bookkeeping_Current_mapObj('Current_norminv_mapObj')=Current_norminv_mapObj;
         Bookkeeping_Current_mapObj('Current_Jocobian_mapObj')=Current_Jocobian_mapObj;
         Bookkeeping_Current_mapObj('Current_R_c_mapObj')=Current_R_c_mapObj;
     end
     
     if mod(it,20)==0
         ['dof']
     Current_dof_mapObj('dof')
     end


     


end