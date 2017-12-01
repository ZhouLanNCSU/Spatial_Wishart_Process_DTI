function [Bookkeeping_Current_mapObj]=MCMC_beta_off_ModifiedModel_NNGP(it,Bookkeeping_Input_mapObj,Bookkeeping_Data_mapObj, Bookkeeping_Priors_mapObj, Bookkeeping_Current_mapObj)

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
     Current_R_b_mapObj=Bookkeeping_Current_mapObj('Current_R_b_mapObj');
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
     sigma_b=0.01;
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
     %R_b
     R_b_B_all=Current_R_b_mapObj('R_b_B_all');
     R_b_F_all=Current_R_b_mapObj('R_b_F_all');
     
     %R_c
     R_c_B_all=Current_R_c_mapObj('R_c_B_all');
     R_c_F_all=Current_R_c_mapObj('R_c_F_all');
     u_t11_central_NNGP=Current_R_c_mapObj('u_t11_central_NNGP');
     u_t22_central_NNGP=Current_R_c_mapObj('u_t22_central_NNGP');
     u_t33_central_NNGP=Current_R_c_mapObj('u_t33_central_NNGP');
     
 
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
     
     
          %% Un Pararell code
     % beta21
     z11_BigX=repmat(sqrt(reshape(z11,[ndim*nsub 1])),[1 ndim*npred]).*BigX;
     %z11_BigX=kron(repmat(1,[1 ndim*npred]),sqrt(reshape(z11,[ndim*nsub 1]))).*BigX;
     
     %z11_tXS21 = arrayfun(@(s) full(transpose(z11_BigX(1+(s-1)*ndim:s*ndim,:))*((sparse(diag(A_inv22(:,s)))*RR_w_inv)*sparse(diag(A_inv22(:,s))))*dof),[1:nsub],'UniformOutput', false);
     z11_tXS21 = arrayfun(@(s) tXS_NNGP(z11_BigX(1+(s-1)*ndim:s*ndim,:),A_inv22(:,s),RR_w_F_all,dof,sigma_w,ndim),[1:nsub],'UniformOutput', false);
     z11_tXS21=(cell2mat(z11_tXS21));
     
     VVV21=inv(z11_tXS21*z11_BigX+kron(diag(1./R_b_F_all)/sigma_b,speye(npred)));
     

     
     %MMM21= VVV21*z11_tXS21*t21(:);
     MMM21= mmtimes(VVV21,z11_tXS21,t21(:));
     
     beta21=MMM21 + sparse(chol(VVV21,'lower'))*normrnd(0,1,[ndim*npred 1]);
     beta21=reshape(beta21,[npred ndim])';
     
     
     % beta31
     %z11_tXS31 = arrayfun(@(s) full(transpose(z11_BigX(1+(s-1)*ndim:s*ndim,:))*((sparse(diag(A_inv33(:,s)))*RR_w_inv)*sparse(diag(A_inv33(:,s))))*dof),[1:nsub],'UniformOutput', false);
     z11_tXS31 = arrayfun(@(s) tXS_NNGP(z11_BigX(1+(s-1)*ndim:s*ndim,:),A_inv33(:,s),RR_w_F_all,dof,sigma_w,ndim),[1:nsub],'UniformOutput', false);
     z11_tXS31=cell2mat(z11_tXS31);
     
     VVV31=inv(z11_tXS31*z11_BigX+kron(diag(1./R_b_F_all)/sigma_b,speye(npred)));
     
     z21=reshape((t21-...
         sqrt(z11).*reshape(sum(X.*repmat(beta21,[nsub,1]),2),[ndim nsub]))...
         ./sqrt(reshape(scale22,[ndim nsub])),[ndim*nsub 1]);
     
%      MMM31= VVV31*z11_tXS31*(t31(:)-z21.*sum(X.*repmat(beta32,[nsub,1]),2));
      MMM31= mmtimes(VVV31,z11_tXS31,(t31(:)-z21.*sum(X.*repmat(beta32,[nsub,1]),2)));
     
     beta31=MMM31 + sparse(chol(VVV31,'lower'))*normrnd(0,1,[ndim*npred 1]);
     beta31=reshape(beta31,[npred ndim])';
     
     
     % beta32
     z21_22_BigX=repmat(sqrt(reshape(z22,[ndim*nsub 1]))+z21,[1 ndim*npred]).*BigX;
     %z21_22_BigX=kron(repmat(1,[1 ndim*npred]),sqrt(reshape(z22,[ndim*nsub 1]))+z21).*BigX;
     
     %z21_22_tXS32 = arrayfun(@(s) full(transpose(z21_22_BigX(1+(s-1)*ndim:s*ndim,:))*((sparse(diag(A_inv33(:,s)))*RR_w_inv)*sparse(diag(A_inv33(:,s))))*dof),[1:nsub],'UniformOutput', false);
     z21_22_tXS32 = arrayfun(@(s) 1/2*tXS_NNGP(z21_22_BigX(1+(s-1)*ndim:s*ndim,:),A_inv33(:,s),RR_w_F_all,dof,sigma_w,ndim),[1:nsub],'UniformOutput', false);
     z21_22_tXS32=cell2mat(z21_22_tXS32);
     
     VVV32=inv(z21_22_tXS32*z21_22_BigX+kron(diag(1./R_b_F_all)/sigma_b,speye(npred)));
     
%      MMM32= VVV21*z21_22_tXS32*(t32(:)+t31(:)-...
%          reshape(sqrt(z11).*reshape(sum(X.*repmat(beta31,[nsub,1]),2),[ndim nsub]),[ndim*nsub 1]));

      MMM32= mmtimes(VVV21,z21_22_tXS32,(t32(:)+t31(:)-...
          reshape(sqrt(z11).*reshape(sum(X.*repmat(beta31,[nsub,1]),2),[ndim nsub]),[ndim*nsub 1])));
     beta32=MMM32 + sparse(chol(VVV32,'lower'))*normrnd(0,1,[ndim*npred 1]);
     beta32=reshape(beta32,[npred ndim])';
     
          %% Finalizing
     t21_central=t21-sqrt(z11).*reshape(sum(X.*repmat(beta21,[nsub,1]),2),[ndim nsub]);
     t31_central=t31-sqrt(z11).*reshape(sum(X.*repmat(beta31,[nsub,1]),2),[ndim nsub])+reshape(z21,[ndim nsub]).*reshape(sum(X.*repmat(beta32,[nsub,1]),2),[ndim nsub]);
     t32_central=t32-sqrt(z22).*reshape(sum(X.*repmat(beta32,[nsub,1]),2),[ndim nsub]);
     
     Current_others_mapObj('t21_central')=t21_central;
     Current_others_mapObj('t31_central')=t31_central;
     Current_others_mapObj('t32_central')=t32_central;
     
     Current_beta_mapObj('beta21')=beta21;
     Current_beta_mapObj('beta31')=beta31;
     Current_beta_mapObj('beta32')=beta32;
     
     Bookkeeping_Current_mapObj('Current_others_mapObj')=Current_others_mapObj;
     Bookkeeping_Current_mapObj('Current_beta_mapObj')=Current_beta_mapObj;
     
     
     
     
 

end