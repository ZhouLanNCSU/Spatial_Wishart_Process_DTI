function [Bookkeeping_Current_mapObj]=Log_Gamma_V1_MCMC_beta_off_ModifiedModel(it,Bookkeeping_Input_mapObj,Bookkeeping_Data_mapObj, Bookkeeping_Priors_mapObj, Bookkeeping_Current_mapObj)
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
     Current_theta_b_mapObj=Bookkeeping_Current_mapObj('Current_theta_b_mapObj');
     Current_theta_w_mapObj=Bookkeeping_Current_mapObj('Current_theta_w_mapObj');
     Current_theta_Xi_mapObj=Bookkeeping_Current_mapObj('Current_theta_Xi_mapObj');
     Current_beta_mapObj=Bookkeeping_Current_mapObj('Current_beta_mapObj');
     Current_scale_mapObj=Bookkeeping_Current_mapObj('Current_scale_mapObj');
     Current_R_Xi_mapObj=Bookkeeping_Current_mapObj('Current_R_Xi_mapObj');
     Current_RR_w_mapObj=Bookkeeping_Current_mapObj('Current_RR_w_mapObj');
     Current_R_b_mapObj=Bookkeeping_Current_mapObj('Current_R_b_mapObj');
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
     %theta_b
     sigma_b=Current_theta_b_mapObj('sigma_b');
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
     R_Xi_Choleksy_Lower=Current_R_Xi_mapObj('R_Xi_Choleksy_Lower');
     R_Xi_Choleksy_Lower_inv=Current_R_Xi_mapObj('R_Xi_Choleksy_Lower_inv');
     R_Xi_det=Current_R_Xi_mapObj('R_Xi_det');
     %RR_w
     RR_w_Choleksy_Lower=Current_RR_w_mapObj('RR_w_Choleksy_Lower');
     RR_w_Choleksy_Lower_inv=Current_RR_w_mapObj('RR_w_Choleksy_Lower_inv');
     RR_w_inv=Current_RR_w_mapObj('RR_w_inv');
     RR_w_det=Current_RR_w_mapObj('RR_w_det');
     %R_b
R_b_inv=Current_R_b_mapObj('R_b_inv');
     %others
     A_inv22=Current_others_mapObj('A_inv22');
     A_inv33=Current_others_mapObj('A_inv33');
     z11=Current_others_mapObj('z11');
     z22=Current_others_mapObj('z22');
     log_t11_central=Current_others_mapObj('log_t11_central');
     log_t22_central=Current_others_mapObj('log_t22_central');
     log_t33_central=Current_others_mapObj('log_t33_central');
     t21_central=Current_others_mapObj('t21_central');
     t31_central=Current_others_mapObj('t31_central');
     t32_central=Current_others_mapObj('t32_central');
     
     
     
     %% Un Pararell code
     % beta21
     z11_BigX=repmat(sqrt(reshape(z11,[ndim*nsub 1])),[1 ndim*npred]).*BigX;
     
     %z11_tXS21 = arrayfun(@(s) full(transpose(z11_BigX(1+(s-1)*ndim:s*ndim,:))*((sparse(diag(A_inv22(:,s)))*RR_w_inv)*sparse(diag(A_inv22(:,s))))*dof),[1:nsub],'UniformOutput', false);
     z11_tXS21 = arrayfun(@(s) tXS(z11_BigX(1+(s-1)*ndim:s*ndim,:),A_inv22(:,s),RR_w_inv,dof,sigma_w,ndim),[1:nsub],'UniformOutput', false);
     z11_tXS21=(cell2mat(z11_tXS21));
     
     VVV21=inv(z11_tXS21*z11_BigX+kron(R_b_inv/sigma_b,speye(npred)));
     
     %MMM21= VVV21*z11_tXS21*t21(:);
     MMM21= mmtimes(VVV21,z11_tXS21,t21(:));
     
     beta21=MMM21 + sparse(chol(VVV21,'lower'))*normrnd(0,1,[ndim*npred 1]);
     beta21=reshape(beta21,[npred ndim])';
     
     
     % beta31
     %z11_tXS31 = arrayfun(@(s) full(transpose(z11_BigX(1+(s-1)*ndim:s*ndim,:))*((sparse(diag(A_inv33(:,s)))*RR_w_inv)*sparse(diag(A_inv33(:,s))))*dof),[1:nsub],'UniformOutput', false);
     z11_tXS31 = arrayfun(@(s) tXS(z11_BigX(1+(s-1)*ndim:s*ndim,:),A_inv33(:,s),RR_w_inv,dof,sigma_w,ndim),[1:nsub],'UniformOutput', false);
     z11_tXS31=cell2mat(z11_tXS31);
     
     VVV31=inv(z11_tXS31*z11_BigX+kron(R_b_inv/sigma_b,speye(npred)));
     
     z21=reshape((t21-...
         sqrt(z11).*reshape(sum(X.*repmat(beta21,[nsub,1]),2),[ndim nsub]))...
         ./sqrt(reshape(scale22,[ndim nsub])),[ndim*nsub 1]);
     
%      MMM31= VVV31*z11_tXS31*(t31(:)-z21.*sum(X.*repmat(beta32,[nsub,1]),2));
      MMM31= mmtimes(VVV31,z11_tXS31,(t31(:)-z21.*sum(X.*repmat(beta32,[nsub,1]),2)));
     
     beta31=MMM31 + sparse(chol(VVV31,'lower'))*normrnd(0,1,[ndim*npred 1]);
     beta31=reshape(beta31,[npred ndim])';
     
     
     % beta32
     z21_22_BigX=repmat(sqrt(reshape(z22,[ndim*nsub 1]))+z21,[1 ndim*npred]).*BigX;
     %z21_22_tXS32 = arrayfun(@(s) full(transpose(z21_22_BigX(1+(s-1)*ndim:s*ndim,:))*((sparse(diag(A_inv33(:,s)))*RR_w_inv)*sparse(diag(A_inv33(:,s))))*dof),[1:nsub],'UniformOutput', false);
     z21_22_tXS32 = arrayfun(@(s) 1/2*tXS(z21_22_BigX(1+(s-1)*ndim:s*ndim,:),A_inv33(:,s),RR_w_inv,dof,sigma_w,ndim),[1:nsub],'UniformOutput', false);
     z21_22_tXS32=cell2mat(z21_22_tXS32);
     
     VVV32=inv(z21_22_tXS32*z21_22_BigX+kron(R_b_inv/sigma_b,speye(npred)));
     
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

