function [Bookkeeping_Current_mapObj]=Log_Gamma_V1_MCMC_beta_diag_ModifiedModel_NNGP(it,Bookkeeping_Input_mapObj,Bookkeeping_Data_mapObj, Bookkeeping_Priors_mapObj, Bookkeeping_Current_mapObj)

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
     sigma_b=sd_beta;
     
     
     %Current
     Current_dof_mapObj=Bookkeeping_Current_mapObj('Current_dof_mapObj');
     Current_theta_w_mapObj=Bookkeeping_Current_mapObj('Current_theta_w_mapObj');
     Current_theta_Xi_mapObj=Bookkeeping_Current_mapObj('Current_theta_Xi_mapObj');
     Current_beta_mapObj=Bookkeeping_Current_mapObj('Current_beta_mapObj');
     Current_scale_mapObj=Bookkeeping_Current_mapObj('Current_scale_mapObj');
     Current_R_Xi_mapObj=Bookkeeping_Current_mapObj('Current_R_Xi_mapObj');
     Current_RR_w_mapObj=Bookkeeping_Current_mapObj('Current_RR_w_mapObj');
     Current_others_mapObj=Bookkeeping_Current_mapObj('Current_others_mapObj');
     Current_R_b_mapObj=Bookkeeping_Current_mapObj('Current_R_b_mapObj');
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

     
          %R_b
     R_b_B_all=Current_R_b_mapObj('R_b_B_all');
     R_b_F_all=Current_R_b_mapObj('R_b_F_all');
     
     %RR_w
     RR_w_B_all=Current_RR_w_mapObj('RR_w_B_all');
     RR_w_F_all=Current_RR_w_mapObj('RR_w_F_all');
     
     t21_central_for_blk=Current_RR_w_mapObj('t21_central_for_blk');
     t31_central_for_blk=Current_RR_w_mapObj('t31_central_for_blk');
     t32_central_for_blk=Current_RR_w_mapObj('t32_central_for_blk');
     
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
     
     
     
     %% beta11
     count=0;
     for dim_index=1:ndim
        for pred_index=1:npred
          
          candidate=normrnd(beta11(dim_index,pred_index),0.1);
          beta11_can=beta11;
          beta11_can(dim_index,pred_index)=candidate;
          
          
          index=dim_index:ndim:(ndim*nsub);
          
          scale11_partial_can=exp(X(index,:)*beta11_can(dim_index,:)');
          scale11_partial=scale11(index);
          
          
          
          %Computing t11
          %aa=R_Xi_inv(:,dim_index)/sigma_Xi11;
          %aa(dim_index)=[];
          %tt=log_t11_central(:,pred_index);
          %tt(dim_index,:)=[];
          log_t11_central_partial_can=log_t11(index)-psi(0,(dof-0)/2)-log(2)+log(dof)-log(scale11_partial_can)';
          log_t11_central_partial=log_t11_central(index);
          
          can11=-0.5*sum(R_Xi_F_all(dim_index)/sigma_Xi11*(log_t11_central_partial_can).^2);
          old11=-0.5*sum(R_Xi_F_all(dim_index)/sigma_Xi11*(log_t11_central_partial).^2);
          
%             can11=0;
%             old11=0;

           
 
          
          %Computing t21
          z11_partial=z11(index);
          z11_partial_can=t11(index)./scale11_partial_can';
          t21_central_partial_can=t21(index)-sqrt(z11_partial_can).*(X(index,:)*beta21(dim_index,:)')';
          t21_central_partial=t21_central(index);
          %[A_inv_cross A_inv_partial] = off_cross_diag( A_inv22,dim_index ,index,ndim,npred);
          [A_inv_partial] = off_cross_diag_NNGP( A_inv22,dim_index ,index,ndim,npred);
          can21=-0.5*sum(1/RR_w_F_all(dim_index)*dof/sigma_w*(t21_central_partial_can.^2.*A_inv_partial.^2));
          old21=-0.5*sum(1/RR_w_F_all(dim_index)*dof/sigma_w*(t21_central_partial.^2.*A_inv_partial.^2));
  
          
          %Computing t31
          t31_central_partial=t31_central(index);
          t31_central_partial_can=t31_central_partial-sqrt(z11_partial_can).*(X(index,:)*beta31(dim_index,:)')'+sqrt(z11_partial).*(X(index,:)*beta31(dim_index,:)')';
  
          %[A_inv_cross A_inv_partial] = off_cross_diag( A_inv33,dim_index ,index,ndim,npred);
          [A_inv_partial] = off_cross_diag_NNGP( A_inv33,dim_index ,index,ndim,npred);
          can31=-0.5*sum(1/RR_w_F_all(dim_index)*dof/sigma_w*(t31_central_partial_can.^2.*A_inv_partial.^2));
          old31=-0.5*sum(1/RR_w_F_all(dim_index)*dof/sigma_w*(t31_central_partial.^2.*A_inv_partial.^2));

          
          %beta prior
 
          bcan11=-0.5*sum(1/R_b_F_all(dim_index)/sigma_b*(candidate)^2);

          bold11=-0.5*sum(1/R_b_F_all(dim_index)/sigma_b*(beta11(dim_index,pred_index))^2);
  
          
          
          
          % Compute prob
          can=can11+can21+can31+bcan11;
          old=old11+old21+old31+bold11;
          prob=min(1,exp(can-old));
          Accept=randsample([0,1],1,true,[1-prob,prob]);
          
          if Accept==1
              beta11=beta11_can;
              scale11(index)=scale11_partial_can;
              log_t11_central(index)=log_t11_central_partial_can;
              z11(index)=z11_partial_can;
              t21_central(index)=t21_central_partial_can;
              t31_central(index)=t31_central_partial_can;
              count=count+1;
          end
          

        end
     end

      if mod(it,20)==0    
 ['beta11 Accpetance Rate'] 
count/(ndim*npred)
      end
 
      
      
      
           %% beta22
     z21=(t21-sqrt(z11).*reshape(sum(X.*repmat(beta21,[nsub,1]),2),[ndim nsub]))./sqrt(reshape(scale22,[ndim nsub]));
     
     count=0;
     for dim_index=1:ndim
        for pred_index=1:npred
          candidate=normrnd(beta22(dim_index,pred_index),0.1);
          beta22_can=beta22;
          beta22_can(dim_index,pred_index)=candidate;
          
 

          index=dim_index:ndim:(ndim*nsub);
          
          scale22_partial_can=exp(X(index,:)*beta22_can(dim_index,:)');
          scale22_partial=scale22(index);
          
          
          %Computing t22
%           aa=R_Xi_inv(:,dim_index)/sigma_Xi22;
%           aa(dim_index)=[];
%           tt=log_t22_central(:,pred_index);
%           tt(dim_index,:)=[];
          log_t22_central_partial_can=log_t22(index)-psi(0,(dof-1)/2)-log(2)+log(dof)-log(scale22_partial_can)';
          log_t22_central_partial=log_t22_central(index);
          
          can22=-0.5*sum(R_Xi_F_all(dim_index)/sigma_Xi22*(log_t22_central_partial_can).^2);
          old22=-0.5*sum(R_Xi_F_all(dim_index)/sigma_Xi22*(log_t22_central_partial).^2);
           
%             can22=0;
%             old22=0;
%  
          
          %Computing t21
          t21_central_partial=t21_central(index);
   
          
          A_inv22_can=A_inv22;
          A_inv22_can(index)=1./sqrt(scale22_partial_can);
          
          %[A_inv_cross A_inv_partial] = off_cross_diag( A_inv22,dim_index ,index,ndim,npred);
          [A_inv_partial] = off_cross_diag_NNGP( A_inv22,dim_index ,index,ndim,npred);
          %[A_inv_cross_can A_inv_partial_can] = off_cross_diag( A_inv22_can,dim_index ,index,ndim,npred);
          [A_inv_partial_can] = off_cross_diag_NNGP( A_inv22_can,dim_index ,index,ndim,npred);
          
          can21=-0.5*sum(1/RR_w_F_all(dim_index)*dof/sigma_w*(t21_central_partial.^2.*A_inv_partial_can.^2));
          old21=-0.5*sum(1/RR_w_F_all(dim_index)*dof/sigma_w*(t21_central_partial.^2.*A_inv_partial.^2));
          
          
          %Computing t32
          z22_partial=z22(index);
          z22_partial_can=t22(index)./scale22_partial_can';
          t32_central_partial=t32_central(index);
          t32_central_partial_can=t32_central_partial-sqrt(z22_partial_can).*(X(index,:)*beta32(dim_index,:)')'+sqrt(z22_partial).*(X(index,:)*beta32(dim_index,:)')';
          
 
          %[A_inv_cross A_inv_partial] = off_cross_diag( A_inv33,dim_index ,index,ndim,npred);
          [A_inv_partial] = off_cross_diag_NNGP( A_inv33,dim_index ,index,ndim,npred);
          can32=-0.5*sum(1/RR_w_F_all(dim_index)*dof/sigma_w*(t32_central_partial_can.^2.*A_inv_partial.^2));
          old32=-0.5*sum(1/RR_w_F_all(dim_index)*dof/sigma_w*(t32_central_partial.^2.*A_inv_partial.^2));
          
          
          %Computing t31
          z21_partial=z21(index);
          z21_partial_can=z21_partial.*sqrt(scale22(index)')./sqrt(scale22_partial_can');
          t31_central_partial=t31_central(index);
          t31_central_partial_can=t31_central_partial-z21_partial_can.*(X(index,:)*beta32(dim_index,:)')'+z21_partial.*(X(index,:)*beta32(dim_index,:)')';

          %[A_inv_cross A_inv_partial] = off_cross_diag( A_inv33,dim_index ,index,ndim,npred);
          [A_inv_partial] = off_cross_diag_NNGP( A_inv33,dim_index ,index,ndim,npred);
          can31=-0.5*sum(1/RR_w_F_all(dim_index)*dof/sigma_w*(t31_central_partial_can.^2.*A_inv_partial.^2));
          old31=-0.5*sum(1/RR_w_F_all(dim_index)*dof/sigma_w*(t31_central_partial.^2.*A_inv_partial.^2));
          
          
          %beta prior
          bcan22=-0.5*sum(1/R_b_F_all(dim_index)/sigma_b*(candidate)^2);

          bold22=-0.5*sum(1/R_b_F_all(dim_index)/sigma_b*(beta11(dim_index,pred_index))^2);
  
           
          
          %bcan22 = 0; bcan22 = bcan22;
          
          

          % Compute prob
          can=can22+can21+can32+can31+bcan22;
          old=old22+old21+old32+old31+bold22;
          prob=min(1,exp(can-old));
          Accept=randsample([0,1],1,true,[1-prob,prob]);
          
          if Accept==1
              beta22=beta22_can;
              scale22(index)=scale22_partial_can;
              log_t22_central(index)=log_t22_central_partial_can;
              z22(index)=z22_partial_can;
              z21(index)=z21_partial_can;
%               t21_central(index)=t21_central_partial_can;
              t32_central(index)=t32_central_partial_can;
              A_inv22=A_inv22_can;
              count=count+1;
          end

        end
     end
if mod(it,20)==0
['beta22 Accpetance Rate'] 
count/(ndim*npred)
end      
     
     
     %% beta33
     count=0;
     for dim_index=1:ndim
        for pred_index=1:npred
          candidate=normrnd(beta33(dim_index,pred_index),0.1);
          beta33_can=beta33;
          beta33_can(dim_index,pred_index)=candidate;
 
          
          index=dim_index:ndim:(ndim*nsub);
          
          scale33_partial_can=exp(X(index,:)*beta33_can(dim_index,:)');
          scale33_partial=scale33(index);
          
          
          
          %Computing t33
%           aa=R_Xi_inv(:,dim_index)/sigma_Xi33;
%           aa(dim_index)=[];
%           tt=log_t33_central(:,pred_index);
%           tt(dim_index,:)=[];
          log_t33_central_partial_can=log_t33(index)-psi(0,(dof-2)/2)-log(2)+log(dof)-log(scale33_partial_can)';
          log_t33_central_partial=log_t33_central(index);
          
          can33=-0.5*sum(R_Xi_F_all(dim_index)/sigma_Xi33*(log_t33_central_partial_can).^2);
          old33=-0.5*sum(R_Xi_F_all(dim_index)/sigma_Xi33*(log_t33_central_partial).^2);
  
%              can33=0;
%              old33=0;
%  
          %Computing t32
          t32_central_partial=t32_central(index);
          A_inv33_can=A_inv33;
          A_inv33_can(index)=1./sqrt(scale33_partial_can);
          
          %[A_inv_cross A_inv_partial] = off_cross_diag( A_inv33,dim_index ,index,ndim,npred);
          [A_inv_partial] = off_cross_diag_NNGP( A_inv33,dim_index ,index,ndim,npred);
          %[A_inv_cross_can A_inv_partial_can] = off_cross_diag( A_inv33_can,dim_index ,index,ndim,npred);
          [A_inv_partial_can] = off_cross_diag_NNGP( A_inv33_can,dim_index ,index,ndim,npred);
          
          can32=-0.5*sum(1/RR_w_F_all(dim_index)*dof/sigma_w*(t32_central_partial.^2.*A_inv_partial_can.^2));
          old32=-0.5*sum(1/RR_w_F_all(dim_index)*dof/sigma_w*(t32_central_partial.^2.*A_inv_partial.^2));
          
          %Computing t31
          t31_central_partial=t31_central(index);
          
 
          
%           A_inv33_can=A_inv33;
%           A_inv33_can(index)=1./sqrt(scale33_partial_can);
%           
%           [A_inv_cross A_inv_partial] = off_cross_diag( A_inv33,dim_index ,index,ndim);
%           [A_inv_cross_can A_inv_partial_can] = off_cross_diag( A_inv33_can,dim_index ,index,ndim);
%           
          
          can31=-0.5*sum(1/RR_w_F_all(dim_index)*dof/sigma_w*(t31_central_partial.^2.*A_inv_partial_can.^2));
          old31=-0.5*sum(1/RR_w_F_all(dim_index)*dof/sigma_w*(t31_central_partial.^2.*A_inv_partial.^2));
          
          
          
          %beta prior
          bcan33=-0.5*sum(1/R_b_F_all(dim_index)/sigma_b*(candidate)^2);

          bold33=-0.5*sum(1/R_b_F_all(dim_index)/sigma_b*(beta11(dim_index,pred_index))^2);
           
          
          %bcan33 = 0; bcan33 = bcan33;
          
          
          % Compute prob
          can=can33+can31+can32+bcan33;
          old=old33+old31+old32+bold33;
          prob=min(1,exp(can-old));
          Accept=randsample([0,1],1,true,[1-prob,prob]);
          
          if Accept==1
              beta33=beta33_can;
              scale33(index)=scale33_partial_can;
              log_t33_central(index)=log_t33_central_partial_can;
   
              A_inv33=A_inv33_can;
              count=count+1;
          end
         
        end
     end
if mod(it,20)==0
['beta33 Accpetance Rate'] 
count/(ndim*npred)
end    

     %% Finalizing
     %beta
     Current_beta_mapObj('beta11')=beta11;
     Current_beta_mapObj('beta22')=beta22;
     Current_beta_mapObj('beta33')=beta33;     
     %scale
     Current_scale_mapObj('scale11')=scale11;
     Current_scale_mapObj('scale22')=scale22;
     Current_scale_mapObj('scale33')=scale33;
     %others
     Current_others_mapObj('A_inv22')=A_inv22;
     Current_others_mapObj('A_inv33')=A_inv33;
     Current_others_mapObj('z11')=z11;
     Current_others_mapObj('z22')=z22;
     Current_others_mapObj('t21_central')=t21_central;
     Current_others_mapObj('t31_central')=t31_central;
     Current_others_mapObj('t32_central')=t32_central;
     Current_others_mapObj('log_t11_central')=log_t11_central;
     Current_others_mapObj('log_t22_central')=log_t22_central;
     Current_others_mapObj('log_t33_central')=log_t33_central;

     
     %Current
     Bookkeeping_Current_mapObj('Current_beta_mapObj')=Current_beta_mapObj;
     Bookkeeping_Current_mapObj('Current_scale_mapObj')=Current_scale_mapObj;
     Bookkeeping_Current_mapObj('Current_others_mapObj')=Current_others_mapObj;





end