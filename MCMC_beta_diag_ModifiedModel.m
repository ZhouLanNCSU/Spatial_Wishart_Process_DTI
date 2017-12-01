function [Bookkeeping_Current_mapObj]=MCMC_beta_diag_ModifiedModel(it,Bookkeeping_Input_mapObj,Bookkeeping_Data_mapObj, Bookkeeping_Priors_mapObj, Bookkeeping_Current_mapObj)
     
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
     Current_theta_b_mapObj=Bookkeeping_Current_mapObj('Current_theta_b_mapObj');
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
     %theta_b
     rho_b=Current_theta_b_mapObj('rho_b');
     nu_b=Current_theta_b_mapObj('nu_b');
     sigma_b=Current_theta_b_mapObj('sigma_b');
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
     %R_b
     R_b_Choleksy_Lower=Current_R_b_mapObj('R_b_Choleksy_Lower');
     R_b_Choleksy_Lower_inv=Current_R_b_mapObj('R_b_Choleksy_Lower_inv');
     R_b_det=Current_R_b_mapObj('R_b_det');
     R_b_inv=Current_R_b_mapObj('R_b_inv');
     %R_c
     R_c_Choleksy_Lower=Current_R_c_mapObj('R_c_Choleksy_Lower');
     R_c_Choleksy_Lower_inv=Current_R_c_mapObj('R_c_Choleksy_Lower_inv');
     R_c_det=Current_R_c_mapObj('R_c_det');
     R_c_inv=Current_R_c_mapObj('R_c_inv');
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
     
     %% beta11
     count=0;
     for dim_index=1:ndim
        for pred_index=1:npred
          candidate=normrnd(beta11(dim_index,pred_index),0.1);
          beta11_can=beta11;
          beta11_can(dim_index,pred_index)=candidate;
 
          %Computing t11
          index=dim_index:ndim:(ndim*nsub);
          
          scale11_partial_can=exp(X(index,:)*beta11_can(dim_index,:)');
          scale11_partial=scale11(index);
          
          q_t11_partial_can=gamcdf(t11(dim_index,:),repmat(dof,[1 nsub])/2,1/dof*2*scale11_partial_can');
          q_t11_partial=q_t11(index);
          
          u_t11_partial_can=norminv(q_t11_partial_can);
          u_t11_partial=u_t11(index);
         
          Jocobian_t11_partial_can=log_gampdf(t11(dim_index,:),repmat(dof,[1 nsub])/2,1/dof*2*scale11_partial_can');
          Jocobian_t11_partial=Jocobian_t11(index);
          
          aa=R_c_inv(:,dim_index);
          aa(dim_index)=[];
          uu=u_t11;
          uu(dim_index,:)=[];
          
          can11=-0.5*sum((R_c_inv(dim_index,dim_index)-1)*u_t11_partial_can.^2)-sum(repmat(aa,[1 nsub]).*uu*u_t11_partial_can')+sum(Jocobian_t11_partial_can);
          old11=-0.5*sum((R_c_inv(dim_index,dim_index)-1)*u_t11_partial.^2)-sum(repmat(aa,[1 nsub]).*uu*u_t11_partial')+sum(Jocobian_t11_partial);
          
          %Computing t21
          z11_partial=z11(index);
          z11_partial_can=t11(index)./scale11_partial_can';
          t21_central_partial_can=t21(index)-sqrt(z11_partial_can).*(X(index,:)*beta21(dim_index,:)')';
          t21_central_partial=t21_central(index);
          aa=RR_w_inv(:,dim_index)*dof/sigma_w;
          aa(dim_index)=[];
          tt=t21_central;
          tt(dim_index,:)=[];
          [A_inv_cross A_inv_partial] = off_cross_diag( A_inv22,dim_index ,index,ndim,npred);
          can21=-0.5*sum(RR_w_inv(dim_index,dim_index)*dof/sigma_w*(t21_central_partial_can.^2.*A_inv_partial.^2))...
              -sum(A_inv_cross.*repmat(aa,[1 nsub]).*tt*t21_central_partial_can');
          old21=-0.5*sum(RR_w_inv(dim_index,dim_index)*dof/sigma_w*(t21_central_partial.^2.*A_inv_partial.^2))...
              -sum(A_inv_cross.*repmat(aa,[1 nsub]).*tt*t21_central_partial');
          
          %Computing t31
          t31_central_partial=t31_central(index);
          t31_central_partial_can=t31_central_partial-sqrt(z11_partial_can).*(X(index,:)*beta31(dim_index,:)')'+sqrt(z11_partial).*(X(index,:)*beta31(dim_index,:)')';
          
          aa=RR_w_inv(:,dim_index)*dof/sigma_w;
          aa(dim_index)=[];
          tt=t31_central;
          tt(dim_index,:)=[];
          [A_inv_cross A_inv_partial] = off_cross_diag( A_inv33,dim_index ,index,ndim,npred);
          can31=-0.5*sum(RR_w_inv(dim_index,dim_index)*dof/sigma_w*(t31_central_partial_can.^2.*A_inv_partial.^2))...
              -sum(A_inv_cross.*repmat(aa,[1 nsub]).*tt*t31_central_partial_can');
          old31=-0.5*sum(RR_w_inv(dim_index,dim_index)*dof/sigma_w*(t31_central_partial.^2.*A_inv_partial.^2))...
              -sum(A_inv_cross.*repmat(aa,[1 nsub]).*tt*t31_central_partial');
           
          %beta prior
          aa=R_b_inv(:,dim_index)/sigma_b;
          aa(dim_index)=[];
          tt=beta11(:,pred_index);
          tt(dim_index,:)=[];
          bcan11=-0.5*sum(R_b_inv(dim_index,dim_index)/sigma_b*(candidate)^2)...
              -sum(aa.*tt*candidate);
          bold11=-0.5*sum(R_b_inv(dim_index,dim_index)/sigma_b*(beta11(dim_index,pred_index))^2)...
              -sum(aa.*tt*beta11(dim_index,pred_index));
           
          
          %bcan11 = 0; bcan11 = bcan11;
          
          
          
          % Compute prob
          can=can11+can21+can31+bcan11;
          old=old11+old21+old31+bold11;
          prob=min(1,exp(can-old));
          Accept=randsample([0,1],1,true,[1-prob,prob]);
          
          if Accept==1
              beta11=beta11_can;
              scale11(index)=scale11_partial_can;
              q_t11(index)=q_t11_partial_can;
              u_t11(index)=u_t11_partial_can;
              Jocobian_t11(index)=Jocobian_t11_partial_can;
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
          
 
          %Computing t22
          index=dim_index:ndim:(ndim*nsub);
          
          scale22_partial_can=exp(X(index,:)*beta22_can(dim_index,:)');
          scale22_partial=scale22(index);
          
          q_t22_partial_can=gamcdf(t22(dim_index,:),repmat(dof-1,[1 nsub])/2,1/(dof)*2*scale22_partial_can');
          q_t22_partial=q_t22(index);
          
          u_t22_partial_can=norminv(q_t22_partial_can);
          u_t22_partial=u_t22(index);
         
          Jocobian_t22_partial_can=log_gampdf(t22(dim_index,:),repmat(dof-1,[1 nsub])/2,1/(dof)*2*scale22_partial_can');
          Jocobian_t22_partial=Jocobian_t22(index);
          
          aa=R_c_inv(:,dim_index);
          aa(dim_index)=[];
          uu=u_t22;
          uu(dim_index,:)=[];
          
          can22=-0.5*sum((R_c_inv(dim_index,dim_index)-1)*u_t22_partial_can.^2)-sum(repmat(aa,[1 nsub]).*uu*u_t22_partial_can')+sum(Jocobian_t22_partial_can);
          old22=-0.5*sum((R_c_inv(dim_index,dim_index)-1)*u_t22_partial.^2)-sum(repmat(aa,[1 nsub]).*uu*u_t22_partial')+sum(Jocobian_t22_partial);
          
          
          %Computing t21
          t21_central_partial=t21_central(index);
          aa=RR_w_inv(:,dim_index)*dof/sigma_w;
          aa(dim_index)=[];
          tt=t21_central;
          tt(dim_index,:)=[];
          
          
          A_inv22_can=A_inv22;
          A_inv22_can(index)=1./sqrt(scale22_partial_can);
          
          [A_inv_cross A_inv_partial] = off_cross_diag( A_inv22,dim_index ,index,ndim,npred);
          [A_inv_cross_can A_inv_partial_can] = off_cross_diag( A_inv22_can,dim_index ,index,ndim,npred);
          
          can21=-0.5*sum(RR_w_inv(dim_index,dim_index)*dof/sigma_w*(t21_central_partial.^2.*A_inv_partial_can.^2))...
              -sum(A_inv_cross_can.*repmat(aa,[1 nsub]).*tt*t21_central_partial');
          old21=-0.5*sum(RR_w_inv(dim_index,dim_index)*dof/sigma_w*(t21_central_partial.^2.*A_inv_partial.^2))...
              -sum(A_inv_cross.*repmat(aa,[1 nsub]).*tt*t21_central_partial');
          
          
          %Computing t32
          z22_partial=z22(index);
          z22_partial_can=t22(index)./scale22_partial_can';
          t32_central_partial=t32_central(index);
          t32_central_partial_can=t32_central_partial-sqrt(z22_partial_can).*(X(index,:)*beta32(dim_index,:)')'+sqrt(z22_partial).*(X(index,:)*beta32(dim_index,:)')';
          
          aa=RR_w_inv(:,dim_index)*dof/sigma_w;
          aa(dim_index)=[];
          tt=t32_central;
          tt(dim_index,:)=[];
          [A_inv_cross A_inv_partial] = off_cross_diag( A_inv33,dim_index ,index,ndim,npred);
          can32=-0.5*sum(RR_w_inv(dim_index,dim_index)*dof/sigma_w*(t32_central_partial_can.^2.*A_inv_partial.^2))...
              -sum(A_inv_cross.*repmat(aa,[1 nsub]).*tt*t32_central_partial_can');
          old32=-0.5*sum(RR_w_inv(dim_index,dim_index)*dof/sigma_w*(t32_central_partial.^2.*A_inv_partial.^2))...
              -sum(A_inv_cross.*repmat(aa,[1 nsub]).*tt*t32_central_partial');
          
          
          %Computing t31
          z21_partial=z21(index);
          z21_partial_can=z21_partial.*sqrt(scale22(index)')./sqrt(scale22_partial_can');
          t31_central_partial=t31_central(index);
          t31_central_partial_can=t31_central_partial-z21_partial_can.*(X(index,:)*beta32(dim_index,:)')'+z21_partial.*(X(index,:)*beta32(dim_index,:)')';
          
          aa=RR_w_inv(:,dim_index)*dof/sigma_w;
          aa(dim_index)=[];
          tt=t31_central;
          tt(dim_index,:)=[];
          [A_inv_cross A_inv_partial] = off_cross_diag( A_inv33,dim_index ,index,ndim,npred);
          can31=-0.5*sum(RR_w_inv(dim_index,dim_index)*dof/sigma_w*(t31_central_partial_can.^2.*A_inv_partial.^2))...
              -sum(A_inv_cross.*repmat(aa,[1 nsub]).*tt*t31_central_partial_can');
          old31=-0.5*sum(RR_w_inv(dim_index,dim_index)*dof/sigma_w*(t31_central_partial.^2.*A_inv_partial.^2))...
              -sum(A_inv_cross.*repmat(aa,[1 nsub]).*tt*t31_central_partial');
          
          
          %beta prior
          aa=R_b_inv(:,dim_index)/sigma_b;
          aa(dim_index)=[];
          tt=beta22(:,pred_index);
          tt(dim_index,:)=[];
          bcan22=-0.5*sum(R_b_inv(dim_index,dim_index)/sigma_b*(candidate)^2)...
              -sum(aa.*tt*candidate);
          bold22=-0.5*sum(R_b_inv(dim_index,dim_index)/sigma_b*(beta22(dim_index,pred_index))^2)...
              -sum(aa.*tt*beta22(dim_index,pred_index));
           
          
          %bcan22 = 0; bcan22 = bcan22;
          
          

          % Compute prob
          can=can22+can21+can32+can31+bcan22;
          old=old22+old21+old32+old31+bold22;
          prob=min(1,exp(can-old));
          Accept=randsample([0,1],1,true,[1-prob,prob]);
          
          if Accept==1
              beta22=beta22_can;
              scale22(index)=scale22_partial_can;
              q_t22(index)=q_t22_partial_can;
              u_t22(index)=u_t22_partial_can;
              Jocobian_t22(index)=Jocobian_t22_partial_can;
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
 
          %Computing t33
          index=dim_index:ndim:(ndim*nsub);
          
          scale33_partial_can=exp(X(index,:)*beta33_can(dim_index,:)');
          scale33_partial=scale33(index);
          
          q_t33_partial_can=gamcdf(t33(dim_index,:),repmat(dof-2,[1 nsub])/2,1/(dof)*2*scale33_partial_can');
          q_t33_partial=q_t33(index);
          
          u_t33_partial_can=norminv(q_t33_partial_can);
          u_t33_partial=u_t33(index);
         
          Jocobian_t33_partial_can=log_gampdf(t33(dim_index,:),repmat(dof-2,[1 nsub])/2,1/(dof)*2*scale33_partial_can');
          Jocobian_t33_partial=Jocobian_t33(index);
          
          aa=R_c_inv(:,dim_index);
          aa(dim_index)=[];
          uu=u_t33;
          uu(dim_index,:)=[];
          
          can33=-0.5*sum((R_c_inv(dim_index,dim_index)-1)*u_t33_partial_can.^2)-sum(repmat(aa,[1 nsub]).*uu*u_t33_partial_can')+sum(Jocobian_t33_partial_can);
          old33=-0.5*sum((R_c_inv(dim_index,dim_index)-1)*u_t33_partial.^2)-sum(repmat(aa,[1 nsub]).*uu*u_t33_partial')+sum(Jocobian_t33_partial);
          
          %Computing t32
          t32_central_partial=t32_central(index);
          
          aa=RR_w_inv(:,dim_index)*dof/sigma_w;
          aa(dim_index)=[];
          tt=t32_central;
          tt(dim_index,:)=[];
          
          A_inv33_can=A_inv33;
          A_inv33_can(index)=1./sqrt(scale33_partial_can);
          
          [A_inv_cross A_inv_partial] = off_cross_diag( A_inv33,dim_index ,index,ndim,npred);
          [A_inv_cross_can A_inv_partial_can] = off_cross_diag( A_inv33_can,dim_index ,index,ndim,npred);
          
          
          can32=-0.5*sum(RR_w_inv(dim_index,dim_index)*dof/sigma_w*(t32_central_partial.^2.*A_inv_partial_can.^2))...
              -sum(A_inv_cross_can.*repmat(aa,[1 nsub]).*tt*t32_central_partial');
          old32=-0.5*sum(RR_w_inv(dim_index,dim_index)*dof/sigma_w*(t32_central_partial.^2.*A_inv_partial.^2))...
              -sum(A_inv_cross.*repmat(aa,[1 nsub]).*tt*t32_central_partial');
          
          %Computing t31
          t31_central_partial=t31_central(index);
          
          aa=RR_w_inv(:,dim_index)*dof/sigma_w;
          aa(dim_index)=[];
          tt=t31_central;
          tt(dim_index,:)=[];
          
%           A_inv33_can=A_inv33;
%           A_inv33_can(index)=1./sqrt(scale33_partial_can);
%           
%           [A_inv_cross A_inv_partial] = off_cross_diag( A_inv33,dim_index ,index,ndim);
%           [A_inv_cross_can A_inv_partial_can] = off_cross_diag( A_inv33_can,dim_index ,index,ndim);
%           
          
          can31=-0.5*sum(RR_w_inv(dim_index,dim_index)*dof/sigma_w*(t31_central_partial.^2.*A_inv_partial_can.^2))...
              -sum(A_inv_cross_can.*repmat(aa,[1 nsub]).*tt*t31_central_partial');
          old31=-0.5*sum(RR_w_inv(dim_index,dim_index)*dof/sigma_w*(t31_central_partial.^2.*A_inv_partial.^2))...
              -sum(A_inv_cross.*repmat(aa,[1 nsub]).*tt*t31_central_partial');
          
          
          
          %beta prior
          aa=R_b_inv(:,dim_index)/sigma_b;
          aa(dim_index)=[];
          tt=beta33(:,pred_index);
          tt(dim_index,:)=[];
          bcan33=-0.5*sum(R_b_inv(dim_index,dim_index)/sigma_b*(candidate)^2)...
              -sum(aa.*tt*candidate);
          bold33=-0.5*sum(R_b_inv(dim_index,dim_index)/sigma_b*(beta33(dim_index,pred_index))^2)...
              -sum(aa.*tt*beta33(dim_index,pred_index));
           
          
          %bcan33 = 0; bcan33 = bcan33;
          
          
          % Compute prob
          can=can33+can31+can32+bcan33;
          old=old33+old31+old32+bold33;
          prob=min(1,exp(can-old));
          Accept=randsample([0,1],1,true,[1-prob,prob]);
          
          if Accept==1
              beta33=beta33_can;
              scale33(index)=scale33_partial_can;
              q_t33(index)=q_t33_partial_can;
              u_t33(index)=u_t33_partial_can;
              Jocobian_t33(index)=Jocobian_t33_partial_can;
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
     %quantile
     Current_quantile_mapObj('q_t11')=q_t11;
     Current_quantile_mapObj('q_t22')=q_t22;
     Current_quantile_mapObj('q_t33')=q_t33;
     %norminv
     Current_norminv_mapObj('u_t11')=u_t11;
     Current_norminv_mapObj('u_t22')=u_t22;
     Current_norminv_mapObj('u_t33')=u_t33;
     %Jocobian 
     Current_Jocobian_mapObj('Jocobian_t11')=Jocobian_t11;
     Current_Jocobian_mapObj('Jocobian_t22')=Jocobian_t22;
     Current_Jocobian_mapObj('Jocobian_t33')=Jocobian_t33;
     %others
     Current_others_mapObj('A_inv22')=A_inv22;
     Current_others_mapObj('A_inv33')=A_inv33;
     Current_others_mapObj('z11')=z11;
     Current_others_mapObj('z22')=z22;
     Current_others_mapObj('t21_central')=t21_central;
     Current_others_mapObj('t31_central')=t31_central;
     Current_others_mapObj('t32_central')=t32_central;

     
     %Current
     Bookkeeping_Current_mapObj('Current_beta_mapObj')=Current_beta_mapObj;
     Bookkeeping_Current_mapObj('Current_scale_mapObj')=Current_scale_mapObj;
     Bookkeeping_Current_mapObj('Current_quantile_mapObj')=Current_quantile_mapObj;
     Bookkeeping_Current_mapObj('Current_norminv_mapObj')=Current_norminv_mapObj;
     Bookkeeping_Current_mapObj('Current_Jocobian_mapObj')=Current_Jocobian_mapObj;
     Bookkeeping_Current_mapObj('Current_others_mapObj')=Current_others_mapObj;
     

  
     
     

end

