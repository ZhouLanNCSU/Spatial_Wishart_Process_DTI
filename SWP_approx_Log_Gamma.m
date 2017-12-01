
clear all
rng(1021)
 [x,y,z] = ndgrid(1:5);
 grid = [x(:),y(:),z(:)];
signal_index=find(3<=grid(:,1) & 4>=grid(:,1) & 3<=grid(:,2) & 4>=grid(:,2) & 3<=grid(:,3) & 4>=grid(:,3))';
grid=grid;


%% Generate Standard Scale Spatial Wishart Process

dof=10;

rho_w=1,nu_w=0.5,sigma_w=1, nugget=0,
rho_Xi=0.00000000000000001,nu_Xi= 1, sigma_Xi11=psi(1,(dof-0)/2), sigma_Xi22=psi(1,(dof-1)/2), sigma_Xi33=psi(1,(dof-2)/2), 
rho_b=0.00000000000000001,nu_b=0.5,sigma_b=0.01
sd_beta=0.01
var_cov_w=correlation_matern(grid, nu_w, rho_w);


ndim=5^3, nsub=20,
npred=1,

mean_nu=-1,sd_nu=1,
mean_range=0,sd_range=1,
mean_r=0,sd_r=10,
a_var=1,
b_var=1


iters=5000
,burn=0


[z11 z22 z33 z21 z31 z32] = Spatial_Wishart_Process_Chol_RND(nsub,dof,var_cov_w,sigma_w*eye(3),true);

%%Design Matrix
X=rand([ndim*nsub,npred]);
X=[repmat([1],[nsub/2*ndim 1]); repmat([0],[nsub/2*ndim 1])];
BigX = X_to_BigX(X,ndim, nsub);

mm1=repmat(0,[ndim 1]);
mm2=repmat(0,[ndim 1]);
mm1(signal_index)=0.50;
mm2(signal_index)=0.25;
%%Generating betas
R_b=correlation_matern(grid, nu_b, rho_b)*1;
 beta11=[ mm1];
 beta22=[ mm1];
 beta33=[ mm1];
 beta21=[ mm2];
 beta31=[ mm2];
 beta32=[ mm2];


R_Xi=correlation_matern(grid, nu_Xi, rho_Xi);
scale11=exp(sum(X.*repmat(beta11,[nsub,1]),2));
scale22=exp(sum(X.*repmat(beta22,[nsub,1]),2));
scale33=exp(sum(X.*repmat(beta33,[nsub,1]),2));
var_cov_w=correlation_matern(grid, nu_w, rho_w)*sigma_w;


 %%Constructing t
scale11=exp(sum(X.*repmat(beta11,[nsub,1]),2));
scale22=exp(sum(X.*repmat(beta22,[nsub,1]),2));
scale33=exp(sum(X.*repmat(beta33,[nsub,1]),2));
t11=z11.*reshape(scale11,[ndim nsub]);
t22=z22.*reshape(scale22,[ndim nsub]);
t33=z33.*reshape(scale33,[ndim nsub]);
log_t11=log(t11);
log_t22=log(t22);
log_t33=log(t33);
t21=sqrt(z11).*reshape(sum(X.*repmat(beta21,[nsub,1]),2),[ndim nsub])+z21.*sqrt(reshape(scale22,[ndim nsub]));
t31=sqrt(z11).*reshape(sum(X.*repmat(beta31,[nsub,1]),2),[ndim nsub])+z21.*reshape(sum(X.*repmat(beta32,[nsub,1]),2),[ndim nsub])+z31.*sqrt(reshape(scale33,[ndim nsub]));
t32=sqrt(z22).*reshape(sum(X.*repmat(beta32,[nsub,1]),2),[ndim nsub])+z32.*sqrt(reshape(scale33,[ndim nsub]));

t21_central=z21.*sqrt(reshape(scale22,[ndim nsub]));
t31_central=z31.*sqrt(reshape(scale33,[ndim nsub]));
t32_central=z32.*sqrt(reshape(scale33,[ndim nsub]));




%%% Bookkeeping

%%Input
Bookkeeping_Input_keySet={'grid','ndim','nsub','npred','X','BigX'};
Bookkeeping_Input_valueSet = {grid,ndim,nsub,npred,X,BigX};
Bookkeeping_Input_mapObj = containers.Map(Bookkeeping_Input_keySet,Bookkeeping_Input_valueSet);

%%Data
Bookkeeping_Data_keySet={'t11','t22','t33','log_t11','log_t22','log_t33','t21','t31','t32'};
Bookkeeping_Data_valueSet = {t11,t22,t33,log_t11,log_t22,log_t33,t21,t31,t32};
Bookkeeping_Data_mapObj = containers.Map(Bookkeeping_Data_keySet,Bookkeeping_Data_valueSet);

%%Priors
Bookkeeping_Priors_keySet={'mean_nu','sd_nu','mean_range','sd_range','mean_r','sd_r','a_var','b_var','sd_beta'};
Bookkeeping_Priors_valueSet = {mean_nu,sd_nu,mean_range,sd_range,mean_r,sd_r,a_var,b_var,sd_beta};
Bookkeeping_Priors_mapObj = containers.Map(Bookkeeping_Priors_keySet,Bookkeeping_Priors_valueSet);

%%Current Values
Bookkeeping_Current_keySet={'Current_dof_mapObj','Current_theta_b_mapObj','Current_theta_w_mapObj','Current_theta_Xi_mapObj',...
    'Current_R_b_mapObj','Current_beta_mapObj','Current_scale_mapObj','Current_R_Xi_mapObj','Current_RR_w_mapObj','Current_others_mapObj'};

%dof 
Current_dof_keySet={'dof'};
Current_dof_valueSet={dof};
Current_dof_mapObj=containers.Map(Current_dof_keySet,Current_dof_valueSet);

%theta_b
Current_theta_b_keySet={'rho_b','nu_b','sigma_b'};
Current_theta_b_valueSet={rho_b,nu_b,sigma_b};
Current_theta_b_mapObj=containers.Map(Current_theta_b_keySet,Current_theta_b_valueSet);

%theta_w
Current_theta_w_keySet={'rho_w','nu_w','sigma_w', 'nugget'};
Current_theta_w_valueSet={rho_w,nu_w,sigma_w, nugget};
Current_theta_w_mapObj=containers.Map(Current_theta_w_keySet,Current_theta_w_valueSet);

%theta_Xi
Current_theta_Xi_keySet={'rho_Xi','nu_Xi','sigma_Xi11','sigma_Xi22','sigma_Xi33'};
Current_theta_Xi_valueSet={rho_Xi,nu_Xi,sigma_Xi11,sigma_Xi22,sigma_Xi33};
Current_theta_Xi_mapObj=containers.Map(Current_theta_Xi_keySet,Current_theta_Xi_valueSet);

%beta
Current_beta_keySet={'beta11','beta22','beta33','beta21','beta31','beta32'};
Current_beta_valueSet={beta11,beta22,beta33,beta21,beta31,beta32};
Current_beta_mapObj=containers.Map(Current_beta_keySet,Current_beta_valueSet);

%scale
scale11=exp(sum(X.*repmat(beta11,[nsub,1]),2));
scale22=exp(sum(X.*repmat(beta22,[nsub,1]),2));
scale33=exp(sum(X.*repmat(beta33,[nsub,1]),2));
Current_scale_keySet={'scale11','scale22','scale33'};
Current_scale_valueSet={scale11,scale22,scale33};
Current_scale_mapObj=containers.Map(Current_scale_keySet,Current_scale_valueSet);

%R_b
R_b_Choleksy_Lower=sparse(chol(R_b,'lower'));
R_b_Choleksy_Lower_inv=inv(R_b_Choleksy_Lower);
R_b_det=prod(diag(R_b_Choleksy_Lower))^2;
R_b_inv=full(R_b_Choleksy_Lower_inv'*R_b_Choleksy_Lower_inv);
Current_R_b_keySet={'R_b_Choleksy_Lower','R_b_Choleksy_Lower_inv','R_b_det','R_b_inv'};
Current_R_b_valueSet={R_b_Choleksy_Lower,R_b_Choleksy_Lower_inv,R_b_det,R_b_inv};
Current_R_b_mapObj=containers.Map(Current_R_b_keySet,Current_R_b_valueSet);


%R_Xi
R_Xi_Choleksy_Lower=sparse(chol(R_Xi,'lower'));
R_Xi_Choleksy_Lower_inv=inv(R_Xi_Choleksy_Lower);
R_Xi_det=prod(diag(R_Xi_Choleksy_Lower))^2;
R_Xi_inv=full(R_Xi_Choleksy_Lower_inv'*R_Xi_Choleksy_Lower_inv);
Current_R_Xi_keySet={'R_Xi_Choleksy_Lower','R_Xi_Choleksy_Lower_inv','R_Xi_det','R_Xi_inv'};
Current_R_Xi_valueSet={R_Xi_Choleksy_Lower,R_Xi_Choleksy_Lower_inv,R_Xi_det,R_Xi_inv};
Current_R_Xi_mapObj=containers.Map(Current_R_Xi_keySet,Current_R_Xi_valueSet);


%RR_w
varcov_w=correlation_matern(grid, nu_w, rho_w);
RR_w_Choleksy_Lower=sparse(chol(varcov_w.*varcov_w,'lower'));
RR_w_Choleksy_Lower_inv=inv(RR_w_Choleksy_Lower);
RR_w_inv=full(RR_w_Choleksy_Lower_inv'*RR_w_Choleksy_Lower_inv);
RR_w_det=prod(diag(RR_w_Choleksy_Lower))^2;
Current_RR_w_keySet={'RR_w_Choleksy_Lower','RR_w_Choleksy_Lower_inv','RR_w_inv','RR_w_det'};
Current_RR_w_valueSet={RR_w_Choleksy_Lower,RR_w_Choleksy_Lower_inv,RR_w_inv,RR_w_det};
Current_RR_w_mapObj=containers.Map(Current_RR_w_keySet,Current_RR_w_valueSet);

%Others
A_inv22=cell2mat(Sigma_Off_z_Generator(X,beta22,nsub,ndim));
A_inv33=cell2mat(Sigma_Off_z_Generator(X,beta33,nsub,ndim));
z11=t11./reshape(scale11,[ndim,nsub]);
z22=t22./reshape(scale22,[ndim,nsub]);
log_t11_central=log_t11-psi(0,(dof-0)/2)-log(2)+log(dof)-log(reshape(scale11,[ndim nsub]));
log_t22_central=log_t22-psi(0,(dof-1)/2)-log(2)+log(dof)-log(reshape(scale22,[ndim nsub]));
log_t33_central=log_t33-psi(0,(dof-2)/2)-log(2)+log(dof)-log(reshape(scale33,[ndim nsub]));

Current_others_keySet={'A_inv22','A_inv33','z11','z22','log_t11_central','log_t22_central','log_t33_central','t21_central','t31_central','t32_central'};
Current_others_valueSet={A_inv22,A_inv33,z11,z22,log_t11_central,log_t22_central,log_t33_central,t21_central,t31_central,t32_central};
Current_others_mapObj=containers.Map(Current_others_keySet,Current_others_valueSet);


Bookkeeping_Current_valueSet={Current_dof_mapObj,Current_theta_b_mapObj,Current_theta_w_mapObj,Current_theta_Xi_mapObj,...
    Current_R_b_mapObj,Current_beta_mapObj,Current_scale_mapObj,Current_R_Xi_mapObj,Current_RR_w_mapObj,Current_others_mapObj};


Bookkeeping_Current_mapObj=containers.Map(Bookkeeping_Current_keySet,Bookkeeping_Current_valueSet);


for it = 1:iters
    if mod(it,20)==0
    it
    end
    %%% Update
      [Bookkeeping_Current_mapObj]=Log_Gamma_V1_MCMC_R_dof_ModifiedModel(it,Bookkeeping_Input_mapObj,Bookkeeping_Data_mapObj, Bookkeeping_Priors_mapObj, Bookkeeping_Current_mapObj);
      %[Bookkeeping_Current_mapObj]=Log_Gamma_V1_MCMC_R_Xi_ModifiedModel(it,Bookkeeping_Input_mapObj,Bookkeeping_Data_mapObj, Bookkeeping_Priors_mapObj, Bookkeeping_Current_mapObj);
      [Bookkeeping_Current_mapObj]=Log_Gamma_V1_MCMC_R_w_ModifiedModel(it,Bookkeeping_Input_mapObj,Bookkeeping_Data_mapObj, Bookkeeping_Priors_mapObj, Bookkeeping_Current_mapObj);
      [Bookkeeping_Current_mapObj]=Log_Gamma_V1_MCMC_beta_off_ModifiedModel(it,Bookkeeping_Input_mapObj,Bookkeeping_Data_mapObj, Bookkeeping_Priors_mapObj, Bookkeeping_Current_mapObj);
      [Bookkeeping_Current_mapObj]=Log_Gamma_V1_MCMC_beta_diag_ModifiedModel(it,Bookkeeping_Input_mapObj,Bookkeeping_Data_mapObj, Bookkeeping_Priors_mapObj, Bookkeeping_Current_mapObj);
% %     

     Bookkeeping_Current_dof_mapObj_copy=copy(Bookkeeping_Current_mapObj('Current_dof_mapObj'));
     Bookkeeping_Current_theta_w_mapObj_copy=copy(Bookkeeping_Current_mapObj('Current_theta_w_mapObj'));
     Bookkeeping_Current_theta_Xi_mapObj_copy=copy(Bookkeeping_Current_mapObj('Current_theta_Xi_mapObj'));
     Bookkeeping_Current_beta_mapObj_copy=copy(Bookkeeping_Current_mapObj('Current_beta_mapObj'));
     Bookkeeping_MCMC(it)={{Bookkeeping_Current_dof_mapObj_copy Bookkeeping_Current_theta_w_mapObj_copy Bookkeeping_Current_theta_Xi_mapObj_copy Bookkeeping_Current_beta_mapObj_copy}};
%     kk=Bookkeeping_MCMC{it};
%     cc=kk{1};
%     cc('dof')
    %     obj=Bookkeeping_MCMC{it};
%     current=obj('Current_theta_w_mapObj');
%     current('rho_w')
end


save('SWP_approx_Log_Gamma.mat')
