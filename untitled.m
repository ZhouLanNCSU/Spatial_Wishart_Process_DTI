for it = 1:1
    if mod(it,1)==0
    it
    end
    %%% Update
      [Bookkeeping_Current_mapObj]=MCMC_dof_ModifiedModel_NNGP(it,Bookkeeping_Input_mapObj,Bookkeeping_Data_mapObj, Bookkeeping_Priors_mapObj, Bookkeeping_Current_mapObj);
      %[Bookkeeping_Current_mapObj]=MCMC_R_b_ModifiedModel_NNGP(it,Bookkeeping_Input_mapObj,Bookkeeping_Data_mapObj, Bookkeeping_Priors_mapObj, Bookkeeping_Current_mapObj);
      [Bookkeeping_Current_mapObj]=MCMC_R_c_ModifiedModel_NNGP(it,Bookkeeping_Input_mapObj,Bookkeeping_Data_mapObj, Bookkeeping_Priors_mapObj, Bookkeeping_Current_mapObj);
      [Bookkeeping_Current_mapObj]=MCMC_R_w_ModifiedModel_NNGP(it,Bookkeeping_Input_mapObj,Bookkeeping_Data_mapObj, Bookkeeping_Priors_mapObj, Bookkeeping_Current_mapObj);
      [Bookkeeping_Current_mapObj]=MCMC_beta_off_ModifiedModel_NNGP(it,Bookkeeping_Input_mapObj,Bookkeeping_Data_mapObj, Bookkeeping_Priors_mapObj, Bookkeeping_Current_mapObj);
      [Bookkeeping_Current_mapObj]=MCMC_beta_diag_ModifiedModel_NNGP(it,Bookkeeping_Input_mapObj,Bookkeeping_Data_mapObj, Bookkeeping_Priors_mapObj, Bookkeeping_Current_mapObj);
%     
    Bookkeeping_Current_dof_mapObj_copy=copy(Bookkeeping_Current_mapObj('Current_dof_mapObj'));
    Bookkeeping_Current_theta_w_mapObj_copy=copy(Bookkeeping_Current_mapObj('Current_theta_w_mapObj'));
    Bookkeeping_Current_theta_c_mapObj_copy=copy(Bookkeeping_Current_mapObj('Current_theta_c_mapObj'));
    Bookkeeping_Current_beta_mapObj_copy=copy(Bookkeeping_Current_mapObj('Current_beta_mapObj'));
    Bookkeeping_MCMC(it)={{Bookkeeping_Current_dof_mapObj_copy Bookkeeping_Current_theta_w_mapObj_copy Bookkeeping_Current_theta_c_mapObj_copy Bookkeeping_Current_beta_mapObj_copy}};
%     kk=Bookkeeping_MCMC{it};
%     cc=kk{1};
%     cc('dof')
    %     obj=Bookkeeping_MCMC{it};
%     current=obj('Current_theta_w_mapObj');
%     current('rho_w')

    if mod(it,50)==0
    save('NNGP_Real_Data_Copula5000.mat')
    end
end