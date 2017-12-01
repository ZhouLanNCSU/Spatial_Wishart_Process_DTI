function [Decision quantile_lower quantile_upper] = MCMC_Hypothesis_Test(Bookkeeping_MCMC, index,burn_in,alpha,p,method)

    iter=length(Bookkeeping_MCMC);
    
    beta_all=cell(6,1);
    beta_label={'beta11' ,'beta22', 'beta33' ,'beta21', 'beta31', 'beta32'};    
    
    for i=1:6       
        mybeta=cell2mat(arrayfun(@(it) Bookkeeping_MCMC{it}{index}(beta_label{i}), [(1+burn_in):iter], 'UniformOutput',0));
        beta_all{i}=mybeta;        
    end
    
    

if isequal(method,'Bonferroni')
    m=6;
    %% Obtaining the quantiles
    quantile_lower=cell(6,1);
    quantile_upper=cell(6,1);


    for i=1:6      
        mybeta=beta_all{i};
        quantile_lower{i}=quantile(mybeta,alpha/m,2);
        quantile_upper{i}=quantile(mybeta,1-alpha/m,2);
    end
    
    
    %% Obtaining indicators and decisions
    Indicator=cell(1,6);
    
    for i=1:6       
        lower=quantile_lower{i}>0;
        upper=quantile_upper{i}<0;        
        Indicator{i}=lower+upper;       
    end
    
    
    Decision=sum(cell2mat(Indicator),2);
    
    Decision(Decision~=0)=1;
    
elseif isequal(method,'Sidak')
     m=6;
    %% Obtaining the quantiles
    quantile_lower=cell(6,1);
    quantile_upper=cell(6,1);


    for i=1:6      
        mybeta=beta_all{i};
        quantile_lower{i}=quantile(mybeta,1-(1-alpha)^m,2);
        quantile_upper{i}=quantile(mybeta,(1-alpha)^m,2);
    end
    
    
    %% Obtaining indicators and decisions
    Indicator=cell(1,6);
    
    for i=1:6       
        lower=quantile_lower{i}>0;
        upper=quantile_upper{i}<0;        
        Indicator{i}=lower+upper;       
    end
    
    
    Decision=sum(cell2mat(Indicator),2);
    
    Decision(Decision~=0)=1;
end
  
        
    
    
    
  





end

