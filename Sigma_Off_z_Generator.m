function [A_inv]= Sigma_Off_z_Generator(X,beta_pp,nsub,ndim)


%z=t./reshape(exp(sum(X.*repmat(beta_kk,[nsub,1]),2)),[ndim nsub]);
A_inv=arrayfun(@(sub) sqrt(exp(-sum(X((sub-1)*ndim+1:sub*ndim,:).*beta_pp,2))), [1:nsub],'UniformOutput', false);


% for sub=1:nsub
%     %B=((sqrt(z(:,sub))*transpose(sqrt(z(:,sub)))).*(RR));
%     A=sqrt(exp(-sum(X((sub-1)*ndim+1:sub*ndim,:).*beta_pp,2)));
%     
%     %varcov_cell(sub)={(diag(sqrt(exp(sum(X((sub-1)*ndim+1:sub*ndim,:).*beta_pp,2)))))*((sqrt(z(:,sub))*transpose(sqrt(z(:,sub)))).*(Rho.*Rho))*(diag(sqrt(exp(sum(X((sub-1)*ndim+1:sub*ndim,:).*beta_pp,2)))))};
%     %varcov_cell(sub)={A*(RR*A)};
% 
% end

end

