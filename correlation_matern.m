function Cov = correlation_matern ( d_matrix, nu, rho)

%[x y] = meshgrid(1:10, 1:10);
%d=[x(:) y(:)];
%d_matrix=distance(d',d');
%d_matrix=pdist2(d,d);
[n1,n2]=size(d_matrix);
Cov=(1*2^(1-nu)/gamma(nu))*(((sqrt(2*nu)*d_matrix./rho).^(nu)).*besselk(nu,sqrt(2*nu)*d_matrix./rho));
Cov(1:n1+1:end)=repmat(1, [n1 1]);


%Cov=reshape(arrayfun( @(i) matern(d_matrix(i),nu, rho, sigma2) ,[1:prod(size(d_matrix))]),size(d_matrix));

%Cov=tril(Cov)+transpose(tril(Cov));

end


  
