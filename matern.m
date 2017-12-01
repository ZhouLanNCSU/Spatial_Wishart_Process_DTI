function [ V ] = matern(d,nu, rho, sigma2)

if d==0
   V=sigma2/2;
else
   V=sigma2*2^(1-nu)/gamma(nu)*((sqrt(2*nu)*d/rho).^(nu).*besselk(nu,sqrt(2*nu)*d/rho));
end

end

