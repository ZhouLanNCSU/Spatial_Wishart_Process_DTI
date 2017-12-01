function FA = fractional_anisotropy(lambdas)

lambda1=lambdas(1);
lambda2=lambdas(2);
lambda3=lambdas(3);


lambda_mean=(lambda1+lambda2+lambda3)/3;

FA=sqrt(3/2)*sqrt(((lambda1-lambda_mean).^2+(lambda2-lambda_mean).^2+(lambda3-lambda_mean).^2)./(lambda1.^2+lambda2.^2+lambda3.^2));

end

