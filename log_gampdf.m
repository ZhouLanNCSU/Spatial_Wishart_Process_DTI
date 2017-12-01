function Y = log_gampdf( X,A,B)

Y=-A.*log(B)-log(gamma(A))+(A-1).*log(X)-X./B;


end

