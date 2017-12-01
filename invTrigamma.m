function X = invTrigamma(Y)

myfun = @(Y,X) (psi(1,X)-Y);  
fun = @(X) myfun(Y,X);    
U = fzero(fun,10)

end
