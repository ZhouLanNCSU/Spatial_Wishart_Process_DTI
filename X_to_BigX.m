function BigX = X_to_BigX(X,ndim, nsub)

    BigX=[];

    for index=1:nsub
        X_sub=X((1+(index-1)*ndim):index*ndim,:);
        BigX_sub=kron(speye(ndim),X_sub);
        BigX_sub=BigX_sub(1:ndim:ndim^2,:);
        BigX=[BigX ;BigX_sub];
    end



end

