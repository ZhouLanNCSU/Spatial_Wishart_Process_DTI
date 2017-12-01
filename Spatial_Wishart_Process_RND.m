function RND = Spatial_Wishart_Process_RND(N,dof,R_M,Scale_M,Scale_Status)

    RND=cell(N,1);
    [n_R n_R]=size(R_M);
    [n_S n_S]=size(Scale_M);
    tensor_R_Scale=kron(R_M,Scale_M);
    for i=1:N
        X=mvnrnd(repmat(0,[n_R*n_S,1]),tensor_R_Scale,dof);
        Data=Wishart_Generator(X,n_S,n_R,dof,Scale_Status);        
        RND(i)={Data};
    end

end

