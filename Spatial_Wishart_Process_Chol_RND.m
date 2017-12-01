function [z11 z22 z33 z21 z31 z32] = Spatial_Wishart_Process_Chol_RND(N,dof,R_M,Scale_M,Scale_Status)


    z11=[];
    z22=[];
    z33=[];
    z21=[];
    z31=[];
    z32=[];
    

    [n_R n_R]=size(R_M);
    [n_S n_S]=size(Scale_M);
    tensor_R_Scale=kron(R_M,Scale_M);
    for i=1:N
        X=mvnrnd(repmat(0,[n_R*n_S,1]),tensor_R_Scale,dof);
        [z11_sub z22_sub z33_sub z21_sub z31_sub z32_sub]=Wishart_Generator_Chol(X,n_S,n_R,dof,Scale_Status);        
        z11=[z11 z11_sub];
        z22=[z22 z22_sub];
        z33=[z33 z33_sub];
        z21=[z21 z21_sub];
        z31=[z31 z31_sub];
        z32=[z32 z32_sub];
    end

end

