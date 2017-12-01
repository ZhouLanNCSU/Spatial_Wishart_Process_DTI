for set=1:20
    nnn=set
rng(1000+set)

[x,y,z] = ndgrid(1:5);
 grid = [x(:),y(:),z(:)];
signal_index=find(3<=grid(:,1) & 4>=grid(:,1) & 3<=grid(:,2) & 4>=grid(:,2) & 3<=grid(:,3) & 4>=grid(:,3))';
grid=pdist2(grid,grid);


%% Generate Standard Scale Spatial Wishart Process



rho_w=1,nu_w=0.5,sigma_w=1, nugget=0,
rho_b=0.0000000001,nu_b=0.5,sigma_b=0.01
rho_c=1,nu_c=0.5,

var_cov_w=correlation_matern(grid, nu_w, rho_w);


ndim=5^3, nsub=20,
npred=1,

mean_nu=-1,sd_nu=1,
mean_range=0,sd_range=1,
mean_r=0,sd_r=10,
a_var=1,
b_var=1


iters=5000
,burn=0

R_c=correlation_matern(grid, nu_c, rho_c);
dof=10;

[z11 z22 z33 z21 z31 z32] = Spatial_Wishart_Process_Chol_RND(nsub,dof,var_cov_w,eye(3),true);



%%Design Matrix
X=rand([ndim*nsub,npred]);
X=[repmat([1],[nsub/2*ndim 1]); repmat([0],[nsub/2*ndim 1])];
BigX = X_to_BigX(X,ndim, nsub);

mm1=repmat(0,[ndim 1]);
mm2=repmat(0,[ndim 1]);
mm1(signal_index)=1;
mm2(signal_index)=0;
%%Generating betas
R_b=correlation_matern(grid, nu_b, rho_b)*1;
 beta11=[ mm1];
 beta22=[ mm2];
 beta33=[ mm2];
 beta21=[ mm2];
 beta31=[ mm2];
 beta32=[ mm2];




 
 %%Constructing t
scale11=exp(sum(X.*repmat(beta11,[nsub,1]),2));
scale22=exp(sum(X.*repmat(beta22,[nsub,1]),2));
scale33=exp(sum(X.*repmat(beta33,[nsub,1]),2));
t11=z11.*reshape(scale11,[ndim nsub]);
t22=z22.*reshape(scale22,[ndim nsub]);
t33=z33.*reshape(scale33,[ndim nsub]);
t21=sqrt(z11).*reshape(sum(X.*repmat(beta21,[nsub,1]),2),[ndim nsub])+z21.*sqrt(reshape(scale22,[ndim nsub]));
t31=sqrt(z11).*reshape(sum(X.*repmat(beta31,[nsub,1]),2),[ndim nsub])+z21.*reshape(sum(X.*repmat(beta32,[nsub,1]),2),[ndim nsub])+z31.*sqrt(reshape(scale33,[ndim nsub]));
t32=sqrt(z22).*reshape(sum(X.*repmat(beta32,[nsub,1]),2),[ndim nsub])+z32.*sqrt(reshape(scale33,[ndim nsub]));


data=repmat(0,[nsub*ndim 1])

for ss=1:nsub*ndim
    
    t=[sqrt(t11(ss)) t21(ss) t31(ss); 0 sqrt(t22(ss)) t32(ss);  0  0 sqrt(t33(ss))];
    matrix=t*t';
    lambdas=eigs(matrix);
    %FA = fractional_anisotropy(lambdas);
    data(ss)=lambdas(3);
    
end

csvwrite(['\\wolftech.ad.ncsu.edu\cos\stat\Redirect\zlan\Desktop\Spatial-Wishart-Process-OPT-FINAL2\Simulation3\FA_DATA\' num2str(set) 'FA_Sim3.dat'],[data X])
end
