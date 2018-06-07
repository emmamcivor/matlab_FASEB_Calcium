
function ss_orai_xyz = fn_ss_PM_j(eta,mu,xj,yj,zj,w,x_orai_i,y_orai_i,Lp,Le)
% boundary solution across plasma membrane (PM) - Orai channel
% returns a matrix of dimensions TODO (things tested)
%
% - eta: eigenvalues
% - mu: TODO

beta_xy=sqrt(kron(eta.^2',ones(length(mu),1)) + kron(ones(length(eta),1),mu'.^2));

bx=sin((xj'-w)*mu)*diag(sin(mu*(xj(x_orai_i)-w)))/w;
by=sin((yj'-w)*eta)*diag(sin(eta*(yj(y_orai_i)-w)))/w;

bz_num= exp(-beta_xy*(Lp-zj)) + exp(-beta_xy*(zj-2*Le+Lp));
bz_den=beta_xy.*(1 - exp(-2*beta_xy*(Lp-Le)));
bz=bz_num./repmat(bz_den,1,length(zj));


bxy_nm=kron(by,bx);
ss_orai_xy_z=bxy_nm*bz;

ss_orai_xyz=reshape(ss_orai_xy_z,length(xj),length(yj),length(zj));


end
