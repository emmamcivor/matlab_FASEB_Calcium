%% boundary solution across ER membrane (ERM) - SERCA pump

function ss_serca_j_xyz = fn_ss_ERM_j(eta,mu,xj,yj,zj,w,x_serca_i,y_serca_i,Lp,Le)

beta_xy=sqrt(kron(eta.^2',ones(length(mu),1)) + kron(ones(length(eta),1),mu'.^2));

bx=sin((xj'-w)*mu)*diag(sin(mu*(xj(x_serca_i)-w)))/w;
by=sin((yj'-w)*eta)*diag(sin(eta*(yj(y_serca_i)-w)))/w;

bz_num= exp(-beta_xy*(2*Lp-zj-Le)) + exp(-beta_xy*(zj-Le));
bz_den=beta_xy.*(exp(-2*beta_xy*(Lp-Le))-1);
bz=bz_num./repmat(bz_den,1,length(zj));

bxy_nm=kron(by,bx);
ss_serca_j_xyz=bxy_nm*bz;

% ss_serca_j_xyz=reshape(bxy_z,length(xj),length(yj),length(zj));

end
