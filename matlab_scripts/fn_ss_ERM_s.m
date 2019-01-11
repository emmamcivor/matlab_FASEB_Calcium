%% boundary solution across ER membrane (ERM) - SERCA pump

function ss_serca_xyz_s = fn_ss_ERM_s(eta_s,mu_s,xs,ys,zs,w,x_serca_i,y_serca_i,Le,Lb)

beta_xy=sqrt(kron(eta_s.^2',ones(length(mu_s),1)) + kron(ones(length(eta_s),1),mu_s'.^2));

% m=n=0 => p=0
ss_serca_xyz_00=zs-Lb;

% m>0, n=0
bx=sin((xs'-w)*mu_s)*diag(sin(mu_s*(xs(x_serca_i)-w)));
by=ones(length(ys),1);
bz_num= exp(-mu_s'*(Le-zs)) - exp(-mu_s'*(zs-2*Lb+Le));
bz_den=mu_s'.*(1 + exp(-2*mu_s'*(Le-Lb)));
bz=bz_num./repmat(bz_den,1,length(zs));
ss_serca_xyz_m0=kron(by,bx)*bz;

% m=0, n>0
by=sin((ys'-w)*eta_s)*diag(sin(eta_s*(ys(y_serca_i)-w)));
bx=ones(length(xs),1);
bz_num= exp(-eta_s'*(Le-zs)) - exp(-eta_s'*(zs-2*Lb+Le));
bz_den=eta_s'.*(1 + exp(-2*eta_s'*(Le-Lb)));
bz=bz_num./repmat(bz_den,1,length(zs));
ss_serca_xyz_0n=kron(by,bx)*bz;

% m>0, n>0
bx=sin((xs'-w)*mu_s)*diag(sin(mu_s*(xs(x_serca_i)-w)));
by=sin((ys'-w)*eta_s)*diag(sin(eta_s*(ys(y_serca_i)-w)));
bz_num= exp(-beta_xy*(Le-zs)) - exp(-beta_xy*(zs-2*Lb+Le));
bz_den=beta_xy.*(1 + exp(-2*beta_xy*(Le-Lb)));
bz=bz_num./repmat(bz_den,1,length(zs));
bxy_nm=kron(by,bx);
ss_serca_xyz_mn=bxy_nm*bz;

ss_serca_xyz_s=(repmat(ss_serca_xyz_00,length(xs)*length(ys),1) + 2*ss_serca_xyz_m0 + 2*ss_serca_xyz_0n + 4*ss_serca_xyz_mn)/(4*w^2);

end
