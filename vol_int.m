%%% create function to do volume integral for f=C or f=Vj (SS soln)
%%% int_{V} G*f dV'

function C_vol_int = vol_int(gxy,gz,x,y,z,dx,dy,dz,C0,ss)
% create matrices to carry out trapezium rule 
%(TODO: pre-compute in script)

trapz_x= dx*[0.5 ones(1,length(x)-2) 0.5];
trapz_y= dy*[0.5 ones(1,length(y)-2) 0.5];
trapz_z= dz*[0.5 ones(1,length(z)-2) 0.5];
trapz_xy= kron(trapz_y,trapz_x);

C_int = gz*diag(trapz_z)*(C0-ss)*diag(trapz_xy)*gxy +ss;

C_vol_int = reshape(C_int',length(x),length(y),length(z));

end
