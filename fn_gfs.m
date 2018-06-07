
function [gzs,gxys] = fn_gfs(Le,Lb,w,xs,ys,zs,mu_s,eta_s,alpha_s,De,dt)
% This function calculates the Green's function for the sub-PM ER
% domain (dt>0)
% returns a matrix of size(|z| x |z|) for the Green's function in the z
% direction
% returns a matrix of size(|x||y| x |x||y|) for the Green's function in the
% x-y plane
% the Greens function at t=0 obeys the convolution property of delta
% functions in the z direction and x-y plane
% inputs: 
% Lb: a number representing the distance to the interface between the sub-PM ER and bulk ER from an internal point
% of ER, L0
% Le: a number representing the distance to the ER membrane from an internal point
% of ER, L0
% w: half width of the sub-PM ER
% xs: vector of length |x| containing all points in x direction
% ys: vector of length |y| containing all points in y direction
% zs: vector of length |z| containing all points in z direction
% mu: eigenvalues of unit 1/nm used to calculate x component of Green's function 
% eta: eigenvalues of unit 1/nm used to calculate y component of Green's function 
% alpha: eigenvalues of unit 1/nm used to calculate z component of Green's function 
% De: a number representing the diffusion coefficient in sub-PM ER (De same for all ER compartments), units nm^2/s
% dt: number representing the time step (seconds)

gxs=(1/(2*w))*(1+2*cos((xs'-w)*mu_s)*diag(exp(-mu_s.^2*De*dt))*cos(mu_s'*(xs-w)));
gys=(1/(2*w))*(1+2*cos((ys'-w)*eta_s)*diag(exp(-eta_s.^2*De*dt))*cos(eta_s'*(ys-w)));
gzs=(2/(Le-Lb))*((sin((zs-Lb)'*alpha_s)*diag(exp(-alpha_s.^2*De*dt))*sin(alpha_s'*(zs-Lb))));

gxys=kron(gys,gxs);

end


