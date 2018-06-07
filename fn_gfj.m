function [gzj,gxyj] = fn_gfj(Lp,Le,w,xj,yj,zj,mu,eta,alpha,Dj,dt)

% This function calculates the Green's function for the ER-PM junction
% domain (dt>0)
% returns a matrix of size(|z| x |z|) for the Green's function in the z
% direction
% returns a matrix of size(|x||y| x |x||y|) for the Green's function in the
% x-y plane
% the Greens function at t=0 obeys the convolution property of delta
% functions in the z direction and x-y plane
% inputs: 
% Lp: a number representing the distance to the PM from an internal point
% of ER, L0
% Le: a number representing the distance to the ER membrane from an internal point
% of ER, L0
% w: half width of the ER-PM junction
% xj: vector of length |x| containing all points in x direction
% yj: vector of length |y| containing all points in y direction
% zj: vector of length |z| containing all points in z direction
% mu: eigenvalues of unit 1/nm used to calculate x component of Green's function 
% eta: eigenvalues of unit 1/nm used to calculate y component of Green's function 
% alpha: eigenvalues of unit 1/nm used to calculate z component of Green's function 
% Dj: a number representing the diffusion coefficient in ER-PM junction, units nm^2/s
% dt: number representing the time step (seconds)

gxj=(1/w)*sin((xj'-w)*mu)*diag(exp(-mu.^2*Dj*dt))*sin(mu'*(xj-w));
gyj=(1/w)*sin((yj'-w)*eta)*diag(exp(-eta.^2*Dj*dt))*sin(eta'*(yj-w));
gzj=(1/(Lp-Le))*(1+2*(cos((zj-Le)'*alpha)*diag(exp(-alpha.^2*Dj*dt))*cos(alpha'*(zj-Le))));


gxyj=kron(gyj,gxj);


end

