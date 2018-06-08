
function [gze,gxye] = fn_gfe(Lb,L0,q,xe,ye,ze,mu_e,eta_e,alpha_e,De,dt)

% This function calculates the Green's function for the bulk ER
% domain (dt>0)
% returns a matrix of size(|z| x |z|) for the Green's function in the z
% direction
% returns a matrix of size(|x||y| x |x||y|) for the Green's function in the
% x-y plane
% the Greens function at t=0 obeys the convolution property of delta
% functions in the z direction and x-y plane
% inputs: 
% Lb: a number representing the distance to the interface between the bulk ER and sub-PM ER from an internal point
% of ER, L0
% L0: a number representing an internal point of the ER. Acts as a
% reference point for distance calculations
% q: half width of the bulk cytoplasm "unit"
% xe: vector of length |x| containing all points in x direction
% ye: vector of length |y| containing all points in y direction
% ze: vector of length |z| containing all points in z direction
% mu: eigenvalues of unit 1/nm used to calculate x component of Green's function 
% eta: eigenvalues of unit 1/nm used to calculate y component of Green's function 
% alpha: eigenvalues of unit 1/nm used to calculate z component of Green's function 
% De: a number representing the diffusion coefficient in ER, units nm^2/s
% dt: number representing the time step (seconds)


gxe=(1/q)*sin((xe'-q)*mu_e)*diag(exp(-mu_e.^2*De*dt))*sin(mu_e'*(xe-q));
gye=(1/q)*sin((ye'-q)*eta_e)*diag(exp(-eta_e.^2*De*dt))*sin(eta_e'*(ye-q));
gze=(2/(Lb-L0))*sin((ze-L0)'*alpha_e)*diag(exp(-alpha_e.^2*De*dt))*sin(alpha_e'*(ze-L0));

gxye=kron(gye,gxe);

end

