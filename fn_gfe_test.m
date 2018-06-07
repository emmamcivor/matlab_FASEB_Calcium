Lp=2500; %distance to PM (nm)
Le=2485; %distance to ER membrane (nm)
Lb=2000; %distance to interface between sub-PM ER and bulk ER (nm)
L0=0; %internal point of bulk ER, we think of this as the origin in z direction

w=100; %half width of cube for ER-PM junction and sub-PM ER (cube runs from -w to +w) (nm)
q=1000; %half width of cube for bulk ER. Chosen so STIM1 can reside anywhere on cube and still localise in ER-PM junction (nm)

Dj=220e6; % Diffusion in ER-PM junction - assume no buffers (nm^2/s)
De=10e6; %Diffusion in sub-PM ER and bulk ER - estimate from Sweitach 2008 and Dayel 1999 (nm^2/s)

dt=1e-6; % time step (s)

dxe=200;
dye=200;
dze=200; %(Lp-Le)/100;
xe=-q:dxe:q;
ye=-q:dye:q;
ze=L0:dze:Lb;

m=1:2*q/dxe;
n=1:2*q/dye;
p=1:(Lb-L0)/dze;

mu_e=m*pi/(2*q);
eta_e=n*pi/(2*q);
alpha_e=(2*p-1)*pi*0.5/(Lb-L0);

trapz_x= dxe*[0.5 ones(1,length(xe)-2) 0.5];
trapz_y= dye*[0.5 ones(1,length(ye)-2) 0.5];
trapz_z= dze*[0.5 ones(1,length(ze)-2) 0.5];
trapz_xy= kron(trapz_y,trapz_x);

tol=1e-8;

gauss_f_z=exp(-((ze'-ze(round(length(ze)/2)))/2).^2);
gauss_f_xy=exp(-((((xe'*ones(1,length(ye))-xe(round(length(xe)/2)))/5000).^2) + (((ones(length(xe),1)*ye-ye(round(length(ye)/2)))/5000).^2)));

%% Test Convolution z direction
dt=0;
[gze,gxye] = fn_gfe(Lb,L0,q,xe,ye,ze,mu_e,eta_e,alpha_e,De,dt);

Iz=gze*diag(trapz_z)*gauss_f_z;

assert(all((Iz-gauss_f_z)<tol))

%% convolution x y direction
dt=0;
[gze,gxye] = fn_gfe(Lb,L0,q,xe,ye,ze,mu_e,eta_e,alpha_e,De,dt);

Ixy=gxye*diag(trapz_xy)*gauss_f_xy(:);

assert(all((Ixy-gauss_f_xy(:))<tol))