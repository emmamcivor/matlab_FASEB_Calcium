% parameters
Lp=2515; %distance to PM (nm)
Le=2500; %distance to ER membrane (nm)
Lb=2000;
Dj=220e6; % Diffusion in ER-PM junction - assume no buffers (nm^2/s)

w=100; %half width of cube for ER-PM junction and sub-PM ER (cube runs from -w to +w) (nm)
dx=5;
dy=5;
dzj=.15; 

x=-w:dx:w;
y=-w:dy:w;
zj=Le:dzj:Lp;

mj=1:2*w/dx; %number of eigenvalues 
nj=1:2*w/dx; %number of eigenvalues 
pj=1:(Lp-Le)/dzj; %number of eigenvalues 

mu_j=mj*pi/(2*w); %eigenvalues
eta_j=nj*pi/(2*w); %eigenvalues
alpha_j=pj*pi/(Lp-Le); %eigenvalues

dt=1e-6; %time step (seconds)

trapz_x= dx*[0.5 ones(1,length(x)-2) 0.5];
trapz_y= dy*[0.5 ones(1,length(y)-2) 0.5];
trapz_z= dzj*[0.5 ones(1,length(zj)-2) 0.5];
trapz_xy= kron(trapz_y,trapz_x);

tol=1e-8;

gauss_f_z=exp(-((zj'-zj(round(length(zj)/2)))/2).^2);
gauss_f_xy=exp(-((((x'*ones(1,length(y))-x(round(length(x)/2)))/10).^2) + (((ones(length(x),1)*y-y(round(length(y)/2)))/10).^2)));
%% Test Convolution z direction
dt=0;
[gzj,gxyj] = fn_gfj(Lp,Le,w,x,y,zj,mu_j,eta_j,alpha_j,Dj,dt);

Iz=gzj*diag(trapz_z)*gauss_f_z;

assert(all((Iz-gauss_f_z)<tol))

%% Test Convolution xy direction
dt=0;
[gzj,gxyj] = fn_gfj(Lp,Le,w,x,y,zj,mu_j,eta_j,alpha_j,Dj,dt);

Ixy=gxyj*diag(trapz_xy)*gauss_f_xy(:);

assert(all((Ixy-gauss_f_xy(:))<tol))