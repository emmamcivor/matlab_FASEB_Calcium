
Lp=2500; %distance to PM (nm)
Le=2485; %distance to ER membrane (nm)
Lb=2000; %distance to interface between sub-PM ER and bulk ER (nm)
L0=0; %internal point of bulk ER, we think of this as the origin in z direction

w=100; %half width of cube for ER-PM junction and sub-PM ER (cube runs from -w to +w) (nm)
q=1000; %half width of cube for bulk ER. Chosen so STIM1 can reside anywhere on cube and still localise in ER-PM junction (nm)

Dj=220e6; % Diffusion in ER-PM junction - assume no buffers (nm^2/s)
De=10e6; %Diffusion in sub-PM ER and bulk ER - estimate from Sweitach 2008 and Dayel 1999 (nm^2/s)

dt=1e-6; % time step (s)

dx=5;
dy=5;
dz=(Le-Lb)/20;
x=-w:dx:w;
y=-w:dy:w;
z=Lb:dz:Le;

ms=1:2*w/dx;
ns=1:2*w/dy;
ps=1:(Le-Lb)/dz;

mu_s=ms*pi/(2*w);
eta_s=ns*pi/(2*w);
alpha_s=(2*ps-1)*pi/(2*(Le-Lb));

trapz_x= dx*[0.5 ones(1,length(x)-2) 0.5];
trapz_y= dy*[0.5 ones(1,length(y)-2) 0.5];
trapz_z= dz*[0.5 ones(1,length(z)-2) 0.5];
trapz_xy= kron(trapz_y,trapz_x);

tol=1e-8;

gauss_f_z=exp(-((z'-z(round(length(z)/2)))/2).^2);
gauss_f_xy=exp(-((((x'*ones(1,length(y))-x(round(length(x)/2)))/20).^2) + (((ones(length(x),1)*y-y(round(length(y)/2)))/20).^2)));

%% Test Convolution z direction
dt=0;
[gzs,gxys] = fn_gfs(Le,Lb,w,x,y,z,mu_s,eta_s,alpha_s,De,dt);

Iz=gzs*diag(trapz_z)*gauss_f_z;

assert(all((Iz-gauss_f_z)<tol))

%% convolution x y direction
dt=0;
[gzs,gxys] = fn_gfs(Le,Lb,w,x,y,z,mu_s,eta_s,alpha_s,De,dt);

Ixy=gxys*diag(trapz_xy)*gauss_f_xy(:);

assert(all((Ixy-gauss_f_xy(:))<tol))