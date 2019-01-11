% parameters
Lp=2500; %distance to PM (nm)
Le=2485; %distance to ER membrane (nm)

w=100; %half width of cube for ER-PM junction and sub-PM ER (cube runs from -w to +w) (nm)
dxj=2;
dyj=2;
dzj=1.5; 

xj=-w:dxj:w;
yj=-w:dyj:w;
zj=Le:dzj:Lp;

mj=1:200; %number of eigenvalues 
nj=1:200; %number of eigenvalues 
pj=1:100; %number of eigenvalues 

mu_j=mj*pi/(2*w); %eigenvalues
eta_j=nj*pi/(2*w); %eigenvalues
alpha_j=pj*pi/(Lp-Le); %eigenvalues

x_serca_i=round(length(xj)/2+30/dxj);
y_serca_i=round(length(yj)/2+30/dyj);


%Volume of cytosol (Berlin, Bassani, Bers 1994)
vol_cyt_L=24.5e-12; %(L) volume of cytosol of ventricular myocyte in litres
vol_cyt_m=vol_cyt_L*1e-3;    %(m^3) volume of cytosol of ventricular myocyte in metres
vol_cyt_nm=vol_cyt_m*1e27;

vol_cyt_nm_L_ratio = vol_cyt_nm/vol_cyt_L;

vol_ER_PM_junction=2*w*2*w*(Lp-Le);

%ER-PM junction
cj0=0.1; %micro moles per litre cytosol
cj=cj0*ones(length(xj),length(yj),length(zj));
cj_DBC=cj0;

%% Test matrix dimensions

ss_serca_xyz = fn_ss_ERM_j(eta_j,mu_j,xj,yj,zj,w,x_serca_i,y_serca_i,Lp,Le); 
ss_serca_xyz=reshape(ss_serca_xyz,length(xj), length(yj), length(zj));

assert(length(xj)==size(ss_serca_xyz,1),'testing x direction')

assert(length(yj)==size(ss_serca_xyz,2),'testing y direction')

assert(length(zj)==size(ss_serca_xyz,3),'testing z direction')

%% Test peak coordinate
ss_serca_xyz = fn_ss_ERM_j(eta_j,mu_j,xj,yj,zj,w,x_serca_i,y_serca_i,Lp,Le); 
ss_serca_xyz=reshape(ss_serca_xyz,length(xj), length(yj), length(zj));

ss_serca_xyz=abs(ss_serca_xyz); %make profile positive for finding peak

a_peak=ss_serca_xyz(x_serca_i,x_serca_i,1); %find peak value
ss_serca_xyz(x_serca_i,y_serca_i,1)=0; %set this value to zero for test

assert(all(ss_serca_xyz(:)<a_peak), 'test function peaks where we prescribe')

