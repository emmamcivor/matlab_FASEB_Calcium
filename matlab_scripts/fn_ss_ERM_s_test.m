% parameters
Lp=2500; %distance to PM (nm)
Le=2485; %distance to ER membrane (nm)
Lb=2000;

w=100; %half width of cube for ER-PM junction and sub-PM ER (cube runs from -w to +w) (nm)
dxs=5;
dys=5;
dzs=5; 

xs=-w:dxs:w;
ys=-w:dys:w;
zs=Lb:dzs:Le;

ms=1:2*w/dxs;
ns=1:2*w/dys;
ps=1:(Le-Lb)/dzs; %number of eigenvalues 

mu_s=ms*pi/(2*w); %eigenvalues
eta_s=ns*pi/(2*w); %eigenvalues
alpha_s=ps*pi/(Lp-Le); %eigenvalues

x_serca_i=round(length(xs)/2+30/dxs);
y_serca_i=round(length(ys)/2+30/dys);

x_serca_i=round(length(xs)/2+30/dxs);
z_ip3_i=round(length(zs)/2);


%Volume of cytosol (Berlin, Bassani, Bers 1994)
vol_cyt_L=24.5e-12; %(L) volume of cytosol of ventricular myocyte in litres
vol_cyt_m=vol_cyt_L*1e-3;    %(m^3) volume of cytosol of ventricular myocyte in metres
vol_cyt_nm=vol_cyt_m*1e27;

vol_cyt_nm_L_ratio = vol_cyt_nm/vol_cyt_L;

vol_ER_PM_junction=2*w*2*w*(Lp-Le);

tol=1e-6;
%% Test matrix dimensions

ss_serca_xyz = fn_ss_ERM_s(eta_s,mu_s,xs,ys,zs,w,x_serca_i,y_serca_i,Le,Lb);
ss_serca_xyz=reshape(ss_serca_xyz,length(xs), length(ys), length(zs));

assert(length(xs)==size(ss_serca_xyz,1),'testing x direction')

assert(length(ys)==size(ss_serca_xyz,2),'testing y direction')

assert(length(zs)==size(ss_serca_xyz,3),'testing z direction')

%% Test peak coordinate
ss_serca_xyz = fn_ss_ERM_s(eta_s,mu_s,xs,ys,zs,w,x_serca_i,y_serca_i,Le,Lb);
ss_serca_xyz=reshape(ss_serca_xyz,length(xs), length(ys), length(zs));

a_peak=ss_serca_xyz(x_serca_i,y_serca_i,end); %find peak value
ss_serca_xyz(x_serca_i,y_serca_i,end)=0; %set this value to zero for test

assert(all(ss_serca_xyz(:)<a_peak), 'test function peaks where we prescribe')

