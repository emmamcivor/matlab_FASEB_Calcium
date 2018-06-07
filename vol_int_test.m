Lp=2500; %distance to PM (nm)
Le=2485; %distance to ER membrane (nm)
Lb=2000; %distance to interface between sub-PM ER and bulk ER (nm)
L0=0; %internal point of bulk ER, we think of this as the origin in z direction

w=100; %half width of cube for ER-PM junction and sub-PM ER (cube runs from -w to +w) (nm)
q=1000; %half width of cube for bulk ER. Chosen so STIM1 can reside anywhere on cube and still localise in ER-PM junction (nm)

Dj=220e6; % Diffusion in ER-PM junction - assume no buffers (nm^2/s)
De=10e6; %Diffusion in sub-PM ER and bulk ER - estimate from Sweitach 2008 and Dayel 1999 (nm^2/s)

dt=1e-4; % time step (s)

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

[gze,gxye] = fn_gfe(Lb,L0,q,xe,ye,ze,mu_e,eta_e,alpha_e,De,dt);
c0=150;
ss=c0*ones(length(xe)*length(ye),length(ze))';  %assuming DBC equal to baseline IC so BCs are satisfied at t=0

%% Test constant IC 
ce0=c0*ones(length(xe)*length(ye),length(ze))';

ce_vol_int = vol_int(gxye,gze,xe,ye,ze,dxe,dye,dze,ce0,ss);    
assert(all(ce_vol_int(:)-ce0(:)==0))

%% Test DBCs satisfied
ce0=c0*ones(length(xe)*length(ye),length(ze))';

ce_vol_int = vol_int(gxye,gze,xe,ye,ze,dxe,dye,dze,ce0,ss);    
ss_xyz=c0*ones(length(xe),length(ye),length(ze));  %assuming DBC equal to baseline IC so BCs are satisfied at t=0

DBC1=ce_vol_int(1,:,:)-ss_xyz(1,:,:);
DBC2=ce_vol_int(end,:,:)-ss_xyz(end,:,:);
DBC3=ce_vol_int(:,1,:)-ss_xyz(:,1,:);
DBC4=ce_vol_int(:,end,:)-ss_xyz(:,end,:);

DBC=[DBC1(:);DBC2(:);DBC3(:);DBC4(:)];

assert(all(DBC(:)==0))

%% Test peak is where we prescribe
ce0_xyz=c0*ones(length(xe),length(ye),length(ze));
ce0_xyz(round(length(xe)/2),round(length(ye)/2),round(length(ze)/2))=c0*2;

ce0=reshape(ce0_xyz,length(xe)*length(ye),length(ze))';

ce_vol_int = vol_int(gxye,gze,xe,ye,ze,dxe,dye,dze,ce0,ss);    

a_peak=ce_vol_int(round(length(xe)/2),round(length(ye)/2),round(length(ze)/2)); %find peak value
ce_vol_int(round(length(xe)/2),round(length(ye)/2),round(length(ze)/2))=0; %set this value to zero for test

assert(all(ce_vol_int(:)<a_peak), 'test function peaks where we prescribe')


%% Test calcium diffuses (peak decreases)
ce0_xyz=c0*ones(length(xe),length(ye),length(ze));
ce0_xyz(round(length(xe)/2),round(length(ye)/2),round(length(ze)/2))=c0*2;

ce0=reshape(ce0_xyz,length(xe)*length(ye),length(ze))';

ce_vol_int = vol_int(gxye,gze,xe,ye,ze,dxe,dye,dze,ce0,ss);    

peak_ce0=ce0_xyz(round(length(xe)/2),round(length(ye)/2),round(length(ze)/2));
ce_vol_int_peak=ce_vol_int(round(length(xe)/2),round(length(ye)/2),round(length(ze)/2)); %find peak value

assert(peak_ce0-ce_vol_int_peak>0)

