%%% simulate SOCE in the ER-PM junction
%% Setup parameters
Lp=2515; %distance to PM (nm)
Le=2500; %distance to ER membrane (nm)
Lb=2000; %distance to interface between sub-PM ER and bulk ER (nm)
L0=0; %internal point of bulk ER, we think of this as the origin in z direction

q=1000; %half width of cube for bulk ER. Chosen so STIM1 can reside anywhere on cube and still localise in ER-PM junction (nm)

Dj=220e6; % Diffusion in ER-PM junction - assume no buffers (nm^2/s)
De=10e6; %Diffusion in sub-PM ER and bulk ER - estimate from Sweitach 2008 and Dayel 1999 (nm^2/s)

%Volume of cytosol (Berlin, Bassani, Bers 1994)
vol_cyt_L=24.5e-12; %(L) volume of cytosol of ventricular myocyte in litres
vol_cyt_m=vol_cyt_L*1e-3;    %(m^3) volume of cytosol of ventricular myocyte in metres
vol_cyt_nm=vol_cyt_m*1e27;

%Volume of Sub-PM ER
vol_subPMER_nm=(Le-Lextd)*4*w^2; %nm^3
vol_subPMER_m=vol_subPMER_nm*1e-27; %m^3
vol_subPMER_L=vol_subPMER_m*1e3; %volume in litres (L)
 
n_timesteps=20;
q=repmat(1:n_timesteps,1,round(length_t/n_timesteps));
P=round(length_t/n_timesteps);  

%% domain parameters ER-PM junction

dzj=0.15; 
zj=Le:dzj:Lp;

mj=1:w; %number of eigenvalues 
nj=1:w; %number of eigenvalues 
pj=1:(Lp-Le)/dzj; %number of eigenvalues 

mu_j=mj*pi/(2*w); %eigenvalues
eta_j=nj*pi/(2*w); %eigenvalues
alpha_j=pj*pi/(Lp-Le); %eigenvalues

%% domain parameters sub-PM ER
dzs=5;
zs=Lextd:dzs:Le;

ms=1:1:w;
ns=1:1:w;
ps=1:(Le-Lextd)/dzs;

mu_s=ms*pi/(2*w);
eta_s=ns*pi/(2*w);
alpha_s=(2*ps-1)*pi/(2*(Le-Lextd));

%% channel/pump flux parameters
%Orai channel flux
F=96485.3365; %Faraday's constant (C/mol)
z_val=2; %Valency of calcium ions (dimensionless)
f_CRAC=I_CRAC/(F*z_val); %moles of calcium entering microdomain per second (mol/s)
f_CRAC_micro_mol=f_CRAC*1e6;    %micro moles per second
J_CRAC=f_CRAC_micro_mol*vol_cyt_nm/vol_cyt_L;
flux_CRAC_magnitude=J_CRAC/Dj;

% SERCA pump parameters (Lytton 1992)
% SERCA2b
Kmf_SERCA2b=0.13;
n_H_SERCA2b=2;
Vmax_ions_SERCA2b=9.1; %Ca2+ per second (according to Hogan 2015) 

%SERCA2a
Kmf_SERCA2a=0.21;
n_H_SERCA2a=2;
Vmax_ions_SERCA2a=20; %Ca2+ per second (according to Hogan 2015) 
Kmr_Bers1998=1700; %micro molar
Q_Shannon2004=2.6; %Temperature factor

%% SERCA2a and SERCA2b pump parameters

if SERCA_choice==1
    disp('SERCA2a')
    Kmf=Kmf_SERCA2a; 
    n_H=n_H_SERCA2a;    %Hill coefficient - I rounded up to 2 for a while, trying Bers again now
    n_avogadro=6.022e23;
    Vmax=(Vmax_ions_SERCA2a/n_avogadro)*1e6; % (36 ca ions/ avogradros number (ions per mol) * 1e6 = X micro moles per second)
elseif SERCA_choice==2
    disp('SERCA2b')
    Kmf=Kmf_SERCA2b; 
    n_H=n_H_SERCA2b;    %Hill coefficient - I rounded up to 2 for a while, trying Bers again now
    n_avogadro=6.022e23;
    Vmax=(Vmax_ions_SERCA2b/n_avogadro)*1e6; % (36 ca ions/ avogradros number (ions per mol) * 1e6 = X micro moles per second)
end

%% initial and boundary conditions
%ER-PM junction
cj0=0.1; %micro moles per litre cytosol
cj=cj0*ones(length(x),length(y),length(zj));
cj_DBC=cj0;

%sub-PM ER
cs0=150; %micro moles per litre cytosol
cs=cs0*ones(length(x),length(y),length(zs));
cs_DBC=cs0;
%% create greens function 
% ER-PM junction
[gzj,gxyj] = fn_gfj(Lp,Le,w,x,y,zj,mu_j,eta_j,alpha_j,Dj,dt);

% sub-PM ER
[gzs,gxys] = fn_gfs(Le,Lextd,w,x,y,zs,mu_s,eta_s,alpha_s,De,dt);

%% create steady state solution for Orai channels on PM
n_orai=length(x_orai);

ss_orai_x_y_z=zeros(length(x),length(y),length(zj));


for j=1:n_orai
x_orai_i=x_orai(j);
y_orai_i=y_orai(j);

ss_orai_x_y_z = ss_orai_x_y_z + activity_level(j)*fn_ss_PM_j(eta_j,mu_j,x,y,zj,w,x_orai_i,y_orai_i,Lp,Le);
end


ss_orai_xyz=reshape(ss_orai_x_y_z,length(x)*length(y),length(zj));
ss_orai_xyz_total=flux_CRAC_magnitude*ss_orai_xyz;

%% create steady state solution for SERCA pumps on ER membrane
% ER-PM junction and sub-PM ER
n_serca=length(x_serca);

ss_serca_xyzn_j=zeros(length(x)*length(y),length(zj),n_serca);
ss_serca_xyzn_s=zeros(length(x)*length(y),length(zs),n_serca);


for j=1:n_serca
x_serca_i=x_serca(j);
y_serca_i=y_serca(j);

ss_serca_xyzn_j(:,:,j) = fn_ss_ERM_j(eta_j,mu_j,x,y,zj,w,x_serca_i,y_serca_i,Lp,Le);
ss_serca_xyzn_s(:,:,j) = fn_ss_ERM_s(eta_s,mu_s,x,y,zs,w,x_serca_i,y_serca_i,Le,Lextd);

end


%% numerically time step solution
cj_n=zeros(length(x)*length(y),length(zj),n_timesteps);
cs_n=zeros(length(x)*length(y),length(zs),n_timesteps);

if d_serca==30
if SERCA_choice==1 
    fn_data_dir=['simulated_data/soce_hek_I_orai_',num2str(I_CRAC),'_',num2str(n_orai),'_Orai_d_orai_',num2str(d_orai),'_',num2str(n_serca),'_SERCA2a_30nm_ring_dt_',num2str(dt),'_dx_',num2str(dx),'_dy_',num2str(dy),'_dzs_',num2str(dzs),'_m_',num2str(length(mj)),'_n_',num2str(length(nj)),'_p_',num2str(length(pj)),'_T_',num2str(length_t),'_Lb_',num2str(Lb),'_Lextd_',num2str(Lextd),'/'];
elseif SERCA_choice==2
    fn_data_dir=['simulated_data/soce_hek_I_orai_',num2str(I_CRAC),'_',num2str(n_orai),'_Orai_d_orai_',num2str(d_orai),'_',num2str(n_serca),'_SERCA2b_30nm_ring_dt_',num2str(dt),'_dx_',num2str(dx),'_dy_',num2str(dy),'_dzs_',num2str(dzs),'_m_',num2str(length(mj)),'_n_',num2str(length(nj)),'_p_',num2str(length(pj)),'_T_',num2str(length_t),'_Lb_',num2str(Lb),'_Lextd_',num2str(Lextd),'/'];
end
else
if SERCA_choice==1 
    fn_data_dir=['simulated_data/soce_hek_I_orai_',num2str(I_CRAC),'_',num2str(n_orai),'_Orai_d_orai_',num2str(d_orai),'_',num2str(n_serca),'_peripheral_SERCA2a_dt_',num2str(dt),'_dx_',num2str(dx),'_dy_',num2str(dy),'_dzs_',num2str(dzs),'_m_',num2str(length(mj)),'_n_',num2str(length(nj)),'_p_',num2str(length(pj)),'_T_',num2str(length_t),'_Lb_',num2str(Lb),'_Lextd_',num2str(Lextd),'/'];
elseif SERCA_choice==2
    fn_data_dir=['simulated_data/soce_hek_I_orai_',num2str(I_CRAC),'_',num2str(n_orai),'_Orai_d_orai_',num2str(d_orai),'_',num2str(n_serca),'_peripheral_SERCA2b_dt_',num2str(dt),'_dx_',num2str(dx),'_dy_',num2str(dy),'_dzs_',num2str(dzs),'_m_',num2str(length(mj)),'_n_',num2str(length(nj)),'_p_',num2str(length(pj)),'_T_',num2str(length_t),'_Lb_',num2str(Lb),'_Lextd_',num2str(Lextd),'/'];
end
end

if exist(fn_data_dir)
else
    mkdir(fn_data_dir);
end

SERCA_activity=zeros(n_serca,length_t);

for p=1:P
    
if d_serca==30
if SERCA_choice==1 
    fn_data=['simulated_data/soce_hek_I_orai_',num2str(I_CRAC),'_',num2str(n_orai),'_Orai_d_orai_',num2str(d_orai),'_',num2str(n_serca),'_SERCA2a_30nm_ring_dt_',num2str(dt),'_dx_',num2str(dx),'_dy_',num2str(dy),'_dzs_',num2str(dzs),'_m_',num2str(length(mj)),'_n_',num2str(length(nj)),'_p_',num2str(length(pj)),'_T_',num2str(length_t),'_Lb_',num2str(Lb),'_Lextd_',num2str(Lextd),'/soce-p_',num2str(p),'.mat'];
elseif SERCA_choice==2
    fn_data=['simulated_data/soce_hek_I_orai_',num2str(I_CRAC),'_',num2str(n_orai),'_Orai_d_orai_',num2str(d_orai),'_',num2str(n_serca),'_SERCA2b_30nm_ring_dt_',num2str(dt),'_dx_',num2str(dx),'_dy_',num2str(dy),'_dzs_',num2str(dzs),'_m_',num2str(length(mj)),'_n_',num2str(length(nj)),'_p_',num2str(length(pj)),'_T_',num2str(length_t),'_Lb_',num2str(Lb),'_Lextd_',num2str(Lextd),'/soce-p_',num2str(p),'.mat'];
end
else
if SERCA_choice==1 
    fn_data=['simulated_data/soce_hek_I_orai_',num2str(I_CRAC),'_',num2str(n_orai),'_Orai_d_orai_',num2str(d_orai),'_',num2str(n_serca),'_peripheral_SERCA2a_dt_',num2str(dt),'_dx_',num2str(dx),'_dy_',num2str(dy),'_dzs_',num2str(dzs),'_m_',num2str(length(mj)),'_n_',num2str(length(nj)),'_p_',num2str(length(pj)),'_T_',num2str(length_t),'_Lb_',num2str(Lb),'_Lextd_',num2str(Lextd),'/soce-p_',num2str(p),'.mat'];
elseif SERCA_choice==2
    fn_data=['simulated_data/soce_hek_I_orai_',num2str(I_CRAC),'_',num2str(n_orai),'_Orai_d_orai_',num2str(d_orai),'_',num2str(n_serca),'_peripheral_SERCA2b_dt_',num2str(dt),'_dx_',num2str(dx),'_dy_',num2str(dy),'_dzs_',num2str(dzs),'_m_',num2str(length(mj)),'_n_',num2str(length(nj)),'_p_',num2str(length(pj)),'_T_',num2str(length_t),'_Lb_',num2str(Lb),'_Lextd_',num2str(Lextd),'/soce-p_',num2str(p),'.mat'];
end
end
     

    for i=(p-1)*n_timesteps+1:p*n_timesteps

ss_serca_xyzn_total_j=zeros(length(x)*length(y),length(zj));
ss_serca_xyzn_total_s=zeros(length(x)*length(y),length(zs));

            for j=1:n_serca
            x_serca_i=x_serca(j);
            y_serca_i=y_serca(j);

        serca_activity = fn_serca_activity(cj,cs,x_serca_i,y_serca_i,Kmf,Kmr_Bers1998,n_H);
        J_SERCA_flux_magnitude_shannon=Q_Shannon2004*Vmax*serca_activity;

        ss_serca_xyzn_total_j=ss_serca_xyzn_total_j+J_SERCA_flux_magnitude_shannon*ss_serca_xyzn_j(:,:,j)*vol_cyt_nm/vol_cyt_L/Dj;
        ss_serca_xyzn_total_s=ss_serca_xyzn_total_s+J_SERCA_flux_magnitude_shannon*ss_serca_xyzn_s(:,:,j)*vol_subPMER_nm/vol_subPMER_L/De;
        SERCA_activity(j,i)=serca_activity;
 
            end        
        cj_z_xy_prev=reshape(cj,length(x)*length(y),length(zj))';
        cs_z_xy_prev=reshape(cs,length(x)*length(y),length(zs))';
        
        ss_all_j=ss_orai_xyz_total' + ss_serca_xyzn_total_j' + cj_DBC;
        ss_all_s = ss_serca_xyzn_total_s' + cs_DBC;
        
        cj = vol_int(gxyj,gzj,x,y,zj,dx,dy,dzj,cj_z_xy_prev,ss_all_j);        
        cj_n(:,:,q(i)) = reshape(cj,length(x)*length(y),length(zj));
        
        cs = vol_int(gxys,gzs,x,y,zs,dx,dy,dzs,cs_z_xy_prev,ss_all_s);        
        cs_n(:,:,q(i)) = reshape(cs,length(x)*length(y),length(zs));
    end
    
if save_data==1
if p==1
    
    save(fn_data,'cj_n','cs_n','-v7.3')
    
    p
    elseif p==(P/2)
    
    save(fn_data,'cj_n','cs_n','-v7.3')
    
    p
    elseif p==P
    save(fn_data,'cj_n','cs_n','-v7.3')
    p
   
end
end

end

if d_serca==30
if SERCA_choice==1 
    fn_serca=['simulated_data/soce_hek_I_orai_',num2str(I_CRAC),'_',num2str(n_orai),'_Orai_d_orai_',num2str(d_orai),'_',num2str(n_serca),'_SERCA2a_30nm_ring_dt_',num2str(dt),'_dx_',num2str(dx),'_dy_',num2str(dy),'_dzs_',num2str(dzs),'_m_',num2str(length(mj)),'_n_',num2str(length(nj)),'_p_',num2str(length(pj)),'_T_',num2str(length_t),'_Lb_',num2str(Lb),'_Lextd_',num2str(Lextd),'/serca_actvity.mat'];
elseif SERCA_choice==2
    fn_serca=['simulated_data/soce_hek_I_orai_',num2str(I_CRAC),'_',num2str(n_orai),'_Orai_d_orai_',num2str(d_orai),'_',num2str(n_serca),'_SERCA2b_30nm_ring_dt_',num2str(dt),'_dx_',num2str(dx),'_dy_',num2str(dy),'_dzs_',num2str(dzs),'_m_',num2str(length(mj)),'_n_',num2str(length(nj)),'_p_',num2str(length(pj)),'_T_',num2str(length_t),'_Lb_',num2str(Lb),'_Lextd_',num2str(Lextd),'/serca_actvity.mat'];
end
else
if SERCA_choice==1 
    fn_serca=['simulated_data/soce_hek_I_orai_',num2str(I_CRAC),'_',num2str(n_orai),'_Orai_d_orai_',num2str(d_orai),'_',num2str(n_serca),'_peripheral_SERCA2a_dt_',num2str(dt),'_dx_',num2str(dx),'_dy_',num2str(dy),'_dzs_',num2str(dzs),'_m_',num2str(length(mj)),'_n_',num2str(length(nj)),'_p_',num2str(length(pj)),'_T_',num2str(length_t),'_Lb_',num2str(Lb),'_Lextd_',num2str(Lextd),'/serca_activity.mat'];
elseif SERCA_choice==2
    fn_serca=['simulated_data/soce_hek_I_orai_',num2str(I_CRAC),'_',num2str(n_orai),'_Orai_d_orai_',num2str(d_orai),'_',num2str(n_serca),'_peripheral_SERCA2b_dt_',num2str(dt),'_dx_',num2str(dx),'_dy_',num2str(dy),'_dzs_',num2str(dzs),'_m_',num2str(length(mj)),'_n_',num2str(length(nj)),'_p_',num2str(length(pj)),'_T_',num2str(length_t),'_Lb_',num2str(Lb),'_Lextd_',num2str(Lextd),'/serca_activity.mat'];
end
end

if save_data==1
save(fn_serca,'SERCA_activity')
end

