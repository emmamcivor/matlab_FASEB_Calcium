% parameters
Lp=2515; %distance to PM (nm)
Le=2500; %distance to ER membrane (nm)
Lb=2000; %distance to interface between sub-PM ER and bulk ER (nm)
L0=0; %internal point of bulk ER, we think of this as the origin in z direction
Lextd=-5000;

w=150; %half width of cube for ER-PM junction and sub-PM ER (cube runs from -w to +w) (nm)
q=1000; %half width of cube for bulk ER. Chosen so STIM1 can reside anywhere on cube and still localise in ER-PM junction (nm)

Dj=220e6; % Diffusion in ER-PM junction - assume no buffers (nm^2/s)
De=10e6; %Diffusion in sub-PM ER and bulk ER - estimate from Sweitach 2008 and Dayel 1999 (nm^2/s)

%Volume of cytosol (Berlin, Bassani, Bers 1994)
vol_cyt_L=24.5e-12; %(L) volume of cytosol of ventricular myocyte in litres
vol_cyt_m=vol_cyt_L*1e-3;    %(m^3) volume of cytosol of ventricular myocyte in metres
vol_cyt_nm=vol_cyt_m*1e27;

dt=1e-5; % time step (s)
length_t=100;   %run time in time steps
T=length_t*dt;
n_timesteps=20;

t=0:dt:T; %s
t=t*1e3; %ms

dx=2;
dy=2;
dzj=0.15; %(Lp-Le)/100;
dzs=5;

x=-w:dx:w;
y=-w:dy:w;
zj=Le:dzj:Lp;
zs=Lb:dzs:Le;
zs_extd=Lextd:dzs:Le;

trapz_x= dx*[0.5 ones(1,length(x)-2) 0.5];
trapz_y= dy*[0.5 ones(1,length(y)-2) 0.5];
trapz_z= dzs*[0.5 ones(1,length(zs)-2) 0.5];
trapz_xy= kron(trapz_y,trapz_x);

fn_d14_ring_ms='/maths/scratch/bhzem/cube/FASEB_sims_20180604/simulated_data/soce_hek_I_orai_9.6e-15_5_Orai_d_orai_14_12_SERCA2b_30nm_ring_dt_1e-05_dx_2_dy_2_dzs_5_m_150_n_150_p_100_T_100_Lb_2000_Lextd_2000/soce-p_5.mat';
load(fn_d14_ring_ms)
cs_d14_ring_ms=reshape(cs_n(:,:,end),length(x)*length(y),length(zs));
int_d14_ring_ms=trapz_z*(cs_d14_ring_ms'-150)*trapz_xy';

fn_d88_ring_ms='/maths/scratch/bhzem/cube/FASEB_sims_20180604/simulated_data/soce_hek_I_orai_9.6e-15_5_Orai_d_orai_14_12_SERCA2b_30nm_ring_dt_1e-05_dx_2_dy_2_dzs_5_m_150_n_150_p_100_T_100_Lb_2000_Lextd_2000/soce-p_5.mat';
load(fn_d88_ring_ms)
cs_d88_ring_ms=reshape(cs_n(:,:,end),length(x)*length(y),length(zs));
int_d88_ring_ms=trapz_z*(cs_d88_ring_ms'-150)*trapz_xy';

fn_d14_peri_ms='/maths/scratch/bhzem/cube/FASEB_sims_20180604/simulated_data/soce_hek_I_orai_9.6e-15_5_Orai_d_orai_14_12_peripheral_SERCA2b_dt_1e-05_dx_2_dy_2_dzs_5_m_150_n_150_p_100_T_100_Lb_2000_Lextd_2000/soce-p_5.mat';
load(fn_d14_peri_ms)
cs_d14_peri_ms=reshape(cs_n(:,:,end),length(x)*length(y),length(zs));
int_d14_peri_ms=trapz_z*(cs_d14_peri_ms'-150)*trapz_xy';

fn_d88_peri_ms='/maths/scratch/bhzem/cube/FASEB_sims_20180604/simulated_data/soce_hek_I_orai_9.6e-15_5_Orai_d_orai_88_12_peripheral_SERCA2b_dt_1e-05_dx_2_dy_2_dzs_5_m_150_n_150_p_100_T_100_Lb_2000_Lextd_2000/soce-p_5.mat';
load(fn_d88_peri_ms)
cs_d88_peri_ms=reshape(cs_n(:,:,end),length(x)*length(y),length(zs));
int_d88_peri_ms=trapz_z*(cs_d88_peri_ms'-150)*trapz_xy';

fn_d14_peri_10s='/maths/scratch/bhzem/cube/FASEB_sims_20180604/simulated_data/soce_hek_I_orai_9.6e-15_5_Orai_d_orai_14_12_peripheral_SERCA2b_dt_0.1_dx_2_dy_2_dzs_5_m_150_n_150_p_100_T_100_Lb_2000_Lextd_-5000/soce-p_5.mat';
load(fn_d14_peri_10s)
cs_d14_peri_10s=reshape(cs_n(:,1+round((Lb-Lextd)/dzs):end,end),length(x)*length(y),length(zs));
int_d14_peri_10s=trapz_z*(cs_d14_peri_10s'-150)*trapz_xy';

fn_d88_peri_10s='/maths/scratch/bhzem/cube/FASEB_sims_20180604/simulated_data/soce_hek_I_orai_9.6e-15_5_Orai_d_orai_88_12_peripheral_SERCA2b_dt_0.1_dx_2_dy_2_dzs_5_m_150_n_150_p_100_T_100_Lb_2000_Lextd_-5000/soce-p_5.mat';
load(fn_d88_peri_10s)
cs_d88_peri_10s=reshape(cs_n(:,1+round((Lb-Lextd)/dzs):end,end),length(x)*length(y),length(zs));
int_d88_peri_10s=trapz_z*(cs_d88_peri_10s'-150)*trapz_xy';

fn_d14_peri_ms_add='/maths/scratch/bhzem/cube/FASEB_sims_20180604/simulated_data/soce_hek_I_orai_9.6e-15_5_Orai_d_orai_14_24_peripheral_SERCA2b_dt_1e-05_dx_2_dy_2_dzs_5_m_150_n_150_p_100_T_100_Lb_2000_Lextd_2000/soce-p_5.mat';
load(fn_d14_peri_ms_add)
cs_d14_peri_ms_add=reshape(cs_n(:,:,end),length(x)*length(y),length(zs));
int_d14_peri_ms_add=trapz_z*(cs_d14_peri_ms_add'-150)*trapz_xy';

fn_d88_peri_ms_add='/maths/scratch/bhzem/cube/FASEB_sims_20180604/simulated_data/soce_hek_I_orai_9.6e-15_5_Orai_d_orai_88_20_peripheral_SERCA2b_dt_1e-05_dx_2_dy_2_dzs_5_m_150_n_150_p_100_T_100_Lb_2000_Lextd_2000/soce-p_5.mat';
load(fn_d88_peri_ms_add)
cs_d88_peri_ms_add=reshape(cs_n(:,:,end),length(x)*length(y),length(zs));
int_d88_peri_ms_add=trapz_z*(cs_d88_peri_ms_add'-150)*trapz_xy';

fn_d14_peri_10s_add='/maths/scratch/bhzem/cube/FASEB_sims_20180604/simulated_data/soce_hek_I_orai_9.6e-15_5_Orai_d_orai_14_24_peripheral_SERCA2b_dt_0.1_dx_2_dy_2_dzs_5_m_150_n_150_p_100_T_100_Lb_2000_Lextd_-5000/soce-p_5.mat';
load(fn_d14_peri_10s_add)
cs_d14_peri_10s_add=reshape(cs_n(:,1+round((Lb-Lextd)/dzs):end,end),length(x)*length(y),length(zs));
int_d14_peri_10s_add=trapz_z*(cs_d14_peri_10s_add'-150)*trapz_xy';

fn_d88_peri_10s_add='/maths/scratch/bhzem/cube/FASEB_sims_20180604/simulated_data/soce_hek_I_orai_9.6e-15_5_Orai_d_orai_88_20_peripheral_SERCA2b_dt_0.1_dx_2_dy_2_dzs_5_m_150_n_150_p_100_T_100_Lb_2000_Lextd_-5000/soce-p_5.mat';
load(fn_d88_peri_10s_add)
cs_d88_peri_10s_add=reshape(cs_n(:,1+round((Lb-Lextd)/dzs):end,end),length(x)*length(y),length(zs));
int_d88_peri_10s_add=trapz_z*(cs_d88_peri_10s_add'-150)*trapz_xy';

ratio_ring_ms=int_d88_ring_ms/int_d14_ring_ms;
% ratio_ring_10s=int_d88_ring_10s/int_d14_ring_10s;
ratio_peri_ms=int_d88_peri_ms/int_d14_peri_ms;
ratio_peri_10s=int_d88_peri_10s/int_d14_peri_10s;
ratio_peri_ms_add=int_d88_peri_ms_add/int_d14_peri_ms_add;
ratio_peri_10s_add=int_d88_peri_10s_add/int_d14_peri_10s_add;

ratios=[ratio_ring_ms; ratio_peri_ms; ratio_peri_10s; ratio_peri_ms_add; ratio_peri_10s_add];
ints_d14=[int_d14_ring_ms; int_d14_peri_ms; int_d14_peri_10s; int_d14_peri_ms_add; int_d14_peri_10s_add];
ints_d88=[int_d88_ring_ms; int_d88_peri_ms; int_d88_peri_10s; int_d88_peri_ms_add; int_d88_peri_10s_add];

fn_ratios_ints='compare_magnitude_ER_refilling.mat';
save(fn_ratios_ints,'ratios','ints_d14','ints_d88');

