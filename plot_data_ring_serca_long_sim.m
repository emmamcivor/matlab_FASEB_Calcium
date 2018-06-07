% plot simulated data

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

dt=1e-2; % time step (s)
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
zs=Lextd:dzs:Le;

[X,Y]=meshgrid(y,x);

[XXj,ZZj]=meshgrid(x,zj);
[XXs,ZZs]=meshgrid(x,zs(1+round((Lb-Lextd)/dzs):end));


fn_data='/maths/scratch/bhzem/cube/FASEB_sims_20180604/simulated_data/soce_hek_I_orai_9.6e-15_5_Orai_d_orai_14_12_peripheral_SERCA2b_dt_0.1_dx_2_dy_2_dzs_5_m_150_n_150_p_100_T_100_Lb_2000_Lextd_-5000/soce-p_5.mat';
fn_serca='/maths/scratch/bhzem/cube/FASEB_sims_20180604/simulated_data/soce_hek_I_orai_9.6e-15_5_Orai_d_orai_14_12_peripheral_SERCA2b_dt_0.1_dx_2_dy_2_dzs_5_m_150_n_150_p_100_T_100_Lb_2000_Lextd_-5000/serca_activity.mat';

load(fn_data)
load(fn_serca)

cj_d14=reshape(cj_n(:,:,end),length(x),length(y),length(zj));
cs_d14=reshape(cs_n(:,:,end),length(x),length(y),length(zs));
serca_d14=SERCA_activity;

fn_data='/maths/scratch/bhzem/cube/FASEB_sims_20180604/simulated_data/soce_hek_I_orai_9.6e-15_5_Orai_d_orai_88_12_peripheral_SERCA2b_dt_0.1_dx_2_dy_2_dzs_5_m_150_n_150_p_100_T_100_Lb_2000_Lextd_-5000/soce-p_5.mat';
fn_serca='/maths/scratch/bhzem/cube/FASEB_sims_20180604/simulated_data/soce_hek_I_orai_9.6e-15_5_Orai_d_orai_88_12_peripheral_SERCA2b_dt_0.1_dx_2_dy_2_dzs_5_m_150_n_150_p_100_T_100_Lb_2000_Lextd_-5000/serca_activity.mat';

load(fn_data)
load(fn_serca)

cj_d88=reshape(cj_n(:,:,end),length(x),length(y),length(zj));
cs_d88=reshape(cs_n(:,:,end),length(x),length(y),length(zs));
serca_d88=SERCA_activity;

fn_plots_dir='/local/bhzem/R-SOCE/cube/matlab/FASEB_sims_20180604/plots_simulated_data_peripheral_serca_long_sim';
if exist(fn_plots_dir)
else
mkdir(fn_plots_dir)
end
cd(fn_plots_dir)

cmax_PM=max([max(max(cj_d14(:,:,end))) max(max(cj_d88(:,:,end)))]);
cmin_PM=min([min(min(cj_d14(:,:,end))) min(min(cj_d88(:,:,end)))]);
cmax_ERM=max([max(max(cj_d14(:,:,1))) max(max(cj_d88(:,:,1)))]);
cmin_ERM=min([min(min(cj_d14(:,:,1))) min(min(cj_d88(:,:,1)))]);

cmax_S=max([max(max(max(cs_d14))) max(max(max(cs_d88)))]);
cmin_S=min([min(min(min(cs_d14))) min(min(min(cs_d88)))]);


figure(1)
pcolor(X,Y,cj_d88(:,:,end))
shading interp
% caxis([cmin_PM cmax_PM])
colormap(jet)
colorbar
c=colorbar;
c.Label.String='C (\mu M)';
xlabel('x')
ylabel('y')
xticks([x(1) x(end)])
yticks([y(1) y(end)])
set(gca,'fontsize',30,'fontweight','bold')
set(gca,'OuterPosition',[0 0 .925 1])
savefig('dorai_88_PM')

figure(2)
pcolor(X,Y,cj_d88(:,:,1))
% caxis([cmin_ERM cmax_ERM])
shading interp
colormap(jet)
colorbar
c=colorbar;
c.Label.String='C (\mu M)';
xlabel('x')
ylabel('y')
xticks([x(1) x(end)])
yticks([y(1) y(end)])
set(gca,'fontsize',30,'fontweight','bold')
set(gca,'OuterPosition',[0 0 .925 1])
savefig('dorai_88_ERM')

figure(3)
pcolor(XXj,ZZj,squeeze(cj_d88(:,round(length(y)/2),:))')
shading interp
colorbar
c=colorbar;
c.Label.String='C (\mu M)';
xlabel('x')
xticks([x(1) x(end)])
yticks([zj(1) zj(end)])
yticklabels({'ERM','PM'})
set(gca,'fontsize',30,'fontweight','bold')
set(gca,'OuterPosition',[0 0 .925 1])
savefig('dorai_88_z')

figure(5)
pcolor(XXs,ZZs,squeeze(cs_d88(:,round(length(y)/2)+round(140/dy),1+round((Lb-Lextd)/dzs):end))')
shading interp
% caxis([cmin_S cmax_S])
colormap(jet)
colorbar
c=colorbar;
c.Label.String='C (\mu M)';
xlabel('x')
xticks([x(1) x(end)])
yticks([zj(1) zj(end)])
yticklabels({'ER_{i}','ERM'})
set(gca,'fontsize',30,'fontweight','bold')
set(gca,'OuterPosition',[0 0 .925 1])
savefig('dorai_88_sub_z')


figure(11)
pcolor(X,Y,cj_d14(:,:,end))
shading interp
% caxis([cmin_PM cmax_PM])
colormap(jet)
colorbar
c=colorbar;
c.Label.String='C (\mu M)';
xlabel('x')
ylabel('y')
xticks([x(1) x(end)])
yticks([y(1) y(end)])
set(gca,'fontsize',30,'fontweight','bold')
set(gca,'OuterPosition',[0 0 .925 1])
savefig('dorai_14_PM')


figure(12)
pcolor(X,Y,cj_d14(:,:,1))
% caxis([cmin_ERM cmax_ERM])
colormap(jet)
colorbar
c=colorbar;
c.Label.String='C (\mu M)';
xlabel('x')
ylabel('y')
xticks([x(1) x(end)])
yticks([y(1) y(end)])
set(gca,'fontsize',30,'fontweight','bold')
set(gca,'OuterPosition',[0 0 .925 1])
savefig('dorai_14_ERM')

figure(13)
pcolor(XXj,ZZj,squeeze(cj_d14(:,round(length(y)/2),:))')
shading interp
% caxis([cmin_PM cmax_PM])
colormap(jet)
colorbar
c=colorbar;
c.Label.String='C (\mu M)';
xlabel('x')
xticks([x(1) x(end)])
yticks([zj(1) zj(end)])
yticklabels({'ERM','PM'})
set(gca,'fontsize',30,'fontweight','bold')
savefig('dorai_14_z')

figure(15)
pcolor(XXs,ZZs,squeeze(cs_d14(:,round(length(y)/2)+round(140/dy),1+round((Lb-Lextd)/dzs):end))')
shading interp
% caxis([cmin_S cmax_S])
colormap(jet)
colorbar
c=colorbar;
c.Label.String='C (\mu M)';
xlabel('x')
xticks([x(1) x(end)])
yticks([zj(1) zj(end)])
yticklabels({'ER_{i}','ERM'})
set(gca,'fontsize',30,'fontweight','bold')
set(gca,'OuterPosition',[0 0 .925 1])
savefig('dorai_14_sub_z')
