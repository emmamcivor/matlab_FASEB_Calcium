% parameters
Lp=2500; %distance to PM (nm)
Le=2485; %distance to ER membrane (nm)
Lb=2000;

w=100; %half width of cube for ER-PM junction and sub-PM ER (cube runs from -w to +w) (nm)
dx=5;
dy=5;
dzj=1.5; 
dzs=48.5; 

x=-w:dx:w;
y=-w:dy:w;
zj=Le:dzj:Lp;
zs=Lb:dzs:Le;

x_serca_i=round(length(x)/2+30/dx);
y_serca_i=round(length(y)/2+30/dy);

%ER-PM junction
cj0=0.1; %micro moles per litre cytosol
cj=cj0*ones(length(x),length(y),length(zj));

%sub-PM ER
cs0=150; %micro moles per litre cytosol
cs=cs0*ones(length(x),length(y),length(zj));

% SERCA pump parameters
Kmf_SERCA2b=0.27;
n_H_SERCA2b=1.7;
Kmr_Bers1998=1700; %micro molar
Kmf=Kmf_SERCA2b; 
n_H=n_H_SERCA2b;    %Hill coefficient - I rounded up to 2 for a while, trying Bers again now


% Assuming that sub-PM ER calcium concentration is constant to check how
% serca activity varies with cytoplasmic calcium concentration

%% Test SERCA activity returns a number

serca_activity = fn_serca_activity(cj,cs,x_serca_i,y_serca_i,Kmf,Kmr_Bers1998,n_H);

assert(isnumeric(serca_activity))

%% Test SERCA activity at baseline calcium returns a number whose magnitude is between 0 and 1 

serca_activity = fn_serca_activity(cj,cs,x_serca_i,y_serca_i,Kmf,Kmr_Bers1998,n_H);

assert(abs(serca_activity)>=0 && abs(serca_activity)<=1)

%% Test SERCA activity at baseline calcium is sensible

serca_activity = fn_serca_activity(cj,cs,x_serca_i,y_serca_i,Kmf,Kmr_Bers1998,n_H);

assert(abs(serca_activity)>= 0.1)

%% Test SERCA activity at high calcium returns a number whose magnitude is between 0 and 1 

cj_high=10; %micro moles per litre cytosol
cj_high=cj_high*ones(length(x),length(y),length(zj));

serca_activity_high = fn_serca_activity(cj_high,cs,x_serca_i,y_serca_i,Kmf,Kmr_Bers1998,n_H);

assert(abs(serca_activity_high)>=0 && abs(serca_activity_high)<=1)

%% Test SERCA activity at high calcium is greater than calcium concentration at baseline (saturating function at high calcium)
cj_high=10; %micro moles per litre cytosol
cj_high=cj_high*ones(length(x),length(y),length(zj));

serca_activity_high = fn_serca_activity(cj_high,cs,x_serca_i,y_serca_i,Kmf,Kmr_Bers1998,n_H);
serca_activity_baseline = fn_serca_activity(cj,cs,x_serca_i,y_serca_i,Kmf,Kmr_Bers1998,n_H);

assert(abs(serca_activity_high)> abs(serca_activity_baseline))
