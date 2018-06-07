
length_t=100;   %run time in time steps
dt=1e-5; % time step (s)
w=150; %half width of cube for ER-PM junction and sub-PM ER (cube runs from -w to +w) (nm)
dx=2;
dy=2;
x=-w:dx:w;
y=-w:dy:w;
Lextd=2000; 

I_CRAC=9.6e-15; %Current through CRAC channel (C/s)

SERCA_choice=2; %Choose which SERCA pump isoform is included (SERCA2a=1 and SERCA2b=2)

 situation=1; % choice of 1-6

if situation==1
%%%%%% Orai channels surrounded by a ring of SERCA pumps
%% Latticed Orai channels with ring of SERCA pumps 30nm away

%%% d_orai=14nm (lattice orai configuration 5 fully active orai channels)
x_orai=round(length(x)/2) + round([0 14 0 -14 0]/dx);
y_orai=round(length(y)/2) + round([0 0 14 0 -14]/dy);
activity_level=[1 1 1 1 1]; %Orai channel activity level
d_orai=14; 

%%% d_orai_serca=30nm when d_orai=14
x_serca=round(length(x)/2)+round([44 36 20 0 -20 -36 -44 -36 -20 0 20 36]/dx);
y_serca=round(length(y)/2)+round([0 20 36 44 36 20 0 -20 -36 -44 -36 -20]/dy);
d_serca=30;

elseif situation==2
%% non-latticed Orai channels (d_orai=88nm) with ring of SERCA pumps 30nm away
 
%%% d_orai=88nm (non-lattice orai configuration 4 fully active orai channels)
x_orai=round(length(x)/2) + round([76 24 -62 -62 24]/dx);
y_orai=round(length(y)/2) + round([0 72 44 -44 -72]/dy);
activity_level=[0.125 0.125 0.125 0.125 0.125]; %Orai channel activity level
d_orai=88; 

%%% d_orai_serca=30nm when d_orai=88
% x_serca=round(length(x)/2)+round([106 90 52 0 -52 -90 -106 -96 -62 -14 40 82 104]/dx);
% 
% y_serca=round(length(y)/2)+round([6 56 92 106 92 56 6 -46 -86 -106 -98 -68 -20]/dy);

x_serca=round(length(x)/2)+round([106 86 46 -6 -56 -92 -106 -92 -56 -6 46 86]/dx);

y_serca=round(length(y)/2)+round([14 62 96 106 90 52 0 -52 -90 -106 -96 -62]/dy);


d_serca=30;

elseif situation==3
%%%% SERCA pumps equidistant along periphery of ER-PM junction
%% Latticed Orai channels with SERCA pumps placed at periphery of ER-PM junction

%%% d_orai=14nm (lattice orai configuration 5 fully active orai channels)
x_orai=round(length(x)/2) + round([0 14 0 -14 0]/dx);
y_orai=round(length(y)/2) + round([0 0 14 0 -14]/dy);
activity_level=[1 1 1 1 1]; %Orai channel activity level
d_orai=14; 

%%% peripheral SERCA pumps
x_serca=round(length(x)/2)+round([140 48 -46 -140 -140 -140 -140 -46 48 140 140 140]/dx);
y_serca=round(length(y)/2)+round([140 140 140 140 48 -46 -140 -140 -140 -140 -46 48]/dy);
d_serca=1;

elseif situation==4
%% non-latticed Orai channels (d_orai=88nm) with SERCA pumps placed at periphery of ER-PM junction

%%% d_orai=88nm (non-lattice orai configuration 4 fully active orai channels)
x_orai=round(length(x)/2) + round([76 24 -62 -62 24]/dx);
y_orai=round(length(y)/2) + round([0 72 44 -44 -72]/dy);
activity_level=[0.125 0.125 0.125 0.125 0.125]; %Orai channel activity level
d_orai=88; 

%%% peripheral SERCA pumps
x_serca=round(length(x)/2)+round([140 48 -46 -140 -140 -140 -140 -46 48 140 140 140]/dx);
y_serca=round(length(y)/2)+round([140 140 140 140 48 -46 -140 -140 -140 -140 -46 48]/dy);
d_serca=1;

elseif situation==5
%%%% SERCA pumps equidistant along periphery of ER-PM junction - extra
%%%% pumps depending on space available within junction
%% Latticed Orai channels with SERCA pumps placed at periphery of ER-PM junction

%%% d_orai=14nm (lattice orai configuration 5 fully active orai channels)
x_orai=round(length(x)/2) + round([0 14 0 -14 0]/dx);
y_orai=round(length(y)/2) + round([0 0 14 0 -14]/dy);
activity_level=[1 1 1 1 1]; %Orai channel activity level
d_orai=14; 

%%% peripheral SERCA pumps
x_serca_1=round(length(x)/2)+round([140 48 -46 -140 -140 -140 -140 -46 48 140 140 140]/dx);
y_serca_1=round(length(y)/2)+round([140 140 140 140 48 -46 -140 -140 -140 -140 -46 48]/dy);

x_serca_2=round(length(x)/2)+round([106 36 -36 -106 -106 -106 -106 -36 36 106 106 106]/dx);
y_serca_2=round(length(y)/2)+round([106 106 106 106 36 -36 -106 -106 -106 -106 -36 36]/dy);

x_serca=[x_serca_1, x_serca_2];
y_serca=[y_serca_1, y_serca_2];

d_serca=1;

elseif situation==6
%% non-latticed Orai channels (d_orai=88nm) with SERCA pumps placed at periphery of ER-PM junction

%%% d_orai=88nm (non-lattice orai configuration 4 fully active orai channels)
x_orai=round(length(x)/2) + round([76 24 -62 -62 24]/dx);
y_orai=round(length(y)/2) + round([0 72 44 -44 -72]/dy);
activity_level=[0.125 0.125 0.125 0.125 0.125]; %Orai channel activity level
d_orai=88; 

%%% peripheral SERCA pumps
x_serca_1=round(length(x)/2)+round([140 48 -46 -140 -140 -140 -140 -46 48 140 140 140]/dx);
y_serca_1=round(length(y)/2)+round([140 140 140 140 48 -46 -140 -140 -140 -140 -46 48]/dy);

x_serca_2=round(length(x)/2)+round([106 0 -106 -106 -106 0 106 106]/dx);
y_serca_2=round(length(y)/2)+round([106 106 106 0 -106 -106 -106 0]/dy);

x_serca=[x_serca_1, x_serca_2];
y_serca=[y_serca_1, y_serca_2];

d_serca=1;

end

run soce_js_HEK.m

