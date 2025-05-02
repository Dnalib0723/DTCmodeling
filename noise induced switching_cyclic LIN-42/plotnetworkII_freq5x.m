clc
clear all
close all
N=25;
N1=250;
N2=500;

figure(1)

load('DAF1AC0_12_11_16_D95_Turn_Lin42_b3_1000xs_2024_12_25_Lin42_gx5.mat')

load('DAF1AC0_15_11_16_D95_t_Lin42_b3_1000xs_2025_01_07_Lin42_gx5.mat')

numbertest=640;

Turn=Unc5_SS_Wt_dist_turn_Noturn
Turn2=Turn(numbertest+1:numbertest+N, N1:N2)
Turn2=Turn(numbertest+1:numbertest+N, N1:N2)

mat=size(Turn2);
y = linspace(1,mat(1),N);
x=t(N1: N2);
z=Turn(numbertest+1:numbertest+N, N1: N2);
[X,Y]=meshgrid(x,y);
hold on;
pcolor(x, y,z);
pcolor(x, y,z);

%colormap summer; 
%colorbar
 
set(gca,'FontSize',70)
x_label=xlabel('Time(hr)'); set(x_label,'FontSize',70);
y_label=ylabel('Trajectories'); set(y_label,'FontSize',70);
