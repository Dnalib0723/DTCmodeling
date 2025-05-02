%clc
%clear all
%close all

figure(4)
bin_size = 0.15;
bin_edge = 14:bin_size:18.5;
load('DAF1AC1_15_11_16_D95_Turningtime_Lin42_b3_1000xs_2025_01_07_Lin42_gx1.mat')
h1 = histc(Unc5_expr_T_Distr_Wt, bin_edge);
color_bar_sto=bar(bin_edge,h1/10 ,'histc'); 
set(color_bar_sto,'FaceColor', [0.12 0.56 1]);
set(color_bar_sto,'EdgeColor', 'none');
set(color_bar_sto,'facealpha',0.5);
hold on

load('DAF1AC1_13_11_16_D95_Turningtime_Lin42_b3_1000xs_2024_12_25_Lin42_gx5.mat')
h2 = histc (Unc5_expr_T_Distr_Wt, bin_edge); 
color_bar_sto=bar(bin_edge,h2/10 ,'histc');
set(color_bar_sto,'FaceColor', [0.98 0.5 0.48]);
set(color_bar_sto,'EdgeColor', 'none');
set(color_bar_sto,'facealpha',0.5);
hold off

%set(h1,'FaceColor',[0 0 1],'EdgeColor',[0 0 1],'facealpha',0.7);
%set(h2,'FaceColor',[1 0 0],'EdgeColor',[1 0 0],'facealpha',0);
%set(h1,'facealpha',0.5);
%set(h2,'facealpha',0.5);

%load('DAF1AC1_13_11_13_D100_mean_Lin42_b2_10xd_2024_03_13_Lin42_gx1.mat')
%mean_det= mean_Unc5_expr_T_Distr_Wt

%load('DafAC2_11_13_D95_mean_100xs_20231004_LIN42b3gx1.mat')
%mean_h1= mean_Unc5_expr_T_Distr_Wt

%load('DafAC_11_13_D95_mean_100xs_20231003_LIN42b15gx1.mat')
%mean_h2= mean_Unc5_expr_T_Distr_Wt

%x_det = [mean_det mean_det];
%x_h1 = [mean_h1 mean_h1];
%x_h2 = [mean_h2 mean_h2];
%y = [0 600];
%line(x_det,y,'Color','black', 'LineWidth',3); hold on;
%line(x_h1,y,'Color','blue','LineWidth',3); hold on;
%line(x_h2,y,'Color','red','LineWidth',3); hold off;
%set(gca,'FontSize',15)

ylim([0,30])
xlim([14 19])
xticks(14:1:19)

lgd =legend('LIN-42 noise freq 1X','LIN-42 noise freq 5X')
set(lgd,'FontSize',16);

x_label=xlabel('Time(hr)'); set(x_label,'FontSize',20);
y_label=ylabel('Fraction of turned worms'); set(y_label,'FontSize',20);
