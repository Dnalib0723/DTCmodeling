figure(1);

load('DAF1AC0_11_16_D95_t_Lin42_1_b3_10xs_2025_01_09_b3gx1.mat')
Time=t;

load('DAF1AC1_11_16_D95_Unc5_Lin42_0_b3_10xd_2025_01_08_b3gx5.mat')
DAF1AC1_Y0= Unc5_sto(1,:);

load('DAF1AC1_11_16_D95_Unc5_Lin42_5.000000e-01_b3_10xd_2025_01_08_b3gx5.mat')
DAF1AC1_Y136= Unc5_sto(1,:); 

load('DAF1AC1_11_16_D95_Unc5_Lin42_1_b3_10xd_2025_01_08_b3gx5.mat')
DAF1AC1_Y272= Unc5_sto(1,:);

subplot(3,3,1); hold on; box on; set(gca,'LineWidth',2.5); set(gca,'FontSize',15)
newcolors = [0 0 0
          0 0 1 
          1 0 0];
colororder(newcolors)
plot(Time, DAF1AC1_Y0,'LineWidth',2.2)
hold on
plot(Time, DAF1AC1_Y136,'LineWidth',2.0)
plot(Time, DAF1AC1_Y272,'LineWidth',1.9)
hold off
axis square
lgd=legend('0', '136' ,'272');
lgd.FontSize=10;
%h_title=title('DAFAC_ det'); set(h_title,'FontSize',15)
%x_label=xlabel(' Time(hr) '); set(x_label,'FontSize',15);
y_label=ylabel(' UNC-5 '); set(y_label,'FontSize',15);
ylim([0,250]);
x = [0 Time(end)];
y = [95 95];
line(x,y,'Color','black','LineStyle',':','LineWidth',1.5)

load('DAF1AC1_11_16_D95_Unc5_Lin42_0_b3_100xs_2025_01_08_b3gx1.mat')
DAF1AC1_Y0= mean(Unc5_sto);

load('DAF1AC1_11_16_D95_Unc5_Lin42_5.000000e-01_b3_100xs_2025_01_08_b3gx1.mat')
DAF1AC1_Y136= mean(Unc5_sto);

load('DAF1AC1_11_16_D95_Unc5_Lin42_1_b3_100xs_2025_01_08_b3gx1.mat')
DAF1AC1_Y272= mean(Unc5_sto);


subplot(3,3,2); hold on; box on; set(gca,'LineWidth',2.5); set(gca,'FontSize',15)
newcolors = [0 0 0
          0 0 1 
          1 0 0];
colororder(newcolors)
plot(Time, DAF1AC1_Y0,'LineWidth',2.2)
hold on
plot(Time, DAF1AC1_Y136,'LineWidth',2.0)
plot(Time, DAF1AC1_Y272,'LineWidth',1.9)
hold off
axis square

h_title=title('LIN-42 noise freq 1X'); set(h_title,'FontSize',15)
%x_label=xlabel('Time(hr)'); set(x_label,'FontSize',15);
%y_label=ylabel(' UNC-5 '); set(y_label,'FontSize',15);
ylim([0,250]);
x = [0 Time(end)];
y = [95 95];
line(x,y,'Color','black','LineStyle',':','LineWidth',1.5)
xlim([0,25]);


load('DAF1AC1_11_16_D95_Unc5_Lin42_0_b3_100xs_2025_01_08_b3gx5.mat')
DAF1AC1_Y0= mean(Unc5_sto);

load('DAF1AC1_11_16_D95_Unc5_Lin42_5.000000e-01_b3_100xs_2025_01_08_b3gx5.mat')
DAF1AC1_Y136= mean(Unc5_sto);

load('DAF1AC1_11_16_D95_Unc5_Lin42_1_b3_100xs_2025_01_08_b3gx5.mat')
DAF1AC1_Y272= mean(Unc5_sto);

subplot(3,3,3); hold on; box on; set(gca,'LineWidth',2.5); set(gca,'FontSize',15)
newcolors = [0 0 0
          0 0 1 
          1 0 0];
colororder(newcolors)
plot(Time, DAF1AC1_Y0,'LineWidth',2.2)
hold on
plot(Time, DAF1AC1_Y136,'LineWidth',2.0)
plot(Time, DAF1AC1_Y272,'LineWidth',1.9)
hold off
axis square
h_title=title('LIN-42 noise freq 5X'); set(h_title,'FontSize',15)
%x_label=xlabel('Time(hr)'); set(x_label,'FontSize',15);
%y_label=ylabel(' UNC-5 '); set(y_label,'FontSize',15);
ylim([0,250]);
x = [0 Time(end)];
y = [95 95];
line(x,y,'Color','black','LineStyle',':','LineWidth',1.5)
xlim([0,25]);


load('DAF1AC0_11_16_D95_Unc5_Lin42_0_b3_10xd_2025_01_08_b3gx5.mat')
DAF1AC0_Y0= Unc5_sto(1,:);

load('DAF1AC0_11_16_D95_Unc5_Lin42_5.000000e-01_b3_10xd_2025_01_08_b3gx5.mat')
DAF1AC0_Y136= Unc5_sto(1,:);

load('DAF1AC0_11_16_D95_Unc5_Lin42_1_b3_10xd_2025_01_08_b3gx5.mat')
DAF1AC0_Y272= Unc5_sto(1,:);


subplot(3,3,4); hold on; box on; set(gca,'LineWidth',2.5); set(gca,'FontSize',15)
newcolors = [0 0 0
          0 0 1 
          1 0 0];
colororder(newcolors)
plot(Time, DAF1AC0_Y0,'LineWidth',2.2)
hold on
plot(Time, DAF1AC0_Y136,'LineWidth',2.0)
plot(Time, DAF1AC0_Y272,'LineWidth',1.9)
hold off
axis square
 
%h_title=title('DAFnoAC_ det'); set(h_title,'FontSize',15)
%x_label=xlabel(' Time(hr) '); set(x_label,'FontSize',15);
y_label=ylabel(' UNC-5 '); set(y_label,'FontSize',15);
ylim([0,250]);
x = [0 Time(end)];
y = [95 95];
line(x,y,'Color','black','LineStyle',':','LineWidth',1.5)

load('DAF1AC0_11_16_D95_Unc5_Lin42_0_b3_100xs_2025_01_08_b3gx1.mat')
DAF1AC0_Y0= mean(Unc5_sto);

load('DAF1AC0_11_16_D95_Unc5_Lin42_5.000000e-01_b3_100xs_2025_01_08_b3gx1.mat')
DAF1AC0_Y136= mean(Unc5_sto);

load('DAF1AC0_11_16_D95_Unc5_Lin42_1_b3_100xs_2025_01_08_b3gx1.mat')
DAF1AC0_Y272= mean(Unc5_sto);


subplot(3,3,5); hold on; box on; set(gca,'LineWidth',2.5); set(gca,'FontSize',15)
newcolors = [0 0 0
          0 0 1 
          1 0 0];
colororder(newcolors)
plot(Time, DAF1AC0_Y0,'LineWidth',2.2)
hold on
plot(Time, DAF1AC0_Y136,'LineWidth',2.0)
plot(Time, DAF1AC0_Y272,'LineWidth',1.9)
hold off
axis square
 
%h_title=title('DAFnoAC _ b1X r1X'); set(h_title,'FontSize',15)
%x_label=xlabel('Time(hr)'); set(x_label,'FontSize',25);
%y_label=ylabel(' UNC-5 '); set(y_label,'FontSize',25);
ylim([0,250]);
x = [0 Time(end)];
y = [95 95];
line(x,y,'Color','black','LineStyle',':','LineWidth',1.5)
xlim([0,25]);


load('DAF1AC0_11_16_D95_Unc5_Lin42_0_b3_100xs_2025_01_08_b3gx5.mat')
DAF1AC0_Y0= mean(Unc5_sto);

load('DAF1AC0_11_16_D95_Unc5_Lin42_5.000000e-01_b3_100xs_2025_01_08_b3gx5.mat')
DAF1AC0_Y136= mean(Unc5_sto);

load('DAF1AC0_11_16_D95_Unc5_Lin42_1_b3_100xs_2025_01_09_b3gx5.mat')
DAF1AC0_Y272= mean(Unc5_sto);


subplot(3,3,6); hold on; box on; set(gca,'LineWidth',2.5); set(gca,'FontSize',15)
newcolors = [0 0 0
          0 0 1 
          1 0 0];
colororder(newcolors)
plot(Time, DAF1AC0_Y0,'LineWidth',2.2)
hold on
plot(Time, DAF1AC0_Y136,'LineWidth',2.0)
plot(Time, DAF1AC0_Y272,'LineWidth',1.9)
hold off
axis square

%h_title=title('DAFnoAC _ b5X r1X'); set(h_title,'FontSize',15)
%x_label=xlabel('Time(hr)'); set(x_label,'FontSize',25);
%y_label=ylabel(' UNC-5 '); set(y_label,'FontSize',25);
ylim([0,250]);
x = [0 Time(end)];
y = [95 95];
line(x,y,'Color','black','LineStyle',':','LineWidth',1.5)
xlim([0,25]);

