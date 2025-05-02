 clc
clear all
tStart = tic;

global PARA_W
global intervals N dt t Time_End T3_End T2_End
global b bm bp a_ON a_OFF
global b_Lin42 a_ON_Lin42 a_OFF_Lin42
global a_Daf12 a_Lin42 a_Dre1
global Z_Daf12 Z_Blmp1 Z_Unc5 Z_Lin42 Z_Dre1 Z_Lin29 Z_geneX 
global M Mdet
global Unc5_det Unc5_sto Lin42_sto
global Daf12 Blmp1 Unc5 Lin29 Lin42 Dre1 geneX blmp1 lin29 unc5 Daf12L Daf12C

PARA_W= grep_mod('parameters_final_testX_ethanfinalm18mod_11_16.txt');


Time_End = 25; %count with hour, total 20 hours in simulation
T2_End = 7.5;
T3_End = 17;

  
DT       = 2.5/60;                % the size of interval every 2.5 min
dt       = 2.5*60;                % dt for stochastic simulation in sec  
intervals= round(Time_End/DT);    % number of intervals(steps)
h        = 3;                     % Hill coefficient
N        = intervals+1;
t        = linspace(0,Time_End,N);% time line

%------- 2STD- 4 times --------------


b         = 3; % burst size of translation for LIN-42, SCF-DRE-1, DAF-12L, DAF-12R and X gene

bm = 3;% burst size of transcription for BLMP-1, LIN-29, and UNC-5
bp = 1;% burst size of translation for BLMP-1, LIN-29, and UNC-5


a_ON	  = 90.6667; % A fixed value to determine burst frequency of transcription and tranlation
a_OFF     = a_ON/5; %basal expressoin of SCF-DRE-1, DAF-12L
b_Lin42   = 3; % burst size of translation for LIN-42 
a_ON_Lin42= 90.6667; % A fixed value to determine burst frequency of LIN-42 tranlation
a_OFF_Lin42=a_ON_Lin42/5; %basal expressoin of LIN-42

 
%% Signal Input  Panel


eval(['signalinput11_24_mod'])

%Gene switch (ON =1 , OFF =0)
Z_Daf12             = 1;            Z_Blmp1     = 1;             Z_Unc5      = 1;        
Z_Lin42             = 1;            Z_Dre1      = 1;             Z_Lin29     = 1;  
Z_geneX             = 1; 

Z_AC = 1;

% Stochastic trojectory for Input signal
Mdet     = 2;
eval(['stochasticsim_mod27det'])

bin_size =0.5;
Det_threshold = 95;


Unc5_expr_T_Distr_Wt = zeros(1,Mdet) +Time_End;Unc5_SS_Wt_dist_turn_Noturn = zeros(Mdet,N);

for j = 1:Mdet; DTC_Turned_orNot =0;
	for i=1:N; 
		if (Unc5_sto(j,i)>= (Det_threshold)); DTC_Turned_orNot = 1; if abs(Unc5_sto(j,i)-Det_threshold) > abs(Unc5_sto(j,i-1)-Det_threshold); Unc5_expr_T_Distr_Wt(j)= t(i-1); Unc5_SS_Wt_dist_turn_Noturn(j,i-1)=0;else Unc5_expr_T_Distr_Wt(j)	= t(i);	end;	end; 
		if DTC_Turned_orNot ~= 1; Unc5_SS_Wt_dist_turn_Noturn(j,i)=1; end; 
		if DTC_Turned_orNot == 1; break; end;
	end; 
end;

Unc5_sto_det=Unc5_sto(1,:);
mean_Unc5_expr_T_Distr_Wt_det = mean(Unc5_expr_T_Distr_Wt); std_Unc5_expr_T_Distr_Wt_det = std(Unc5_expr_T_Distr_Wt); 

bin_size 	= 1.0;
bin_edge = 1:bin_size:Time_End;
elements_U5_T_Distr_det = histc(Unc5_expr_T_Distr_Wt,bin_edge);

s=0;

 
M        = 10;                %Number of trials
eval(['stochasticsim_mod27sc'])


Unc5_expr_T_Distr_Wt = zeros(1,M) +Time_End;Unc5_SS_Wt_dist_turn_Noturn = zeros(M,N);

for j = 1:M; DTC_Turned_orNot =0;
	for i=1:N; 
		if (Unc5_sto(j,i)>= (Det_threshold)); DTC_Turned_orNot = 1; if abs(Unc5_sto(j,i)-Det_threshold) > abs(Unc5_sto(j,i-1)-Det_threshold); Unc5_expr_T_Distr_Wt(j)= t(i-1); Unc5_SS_Wt_dist_turn_Noturn(j,i-1)=0;else Unc5_expr_T_Distr_Wt(j)	= t(i);	end;	end; 
		if DTC_Turned_orNot ~= 1; Unc5_SS_Wt_dist_turn_Noturn(j,i)=1; end; 
		if DTC_Turned_orNot == 1; break; end;
	end; 
end;
mean_Unc5_sto = mean(Unc5_sto); std_Unc5_sto = std(Unc5_sto); 
mean_Unc5_expr_T_Distr_Wt = mean(Unc5_expr_T_Distr_Wt); std_Unc5_expr_T_Distr_Wt = std(Unc5_expr_T_Distr_Wt); 

bin_size 	= 1.0;
bin_edge = 1:bin_size:Time_End;
elements_U5_T_Distr = histc(Unc5_expr_T_Distr_Wt,bin_edge);




Unc5_chop = Unc5_SS_Wt_dist_turn_Noturn.*Unc5_sto;
tmax=t.*Unc5_SS_Wt_dist_turn_Noturn;
tplot=tmax;
tplot(tplot==0)=nan;
tplot(:,1)=0;
yplot=Unc5_chop;
Unc5_chop(Unc5_chop==0)=nan;

bin_size 	= 1.0;
bin_edge = 1:bin_size:Time_End;
elements_U5_T_Distr = histc(Unc5_expr_T_Distr_Wt,bin_edge);

figure(1); 
subplot(1,2,1); hold on; box on; set(gca,'LineWidth',2.5); set(gca,'FontSize',16);
axis square;
t_err = [t, fliplr(t)];
Unc5_sto_err = [mean_Unc5_sto + std_Unc5_sto, fliplr(mean_Unc5_sto - std_Unc5_sto)];


patch([t fliplr(t)], [mean_Unc5_sto - std_Unc5_sto, fliplr(mean_Unc5_sto + std_Unc5_sto)],[0.65 0.81 0.94],'FaceAlpha',.6,'EdgeColor','none'); hold on;

for l=1:1:10
    plot(tplot(l,:),yplot(l,:),'-g','LineWidth',1.0,'MarkerEdgeColor','g','MarkerSize',5); hold on;
end


plot(t,Unc5_sto_det,'-r','LineWidth',2.5,'MarkerEdgeColor','r','MarkerSize',5); hold on;
plot(Unc5_expr_T_Distr_Wt(1:10),Det_threshold,'o','MarkerSize',2,'color',[0 0.5 0]);hold off;
ylim([0 400]);
x = [0 Time_End];
y = [Det_threshold Det_threshold];
line(x,y,'Color','red','LineStyle','--','LineWidth',1.5);
%ylim([0 500])
h_title=title('WT'); set(h_title,'FontSize',14);
x_label=xlabel('Time(hr)'); set(x_label,'FontSize',14);
y_label=ylabel('UNC-5 protein'); set(y_label,'FontSize',14);
v1 = [0 0; 0 400; T2_End 400; T2_End 0];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'EdgeColor','none','FaceColor',[0.7 0.7 0.7],'FaceAlpha',.4);

v2 = [T2_End 0; T2_End 400; T3_End 400; T3_End 0];
f2 = [1 2 3 4];
patch('Faces',f2,'Vertices',v2,'EdgeColor','none','FaceColor',[0.4 0.4 0.4],'FaceAlpha',.4);

v3 = [T3_End 0; T3_End 400; Time_End 400; Time_End 0];
f3 = [1 2 3 4];
patch('Faces',f3,'Vertices',v3,'EdgeColor','none','FaceColor',[0.1 0.1 0.1],'FaceAlpha',.4);

subplot(1,2,2); hold on; box on; set(gca,'LineWidth',2.5); set(gca,'FontSize',16);
%axis('equal')
color_bar_det=bar(bin_edge,elements_U5_T_Distr_det*50); set(color_bar_det,'FaceColor',[0.9290, 0.6940, 0.1250]); hold on;
color_bar_sto=bar(bin_edge,elements_U5_T_Distr); set(color_bar_sto,'FaceColor',[0, 0.4470, 0.7410]);hold off;
axis square; h_title=title('WT'); set(h_title,'FontSize',14);
x_label=xlabel('Time(hr)'); set(x_label,'FontSize',14);
y_label=ylabel('fraction of turned worms'); set(y_label,'FontSize',14);
xlim([0,Time_End]);
ylim([0,100]);
v1 = [0 0; 0 100; T2_End 100; T2_End 0];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'EdgeColor','none','FaceColor',[0.7 0.7 0.7],'FaceAlpha',.4);

v2 = [T2_End 0; T2_End 100; T3_End 100; T3_End 0];
f2 = [1 2 3 4];
patch('Faces',f2,'Vertices',v2,'EdgeColor','none','FaceColor',[0.4 0.4 0.4],'FaceAlpha',.4);

v3 = [T3_End 0; T3_End 100; Time_End 100; Time_End 0];
f3 = [1 2 3 4];
patch('Faces',f3,'Vertices',v3,'EdgeColor','none','FaceColor',[0.1 0.1 0.1],'FaceAlpha',.4);

trial=15;

time = datestr(now, 'yyyy_mm_dd');
Filename_Lin42 = "DAF%dAC%d_%d_11_16_D%d_Lin42_Lin42_b%d_%dxs_%s_Lin42_gx1.mat";
Filename_Unc5 = "DAF%dAC%d_%d_11_16_D%d_Unc5_Lin42_b%d_%dxs_%s_Lin42_gx1.mat";
Filename_Turn = "DAF%dAC%d_%d_11_16_D%d_Turn_Lin42_b%d_%dxs_%s_Lin42_gx1.mat";
Filename_mean = "DAF%dAC%d_%d_11_16_D%d_mean_Lin42_b%d_%dxs_%s_Lin42_gx1.mat";
Filename_std = "DAF%dAC%d_%d_11_16_D%d_std_Lin42_b%d_%dxs_%s_Lin42_gx1.mat";
Filename_Turningtime = "DAF%dAC%d_%d_11_16_D%d_Turningtime_Lin42_b%d_%dxs_%s_Lin42_gx1.mat";
Filename_t = "DAF%dAC%d_%d_11_16_D%d_t_Lin42_b%d_%dxs_%s_Lin42_gx1.mat";

%format bank;L42_scale;

filename_Lin42 = sprintf(Filename_Lin42,Z_Daf12,Z_AC,trial,Det_threshold,b_Lin42,M,time);
filename_Unc5 = sprintf(Filename_Unc5,Z_Daf12,Z_AC,trial,Det_threshold,b_Lin42,M,time);
filename_Turn = sprintf(Filename_Turn,Z_Daf12,Z_AC,trial,Det_threshold,b_Lin42,M,time);
filename_mean = sprintf(Filename_mean,Z_Daf12,Z_AC,trial,Det_threshold,b_Lin42,M,time);
filename_std = sprintf(Filename_std,Z_Daf12,Z_AC,trial,Det_threshold,b_Lin42,M,time);
filename_Turningtime = sprintf(Filename_Turningtime,Z_Daf12,Z_AC,trial,Det_threshold,b_Lin42,M,time);
filename_t = sprintf(Filename_t,Z_Daf12,Z_AC,trial,Det_threshold,b_Lin42,M,time);

save(filename_Lin42 ,'Lin42_sto');
save(filename_Unc5 ,'Unc5_sto');
save(filename_Turn ,'Unc5_SS_Wt_dist_turn_Noturn');
save(filename_mean ,'mean_Unc5_expr_T_Distr_Wt');
save(filename_std ,'std_Unc5_expr_T_Distr_Wt');
save(filename_Turningtime ,'Unc5_expr_T_Distr_Wt');
save(filename_t ,'t');

%Det_threshold = y_det_Wt(turning_time_index);

tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));

