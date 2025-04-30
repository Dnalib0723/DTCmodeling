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

Time_End = 20; %count with hour, total 20 hours in simulation
T2_End = 7.5; %end of L2 period
T3_End = 17; %end of L3 period
 
DT       = 2.5/60;                % the size of interval every 2.5 min
dt       = 2.5*60;                % dt for stochastic simulation in sec  
intervals= round(Time_End/DT);    % number of intervals(steps)
h        = 3;                     % Hill coefficient
N        = intervals+1;
t        = linspace(0,Time_End,N);% time line

%------- 2STD- 4 times --------------


b = 3; % burst size of translation for LIN-42, SCF-DRE-1, DAF-12L, DAF-12R and X gene
bm = 3; % burst size of transcription for BLMP-1, LIN-29, and UNC-5
bp = 1; % burst size of translation for BLMP-1, LIN-29, and UNC-5

a_ON	  = 90.6667; % A fixed value to determine burst frequency of transcription and tranlation
a_OFF     = a_ON/5; %basal expressoin of SCF-DRE-1, DAF-12L

b_Lin42   = 3; % burst size of translation for LIN-42 
a_ON_Lin42= 90.6667; % A fixed value to determine burst frequency of LIN-42 tranlation
a_OFF_Lin42=a_ON_Lin42/5; %basal expressoin of LIN-42

eval(['signalinput11_24_mod'])  %define the gene expression of SCF-DRE-1, DAF-12L and LIN-42

%Gene switch (ON =1 , OFF =0)
Z_Daf12             = 1;            Z_Blmp1     = 1;             Z_Unc5      = 1;        
Z_Lin42             = 1;            Z_Dre1      = 1;             Z_Lin29     = 1;  
Z_geneX             = 1; 

% Stochastic trojectory for Input signal
Mdet     = 2;
eval(['stochasticsim_mod27det']) % detemonsitic simulation 

bin_size =0.5; % bin size for distribution of turning time
Det_threshold = 95; % UNC-5 level require for DTC turning


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

M        = 10;      %Number of trials
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

Unc5_chop = Unc5_SS_Wt_dist_turn_Noturn.*Unc5_sto;
Unc5_chop(Unc5_chop==0)=nan;

Logic_Unc5_chop=isnan(Unc5_chop);
Logic_Unc5_chop_T=Logic_Unc5_chop';
SUM_chop=length(Unc5_chop)-sum(Logic_Unc5_chop_T)+1;
for itt=1:M
Unc5_chop(itt,SUM_chop(itt))= Det_threshold;
end
if length(Unc5_chop)~=intervals+1;
Unc5_chop(:,length(Unc5_chop))=[];
end

SUM_Noturn=sum(Unc5_SS_Wt_dist_turn_Noturn')+1;
for ittt=1:M
Unc5_SS_Wt_dist_turn_Noturn(ittt, SUM_Noturn(ittt))= 1;
end
if length(Unc5_SS_Wt_dist_turn_Noturn)~=intervals+1;
Unc5_SS_Wt_dist_turn_Noturn(:,length(Unc5_SS_Wt_dist_turn_Noturn))=[];
end

tmax=t.*Unc5_SS_Wt_dist_turn_Noturn;
tplot=tmax;
tplot(tplot==0)=nan;
tplot(:,1)=0;
yplot=Unc5_chop;

bin_size 	= 1.0;
bin_edge = 1:bin_size:Time_End;
elements_U5_T_Distr = histc(Unc5_expr_T_Distr_Wt,bin_edge);

handle=figure(1); 
subplot(1,2,1); hold on; box on; set(gca,'LineWidth',2.5); set(gca,'FontSize',16)
axis square;
t_err = [t, fliplr(t)];
Unc5_sto_err = [mean_Unc5_sto + std_Unc5_sto, fliplr(mean_Unc5_sto - std_Unc5_sto)];


patch([t fliplr(t)], [mean_Unc5_sto - std_Unc5_sto, fliplr(mean_Unc5_sto + std_Unc5_sto)],[0.63 0.80 0.19],'FaceAlpha',.6,'EdgeColor','none'); hold on;

for l=1:1:10
    plot(tplot(l,:),yplot(l,:),'-g','LineWidth',1.0,'color',[0.13, 0.54, 0.13],'MarkerSize',5); hold on;
end

title_name= 'WT';
plot(t,Unc5_sto_det,'-r','LineWidth',2.5,'color',[1 0.55 0],'MarkerSize',5); hold on;
plot(Unc5_expr_T_Distr_Wt(1:10),Det_threshold,'o','MarkerSize',5,'color',[1 0 0]);hold off;
ylim([0 300]);
x = [0 Time_End];
y = [Det_threshold Det_threshold];
line(x,y,'Color','black','LineStyle','--','LineWidth',1.5)
%ylim([0 500])
h_title=title(title_name); set(h_title,'FontSize',14);
x_label=xlabel('Time(hr)'); set(x_label,'FontSize',14);
y_label=ylabel('UNC-5 protein'); set(y_label,'FontSize',14);
v1 = [0 0; 0 300; T2_End 300; T2_End 0];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'EdgeColor','none','FaceColor',[0.7 0.7 0.7],'FaceAlpha',.4);

v2 = [T2_End 0; T2_End 300; T3_End 300; T3_End 0];
f2 = [1 2 3 4];
patch('Faces',f2,'Vertices',v2,'EdgeColor','none','FaceColor',[0.4 0.4 0.4],'FaceAlpha',.4);

v3 = [T3_End 0; T3_End 300; Time_End 300; Time_End 0];
f3 = [1 2 3 4];
patch('Faces',f3,'Vertices',v3,'EdgeColor','none','FaceColor',[0.1 0.1 0.1],'FaceAlpha',.4);

subplot(1,2,2); hold on; box on; set(gca,'LineWidth',2.5); set(gca,'FontSize',16)
color_bar_det=bar(bin_edge,elements_U5_T_Distr_det*50); set(color_bar_det,'FaceColor',[1 0.55 0]); hold on;
color_bar_sto=bar(bin_edge,elements_U5_T_Distr); set(color_bar_sto,'FaceColor',[0.13, 0.54, 0.13]);hold off;
axis square; h_title=title(title_name); set(h_title,'FontSize',14);
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


tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));



