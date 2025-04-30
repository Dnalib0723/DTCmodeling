clc
clear all
tStart = tic;

global PARA_W
global intervals N dt t
global b bm bp a_ON a_OFF
global b_Lin42 a_ON_Lin42 a_OFF_Lin42
global a_Daf12 a_Lin42 a_Dre1
global Z_Daf12 Z_Blmp1 Z_Unc5 Z_Lin42 Z_Dre1 Z_Lin29 Z_geneX 
global M
global Unc5_det
global Lin42_k blmp1_k Blmp1_k lin29_k Lin29_k unc5_k Unc5_k blmp1_km Blmp1_km lin29_km Lin29_km unc5_km Unc5_km geneX_k geneX_km geneX_km  
global Daf12 Blmp1 Unc5 Lin29 Lin42 Dre1 geneX blmp1 lin29 unc5 Daf12L Daf12C

PARA_W= grep_mod('parameters_final_testX_ethanfinalm18mod_11_16.txt');

%Time_End = 40; %count with hour, total 40 hours in simulation
Time_End = 10; %count with hour, total 10 hours in simulation
DT       = 0.75/60;                % the size of interval every 0.75 min
dt       = 0.75*60;                % dt for stochastic simulation in sex  
intervals= round(Time_End/DT);    % number of intervals(steps)
h        = 3;                     % Hill coefficient
N        = intervals+1;
t        = linspace(0,Time_End,N);% time line

%------- 2STD- 4 times --------------

b         = 3;% burst size of translation for LIN-42, SCF-DRE-1, DAF-12L, DAF-12R and X gene
bm = 3;% burst size of transcription for BLMP-1, LIN-29, and UNC-5
bp = 1; % burst size of translation for BLMP-1, LIN-29, and UNC-5

a_ON	  = 90.6667; % A fixed value to determine burst frequency of transcription and tranlation
a_OFF     = a_ON/5;%basal expressoin of SCF-DRE-1, DAF-12L
b_Lin42   = 3; % burst size of translation for LIN-42 
a_ON_Lin42= 90.6667; % A fixed value to determine burst frequency of LIN-42 tranlation

a_OFF_Lin42=a_ON_Lin42/5;  %basal expressoin of LIN-42

%% Signal Input  Panel

eval(['signalinput18_mod']) %define the gene expression of SCF-DRE-1, DAF-12L and LIN-42

%Gene switch (ON =1 , OFF =0)
Z_Daf12             = 1;            Z_Blmp1     = 1;            Z_Unc5      = 1;        
Z_Lin42             = 1;            Z_Dre1      = 1;            Z_Lin29     = 1;  
Z_geneX             = 1; 

% Stochastic trojectory for Input signal
M        = 2;                 %Number of trials

%select circuit condition for analysis(mod30 = circuitI; mod30_noAC =circuit II)
%simulation with a high UNC-5 steady state level
eval(['deterministic_mod30'])
%eval(['deterministic_mod30_noAC'])

%Gene switch (ON =1 , OFF =0)
Z_Daf12             = 1;            Z_Blmp1     = 1;            Z_Unc5      = 1;        
Z_Lin42             = 1;            Z_Dre1      = 1;            Z_Lin29     = 1;  
Z_geneX             = 1; 

% Stochastic trojectory for Input signal
M        = 2;                 %Number of trials

%select circuit condition for analysis(mod30m = circuitI; mod30m_noAC =circuit II)
%simulation with a low UNC-5 steady state level
eval(['deterministic_mod30m'])
%eval(['deterministic_mod30m_noAC'])


handle=figure(1);

subplot(3,1,1); hold on; box on; set(gca,'LineWidth',2.5); set(gca,'FontSize',10)
axis square;
plot(Lin42_k(1:length(Lin42_k)),Blmp1_km,'b-o','MarkerSize',2.5);
hold on;
plot(Lin42_k(1:length(Lin42_k)),Blmp1_k,':o','MarkerSize',2,'color',[0.93 0.46 0.13]);
hold off;

xlim([0 300]);
ylim([0 300]);
xlabel('LIN-42');
ylabel('BLMP-1');

subplot(3,1,2);  hold on; box on; set(gca,'LineWidth',2.5); set(gca,'FontSize',10)
axis square;

plot(Lin42_k(1:length(Lin42_k)),Lin29_km,'b-o','MarkerSize',2.5); 
hold on;
plot(Lin42_k(1:length(Lin42_k)),Lin29_k,':o','MarkerSize',2,'color',[0.93 0.46 0.13]);
hold off;

xlim([0 300]);
ylim([0 300]);
xlabel('LIN-42');
ylabel('LIN-29');

subplot(3,1,3);  hold on; box on; set(gca,'LineWidth',2.5); set(gca,'FontSize',10)
axis square;

plot(Lin42_k(1:length(Lin42_k)),Unc5_km,'b-o','MarkerSize',2.5); 
hold on;
plot(Lin42_k(1:length(Lin42_k)),Unc5_k,':o','MarkerSize',2,'color',[0.93 0.46 0.13]);
hold off;

xlim([0 300]);
ylim([0 300]);
xlabel('LIN-42');
ylabel('UNC-5');

%saving file for marking with for marking circuit
%circuitI: Z_AC=1; circuit II:  Z_AC=0
Z_AC = 1;

% DRE-1 level is set at 54 protein number/per cell
Z_Dre1 = 54;

%save the simulation result
time = datestr(now, 'yyyy_mm_dd');
Filename = "DRE%dDAF%dAC%d_%s";
Filename_Lin42_k = "DRE%dDAF%dAC%d_m18mod_11_16_Lin42_k_%s.mat";
Filename_Unc5_k = "DRE%dDAF%dAC%d_m18mod_11_16_Unc5_k_%s.mat";
Filename_Unc5_km = "DRE%dDAF%dAC%d_m18mod_11_16_Unc5_km_%s.mat";
Filename_Blmp1_k = "DRE%dDAF%dAC%d_m18mod_11_16_Blmp1_k_%s.mat";
Filename_Blmp1_km = "DRE%dDAF%dAC%d_m18mod_11_16_Blmp1_km_%s.mat";
Filename_Lin29_k = "DRE%dDAF%dAC%d_m18mod_11_16_Lin29_k_%s.mat";
Filename_Lin29_km = "DRE%dDAF%dAC%d_m18mod_11_16_Lin29_km_%s.mat";

filename = sprintf(Filename,Z_Dre1,Z_Daf12,Z_AC,time);
filename_Lin42_k = sprintf(Filename_Unc5_k,Z_Dre1,Z_Daf12,Z_AC,time);
filename_Unc5_k = sprintf(Filename_Unc5_k,Z_Dre1,Z_Daf12,Z_AC,time);
filename_Unc5_km = sprintf(Filename_Unc5_km,Z_Dre1,Z_Daf12,Z_AC,time);
filename_Blmp1_k = sprintf(Filename_Blmp1_k,Z_Dre1,Z_Daf12,Z_AC,time);
filename_Blmp1_km = sprintf(Filename_Blmp1_km,Z_Dre1,Z_Daf12,Z_AC,time);
filename_Lin29_k = sprintf(Filename_Lin29_k,Z_Dre1,Z_Daf12,Z_AC,time);
filename_Lin29_km = sprintf(Filename_Lin29_km,Z_Dre1,Z_Daf12,Z_AC,time);

saveas(handle,filename,'fig');
save(filename_Lin42_k ,'Lin42_k');
save(filename_Unc5_k ,'Unc5_k');
save(filename_Unc5_km ,'Unc5_km');
save(filename_Blmp1_k ,'Blmp1_k');
save(filename_Blmp1_km ,'Blmp1_km');
save(filename_Lin29_k ,'Lin29_k');
save(filename_Lin29_km ,'Lin29_km');

tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));