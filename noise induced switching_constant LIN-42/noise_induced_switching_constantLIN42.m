clc
clear all
tStart = tic;

global PARA_W
global intervals N dt t Time_End T3_End T2_End
global b bm bp a_ON a_OFF L42_scale
global b_Lin42 a_ON_Lin42 a_OFF_Lin42
global a_Daf12 a_Lin42 a_Dre1
global Z_Daf12 Z_Blmp1 Z_Unc5 Z_Lin42 Z_Dre1 Z_Lin29 Z_geneX 
global M Mdet
global Unc5_det Unc5_sto Lin42_sto
global Daf12 Blmp1 Unc5 Lin29 Lin42 Dre1 geneX blmp1 lin29 unc5 Daf12L Daf12C



PARA_W= grep_mod('parameters_final_testX_ethanfinalm18mod_11_16.txt')

Time_End = 25; %count with hour, total 20 hours in simulation

  
DT       = 2.5/60;                % the size of interval every 2.5 min
dt       = 2.5*60;                % dt for stochastic simulation in sec  
intervals= round(Time_End/DT);    % number of intervals(steps)
h        = 3;                     % Hill coefficient
N        = intervals+1;
t        = linspace(0,Time_End,N);% time line

%------- 2STD- 4 times --------------


b         = 3;% burst size of translation for LIN-42, SCF-DRE-1, DAF-12L, DAF-12R and X gene
bm = 3;% burst size of transcription for BLMP-1, LIN-29, and UNC-5
bp = 1;% burst size of translation for BLMP-1, LIN-29, and UNC-5


a_ON	  = 90.6667;% A fixed value to determine burst frequency of transcription and tranlation
a_OFF     = a_ON/5;%basal expressoin of SCF-DRE-1, DAF-12L

b_Lin42   = 3; % burst size of translation for LIN-42 
a_ON_Lin42= 90.6667; % A fixed value to determine burst frequency of LIN-42 tranlation
a_OFF_Lin42=a_ON_Lin42/5; %basal expressoin of LIN-42


%eval(['signalinput11_24_mod'])

%Gene switch (ON =1 , OFF =0)
Z_Daf12             = 1;            Z_Blmp1     = 1;             Z_Unc5      = 1;        
Z_Lin42             = 1;            Z_Dre1      = 1;             Z_Lin29     = 1;  
Z_geneX             = 1; 


% index for BLMP-1 activation of lin-29
Z_AC = 0;

% Stochastic trojectory for Input signal
Mdet     = 2;

% tuning a constant LIN-42 level for simulation 
L42_sc=1;

bin_size =1.0;

Det_threshold = 95;


for rep = 1:length(L42_sc)
L42_scale    = L42_sc(rep)

eval(['signalinput11_24_mod'])

eval(['stochasticsim_mod27det'])

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

%M        = 100;                 %Number of trials
M        = 10;                 %Number of trials
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


time = datestr(now, 'yyyy_mm_dd');
Filename_Lin42 = "DAF%dAC%d_11_16_D%d_Lin42_Lin42_%d_b%d_%dxs_%s_b3gx5.mat";
Filename_Unc5 = "DAF%dAC%d_11_16_D%d_Unc5_Lin42_%d_b%d_%dxs_%s_b3gx5.mat";
Filename_Turn = "DAF%dAC%d_11_16_D%d_Turn_Lin42_%d_b%d_%dxs_%s_b3gx5.mat";
Filename_mean = "DAF%dAC%d_11_16_D%d_mean_Lin42_%d_b%d_%dxs_%s_b3gx5.mat";
Filename_std = "DAF%dAC%d_11_16_D%d_std_Lin42_%d_b%d_%dxs_%s_b3gx5.mat";
Filename_t = "DAF%dAC%d_11_16_D%d_t_Lin42_%d_b%d_%dxs_%s_b3gx5.mat";

format bank;L42_scale;

filename_Lin42 = sprintf(Filename_Lin42,Z_Daf12,Z_AC,Det_threshold,L42_scale,b_Lin42,M,time);
filename_Unc5 = sprintf(Filename_Unc5,Z_Daf12,Z_AC,Det_threshold,L42_scale,b_Lin42,M,time);
filename_Turn = sprintf(Filename_Turn,Z_Daf12,Z_AC,Det_threshold,L42_scale,b_Lin42,M,time);
filename_mean = sprintf(Filename_mean,Z_Daf12,Z_AC,Det_threshold,L42_scale,b_Lin42,M,time);
filename_std = sprintf(Filename_std,Z_Daf12,Z_AC,Det_threshold,L42_scale,b_Lin42,M,time);
filename_t = sprintf(Filename_t,Z_Daf12,Z_AC,Det_threshold,L42_scale,b_Lin42,M,time);

save(filename_Lin42 ,'Lin42');
save(filename_Unc5 ,'Unc5_sto');
save(filename_Turn ,'Unc5_SS_Wt_dist_turn_Noturn');
save(filename_mean ,'mean_Unc5_expr_T_Distr_Wt');
save(filename_std ,'std_Unc5_expr_T_Distr_Wt');
save(filename_t ,'t');

end



tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));

