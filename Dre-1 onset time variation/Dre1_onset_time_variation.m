clc
clear all
tStart = tic;

global PARA_W D_scale 
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
%T2_End = 7;
%T3_End = 16;
 
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

a_ON	  = 90.6667; % A fixed value to determine burst frequency of transcription and tranlation
a_OFF     = a_ON/5;%basal expressoin of SCF-DRE-1, DAF-12L
b_Lin42   = 3;% burst size of translation for LIN-42 
a_ON_Lin42= 90.6667; % A fixed value to determine burst frequency of LIN-42 tranlation
a_OFF_Lin42=a_ON_Lin42/5; %basal expressoin of LIN-42

%Gene switch (ON =1 , OFF =0)
Z_Daf12             = 1;            Z_Blmp1     = 1;             Z_Unc5      = 1;        
Z_Lin42             = 1;            Z_Dre1      = 1;             Z_Lin29     = 1;  
Z_geneX             = 1; 

% Stochastic trojectory for Input signal
Mdet     = 2;

%parameter for onsettime variation
d_sc=-25:10:25;

output=zeros(1,length(d_sc));
output_s=zeros(4,length(d_sc));

for rep = 1:length(d_sc)
D_scale    = d_sc(rep);
eval(['signalinput11_31_mod_scale'])

eval(['stochasticsim_mod29det'])

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

output(rep)=mean_Unc5_expr_T_Distr_Wt_det(1);

end
output_s(1,:)=output;
Dre1_ON_time=t(251:10:301);

figure(1)
plot(Dre1_ON_time,output,'LineWidth',1.5, 'color', [0 0.5 0.5]); hold on;

for rep = 1:length(d_sc)

    D_scale    = d_sc(rep);
eval(['signalinput11_31_mod_scale'])

eval(['stochasticsim_mod29det_noAC'])

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

output(rep)=mean_Unc5_expr_T_Distr_Wt_det(1);

end
output_s(2,:)=output;
figure(1)
plot(Dre1_ON_time,output,'LineWidth',1.5, 'color', [0.5 0 0.5]); hold on;

%Gene switch (ON =1 , OFF =0)
Z_Daf12             = 0;            Z_Blmp1     = 1;             Z_Unc5      = 1;        
Z_Lin42             = 1;            Z_Dre1      = 1;             Z_Lin29     = 1;  
Z_geneX             = 1; 

% Stochastic trojectory for Input signal
Mdet     = 2;


for rep = 1:length(d_sc)
    
D_scale    = d_sc(rep);
eval(['signalinput11_31_mod_scale'])

eval(['stochasticsim_mod29det'])

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

output(rep)=mean_Unc5_expr_T_Distr_Wt_det(1);

end

output_s(3,:)=output;
figure(1)
plot(Dre1_ON_time,output,'LineWidth',1.5, 'color', [0.9100 0.4100 0.1700]); hold on;

%Gene switch (ON =1 , OFF =0)
Z_Daf12             = 0;            Z_Blmp1     = 1;             Z_Unc5      = 1;        
Z_Lin42             = 1;            Z_Dre1      = 1;             Z_Lin29     = 1;  
Z_geneX             = 1; 

% Stochastic trojectory for Input signal
Mdet     = 2;


for rep = 1:length(d_sc)
    
D_scale    = d_sc(rep);
eval(['signalinput11_31_mod_scale'])

eval(['stochasticsim_mod29det_noAC'])

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

output(rep)=mean_Unc5_expr_T_Distr_Wt_det(1);

end

output_s(4,:)=output;
figure(1)
plot(Dre1_ON_time,output,'LineWidth',1.5, 'color', [0 0 0]); hold off;
legend('complete','noAC','noDAF-12','noDAf-12noAC')

x_label=xlabel('DRE-1 onset time(hr)'); set(x_label,'FontSize',14);
y_label=ylabel('DTC turning time(hr)'); set(y_label,'FontSize',14);

tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));
