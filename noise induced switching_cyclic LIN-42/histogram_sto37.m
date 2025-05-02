clc; clear all;

Threshold= 10;
%Threshold= 40;
%Threshold= 100;

load('DafAC2_18119_D95_t_100x_1220.mat')%DafAC_623_D110_t_100xs_0406.mat
load('DafAC2_18119_D95_Lin42_100x_1220.mat')%DafAC_623_D110_Lin42_100xs_0406.mat

t_size=size(t);
percent_complete=zeros(1,t_size(2));
percent_DafnoAC=zeros(1,t_size(2));
percent_noDafAC=zeros(1,t_size(2));
percent_noDafnoAC=zeros(1,t_size(2));


for n=2: t_size(2)
    
load('DafAC2_18119_D95_Unc5_100x_1220.mat')%DafAC_623_D110_Unc5_100xs_0406.mat

%time=280;
%time=320;
time=n;
%% 
%Lin42_sto=Lin42_sto(:,time);
Unc5_sto=Unc5_sto(:,time);
%Lin42.edges = [0:25:500];
%Unc5.edges = [0:10:300];
Unc5.edges = [0:10:350];
%subplot(5,1,1)
%hLin42 = histogram(Lin42_sto,Lin42.edges,'Normalization','probability');
%ylabel('LIN-42 %')
%ylim([0 0.5]);
%subplot(5,1,2)
hUnc5 = histogram(Unc5_sto,Unc5.edges,'Normalization','probability');
%ylabel('UNC-5 %')
%ylim([0 1.0]);
%hUnc5.FaceColor=[0 0.5 0.5];
hUnc5_cut=hUnc5.Values(Threshold:end);
percent_complete(n)=sum(hUnc5_cut);
percent_complete(n)=(percent_complete(n)>=percent_complete(n-1)).*percent_complete(n)+ ~(percent_complete(n)>=percent_complete(n-1)).*percent_complete(n-1);

%load('noDafAC_Lin42_100x_0225.mat')
load('DafnoAC3_18119_D95_Unc5_100x_1220.mat')%DafnoAC_623_D110_Unc5_100xs_0406.mat
%Lin42_sto=Lin42_sto(:,time);
Unc5_sto=Unc5_sto(:,time);
%Lin42.edges = [0:25:500];
%Unc5.edges = [0:10:300];
Unc5.edges = [0:10:350];
%subplot(2,2,2)
%hLin42 = histogram(Lin42_sto,Lin42.edges,'Normalization','probability');
%ylim([0 0.5]);
%subplot(5,1,3)
hUnc5 = histogram(Unc5_sto,Unc5.edges,'Normalization','probability');
%ylabel('UNC-5 %')
%ylim([0 1.0]);
%hUnc5.FaceColor=[0.5 0 0.5];
hUnc5_cut=hUnc5.Values(Threshold:end);
percent_DafnoAC(n)=sum(hUnc5_cut);
percent_DafnoAC(n)=(percent_DafnoAC(n)>=percent_DafnoAC(n-1)).*percent_DafnoAC(n)+ ~(percent_DafnoAC(n)>=percent_DafnoAC(n-1)).*percent_DafnoAC(n-1);


load('noDafAC2_18119_D95_Unc5_100x_1220.mat')%noDafAC_623_D110_Unc5_100xs_0406.mat
%Lin42_sto=Lin42_sto(:,time);
Unc5_sto=Unc5_sto(:,time);
%Lin42.edges = [0:25:500];
%Unc5.edges = [0:10:300];
Unc5.edges = [0:10:350];
%subplot(2,2,2)
%hLin42 = histogram(Lin42_sto,Lin42.edges,'Normalization','probability');
%ylim([0 0.5]);
%subplot(5,1,4)
hUnc5 = histogram(Unc5_sto,Unc5.edges,'Normalization','probability');
%ylim([0 1.0]);
%ylabel('UNC-5 %')
%xlabel('Protein level')
%hUnc5.FaceColor=[0.8 0.8 0];
hUnc5_cut=hUnc5.Values(Threshold:end);
percent_noDafAC(n)=sum(hUnc5_cut);
percent_noDafAC(n)=(percent_noDafAC(n)>=percent_noDafAC(n-1)).*percent_noDafAC(n)+ ~(percent_noDafAC(n)>=percent_noDafAC(n-1)).*percent_noDafAC(n-1);


load('noDafnoAC_18119_D95_Unc5_100x_1220.mat')%noDafnoAC_623_D110_Unc5_100xs_0406.mat
%Lin42_sto=Lin42_sto(:,time);
Unc5_sto=Unc5_sto(:,time);
%Lin42.edges = [0:25:500];
%Unc5.edges = [0:10:300];
Unc5.edges = [0:10:350];
%subplot(2,2,2)
%hLin42 = histogram(Lin42_sto,Lin42.edges,'Normalization','probability');
%ylim([0 0.5]);
%subplot(5,1,5)
hUnc5 = histogram(Unc5_sto,Unc5.edges,'Normalization','probability');
%ylim([0 1.0]);
%ylabel('UNC-5 %')
%xlabel('Protein level')
%hUnc5.FaceColor=[0.6 0.6 0.6];
hUnc5_cut=hUnc5.Values(Threshold:end);
percent_noDafnoAC(n)=sum(hUnc5_cut);
percent_noDafnoAC(n)=(percent_noDafnoAC(n)>=percent_noDafnoAC(n-1)).*percent_noDafnoAC(n)+ ~(percent_noDafnoAC(n)>=percent_noDafnoAC(n-1)).*percent_noDafnoAC(n-1);

    
    
n
end

figure(1)
plot(t, percent_complete,'LineWidth',1.5, 'color', [0 0.5 0.5]); xlim([10, 20]); ylim([0 1]); hold on;
plot(t, percent_DafnoAC,'LineWidth',1.5, 'color', [0.5 0 0.5]); xlim([10, 20]); ylim([0 1]); hold on;
plot(t, percent_noDafAC,'LineWidth',1.5, 'color', [0.9100 0.4100 0.1700]); xlim([10, 20]); ylim([0 1]); hold on;
plot(t, percent_noDafnoAC,'LineWidth',1.5, 'color', [0 0 0]); xlim([10, 20]); ylim([0 1]); hold off;
xlabel('Time(hr)'); ylabel('% of turned DTC ');


v2 = [10 0; 10 1; 16 1; 16 0];
f2 = [1 2 3 4];
patch('Faces',f2,'Vertices',v2,'EdgeColor','none','FaceColor',[0.4 0.4 0.4],'FaceAlpha',.4);

v3 = [16 0; 16 1; 20 1; 20 0];
f3 = [1 2 3 4];
patch('Faces',f3,'Vertices',v3,'EdgeColor','none','FaceColor',[0.1 0.1 0.1],'FaceAlpha',.4);


percent_complete2=percent_complete(2:end);
percent_complete1=percent_complete(1:n-1);
percent_complete3=percent_complete2-percent_complete1;

percent_DafnoAC2=percent_DafnoAC(2:end);
percent_DafnoAC1=percent_DafnoAC(1:n-1);
percent_DafnoAC3=percent_DafnoAC2-percent_DafnoAC1;

percent_noDafAC2=percent_noDafAC(2:end);
percent_noDafAC1=percent_noDafAC(1:n-1);
percent_noDafAC3=percent_noDafAC2-percent_noDafAC1;

percent_noDafnoAC2=percent_noDafnoAC(2:end);
percent_noDafnoAC1=percent_noDafnoAC(1:n-1);
percent_noDafnoAC3=percent_noDafnoAC2-percent_noDafnoAC1;

figure(2);
%bar(percent_complete3,'FaceColor',[0 0.5 0.5]);
plot(t(2:end),percent_complete3,'LineWidth',1.5, 'color', [0 0.5 0.5]);
ylim([0,0.2]);
xlim([10, 20]);
hold on
%bar(percent_DafnoAC3,'FaceColor',[0.5 0 0.5]);
plot(t(2:end),percent_DafnoAC3, 'LineWidth',1.5, 'color', [0.5 0 0.5]);
ylim([0,0.2]);
xlim([10, 20]);
hold on
%bar(percent_noDafAC3,'FaceColor',[0.8 0.8 0]);
plot(t(2:end),percent_noDafAC3, 'LineWidth',1.5, 'color', [0.9100 0.4100 0.1700]);
ylim([0,0.2]);
xlim([10, 20]);
hold on
%bar(percent_noDafnoAC3,'FaceColor',[0.6 0.6 0.6]);
plot(t(2:end),percent_noDafnoAC3,'LineWidth',1.5, 'color', [0 0 0]);
ylim([0,0.2]);
xlim([10, 20]);

v2 = [10 0; 10 0.2; 16 0.2; 16 0];
f2 = [1 2 3 4];
patch('Faces',f2,'Vertices',v2,'EdgeColor','none','FaceColor',[0.4 0.4 0.4],'FaceAlpha',.4);

v3 = [16 0; 16 0.2; 20 0.2; 20 0];
f3 = [1 2 3 4];
patch('Faces',f3,'Vertices',v3,'EdgeColor','none','FaceColor',[0.1 0.1 0.1],'FaceAlpha',.4);

xlabel('Time(hr)'); ylabel('% DTC with dorsal turn');
