close all
clear

gamma_commom=0.002154; % common protein degradation rate
a_ON=90.67; % A fixed value to determine burst frequency of transcription and tranlation
b=3; % burst size of translation

% threshold value of gene regualtion function

KM_LIN42_Rp_DAF12cep = 100;
KM_LIN42_Rp_lin29 = 100;
KM_LIN42_Ac_blmp1 = 50;
KM_DRE1_Rp_BLMP1 = 130;
KM_LIN29rpblmp1 = 50;
KM_blmp1_pos = 105;
KM_BLMP1rplin29 = 85;
KM_BLMP1rpgeneX = 15;
KM_geneXrplin29 = 200;
KM_DAF12cox_Rp_blmp1 = 180;

n=3; % Hill coefficient of gene regulation functionHill coefficient of gene regulation function

    
%%steady state of upstream signals
DAF12lg=0*((gamma_commom*a_ON*b)/gamma_commom);
DRE1=0.2*((gamma_commom*a_ON*b)/gamma_commom);
LIN42=0.1*0.5*((gamma_commom*a_ON*b)/gamma_commom);


%% regulatory function by the upstream signals
neg_LIN42rplin29 = (KM_LIN42_Rp_lin29^n)/(LIN42^n+KM_LIN42_Rp_lin29^n);
pos_LIN42acblmp1 = LIN42^n/(LIN42^n+KM_LIN42_Ac_blmp1^n);

neg_LIN42rpDAF12cep = (KM_LIN42_Rp_DAF12cep^n)/(LIN42^n+KM_LIN42_Rp_DAF12cep^n);
DAF12cep=0.1*((gamma_commom*a_ON*b)/gamma_commom)+0.9*neg_LIN42rpDAF12cep*((gamma_commom*a_ON*b)/gamma_commom);
DAF12cox=(DAF12lg/(a_ON*b))*(DAF12cep/(a_ON*b))*((gamma_commom*a_ON*b)/gamma_commom);
neg_DAF12coxrpblmp1 = (KM_DAF12cox_Rp_blmp1^n)/(DAF12cox^n+KM_DAF12cox_Rp_blmp1^n);

neg_DRE1rpBLMP1 = (DRE1^n)/(DRE1^n+KM_DRE1_Rp_BLMP1^n);


%% set parameters

k_BLMP1_deg = gamma_commom;
k_BLMP1_prod = 5*k_BLMP1_deg;
k_BLMP1_DRE1 = 0.0020;
k_blmp1_prod = k_BLMP1_deg*a_ON*b*0.9;
k_blmp1_prod_base = k_BLMP1_deg*a_ON*b*0.1;
k_blmp1_deg = 5*k_BLMP1_deg;
k_BLMP1_prod_base = k_BLMP1_prod*0.1*(k_blmp1_prod/k_blmp1_deg);

Total = 300; % setting of the mac steady state

k_LIN29_deg = gamma_commom*0.25;
k_LIN29_prod = 5*k_LIN29_deg;
k_lin29_prod = k_LIN29_deg*a_ON*b*0.9;
k_lin29_prod_base= k_LIN29_deg*a_ON*b*0.1;
k_lin29_deg = 5*k_LIN29_deg;
k_LIN29_prod_base = k_LIN29_prod*0.1*(k_lin29_prod/k_lin29_deg);

%% tuning for the newly established BLMP-1 activation of lin-29(0=ON; 1=OFF)

neg_strength_BLMP1rpgeneX = 0;

%% find the lin29 roots at steady state
for ii = 0:30
    BLMP1(ii+1) = ii*10;
    neg_BLMP1rpgeneX = (neg_strength_BLMP1rpgeneX*BLMP1(ii+1)^n+KM_BLMP1rpgeneX^n)/(BLMP1(ii+1)^n+KM_BLMP1rpgeneX^n);
    geneX = 0.1*(gamma_commom*a_ON*b/gamma_commom)+0.9*neg_BLMP1rpgeneX*(gamma_commom*a_ON*b/gamma_commom);
    neg_geneXrplin29 = KM_geneXrplin29^n/(geneX^n+KM_geneXrplin29^n);
    neg_BLMP1rplin29 = KM_BLMP1rplin29^n/(BLMP1(ii+1)^n+KM_BLMP1rplin29^n);
    lin29_reg_func = neg_LIN42rplin29*neg_geneXrplin29*neg_BLMP1rplin29; 
    LIN29_roots(ii+1,:) = roots([(-k_LIN29_deg) k_LIN29_prod*((k_lin29_prod_base+k_lin29_prod*lin29_reg_func)/k_lin29_deg)]);

end

figure(1)
set(gca,'FontSize',30)
%idxx=10-iii;
%subplot(3,3,idxx)

%% get rid of imaginary roots
image_roots = abs(imag(LIN29_roots))>20;
LIN29_roots2 = LIN29_roots;
LIN29_roots2(image_roots) = nan;
LIN29_line = real(LIN29_roots2);
BLMP1_line = BLMP1;
plot(BLMP1_line, LIN29_line, 'b', 'LineWidth', 3)
ylim([0 Total])
ylim([-100 500])
hold on

%% find the blmp1 roots at steady state
for ii = 0:Total
    LIN29(ii+1) = ii;
    neg_LIN29rpblmp1 = (KM_LIN29rpblmp1^n)/(LIN29(ii+1)^n+KM_LIN29rpblmp1^n);
    blmp1_reg_func = neg_LIN29rpblmp1*pos_LIN42acblmp1*neg_DAF12coxrpblmp1; 
    BLMP1_roots(ii+1,:) = roots([(-k_BLMP1_deg-k_BLMP1_DRE1*neg_DRE1rpBLMP1) ((k_BLMP1_prod/k_blmp1_deg)*(k_blmp1_prod+k_blmp1_prod_base))  0 -1*(k_BLMP1_DRE1*neg_DRE1rpBLMP1+k_BLMP1_deg)*KM_blmp1_pos^n  KM_blmp1_pos^n*(k_BLMP1_prod/k_blmp1_deg)*(k_blmp1_prod_base+k_blmp1_prod*blmp1_reg_func)]);
end
%%
%% get rid of imaginary roots
image_roots = abs(imag(BLMP1_roots))>2;
BLMP1_roots2 = BLMP1_roots;
BLMP1_roots2(image_roots) = nan;
BLMP1_line2 = real(BLMP1_roots2);
BLMP1_line2_trans = transpose(BLMP1_line2);
BLMP1_line2_sort = sort(BLMP1_line2_trans);
BLMP1_line2 = transpose(BLMP1_line2_sort);
BLMP1_line2(BLMP1_line2<=0) = nan;
clear size
size=size(BLMP1_line2);
idx1 = 2;
idx2 = 4;
for i = 1: size(1)
    if isnan(BLMP1_line2(i, 3)) &  (BLMP1_line2(i, 2)>50)  
        BLMP1_line2(i, [idx1,idx2]) = BLMP1_line2(i, [idx2,idx1])
    end
end
LIN29_line2 = LIN29;
%%


plot(BLMP1_line2,LIN29_line2, 'r', 'LineWidth', 2)
axis square
xticks(0:50:Total)
yticks(0:50:Total)
set(gca,'FontSize',20)
%% get the vector for quiver plot, and save data for later plotting with python 
BLMP1_2 = 0:30:Total;
LIN29_2 = 0:30:Total;
[mBLMP1, mLIN29] = meshgrid(BLMP1_2,LIN29_2);

neg_BLMP1rpgeneX = (neg_strength_BLMP1rpgeneX*mBLMP1.^n+KM_BLMP1rpgeneX^n)./(mBLMP1.^n+KM_BLMP1rpgeneX^n);
geneX =0.1*(gamma_commom*a_ON*b/gamma_commom)+0.9*neg_BLMP1rpgeneX*(gamma_commom*a_ON*b/gamma_commom);
neg_geneXrplin29 = (KM_geneXrplin29^n)./(geneX.^n+KM_geneXrplin29^n);
neg_BLMP1rplin29 = (KM_BLMP1rplin29^n)./(mBLMP1.^n+KM_BLMP1rplin29^n);
lin29_reg_func = neg_LIN42rplin29*neg_geneXrplin29.*neg_BLMP1rplin29; 

dlin29 = k_LIN29_prod*((k_lin29_prod_base+k_lin29_prod*lin29_reg_func)/k_lin29_deg) - k_LIN29_deg*mLIN29;

neg_LIN29rpblmp1 = (KM_LIN29rpblmp1^n)./(mLIN29.^n+KM_LIN29rpblmp1^n);
blmp1_reg_func_base = pos_LIN42acblmp1*neg_DAF12coxrpblmp1*neg_LIN29rpblmp1; 
blmp1_reg_func= blmp1_reg_func_base + (mBLMP1.^n./(mBLMP1.^n+KM_blmp1_pos^n))- blmp1_reg_func_base.*(mBLMP1.^n./(mBLMP1.^n+KM_blmp1_pos^n))

dblmp1 = k_BLMP1_prod*((k_blmp1_prod_base+k_blmp1_prod*blmp1_reg_func)/k_blmp1_deg)- (k_BLMP1_deg+k_BLMP1_DRE1*neg_DRE1rpBLMP1)*mBLMP1;
% 

hq = quiver(mBLMP1,mLIN29,dblmp1,dlin29,0.5,'color',[0.8500 0.3250 0.0980]);
set(hq,'MaxheadSize',2)
set(hq,'autoscalefactor',2) ;
set(hq,'ShowArrowhead','on') 
hold on;

xlabel('BLMP-1','FontSize',20);
ylabel('LIN-29','FontSize',20)
xlim([0 Total])
xlim([0 Total])
ylim([0 Total])

% production rate of unc-5 mRNA and protein 

k_UNC5_deg = gamma_commom*0.25;
k_UNC5_prod = 5*k_UNC5_deg;
k_UNC5_prod_base = k_UNC5_prod*0.1;
k_unc5_prod = k_UNC5_deg*a_ON*b;
k_unc5_deg = 5*k_UNC5_deg;

KM_BLMP1_Rp_unc5 = 125;
KM_DAF12cox_Ac_unc5 = 130;
KM_LIN29_Ac_unc5 = 110;

neg_BLMP1_rp_unc5 = (KM_BLMP1_Rp_unc5^n)./(mBLMP1.^n+KM_BLMP1_Rp_unc5^n);

pos_DAF12cox_ac_unc5 = DAF12cox^n/(DAF12cox^n+KM_DAF12cox_Ac_unc5^n);
pos_LIN29_ac_unc5 = mLIN29.^n./(mLIN29.^n+KM_LIN29_Ac_unc5^n);
unc5_ac_func = pos_DAF12cox_ac_unc5 + pos_LIN29_ac_unc5 - pos_DAF12cox_ac_unc5*pos_LIN29_ac_unc5; 

UNC5 = (k_unc5_prod/k_unc5_deg)*(k_UNC5_prod*0.1 + 0.9*k_UNC5_prod*neg_BLMP1_rp_unc5.*unc5_ac_func)*(1/k_UNC5_deg);
[BLMP1_i,LIN29_i] = meshgrid(0:1:Total, 0:1:Total);

%// Interpolate the data and show the output
outData = interp2(mBLMP1, mLIN29, UNC5, BLMP1_i, LIN29_i, 'linear');

I=imagesc(outData);
I.AlphaData=0.4;
hold off;


%// Add colour bar
xlabel('BLMP-1','FontSize',20);
ylabel('LIN-29','FontSize',20);
%title(['\theta_{LIN-42}=' num2str(0.01*iii)])

a=colorbar;
set(a, 'ylim', [30 270])
set(a,'Visible','on');
c = colorbar();
drawnow;
alphaVal = 0.5;

cdata = c.Face.Texture.CData;
cdata(end,:) = uint8(alphaVal * cdata(end,:));
c.Face.Texture.ColorType = 'truecoloralpha';
c.Face.Texture.CData = cdata


%'FaceAlpha',.5
%set(a, 'FaceAlpha', 0.5)
label_unc5=ylabel(c,'UNC-5','FontSize',20,'Rotation',270);
label_unc5.Position(1) = 3.4;



 