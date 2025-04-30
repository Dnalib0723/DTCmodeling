global PARA_W scale

SetNum                = PARA_W(1);
PheTypes              = PARA_W(2);
Integr_Time           = PARA_W(3);
stdNum                = PARA_W(4);
DTC_Turning_Time      = PARA_W(5);
gama_common           = PARA_W(6);
gama_Blmp1            = PARA_W(7);

%gama_Blmp1_Dre1       = PARA_W(8)*scale;
gama_Blmp1_Dre1       = PARA_W(8);

K_Daf12_ac_Unc5       = PARA_W(9);
K_Daf12_re_Blmp1      = PARA_W(10);
K_Lin42_re_Lin29      = PARA_W(11);
K_Lin42_ac_Blmp1      = PARA_W(12);

K_Dre1_re_Blmp1       = PARA_W(13)*scale;
%K_Dre1_re_Blmp1       = PARA_W(13);

K_Lin29_re_Blmp1      = PARA_W(14);
K_Lin29_ac_Unc5       = PARA_W(15);
K_Blmp1_ac_Blmp1      = PARA_W(16); 
K_Blmp1_re_Lin29      = PARA_W(17);
K_Blmp1_re_Unc5       = PARA_W(18);
K_Blmp1_re_geneX      = PARA_W(19);
K_geneX_re_Lin29      = PARA_W(20);

K_Lin42_re_Daf12      =  100 ;
K_Daf12_ac_Daf12C     =  130  ; 
K_Daf12L_ac_Daf12C    =  100  ;