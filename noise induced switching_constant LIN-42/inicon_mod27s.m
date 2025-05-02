function inicon_mod27s()

eval(['para'])

global b bm bp a_ON a_OFF N b_Lin42
global a_Daf12 a_Lin42 a_Dre1
global Z_Daf12 Z_Blmp1 Z_Lin42 Z_Dre1
global basal_Lin29 basal_geneX basal_Unc5 basal_Blmp1 basal_Daf12C basal_Daf12L
global Daf12 Blmp1 Unc5 Lin29 Lin42 Dre1 geneX blmp1 lin29 unc5 Daf12L Daf12C

Daf12(1)  = a_Daf12(1)*b*Z_Daf12    ;
Daf12L(1) = a_ON*b*basal_Daf12L;
Daf12C(1) = a_ON*b*basal_Daf12C;

Unc5(1)	 = a_ON*b*basal_Unc5   ;
Lin29(1) = a_ON*b*basal_Lin29	  ;
Lin42(1) = a_Lin42(N)*b_Lin42*Z_Lin42  ;
Dre1(1)	 = a_Dre1(1)*b*Z_Dre1;
geneX(1) = a_ON*b;


unc5(1)	 = a_ON*basal_Unc5*(1/5)*bm;
lin29(1) = a_ON*basal_Lin29*(1/5)*bm;
blmp1(1) = a_ON*bm*(1/5);


LogBlmp1  = LogicBlmp1(Daf12(1), Lin42(1), Lin29(1), Blmp1(1), Dre1(1));
dre1ACblmp    = H_dre1ACblmp(Dre1(1));
Blmp1(1)    = a_ON*b;


end
