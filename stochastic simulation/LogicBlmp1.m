function Blmp1Logic = LogicBlmp1(Daf12i,Lin42i,Lin29i,Blmp1i,Dre1i)
       global h SetNum 
	   h=3;
	   eval(['para'])		
		
		H_dafREblmp     = (K_Daf12_re_Blmp1^h)/(K_Daf12_re_Blmp1^h +Daf12i^h);
        H_lin42ACblmp   = (Lin42i^h)/(K_Lin42_ac_Blmp1^h +Lin42i^h);
        H_lin29REblmp   = (K_Lin29_re_Blmp1^h)/(K_Lin29_re_Blmp1^h +Lin29i^h);
        H_blmpACblmp    = (Blmp1i^h)/(K_Blmp1_ac_Blmp1^h +Blmp1i^h);
        H_dre1ACblmp    = (Dre1i^h)/(K_Dre1_re_Blmp1^h +Dre1i^h);
        H_dafACblmp     = (Daf12i^h)/(K_Daf12_re_Blmp1^h +Daf12i^h);   
        % pFB- AND gate
        Blmp1Logic     = Blmp1_logic(H_dafREblmp,H_lin42ACblmp,H_lin29REblmp,H_blmpACblmp,H_dre1ACblmp,H_dafACblmp,SetNum);


end

