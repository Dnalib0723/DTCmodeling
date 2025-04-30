function Unc5Logic  = LogicUnc5(Daf12i,Lin29i,Blmp1i)
       global h  
	   h=3;
	   eval(['para28'])		
		
		H_dafACunc        = (Daf12i^h)/(K_Daf12_ac_Unc5^h +Daf12i^h);
        H_lin29ACunc      = (Lin29i^h)/(K_Lin29_ac_Unc5^h +Lin29i^h); 
        H_blmpREunc       = (K_Blmp1_re_Unc5^h)/(K_Blmp1_re_Unc5^h +Blmp1i^h);
		
        Unc5Logic        = ((H_lin29ACunc + H_dafACunc)-(H_lin29ACunc *H_dafACunc))*H_blmpREunc; 


end

