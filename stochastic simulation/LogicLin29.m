function Lin29Logic = LogicLin29(Lin42i,Blmp1i,geneXi)
       global h
	   h=3;
	   eval(['para'])		
		
		H_lin42RElin29    = (K_Lin42_re_Lin29^h)/(K_Lin42_re_Lin29^h +Lin42i^h);
		H_blmpRElin29     = (K_Blmp1_re_Lin29^h)/(K_Blmp1_re_Lin29^h +Blmp1i^h);
		H_geneXRElin29    = (K_geneX_re_Lin29^h)/(K_geneX_re_Lin29^h +geneXi^h);	
		
		Lin29Logic       = H_lin42RElin29 * H_blmpRElin29 * H_geneXRElin29;
        
        
end

