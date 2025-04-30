function geneXLogic = LogicgeneX_noAC(Blmp1i)
       global h
	   h=3;
	   eval(['para28'])		
		
        H_blmpREgeneX  	  = (K_Blmp1_re_geneX^h)/(K_Blmp1_re_geneX^h +Blmp1i^h);  
		
		%geneXLogic        = H_blmpREgeneX;
        geneXLogic        = 1;  

end

