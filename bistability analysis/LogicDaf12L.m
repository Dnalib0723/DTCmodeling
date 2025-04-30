function Daf12LLogic = LogicDaf12L(Lin42i)
       global h
	   h=3;
	   eval(['para'])		
		
        H_Lin42REDaf12L  	  = (K_Lin42_re_Daf12^h)/(K_Lin42_re_Daf12^h +Lin42i^h);  
		
		Daf12LLogic        =  H_Lin42REDaf12L;
        %geneXLogic        = 1;  

end

