function Daf12CLogic = LogicDaf12C(Daf12i,Daf12Li)
       global h a_ON b
	   h=3;
	   eval(['para'])		
		
        H_Daf12ACDaf12C  	  = (Daf12i)/(a_ON*b);  
        H_Daf12LACDaf12C  	  = (Daf12Li)/(a_ON*b); 
        %H_Daf12ACDaf12C  	  = (Daf12i^h)/(K_Daf12_ac_Daf12C^h +Daf12i^h);  
        %H_Daf12LACDaf12C  	  = (Daf12Li^h)/(K_Daf12L_ac_Daf12C^h +Daf12Li^h); 
		
		Daf12CLogic        = H_Daf12ACDaf12C*H_Daf12LACDaf12C; 

end

