function Daf12CLogic = LogicDaf12C(Daf12i,Daf12Li)
       global h a_ON b
	   h=3;
	   eval(['para28'])		
		
        H_Daf12ACDaf12C  	  = (Daf12i)/(a_ON*b);  
        H_Daf12LACDaf12C  	  = (Daf12Li)/(a_ON*b); 
		
		Daf12CLogic        = H_Daf12ACDaf12C*H_Daf12LACDaf12C; 

end

