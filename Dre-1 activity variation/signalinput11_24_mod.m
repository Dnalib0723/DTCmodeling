function signalinput11_24_mod()

global intervals N t Time_End T3_End T2_End
global a_ON a_OFF 
global L2_ON L2_OFF L3_ON L3_OFF L4_ON a_ON_Lin42
global a_Daf12 a_Lin42 a_Dre1


% *** Daf 12***
a_Daf12        = zeros(1,N);
for i = 1:((intervals/Time_End*4.0)+1) 
    a_Daf12(i) = a_OFF; 
end
for i = ((intervals/Time_End*4.0)+2):N 
    a_Daf12(i) = a_ON; 
end

% *** Lin-42***
a_Lin42        = zeros(1,N);
			

for i = 1:((intervals/Time_End*7)+1)
    a_Lin42(i) = a_ON_Lin42*sin(2*pi/15.*(t(i)+0.5));	
end
     
	Lin42_min = a_Lin42(i);
	
for i = ((intervals/Time_End*6.5)+2):((intervals/Time_End*8)+1)
    a_Lin42(i) = a_ON_Lin42*1/5;
end	

for i = ((intervals/Time_End*8)+2):((intervals/Time_End*15)+1)
    a_Lin42(i) = a_ON_Lin42*sin(2*pi/15.*(t(i)-7.5));
end
     
	Lin42_min = a_Lin42(i);
	
for i = ((intervals/Time_End*14.5)+2):N
    a_Lin42(i) = a_ON_Lin42*1/5;
end	



% *** Dre1 modified***
a_Dre1        = zeros(1,N);
for i = 1:((intervals/Time_End*11.5)+1) %0:0.35
    a_Dre1(i) = a_OFF; 
end
for i = ((intervals/Time_End*11.5)+2):N %0.35:1
    a_Dre1(i) = a_ON; 
end



end

