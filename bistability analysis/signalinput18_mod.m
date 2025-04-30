function signalinput18_mod()

global intervals N t
global a_ON a_OFF
global L2_ON L2_OFF L3_ON L3_OFF L4_ON a_ON_Lin42
global a_Daf12 a_Lin42 a_Dre1


% *** Daf 12***
a_Daf12        = zeros(1,N);
for i = 1:((intervals/20*L3_ON)+1) 
    a_Daf12(i) = a_ON; 
end
for i = ((intervals/20*L3_ON)+2):N 
    a_Daf12(i) = a_ON; 
end

% *** Lin-42***
a_Lin42        = zeros(1,N);
			

for i = 1:((intervals/20*L4_ON)+1)
    a_Lin42(i) = a_ON_Lin42*(0.4*(sin((2*pi/8).*(t(i)-2)))+0.4+0.2);	
end
     
	Lin42_min = a_Lin42(i);
	
for i = ((intervals/20*L4_ON)+2):N
    a_Lin42(i) = Lin42_min;
end	

% *** Dre1 modified***
a_Dre1        = zeros(1,N);
for i = 1:((intervals/20*L3_OFF)+1) %0:0.35
%    a_Dre1(i) = a_OFF;   
    a_Dre1(i) = a_ON*0.2; 
end
for i = ((intervals/20*L3_OFF)+2):N %0.35:1
    a_Dre1(i) = a_ON*0.2; 
end



end

