function PARA_W = grep(filename)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global PARA_W

fid = fopen(filename,'r');
current_lineNo = 0;
Line_Number   = 4;

reduction_factor = 5; % to scale the threshold values

while ~feof(fid)
    current_lineNo = current_lineNo+1;
    param_line=fgetl(fid);%# read line by line
   
    if current_lineNo == Line_Number
	    fprintf('current_lineNo    = %2.0f \n',current_lineNo)
        A = sscanf(param_line,'%f');
		
        if length(A) >= 42; %160
            n=1;
            SetNum         	 = A(n,1); n = n+1; %A(1,1); 
            PheTypes         = A(n,1); n = n+1; %A(2,1); 
            Integr_Time  	 = A(n,1); n = n+1; %A(3,1);
			stdNum			 = 1;
            DTC_Turning_Time = A(n,1); n = n+1; %A(4,1); 
            gama_common      = A(n,1); n = n+1; %A(5,1); 
            gama_Blmp1       = A(n,1); n = n+1; %A(6,1); 
            gama_Blmp1_Dre1  = A(n,1); n = n+1; %A(7,1); 
            
			K_Daf12_ac_Unc5  = A(n,1)/reduction_factor; n = n+1; %A(8,1); 
            K_Daf12_re_Blmp1 = A(n,1)/reduction_factor; n = n+1; %A(9,1); 
            K_Lin42_re_Lin29 = A(n,1)/reduction_factor; n = n+1; %A(10,1); 
            K_Lin42_ac_Blmp1 = A(n,1)/reduction_factor; n = n+1; %A(11,1); 
            K_Dre1_re_Blmp1  = A(n,1)/reduction_factor; n = n+1; %A(12,1);  
            K_Lin29_re_Blmp1 = A(n,1)/reduction_factor; n = n+1; %A(13,1);  
            K_Lin29_ac_Unc5  = A(n,1)/reduction_factor; n = n+1; %A(14,1);  
            K_Blmp1_ac_Blmp1 = A(n,1)/reduction_factor; n = n+1; %A(15,1);  
            K_Blmp1_re_Lin29 = A(n,1)/reduction_factor; n = n+1; %A(16,1);  
            K_Blmp1_re_Unc5  = A(n,1)/reduction_factor; n = n+1; %A(17,1); 
			K_Blmp1_re_geneX = A(n,1)/reduction_factor; n = n+1; %A(18,1);
			K_geneX_re_Lin29 = A(n,1)/reduction_factor; n = n+1; %A(19,1);
			
			PARA_W = zeros(20,1);
            PARA_W(1) = SetNum;
            PARA_W(2) = PheTypes;
            PARA_W(3) = Integr_Time;
			PARA_W(4) = stdNum;
			PARA_W(5) = DTC_Turning_Time;
			PARA_W(6) = gama_common;
			PARA_W(7) = gama_Blmp1;
			PARA_W(8) = gama_Blmp1_Dre1;
			PARA_W(9) = K_Daf12_ac_Unc5;
			PARA_W(10) = K_Daf12_re_Blmp1;
			PARA_W(11) = K_Lin42_re_Lin29;
		    PARA_W(12) = K_Lin42_ac_Blmp1;
			PARA_W(13) = K_Dre1_re_Blmp1;
			PARA_W(14) = K_Lin29_re_Blmp1;
			PARA_W(15) = K_Lin29_ac_Unc5;
			PARA_W(16) = K_Blmp1_ac_Blmp1; 
			PARA_W(17) = K_Blmp1_re_Lin29;
			PARA_W(18) = K_Blmp1_re_Unc5;
			PARA_W(19) = K_Blmp1_re_geneX;
			PARA_W(20) = K_geneX_re_Lin29;
        end
		%K_Daf12_ac_Unc5  = 140;

		fprintf('SetNum           = %2.0f \n',SetNum)
        fprintf('PheTypes         = %2.0f \n',PheTypes)
        fprintf('Integr_Time      = %2.1f \n',Integr_Time)
        fprintf('stdNum           = %2.0f \n',stdNum)
        fprintf('DTC_Turning_Time = %2.0f \n',DTC_Turning_Time)
        fprintf('gama_common      = %2.8f \n',gama_common)
        fprintf('gama_Blmp1       = %2.10f \n',gama_Blmp1)
        fprintf('gama_Blmp1_Dre1  = %2.8f \n\n',gama_Blmp1_Dre1)
        
        fprintf('K_Daf12_ac_Unc5  = %2.2f \n',K_Daf12_ac_Unc5)
        fprintf('K_Daf12_re_Blmp1 = %2.2f \n',K_Daf12_re_Blmp1)
        fprintf('K_Lin42_re_Lin29 = %2.2f \n',K_Lin42_re_Lin29)
        fprintf('K_Lin42_ac_Blmp1 = %2.2f \n',K_Lin42_ac_Blmp1)
        fprintf('K_Dre1_re_Blmp1  = %2.2f \n',K_Dre1_re_Blmp1)
        fprintf('K_Lin29_re_Blmp1 = %2.2f \n',K_Lin29_re_Blmp1)
        fprintf('K_Lin29_ac_Unc5  = %2.2f \n',K_Lin29_ac_Unc5)
        fprintf('K_Blmp1_ac_Blmp1 = %2.2f \n',K_Blmp1_ac_Blmp1)
        fprintf('K_Blmp1_re_Lin29 = %2.2f \n',K_Blmp1_re_Lin29)
        fprintf('K_Blmp1_re_Unc5  = %2.2f \n',K_Blmp1_re_Unc5)
		fprintf('K_Blmp1_re_geneX = %2.2f \n',K_Blmp1_re_geneX)
		fprintf('K_geneX_re_Lin29  = %2.2f\n\n',K_geneX_re_Lin29)
%% Parameters panel

%% Parameters panel

    end %end of if length
end % end of while

fclose(fid);
end % function end

