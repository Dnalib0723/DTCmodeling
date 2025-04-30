function stochasticsim_mod29det()

global PARA_W scale
global Z_Daf12 Z_Blmp1 Z_Unc5 Z_Lin42 Z_Dre1 Z_Lin29 Z_geneX
global b bm bp a_ON a_OFF
global b_Lin42 a_ON_Lin42 a_OFF_Lin42
global M Mdet N dt t Time_End T3_End T2_End
global basal_Lin29 basal_geneX basal_Unc5 basal_Blmp1 basal_Daf12L basal_Daf12C
global a_Daf12 a_Lin42 a_Dre1
global Daf12 Blmp1 Unc5 Lin29 Lin42 Dre1 geneX blmp1 lin29 unc5 Daf12L Daf12C
global Unc5_sto Lin42_sto

eval(['para'])

% basal production rate
basal_Daf12         = Z_Daf12*0.1;	basal_Blmp1	= Z_Blmp1*0.1;	basal_blmp1	= Z_Blmp1*0.1;    basal_Unc5  = Z_Unc5*0.1;   basal_unc5  = Z_Unc5*0.1;
basal_Lin42         = Z_Lin42*0.1;	basal_Dre1  = Z_Dre1 *0.1; 	basal_Lin29 = Z_Lin29*0.1;    basal_lin29 = Z_Lin29*0.1;  basal_geneX = Z_geneX*0.1;
basal_Daf12L         = Z_Daf12*0.1;  basal_Daf12C = Z_Daf12*0;

%randn('state',100) % set the random number

for j = 1:Mdet
    %Preallocation
    Daf12                 = zeros(1,N); blmp1         = zeros(1,N);   Blmp1         = zeros(1,N);  unc5          = zeros(1,N); Unc5          = zeros(1,N); lin29       = zeros(1,N) ; Lin29       = zeros(1,N);   Lin42      = zeros(1,N);    Dre1      	= zeros(1,N);    geneX       = zeros(1,N);  Daf12L       = zeros(1,N);  Daf12C       = zeros(1,N);  
    Daf12_BP              = zeros(1,N); blmp1_BP      = zeros(1,N);   Blmp1_BP      = zeros(1,N);  unc5_BP       = zeros(1,N); Unc5_BP       = zeros(1,N); lin29_BP	    = zeros(1,N);Lin29_BP	= zeros(1,N);   Lin42_BP   = zeros(1,N);    Dre1_BP   	= zeros(1,N);    geneX_BP    = zeros(1,N);  Daf12L_BP    = zeros(1,N);  Daf12C_BP    = zeros(1,N);        
    Daf12_DP              = zeros(1,N); blmp1_DP      = zeros(1,N);   Blmp1_DP      = zeros(1,N);  unc5_DP       = zeros(1,N); Unc5_DP       = zeros(1,N); lin29_DP 	= zeros(1,N); Lin29_DP 	= zeros(1,N);   Lin42_DP   = zeros(1,N);    Dre1_DP  	= zeros(1,N);    geneX_DP    = zeros(1,N);  Daf12L_DP    = zeros(1,N);  Daf12C_DP    = zeros(1,N);
    Daf12_GN              = zeros(1,N); blmp1_GN      = zeros(1,N);   Blmp1_GN      = zeros(1,N);  unc5_GN       = zeros(1,N); Unc5_GN       = zeros(1,N); lin29_GN    = zeros(1,N) ; Lin29_GN    = zeros(1,N);   Lin42_GN   = zeros(1,N);    Dre1_GN     = zeros(1,N);    geneX_GN    = zeros(1,N);  Daf12L_GN    = zeros(1,N);  Daf12C_GN    = zeros(1,N);
    Daf12_BPM             = zeros(1,N); blmp1_BPM     = zeros(1,N);   Blmp1_BPM     = zeros(1,N);  unc5_BPM      = zeros(1,N); Unc5_BPM      = zeros(1,N); lin29_BPM   = zeros(1,N) ; Lin29_BPM   = zeros(1,N);   Lin42_BPM  = zeros(1,N);    Dre1_BPM  	= zeros(1,N);    geneX_BPM   = zeros(1,N);  Daf12L_BPM   = zeros(1,N);  Daf12C_BPM   = zeros(1,N);
    Daf12_BPN             = zeros(1,N); blmp1_BPN     = zeros(1,N);   Blmp1_BPN     = zeros(1,N);  unc5_BPN      = zeros(1,N); Unc5_BPN      = zeros(1,N); lin29_BPN   = zeros(1,N) ;Lin29_BPN   = zeros(1,N);   Lin42_BPN  = zeros(1,N);    Dre1_BPN	= zeros(1,N);    geneX_BPN   = zeros(1,N);  Daf12L_BPN   = zeros(1,N);  Daf12C_BPN   = zeros(1,N);
    Daf12_DPM             = zeros(1,N); blmp1_DPM     = zeros(1,N);   Blmp1_DPM     = zeros(1,N);  unc5_DPM      = zeros(1,N); Unc5_DPM      = zeros(1,N); lin29_DPM   = zeros(1,N) ;Lin29_DPM   = zeros(1,N);   Lin42_DPM  = zeros(1,N);    Dre1_DPM	= zeros(1,N);    geneX_DPM   = zeros(1,N);  Daf12L_DPM   = zeros(1,N);  Daf12C_DPM   = zeros(1,N);
    Daf12_DPN             = zeros(1,N); blmp1_DPN     = zeros(1,N);   Blmp1_DPN     = zeros(1,N);  unc5_DPN      = zeros(1,N); Unc5_DPN      = zeros(1,N); lin29_DPN   = zeros(1,N) ;Lin29_DPN   = zeros(1,N);   Lin42_DPN  = zeros(1,N);    Dre1_DPN    = zeros(1,N);    geneX_DPN   = zeros(1,N);  Daf12L_DPN   = zeros(1,N);  Daf12C_DPN   = zeros(1,N);
    Logic_Blmp1           = zeros(1,N);   Logic_Unc5    = zeros(1,N);   Logic_Lin29   = zeros(1,N); Logic_geneX  = zeros(1,N); Logic_Daf12L  = zeros(1,N);   Logic_Daf12C  = zeros(1,N);   
    
    %Initial Concentration

	eval(['inicon_mod27s'])
	
	%buffering set
	
    Daf12_BP_buff         = 0;            Daf12_DP_buff = 0;            Blmp1_BP_buff = 0;          Blmp1_DP_buff = 0;	   Unc5_BP_buff  = 0;	  Unc5_DP_buff  = 0;
    Lin29_BP_buff         = 0;            Lin29_DP_buff = 0;            Lin42_BP_buff = 0;          Lin42_DP_buff = 0;     Dre1_BP_buff  = 0;     Dre1_DP_buff  = 0;
	geneX_BP_buff	      = 0;            geneX_DP_buff = 0;            blmp1_BP_buff = 0;          blmp1_DP_buff = 0;	   unc5_BP_buff  = 0;	  unc5_DP_buff  = 0;
    lin29_BP_buff         = 0;            lin29_DP_buff = 0;            Daf12L_BP_buff = 0;         Daf12L_DP_buff = 0;    Daf12C_BP_buff = 0;    Daf12C_DP_buff = 0;
	
    % Noise tuning 
    noise_Daf12           = 0;            noise_Lin42	= 0;            noise_Dre1	  = 0;          noise_geneX    = 0;   
    noise_Lin29           = 0;            noise_Blmp1	= 0;            noise_Unc5	  = 0;          noise_Daf12L   = 0;  
	noise_lin29           = 0;            noise_blmp1	= 0;            noise_unc5	  = 0;          noise_Daf12C   = 0;
	
	% degration rate %/sec
    gama_Daf12      = gama_common*0.03;      gama_Unc5   = gama_common*0.25;      gama_Blmp1 = gama_Blmp1;     gama_geneX      = gama_common;
    gama_Lin29      = gama_common*0.25;      gama_Lin42  = gama_common*1;      gama_Dre1    = gama_common;  gama_Daf12L      = gama_common*0.3;
	gama_unc5       = 5*gama_Unc5;      gama_blmp1 = 5*gama_Blmp1;        gama_lin29  = 5*gama_Lin29;    gama_Daf12C     = gama_common;
	
	Daf12_SS_Wt_dist(j,1)	= Daf12(1);  Lin42_SS_Wt_dist(j,1) = Lin42(1);  Dre1_SS_Wt_dist(j,1) = Dre1(1);  Blmp1_SS_Wt_dist(j,1) = Blmp1(1);  Lin29_SS_Wt_dist(j,1) = Lin29(1);    Unc5_SS_Wt_dist(j,1) = Unc5(1);
	geneX_SS_Wt_dist(j,1)   = geneX(1);  blmp1_SS_Wt_dist(j,1) = blmp1(1);  lin29_SS_Wt_dist(j,1) = lin29(1);  unc5_SS_Wt_dist(j,1) = unc5(1);  Daf12L_SS_Wt_dist(j,1)	= Daf12L(1);  Daf12C_SS_Wt_dist(j,1)	= Daf12C(1);
	
	for i=2:N
	    %randn
	    eval(['num_mod3'])
% * Daf12 *
        % BPM = Burst Production Mean DPM = Degredation of mRNA Mean
		a_Daf12_i = a_Daf12(i);
        Daf12_BPM(i)      = (eval(Daf12_num{1})+eval(Daf12_num{2}))*b; 
        Daf12_BPN(i)      = noise_Daf12*sqrt(2*(eval(Daf12_num{1})+eval(Daf12_num{2})))*b*randn;
        Daf12_BP(i)       = Daf12_BPM(i)+Daf12_BPN(i);
		
		if noise_Daf12~=0;
		Daf12_BP(i)=round(Daf12_BP(i));
		end
		
        if (Daf12_BP_buff + Daf12_BP(i)) < 0 ; Daf12_BP_buff  = Daf12_BP_buff + Daf12_BP(i);    Daf12_BP(i)   = 0;
        else                                    Daf12_BP(i)   = Daf12_BP_buff + Daf12_BP(i);    Daf12_BP_buff = 0; end
		Daf12_ir = Daf12(i-1);
        Daf12_DPM(i)      = eval(Daf12_num{3}) ;
        Daf12_DPN(i)      = noise_Daf12*sqrt(eval(Daf12_num{3}))*randn ;
        Daf12_DP(i)       = Daf12_DPM(i) +Daf12_DPN(i);
		
		if noise_Daf12~=0;
		Daf12_DP(i)=round(Daf12_DP(i));
		end
		
        if (Daf12_DP_buff + Daf12_DP(i)) < 0 ; Daf12_DP_buff  = Daf12_DP_buff + Daf12_DP(i);    Daf12_DP(i)   = 0;
        else                                   Daf12_DP(i)    = Daf12_DP_buff + Daf12_DP(i);    Daf12_DP_buff = 0; end
        Daf12(i)          = Daf12(i-1) +Daf12_BP(i) -Daf12_DP(i);
        if Daf12(i)       < 0 ;  Daf12(i) = 0; end
        
        % * Lin42 *
        a_Lin42_i = a_Lin42(i);
        Lin42_BPM(i)      = (eval(Lin42_num{1})+eval(Lin42_num{2}))*b_Lin42; 
        Lin42_BPN(i)      = noise_Lin42*sqrt(2*(eval(Lin42_num{1})+eval(Lin42_num{2})))*b_Lin42*randn;
        Lin42_BP(i)       = Lin42_BPM(i)+Lin42_BPN(i);
		
		if noise_Lin42~=0;
		Lin42_BP(i)=round(Lin42_BP(i));
		end
		
        if (Lin42_BP_buff + Lin42_BP(i)) < 0 ; Lin42_BP_buff  = Lin42_BP_buff + Lin42_BP(i);    Lin42_BP(i)   = 0;
        else                                    Lin42_BP(i)   = Lin42_BP_buff + Lin42_BP(i);    Lin42_BP_buff = 0; end
		Lin42_ir = Lin42(i-1);
        Lin42_DPM(i)      = eval(Lin42_num{3}) ;
        Lin42_DPN(i)      = noise_Lin42*sqrt(eval(Lin42_num{3}))*randn;
        Lin42_DP(i)       = Lin42_DPM(i) +Lin42_DPN(i);
		
		if noise_Lin42~=0;
		Lin42_DP(i)=round(Lin42_DP(i));
		end
		
        if (Lin42_DP_buff + Lin42_DP(i)) < 0 ; Lin42_DP_buff  = Lin42_DP_buff + Lin42_DP(i);    Lin42_DP(i)   = 0;
        else                                   Lin42_DP(i)    = Lin42_DP_buff + Lin42_DP(i);    Lin42_DP_buff = 0; end
        Lin42(i)          = Lin42(i-1) +Lin42_BP(i) -Lin42_DP(i) ;
        if Lin42(i)       < 0 ;  Lin42(i) = 0; end
        
        % * Daf12L *
		Logic_Daf12L_i     = LogicDaf12L(Lin42(i-1));
		Daf12L_BPM(i)      = (eval(Daf12L_num{1})+eval(Daf12L_num{2}))*b;
		Daf12L_BPN(i)      = noise_Daf12L*sqrt(2*(eval(Daf12L_num{1})+eval(Daf12L_num{2})))*b*randn;
		Daf12L_BP(i)       = Daf12L_BPM(i)+Daf12L_BPN(i);
		
	    if noise_Daf12L~=0;
		Daf12L_BP(i)=round(Daf12L_BP(i));
		end
		
		if (Daf12L_BP_buff + Daf12L_BP(i)) < 0 ; Daf12L_BP_buff  = Daf12L_BP_buff + Daf12L_BP(i);    Daf12L_BP(i)   = 0;
		else                                    Daf12L_BP(i)   = Daf12L_BP_buff + Daf12L_BP(i);    Daf12L_BP_buff = 0; end
		Daf12L_ir          = Daf12L(i-1);
		Daf12L_DPM(i)      = eval(Daf12L_num{3});
		Daf12L_DPN(i)      = noise_Daf12L*sqrt(eval(Daf12L_num{3}))*randn ;
		Daf12L_DP(i)       = Daf12L_DPM(i) +Daf12L_DPN(i);
		
	    if noise_Daf12L~=0;
		Daf12L_DP(i)=round(Daf12L_DP(i));
		end
		
		if (Daf12L_DP_buff + Daf12L_DP(i)) < 0 ; Daf12L_DP_buff  = Daf12L_DP_buff + geneX_DP(i);    Daf12L_DP(i)   = 0;
		else                                   Daf12L_DP(i)    = Daf12L_DP_buff + Daf12L_DP(i);    Daf12L_DP_buff = 0; end
		Daf12L(i)          = Daf12L(i-1) +Daf12L_BP(i) -Daf12L_DP(i);
		if Daf12L(i)       < 0 ; Daf12L(i) = 0; end
        
        % * Daf12C *
		Logic_Daf12C_i     = LogicDaf12C(Daf12(i-1),Daf12L(i-1));
		Daf12C_BPM(i)      = (eval(Daf12C_num{1})+eval(Daf12C_num{2}))*b;
		Daf12C_BPN(i)      = noise_Daf12C*sqrt(2*(eval(Daf12C_num{1})+eval(Daf12C_num{2})))*b*randn;
		Daf12C_BP(i)       = Daf12C_BPM(i)+Daf12C_BPN(i);
		
	    if noise_Daf12C~=0;
		Daf12C_BP(i)=round(Daf12C_BP(i));
		end
		
		if (Daf12C_BP_buff + Daf12C_BP(i)) < 0 ; Daf12C_BP_buff  = Daf12C_BP_buff + Daf12C_BP(i);    Daf12C_BP(i)   = 0;
		else                                    Daf12C_BP(i)   = Daf12C_BP_buff + Daf12C_BP(i);    Daf12C_BP_buff = 0; end
		Daf12C_ir          = Daf12C(i-1);
		Daf12C_DPM(i)      = eval(Daf12C_num{3});
		Daf12C_DPN(i)      = noise_Daf12C*sqrt(eval(Daf12C_num{3}))*randn ;
		Daf12C_DP(i)       = Daf12C_DPM(i) +Daf12C_DPN(i);
		
	    if noise_Daf12C~=0;
		Daf12C_DP(i)=round(Daf12C_DP(i));
		end
		
		if (Daf12C_DP_buff + Daf12C_DP(i)) < 0 ; Daf12C_DP_buff  = Daf12C_DP_buff + Daf12C_DP(i);    Daf12C_DP(i)   = 0;
		else                                   Daf12C_DP(i)    = Daf12C_DP_buff + Daf12C_DP(i);    Daf12C_DP_buff = 0; end
		Daf12C(i)          = Daf12C(i-1) +Daf12C_BP(i) -Daf12C_DP(i);
		if Daf12C(i)       < 0 ; Daf12C(i) = 0; end
        
        % * Dre1 *
        a_Dre1_i = a_Dre1(i); 
        Dre1_BPM(i)       = (eval(Dre1_num{1})+eval(Dre1_num{2}))*b; 
        Dre1_BPN(i)       = noise_Dre1*sqrt(2*(eval(Dre1_num{1})+eval(Dre1_num{2})))*b*randn;
        Dre1_BP(i)        = Dre1_BPM(i)+Dre1_BPN(i);
		
		if noise_Dre1~=0;
		Dre1_BP(i)=round(Dre1_BP(i));
		end
		
        if (Dre1_BP_buff  + Dre1_BP(i)) < 0 ; Dre1_BP_buff  = Dre1_BP_buff + Dre1_BP(i);    Dre1_BP(i)   = 0;
        else                                    Dre1_BP(i)   = Dre1_BP_buff + Dre1_BP(i);    Dre1_BP_buff = 0; end
		Dre1_ir = Dre1(i-1);
        Dre1_DPM(i)       = eval(Dre1_num{3}) ;
        Dre1_DPN(i)       = noise_Dre1*sqrt(eval(Dre1_num{3}))*randn ;
        Dre1_DP(i)        = Dre1_DPM(i) +Dre1_DPN(i);
		
		if noise_Dre1~=0;
		Dre1_DP(i)=round(Dre1_DP(i));
		end
		
        if (Dre1_DP_buff  + Dre1_DP(i)) < 0 ; Dre1_DP_buff  = Dre1_DP_buff + Dre1_DP(i);    Dre1_DP(i)   = 0;
        else                                   Dre1_DP(i)    = Dre1_DP_buff + Dre1_DP(i);    Dre1_DP_buff = 0; end
        Dre1(i)           = Dre1(i-1) +Dre1_BP(i) -Dre1_DP(i);
        if Dre1(i)        < 0 ;  Dre1(i) = 0; end
        
		% * blmp-1 *   
        Logic_Blmp1_i     = LogicBlmp1(Daf12C(i-1),Lin42(i-1),Lin29(i-1),Blmp1(i-1),Dre1(i-1)); 
        blmp1_BPM(i)      = (eval(blmp1_num{1})+eval(blmp1_num{2}))*bm; 
        blmp1_BPN(i)      = noise_blmp1*sqrt(2*(eval(blmp1_num{1})+eval(blmp1_num{2})))*bm*randn;
        blmp1_BP(i)       = blmp1_BPM(i)+blmp1_BPN(i);
		
		if noise_blmp1~=0;
		blmp1_BP(i)=round(blmp1_BP(i));
		end
		
        if (blmp1_BP_buff + blmp1_BP(i)) < 0 ; blmp1_BP_buff  = blmp1_BP_buff + blmp1_BP(i);    blmp1_BP(i)   = 0;
        else                                    blmp1_BP(i)   = blmp1_BP_buff + blmp1_BP(i);    blmp1_BP_buff = 0; end
		%Hdre1ACblmp       = H_dre1ACblmp(Dre1(i-1));
		blmp1_ir          = blmp1(i-1);
        blmp1_DPM(i)      = eval(blmp1_num{3});
        blmp1_DPN(i)      = noise_blmp1*sqrt(eval(blmp1_num{3}))*randn ; %(gama_Blmp1+(gama_Blmp1_Dre1*H_dre1ACblmp)
        blmp1_DP(i)       = blmp1_DPM(i) +blmp1_DPN(i);
		
		if noise_blmp1~=0;
		blmp1_DP(i)=round(blmp1_DP(i));
		end
        
        if (blmp1_DP_buff + blmp1_DP(i)) < 0 ; blmp1_DP_buff  = blmp1_DP_buff + blmp1_DP(i);    blmp1_DP(i)   = 0;
        else                                   blmp1_DP(i)    = blmp1_DP_buff + blmp1_DP(i);    blmp1_DP_buff = 0; end
                
        blmp1(i)          = blmp1(i-1) +blmp1_BP(i) -blmp1_DP(i) ;
        if blmp1(i)       < 0 ; blmp1(i) = 0; end
		
        % * Blmp-1 *   
        %Logic_Blmp1_i     = LogicBlmp1(Daf12(i-1),Lin42(i-1),Lin29(i-1),Blmp1(i-1),Dre1(i-1)); 
        Blmp1_BPM(i)      = (eval(Blmp1_num{1})+eval(Blmp1_num{2}))*bp; 
        Blmp1_BPN(i)      = noise_Blmp1*sqrt(1*(eval(Blmp1_num{1})+eval(Blmp1_num{2})))*bp*randn;
        Blmp1_BP(i)       = Blmp1_BPM(i)+Blmp1_BPN(i);
		
		if noise_Blmp1~=0;
		Blmp1_BP(i)=round(Blmp1_BP(i));
		end
		
        if (Blmp1_BP_buff + Blmp1_BP(i)) < 0 ; Blmp1_BP_buff  = Blmp1_BP_buff + Blmp1_BP(i);    Blmp1_BP(i)   = 0;
        else                                    Blmp1_BP(i)   = Blmp1_BP_buff + Blmp1_BP(i);    Blmp1_BP_buff = 0; end
		Hdre1ACblmp       = H_dre1ACblmp(Dre1(i-1));
		Blmp1_ir          = Blmp1(i-1);
        Blmp1_DPM(i)      = eval(Blmp1_num{3});
        Blmp1_DPN(i)      = noise_Blmp1*sqrt(eval(Blmp1_num{3}))*randn ; %(gama_Blmp1+(gama_Blmp1_Dre1*H_dre1ACblmp)
        Blmp1_DP(i)       = Blmp1_DPM(i) +Blmp1_DPN(i);
		
		if noise_Blmp1~=0;
		Blmp1_DP(i)=round(Blmp1_DP(i));
		end
        
        if (Blmp1_DP_buff + Blmp1_DP(i)) < 0 ; Blmp1_DP_buff  = Blmp1_DP_buff + Blmp1_DP(i);    Blmp1_DP(i)   = 0;
        else                                   Blmp1_DP(i)    = Blmp1_DP_buff + Blmp1_DP(i);    Blmp1_DP_buff = 0; end
                
        Blmp1(i)          = Blmp1(i-1) +Blmp1_BP(i) -Blmp1_DP(i) ;
        if Blmp1(i)       < 0 ; Blmp1(i) = 0; end
		
        % * geneX *
		Logic_geneX_i     = LogicgeneX(Blmp1(i-1));
		%Logic_geneX_i     = LogicgeneX_noAC(Blmp1(i-1));
		geneX_BPM(i)      = (eval(geneX_num{1})+eval(geneX_num{2}))*b;
		geneX_BPN(i)      = noise_geneX*sqrt(2*(eval(geneX_num{1})+eval(geneX_num{2})))*b*randn;
		geneX_BP(i)       = geneX_BPM(i)+geneX_BPN(i);
		
	    if noise_geneX~=0;
		geneX_BP(i)=round(geneX_BP(i));
		end
		
		if (geneX_BP_buff + geneX_BP(i)) < 0 ; geneX_BP_buff  = geneX_BP_buff + geneX_BP(i);    geneX_BP(i)   = 0;
		else                                    geneX_BP(i)   = geneX_BP_buff + geneX_BP(i);    geneX_BP_buff = 0; end
		geneX_ir          = geneX(i-1);
		geneX_DPM(i)      = eval(geneX_num{3});
		geneX_DPN(i)      = noise_geneX*sqrt(eval(geneX_num{3}))*randn ;
		geneX_DP(i)       = geneX_DPM(i) +geneX_DPN(i);
		
	    if noise_geneX~=0;
		geneX_DP(i)=round(geneX_DP(i));
		end
		
		if (geneX_DP_buff + geneX_DP(i)) < 0 ; geneX_DP_buff  = geneX_DP_buff + geneX_DP(i);    geneX_DP(i)   = 0;
		else                                   geneX_DP(i)    = geneX_DP_buff + geneX_DP(i);    geneX_DP_buff = 0; end
		geneX(i)          = geneX(i-1) +geneX_BP(i) -geneX_DP(i);
		if geneX(i)       < 0 ; geneX(i) = 0; end
	
	
	    % * lin-29 *
		Logic_Lin29_i     = LogicLin29(Lin42(i-1),Blmp1(i-1),geneX(i-1));
		lin29_BPM(i)      = (eval(lin29_num{1})+eval(lin29_num{2}))*bm; 
		lin29_BPN(i)      = noise_lin29*sqrt(2*(eval(lin29_num{1})+eval(lin29_num{2})))*bm*randn;
		lin29_BP(i)       = lin29_BPM(i)+lin29_BPN(i);
		
	    if noise_lin29~=0;
		lin29_BP(i)=round(lin29_BP(i));
		end
		
		if (lin29_BP_buff + lin29_BP(i)) < 0 ; lin29_BP_buff  = lin29_BP_buff + lin29_BP(i);    lin29_BP(i)   = 0;
		else                                   lin29_BP(i)   = lin29_BP_buff + lin29_BP(i);    lin29_BP_buff = 0; end
		lin29_ir          = lin29(i-1);
		lin29_DPM(i)      = eval(lin29_num{3});
		lin29_DPN(i)      = noise_lin29*sqrt(eval(lin29_num{3}))*randn ;
		lin29_DP(i)       = lin29_DPM(i) +lin29_DPN(i);
		
		if noise_lin29~=0;
		lin29_DP(i)=round(lin29_DP(i));
		end
		
		if (lin29_DP_buff + lin29_DP(i)) < 0 ; lin29_DP_buff  = lin29_DP_buff + lin29_DP(i);    lin29_DP(i)   = 0;
		else                                   lin29_DP(i)    = lin29_DP_buff + lin29_DP(i);    lin29_DP_buff = 0; end
		lin29(i)          = lin29(i-1) +lin29_BP(i) -lin29_DP(i) ;
		if lin29(i)       < 0 ; lin29(i) = 0; end
		
		% * Lin-29 *
		%Logic_Lin29_i     = LogicLin29(Lin42(i-1),Blmp1(i-1),geneX(i-1));
		Lin29_BPM(i)      = (eval(Lin29_num{1})+eval(Lin29_num{2}))*bp; 
		Lin29_BPN(i)      = noise_Lin29*sqrt(1*(eval(Lin29_num{1})+eval(Lin29_num{2})))*bp*randn;
		Lin29_BP(i)       = Lin29_BPM(i)+Lin29_BPN(i);
		
		if noise_Lin29~=0;
		Lin29_BP(i)=round(Lin29_BP(i));
		end
		
		if (Lin29_BP_buff + Lin29_BP(i)) < 0 ; Lin29_BP_buff  = Lin29_BP_buff + Lin29_BP(i);    Lin29_BP(i)   = 0;
		else                                    Lin29_BP(i)   = Lin29_BP_buff + Lin29_BP(i);    Lin29_BP_buff = 0; end
		Lin29_ir          = Lin29(i-1);
		Lin29_DPM(i)      = eval(Lin29_num{3});
		Lin29_DPN(i)      = noise_Lin29*sqrt(eval(Lin29_num{3}))*randn ;
		Lin29_DP(i)       = Lin29_DPM(i) +Lin29_DPN(i);
		
		if noise_Lin29~=0;
		Lin29_DP(i)=round(Lin29_DP(i));
		end
		
		if (Lin29_DP_buff + Lin29_DP(i)) < 0 ; Lin29_DP_buff  = Lin29_DP_buff + Lin29_DP(i);    Lin29_DP(i)   = 0;
		else                                   Lin29_DP(i)    = Lin29_DP_buff + Lin29_DP(i);    Lin29_DP_buff = 0; end
		Lin29(i)          = Lin29(i-1) +Lin29_BP(i) -Lin29_DP(i) ;
		if Lin29(i)       < 0 ; Lin29(i) = 0; end
        
		% * unc-5 *
        Logic_Unc5_i	  = LogicUnc5(Daf12C(i-1),Lin29(i-1),Blmp1(i-1));
        unc5_BPM(i)       = (eval(unc5_num{1})+eval(unc5_num{2}))*bm; 
        unc5_BPN(i)       = noise_unc5*sqrt(2*(eval(unc5_num{1})+eval(unc5_num{2})))*bm*randn;
        unc5_BP(i)        = unc5_BPM(i)+unc5_BPN(i);
		
		if noise_unc5~=0;
		unc5_BP(i)=round(unc5_BP(i));
		end
		
        if (unc5_BP_buff  + unc5_BP(i)) < 0 ; unc5_BP_buff  = unc5_BP_buff + unc5_BP(i);    unc5_BP(i)   = 0;
        else                                    unc5_BP(i)  = unc5_BP_buff + unc5_BP(i);    unc5_BP_buff = 0; end
		unc5_ir           = unc5(i-1);
        unc5_DPM(i)       = eval(unc5_num{3});
        unc5_DPN(i)       = noise_unc5*sqrt(eval(unc5_num{3}))*randn ;
        unc5_DP(i)        = unc5_DPM(i) +unc5_DPN(i);
		
		if noise_unc5~=0;
		unc5_DP(i)=round(unc5_DP(i));
		end
		
        if (unc5_DP_buff  + unc5_DP(i)) < 0 ; unc5_DP_buff  = unc5_DP_buff + unc5_DP(i);    unc5_DP(i)   = 0;
        else                                  unc5_DP(i)    = unc5_DP_buff + unc5_DP(i);    unc5_DP_buff = 0; end
        unc5(i)           = unc5(i-1) +unc5_BP(i) -unc5_DP(i) ;
        if unc5(i)        < 0 ;  unc5(i) = 0; end
        
        % * Unc-5 *
        %Logic_Unc5_i	  = LogicUnc5(Daf12(i-1),Lin29(i-1),Blmp1(i-1));
        Unc5_BPM(i)       = (eval(Unc5_num{1})+eval(Unc5_num{2}))*bp*0.7; 
        Unc5_BPN(i)       = noise_Unc5*sqrt(1*(eval(Unc5_num{1})+eval(Unc5_num{2})))*bp*0.7*randn;
        Unc5_BP(i)        = Unc5_BPM(i)+Unc5_BPN(i);
		
		if noise_Unc5~=0;
		Unc5_BP(i)=round(Unc5_BP(i));
		end
		
        if (Unc5_BP_buff  + Unc5_BP(i)) < 0 ; Unc5_BP_buff  = Unc5_BP_buff + Unc5_BP(i);    Unc5_BP(i)   = 0;
        else                                    Unc5_BP(i)   = Unc5_BP_buff + Unc5_BP(i);    Unc5_BP_buff = 0; end
		Unc5_ir           = Unc5(i-1);
        Unc5_DPM(i)       = eval(Unc5_num{3});
        Unc5_DPN(i)       = noise_Unc5*sqrt(eval(Unc5_num{3}))*randn ;
        Unc5_DP(i)        = Unc5_DPM(i) +Unc5_DPN(i);
		
		if noise_Unc5~=0;
		Unc5_DP(i)=round(Unc5_DP(i));
		end
		
        if (Unc5_DP_buff  + Unc5_DP(i)) < 0 ; Unc5_DP_buff  = Unc5_DP_buff + Unc5_DP(i);    Unc5_DP(i)   = 0;
        else                                   Unc5_DP(i)    = Unc5_DP_buff + Unc5_DP(i);    Unc5_DP_buff = 0; end
        Unc5(i)           = Unc5(i-1) +Unc5_BP(i) -Unc5_DP(i) ;
        if Unc5(i)        < 0 ;  Unc5(i) = 0; end
        
        % [Steady state data collecting]
       Daf12_SS_Wt_dist(j,i)	= Daf12(i);  Lin42_SS_Wt_dist(j,i) = Lin42(i);  Dre1_SS_Wt_dist(j,i) = Dre1(i);  Blmp1_SS_Wt_dist(j,i) = Blmp1(i);  Lin29_SS_Wt_dist(j,i) = Lin29(i);  Unc5_SS_Wt_dist(j,i) = Unc5(i);
       geneX_SS_Wt_dist(j,i)	= geneX(i);  blmp1_SS_Wt_dist(j,i) = blmp1(i);  lin29_SS_Wt_dist(j,i) = lin29(i);  unc5_SS_Wt_dist(j,i) = unc5(i);  Daf12L_SS_Wt_dist(j,i)	= Daf12L(i);  Daf12C_SS_Wt_dist(j,i)	= Daf12C(i);
	end
end

Unc5_sto      = zeros(M,N); Unc5_sto = Unc5_SS_Wt_dist;
Lin42_sto      = zeros(M,N); Lin42_sto = Lin42_SS_Wt_dist;


end