global Daf12_num Lin42_num Dre1_num Blmp1_num Lin29_num geneX_num Unc5_num lin29_num unc5_num blmp1_num Daf12L_num Daf12C_num

Daf12_num = {'basal_Daf12*gama_Daf12*a_Daf12_i*dt','(1-basal_Daf12)*gama_Daf12*a_Daf12_i*Z_Daf12*dt','gama_Daf12*Daf12_ir*dt'};

Daf12L_num = {'basal_Daf12L*gama_Daf12L*a_ON*dt','Logic_Daf12L_i*(1-basal_Daf12L)*gama_Daf12L*a_ON*Z_Daf12*dt','gama_Daf12L*Daf12L_ir*dt'};

Daf12C_num = {'basal_Daf12C*gama_Daf12C*a_ON*dt','Logic_Daf12C_i*(1-basal_Daf12C)*gama_Daf12C*a_ON*Z_Daf12*dt','gama_Daf12C*Daf12C_ir*dt'};

Lin42_num = {'basal_Lin42*gama_Lin42*a_Lin42_i*dt','(1-basal_Lin42)*gama_Lin42*a_Lin42_i*Z_Lin42*dt','gama_Lin42*Lin42_ir*dt'};

Dre1_num = {'basal_Dre1*gama_Dre1*a_Dre1_i*dt','(1-basal_Dre1)*gama_Dre1*a_Dre1_i*Z_Dre1*dt','gama_Dre1*Dre1_ir*dt'};

Blmp1_num = {'basal_Blmp1*gama_blmp1*blmp1_ir*dt','(1-basal_Blmp1)*gama_blmp1*blmp1_ir*dt','(gama_Blmp1+(gama_Blmp1_Dre1*Hdre1ACblmp))*Blmp1_ir*dt'};

Lin29_num = {'basal_Lin29*gama_lin29*lin29_ir*dt','(1-basal_Lin29)*gama_lin29*lin29_ir*dt','gama_Lin29*Lin29_ir*dt'};

geneX_num = {'basal_geneX*gama_geneX*a_ON*dt','Logic_geneX_i*(1-basal_geneX)*gama_geneX*a_ON*Z_geneX*dt','gama_geneX*geneX_ir*dt'};

Unc5_num = {'basal_Unc5*gama_unc5*unc5_ir*dt','(1-basal_Unc5)*gama_unc5*unc5_ir*dt','gama_Unc5*Unc5_ir*dt'};

lin29_num = {'basal_lin29*gama_Lin29*a_ON*dt','Logic_Lin29_i*(1-basal_lin29)*gama_Lin29*a_ON*Z_Lin29*dt','gama_lin29*lin29_ir*dt'};

unc5_num = {'basal_unc5*gama_Unc5*a_ON*dt','Logic_Unc5_i*(1-basal_unc5)*gama_Unc5*a_ON*Z_Unc5*dt','gama_unc5*unc5_ir*dt'};

blmp1_num = {'basal_blmp1*gama_Blmp1*a_ON*dt','Logic_Blmp1_i*(1-basal_blmp1)*gama_Blmp1*a_ON*Z_Blmp1*dt','gama_blmp1*blmp1_ir*dt'};