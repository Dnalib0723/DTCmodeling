## Original code for "Cyclic LIN-42/PERIOD Precisely Times Stage-Specific Cell Migration Through Gene Circuit Dynamics" 
[![DOI](https://zenodo.org/badge/975369832.svg)](https://doi.org/10.5281/zenodo.15621826)

Here you can find the codes for the paper : Cyclic LIN-42/PERIOD precisely times stage-specific cell migration through gene circuit dynamics
## Package Installation
• Matlab 2021b
• Python 3.9.7

# Instructions for Running Simulations and Analyses

This document outlines how to use the provided MATLAB codebase for simulating and analyzing a gene regulatory network. It details procedures for deterministic and stochastic simulations, nullcline analysis, and investigations into SCF-DRE-1 activity and its onset time variations.

---

## 1. Deterministic and Stochastic Simulation of WT and Mutant Phenotype

This section covers running deterministic and stochastic simulations for both wild-type (WT) and mutant phenotypes, corresponding to Figures 2B-G, 3B, 3C, 3G, S3B, S4B, and S9.

* To run the simulation, execute **`stochastic simulation.m`** located in the **`./stochastic simulation/`** folder.
    * **1.1 Gene Mutation:** To simulate a gene mutation, set the **`gene switch`** variable to **`0`** within **`stochastic simulation.m`**.
    * **1.2 Deterministic Simulation:** For deterministic simulations, set the **`noise tuning`** variable to **`0`** in both **`stochasticsim_mod27det.m`** and **`stochasticsim_mod27sc.m`**.
    * **1.3 Excluding BLMP-1 Activation:** To exclude BLMP-1 activation of lin-29 through regulating gene X, switch the regulatory function of gene X from **`LogicgeneX`** to **`LogicgeneX_noAC`** in **`stochasticsim_mod27det.m`** and **`stochasticsim_mod27sc.m`**.

---

## 2. Nullcline Analysis

This section explains how to perform nullcline analysis, corresponding to Figure 4.

* To run the analysis, execute **`nullcline_analysis.m`**.
    * **2.1 Phase Map Generation:** Generate the phase map by setting the ratio of maximum LIN-42, SCF-DRE-1, and DAF-12 ligand production.
    * **2.2 BLMP-1 Activation:** To include BLMP-1 activation of lin-29, set the parameter **`neg_strength_BLMP1rpgeneX`** to **`0`**. To exclude it, set it to **`1`**.

---

## 3. SCF-DRE-1 Activity Variation

This section describes how to simulate SCF-DRE-1 activity variations, corresponding to Figures 5Bii and 5Biii.

* To run the simulation, execute **`Dre1_activity_variation.m`** located in the **`./Dre-1 activity variation/`** folder.
    * **3.1 Plotting Output:** In **`Dre1_activity_variation.m`**, choose either **`K_Dre1_re_Blmp1_sc`** or **`gama_Dre1_re_Blmp1_sc`** to plot the output when tuning the threshold of BLMP-1 degradation or the degradation rate of BLMP-1 regulated by SCF-DRE-1, respectively.
    * **3.2 Tuning Parameters:** In **`para28.m`**, multiply the scaling factor for either parameter **`K_Dre1_re_Blmp1_sc`** or **`gama_Dre1_re_Blmp1_sc`** to tune the threshold of BLMP-1 degradation or the degradation rate of BLMP-1 by SCF-DRE-1, respectively.

---

## 4. Dre-1 Onset Time Variation

This section details how to simulate Dre-1 onset time variations, corresponding to Figure 5Bii.

* To run the simulation, execute **`Dre1_onset_time_variation.m`** located in the **`./Dre-1_onset_time_variation/`** folder.

---

## 5. Noise-Induced Switching Test (Constant LIN-42 Level Input)

This section covers the noise-induced switching test using a constant LIN-42 level as input, corresponding to Figures 6Ai and 6Bi.

* To run the simulation, execute **`noise induced switching_constant LIN-42.m`** located in the **`./noise induced switching constant LIN-42/`** folder.
    * **5.1 LIN-42 Production Ratio:** Use the parameter **`L42_sc`** to set the ratio of maximum LIN-42 production within **`noise induced switching_constant LIN-42.m`**.
    * **5.2 Gene Mutation:** To set a gene mutation, switch the **`gene switch`** variable to **`0`** in **`noise induced switching_constant LIN-42.m`**.
    * **5.3 Deterministic Simulation:** For deterministic simulations, switch the **`noise tuning`** variable to **`0`** in both **`stochasticsim_mod27det.m`** and **`stochasticsim_mod27sc.m`**.
    * **5.4 Compare Results:** Run **`swithingplot.m`** to compare the simulation results.

---

## 6. Noise-Induced Switching Test (Cyclic LIN-42 Level Input)

This section covers the noise-induced switching test using a cyclic LIN-42 level as input, corresponding to Figures 6Aii-iii and 6Bii-iii.

* To run the simulation, execute **`noise_induced_switching_cyclicLIN42.m`** located in the **`./noise induced switching cyclic LIN-42/`** folder.
    * **6.1 Gene Mutation:** To set a gene mutation, switch the **`gene switch`** variable to **`0`** in **`noise induced switching_constant LIN-42.m`**.
    * **6.2 Deterministic Simulation:** For deterministic simulations, switch the **`noise tuning`** variable to **`0`** in both **`stochasticsim_mod27det.m`** and **`stochasticsim_mod27sc.m`**.
    * **6.3 LIN-42 Noise Frequency:** Set **`gamma_L42`** to 5 times its default value to increase the LIN-42 noise frequency by fivefold in **`stochasticsim_mod27det.m`** and **`stochasticsim_mod27sc.m`**.
    * **6.4 Compare DTC Turning Times (Distribution):** Run **`plotnetworkI.m`** to compare the distribution of DTC turning times simulated with either 1X (default) or 5X LIN-42 noise frequency.
    * **6.5 Compare DTC Turning Times (Heatmap):** Run **`plotnetworkII.m`** to compare the heatmap of DTC turning time simulated with 1X (default) or 5X LIN-42 noise frequency.

---

## 7. Bistability Analysis

This section outlines how to perform bistability analysis, corresponding to Figures S7 and S8.

* To run the analysis, execute **`bistability_analysis.m`** located in the **`./bistability analysis/`** folder.
    * **7.1 Gene Mutation:** To set a gene mutation, switch the **`gene switch`** variable to **`0`** in **`bistability_analysis.m`**.
    * **7.2 High UNC-5 Simulation:** Run **`deterministic_mod30.m`** and **`deterministic_mod30_noAC.m`** for simulating from a high UNC-5 level. Files with **`_noAC`** indicate simulations excluding BLMP-1 activation of lin-29.
    * **7.3 Low UNC-5 Simulation:** Run **`deterministic_mod30m.m`** and **`deterministic_mod30m_noAC.m`** for simulating from a low UNC-5 level. Files with **`_noAC`** indicate simulations excluding BLMP-1 activation of lin-29.

---

## 8. Simple Regulation Test

This section describes how to perform a simple regulation test, corresponding to Figure S10.

* To run the analysis, execute **`simple model.py`**.
    * **8.1 Stochastic Simulation (Including BLMP-1 Activation):**
        * For a 'low to high' LIN-42 transition, use the regulatory function **`f_as_g`**.
        * For a 'high to low' LIN-42 transition, use the regulatory function **`f_de_r`**.
    * **8.2 Stochastic Simulation (Excluding BLMP-1 Activation):**
        * For a 'low to high' LIN-42 transition, use the regulatory function **`f_noAC_as_g`**.
        * For a 'high to low' LIN-42 transition, use the regulatory function **`f_noAC_de_r`**.

---

This document provides a comprehensive guide to utilizing the provided code for various biological simulations and analyses. By following these instructions, you can reproduce the figures mentioned and explore the underlying regulatory mechanisms.
