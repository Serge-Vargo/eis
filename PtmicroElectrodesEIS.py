''' @verbatim
EIS circuit model for Pt, Pt-black and TiN bio-micro-electrodes for neural
recording applications is developed.

This script replicates the results demonstrated in: 
W. Franks, I. Schenker, P. Schmutz and A. Hierlemann, 
"Impedance characterization and modeling of electrodes for biomedical applications,"
in IEEE Transactions on Biomedical Engineering, vol. 52, no. 7, pp. 1295-1302,
July 2005, doi: 10.1109/TBME.2005.847523.

@endverbatim '''

import numpy as np
import matplotlib.pyplot as plt
import eismodels as eis

''' Workspace variables '''
## Constants
R_GAS = 8.3145      # (J/mol.K) gas constant 
ROOM_TEMP = 298     # (K)
FARADAY = 96485     # (s.A/mol) Faraday's constant 
EPS_0 = 8.854e-12   # (F/m) permittibity of free space
Q_ELEM = 1.602e-19  # (C) elementary charge (of electron)
K_BOLTZ = 1.381e-23 # (J/K) boltazmann constant

## Electrolyte parameters
eps_r = 78        # () relative permittivity of double layer - dilute aqueous solution at 25 C
d_OHP = 5e-10     # (m) thickness of double layer - physiological saline at 25 C
z = 4             # ()  charge of ion in solution - O_2 (oxigen) reduction reaction
n0 = 9.3e25       # (ions/m^3) bulk concentration - physiological saline at 25 C 
rho = 72          # (ohms.cm) resistivity - physiological saline at 25 C 
Ut = K_BOLTZ*ROOM_TEMP/Q_ELEM # (V) thermal voltage
#Ut = 0.0259       # (V) thermal voltage


## Frequency range
freq_low  = 1e-2     # (Hz) lowest freq
freq_high = 1e5      # (Hz) highest freq 1 MHz
freq_range = np.linspace(freq_low,freq_high,100000+1)  # Frequency range

## Experimental Values 
Vo = 5e-3         # (V) applied electrode potential (perturbation signal)
V_step = np.linspace(Vo,Vo + 1.5, 11)  # step potential
n = 0.8

''' Replicate results of paper.
The circuit model parameters:
    --> Constant-Phase-Element Impedance (Z_CPE), which represents interface capacitance impedance
    --> Charge Transfer Resistance R_ct
    * these two impedances above are shunted (in parallel)
    --> Solution (electrolyte) resistance R_s
    * in series with the other two elements.
    
    ** The Warburgh impedance due to diffusion of the chemical reactants in 
       solution is not included. Based on materials and frequency, not significant
 '''
## Constant-Phase-Element Impedance
# Interferance Capcitance
cap_H = eis.cap_stern_layer(d_OHP, EPS_0, eps_r)    # Stern Layer capacitance
diffu_length = eis.length_debye(EPS_0, eps_r, Ut, n0, z, Q_ELEM)   # diffusion layer length
cap_G = eis.cap_diff_layer(diffu_length, EPS_0, eps_r, z, Vo, Ut)  # Diffusion layer capacitance
cap_dl = eis.cap_double_layer(cap_H, cap_G)   # double layer capacitance

# CPE impedance
z_cpe = eis.imp_cpe(freq_range, n, cap_dl)    #CPE impedance 

## Charge Transfer Resistance


# =============================================================================
# 
# # Equilibrium Exchange Current Density and R_ct
# eq_current_den = eis.Jo_eq(Faraday, k_c, C_a, beta, V_eq, R, T):
# 
# 
# 
# # Randles with Warburg impedance
# R_solu = 20         # (Ohms) electrolyte solution resistance 
# warburg_coeff = 150 # from text (given 100uM concentration & 1.6e-5cm^2/c diff coef)
# Z_wab = eis.Zwarburg(warburg_coeff, freq_range)
# R_ct = R_polariz
# 
# Z_rand_cell = eis.ZinSeries(R_solu, eis.ZinParallel(Z_Cdl,eis.ZinSeries(R_ct,Z_wab)))
# 
# eis.bodePlot(Z_rand_cell, freq_range, "Bode Plot - Randles w/ Warburg Imp. Cell")
# eis.nyquistPlot(Z_rand_cell, "Re", "-Im", "Nyquist Plot - Randles w/ Warburg Imp. Cell") 
# 
# =============================================================================
plt.show()

