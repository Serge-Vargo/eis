'''
This script serves as an example of how to use the module eismodels.py
The examples used are taken from the Application Note by GAMRY Instruments
introducing the basics of Electrochemical Impedance Spectroscopy EIS
'''

import numpy as np
import matplotlib.pyplot as plt
import eismodels as eis

''' Workspace variables '''
## Constants
R_GAS = 8.3145      # gas constant (J/mol.K)
ROOM_TEMP = 298     # K
FARADAY = 96485     # FARADAY's constant (s.A/mol)
EPS_0 = 8.854e-12   # permittibity of free space

# Frequency range
f_lo  = 0.1     # (Hz) lowest freq
f_hi = 1e6      # (Hz) highest freq 1 MHz
f_range = np.linspace(f_lo, f_hi, int(f_hi/f_lo) + 1)  # Frequency range

''' Replicate results of App note "Basics of EIS" '''

# =============================================================================
# # --- Lissajous Figure --- #
# t = np.arange(0, 4*np.pi, 0.1)
# x = np.sin(t + np.pi/5) + 5
# y = 2*np.cos(t) + 5
#
# eis.lissajousFig(x, y, r'$sin(t + \frac{\pi}{5}) + 5$', r'$2cos(t)$', \
#                   'Lissajous Plot')
# =============================================================================

# ----------------------------------------------------------------------------
# --- Purely Capacitive Coating (Excelent Coating) --- #
r_electrolyte = 500     # (Ohms) poorly conductive solution
eps_r_coat = 6          # () relative permittivity of the coating layer
surf_a_coat = 1 * 1e-4  # (cm^2 = 1e-4 m^2) surface area, assume 1cm^2
depth_coat = 25 * 1e-6  # (um = 1 * 1e-6 m) depth of the capacitive coating (25 um)
c_coat = eis.cap_double_plate(eps_r_coat, EPS_0, surf_a_coat, depth_coat) # (F) capacitance of coating
#c_coat_approx = 200e-12     # (F) realistic for a 1 cm^2 sample, a 25um coating and eps_r = 6

z_coating = eis.imp_cap(c_coat, f_range)      #Impedance of capacitive coating
z_sys_pure_cap_coating = eis.imp_series(r_electrolyte, z_coating)  #Z of interface 

eis.bode_plot(z_sys_pure_cap_coating, f_range, "Bode Plot - Capacitive Coating")
eis.nyquist_plot(z_sys_pure_cap_coating, "Re", "-Im", "Nyquist Plot - Capacitive Coating")

# ----------------------------------------------------------------------------
# --- Warburg Impedance ---#
# Electrolyte parameters
z_ion = 2                   # ()  charge of ion in solution
c_bulk = 100 * 1e-3         # (uM = 1e-3 mol/m^3) bulk concentration
diff_coeff = 1.6e-5 * 1e-4  # (cm^2/s = 1e-4 m^2/s) diffusion coefficient
surf_area = 1 * 1e-4        # (cm^2 = 1e-4 m^2) surface area, assume 1cm^2

sigma_coef = eis.warburg_coeff(R_GAS, ROOM_TEMP, z_ion, FARADAY, \
                               surf_area, c_bulk*2, diff_coeff,  \
                               c_bulk*2, diff_coeff)

warb_imp = eis.imp_warburg(sigma_coef, f_range)


# ----------------------------------------------------------------------------
# ---  Simplified Randles Cell --- #
r_solu_simp_randle = 20 # (Ohms) electrolyte solution resistance
r_polariz = 250         # (Ohms) polarization
c_dob_layer = 40e-6     # (F) 1cm^2 sample with 40uF/cm^2 cap per surface

z_cap_dl = eis.imp_cap(c_dob_layer, f_range)
z_paralel_simp_randle = eis.imp_parallel(r_polariz, z_cap_dl)
z_simp_randle = eis.imp_series(r_solu_simp_randle, z_paralel_simp_randle)

eis.bode_plot(z_simp_randle, f_range, "Bode Plot - Randles Cell")
eis.nyquist_plot(z_simp_randle, "Re", "-Im", "Nyquist Plot - Randles Cell")

# ----------------------------------------------------------------------------
# --- Randles with Warburg impedance --- #
r_solu = 20         # (Ohms) electrolyte solution resistance
warburg_coeff = 150 # from text (given 100uM concentration & 1.6e-5cm^2/c diff coef)
z_warb = eis.imp_warburg(warburg_coeff, f_range)
r_ct = r_polariz

z_rand_cell = eis.imp_series(r_solu, eis.imp_parallel(z_cap_dl,eis.imp_series(r_ct,z_warb)))

eis.bode_plot(z_rand_cell, f_range, "Bode Plot - Randles w/ Warburg Imp. Cell")
eis.nyquist_plot(z_rand_cell, "Re", "-Im", "Nyquist Plot - Randles w/ Warburg Imp. Cell")

plt.show()
