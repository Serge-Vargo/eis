# @verbatim
# This file contains functions that can be used to build circuit models
# representing the electrical behaviour of neural electrodes.
# The models used are
# @endverbatim

import numpy as np
import matplotlib.pyplot as plt

## Impedances
def imp_cap(cap, f):
    # Equivalent impedance of a capacitor
    # cap -> capacitance
    # f -> frequency
    return 1/(1j*2*np.pi*f*cap)

def imp_series(z1, z2):
    # Equivalent series impedance of two elements z1 & z2
    return (z1 + z2)

def imp_parallel(z1, z2):
    # Equivalent parallel impedance of two elements z1 & z2 
    return (z1*z2)/(z1 + z2)

def res_elyte(k, l, cross_a):    
    # Electrolyte Resistance (depends of electrolyte properties (k) and geometry)
    # k -> conductivity of the solution
    # l -> length of electrolyte
    # cross_a -> crossectional area
    return l/(k*cross_a) 
    
def imp_warburg(sigma, freq): 
    # Impedance due to Diffusion of infinite thickness
    # sigma -> Warburg coefficient
    # freq  -> frequency
    return (sigma / np.sqrt(2*np.pi*freq)) * (1-1j)

def imp_warburg_nonideal(sigma, freq, thck, diff_coef):
    # Impedance due to Diffusion for a finite thickness
    # sigma     -> Warburg coefficient
    # freq      -> frequency
    # thck      -> Nerst Diffusion layer thickness
    # diff_coef -> some avg value of diffusion coefficients (anode & cathode)
    return (sigma / np.sqrt(2*np.pi*freq)) * (1 - 1j)  \
            * np.tanh(thck*np.sqrt(1j*2*np.pi*freq / diff_coef))
            
def warburg_coeff(R, T, n, F, surf_a, Cob, Do, Crb, Dr):
    # Warburg coefficient
    # R   -> gas constant 
    # T   -> temperature
    # n   -> number of electrons involved
    # F   -> Faraday's constant
    # surf_a   -> Surface area of electrode
    # Cob -> Bulk Concentration of oxidant
    # Do  -> Diffusion coefficient of oxidant
    # Crb -> Bulk Concentration of reductant
    # Dr  -> Diffusion coefficient of reductant
    return (R*T/(np.sqrt(2)*n**2*F**2*surf_a)) * (1/(Cob*np.sqrt(Do)) + 1/(Crb*np.sqrt(Dr)))

def cap_double_plate(eps_r, eps_o, surf_a, d):
    # Capacitance of a double plate capacitor
    # eps_r -> relative permitivity of the electrolyte (and/or other materials)
    # eps_o -> permittivity of free space
    # surf_a-> surface area of plates
    # d     -> thickness between plates
    return eps_r*eps_o*surf_a/d

def imp_cpe(f, n, c_cpe):
    # Impedance of a constant phase element (relevant for double layer capacitors)
    # f -> frequency 
    # n -> non-ideality constant, represents inhomogeneities. n <= 1 Generally 0.9 - 1.0   (1.0 -> ideal cap)
    # c_cpe -> Interface Capacitance 
    return 1/((1j*2*np.pi*f * c_cpe)**n)

def cap_double_layer(cap_stern, cap_diff):
    # Capacitance of a constant phase element based on double layer (diff + boundry in series)
    # cap_stern --> capacitance of stern layer
    # cap_diff  --> capacitance of diffuse layer
    return  1/( 1/cap_stern + 1/cap_diff )

def cap_stern_layer(d_OHP, eps_o, eps_r):
    # Capacitance of the stern layer (part of the double layer)
    # d_OHP --> thickness of layer
    # eps_o --> permitivitty of free space
    # eps_r --> permitivitty of the double layer
    return eps_o*eps_r/d_OHP

def cap_diff_layer(len_deb, eps_o, eps_r, z, Vo, v_t):
    # Capacitance of the diffusive layer (part of the double layer)
    # eps_o   --> permitivitty of free space
    # eps_r   --> permitivitty of the double layer
    # len_deb --> Debye length (diffusion length)
    # z       --> charge on the ion in solution
    # Vo      --> applied electrode potential
    # v_t     --> thermal voltage
    return (eps_o*eps_r*np.cosh(0.5*z*Vo/v_t))/len_deb

def length_debye(eps_0, eps_r, v_t, n0, z, q):
    # Debye Length, representing the length of diffusion concentration gradient
    # eps_o --> permitivitty of free space
    # eps_r --> permitivitty of the double layer
    # v_t    --> thermal voltage
    # n0    --> bulk number concentration
    # z     --> charge on the ion in solution
    # q     --> elementary charge
    return np.sqrt(eps_0*eps_r*v_t/(2*n0*z**2*q))

def jo_eq(F, k_c, C_a, beta, V_eq, R, T):
    # At equilibrium equal, and oposite reduction and oxidation currents flow 
    # accross the electrode-electrolyte interface. The mag of this current -> Jo_eq
    # The equilibrium current is a measure of the eletrode's ability to participate
    # in a exchange current reactions. 
    # For ideally polarizable electrode Jo_eq = 0
    # For ideally unpolarizable electrode Jo_eq = infinity
    # F    --> Faraday's constant
    # k_c  --> reduction reaction rate constant
    # c_a  --> concentration of electron-acceptor ions "a" in solution plane of the interface
    # beta --> symetry factor (between oxidation and reduction?)
    # V_eq --> equilibrium potential
    # R    --> gas constant
    # T    --> temperature
    return F*k_c*C_a*np.exp(-beta*F*V_eq/(R*T))

def j_low_field_approx(J_o, F, eta, R, T):
    # Current density at equilibrium using the low-field approximation to the 
    # Butler-Volmer equation (eta < 0.005/z) 
    # J_o  --> equilibrium exchange current density
    # F    --> Faraday's constant
    # eta  --> applied overpotential 
    # R    --> gas constant 
    # T    --> temperature
    return (J_o*F*eta)/(R*T)
    
def res_charge_transf(R, T, z, F, J_o): 
    # Charge Transfer Resistance (based on speed of reaction) at equilibrium 
    # exchange current
    # R    --> gas constant 
    # T    --> temperature
    # z    --> number of electrons involved
    # F    --> Faraday's constant
    # J_o  --> equilibrium exchange current density
    return (R*T)/(z*F*J_o)


## Plots
def lissajous_fig(x, y, xlabel='', ylabel='', tle=''):
    # Plot Lissajous Graph
    # x -> x values
    # y -> y values
    # xlabel -> x-axis label 
    # ylabel -> y-axis label
    # tle -> title of plot
    plt.figure()
    plt.plot(x, y)
    plt.xlabel(xlabel)
    plt.xlim(xmin=0)
    plt.ylabel(ylabel)
    plt.ylim(ymin=0)
    plt.title(tle)
    plt.grid(True, which='both')
    plt.show(block = False)
    return

def bode_plot(Z, f, title=''):
    # Bode Plot for magnitude and phase in same plot
    # Z -> impedance
    # f -> frequency
    # title -> title of plot
    mag = np.log10(np.abs(Z))
    phase = 180 / np.pi * np.angle(Z)
    fig, ax1 = plt.subplots()

    ax1.plot(np.log10(f), mag, color='tab:red')
    ax1.set_xlabel(r'$log_{10}$ of Freq (Hz)') 
    ax1.set_ylabel(r'$log_{10}$ of |Z|', color='tab:red')
    ax1.set_yticks(np.linspace(ax1.get_yticks()[0], ax1.get_yticks()[-1], \
                               len(ax1.get_yticks())))
    ax1.tick_params(axis='y', labelcolor='tab:red')
    ax1.grid()

    ax2 = ax1.twinx()  # second y-axis with same x-axis values
    ax2.plot(np.log10(f), phase, color='tab:blue')
    ax2.set_ylabel(r'Phase ($\phi$)', color='tab:blue')
    ax2.set_yticks(np.linspace(ax2.get_yticks()[0], ax2.get_yticks()[-1], \
                               len(ax1.get_yticks())))
    ax2.tick_params(axis='y', labelcolor='tab:blue')
    ax2.grid()
    
    plt.title(title)
    #plt.grid(True, which='both')
    plt.show(block = False)
    return

def nyquist_plot(Z, xlabel='', ylabel='', tle=''): 
    # Nyquist Plot for given impedance 
    # Z -> impedance
    # xlabel -> x-axis label 
    # ylabel -> y-axis label
    # tle    -> title of plot
    Re = Z.real
    negIm = -Z.imag
    plt.figure()
    plt.plot(Re, negIm, color='tab:red')
    plt.xlabel(xlabel)
    plt.xlim(xmin=0)
    plt.ylabel(ylabel)
    plt.ylim(ymin=0)
    plt.title(tle)
    plt.grid(True, which='both')
    plt.show(block = False)
    return