# all model parameters
# create dictionary based on this
import math
import numpy as np

"""
% Jordan E. Massad
% Oct 18, 2004
% Define Model Material Parameters
%MATPARAM
%[mparams, cparams] = matparam(cs)   
%Input  
% n - null input.
%Output
% mparams - material parameters (see below)
function mparams = matparam(n)    

mparams = zeros(19,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Material Parameters
mparams(1) = 37.2e3;                %EA [MPa]
mparams(2) = 25.5e3;                %EM [MPa]
mparams(3) = 0.0300;                %eT []
mparams(4) = 100.00;                %dsig (Mean Gap Stress) [MPa]
h = 1e0;                            %Convective Heat Transfer Coeff [W/m^2/K]
mparams(5) = 2e-3*(1e3/5e2+1/7)*h;  %Geometry-dependent HTC [MW/m^3/K]
                                	%   Formula: surf_area/volume * h.
rho = 6.45e3;                   	%SMA density [kg/m^3]
cv = 2.30e2;	                    %Specific heat (constant volume) [J/(kg*K)]
mparams(6) = rho*cv*1e-6;           %pcv [MJ/m^3/K]
mparams(7) = 11.0;                  %lamA [1e-6/K]  Austenite Thermal Expansion
mparams(8) = 6.6;               	%lamM [1e-6/K]  Martensite Thermal Expansion
mparams(9) = 100.0;             	%rhoA [1e-8 Ohm*m]  Austenite Resistivity
mparams(10) = 80.0;             	%rhoM [1e-8 Ohm*m]  Martensite Resistivity
mparams(11) = 248.56;           	%T_M [K]  Measured Mf<T_M<Ms
mparams(12) = 260;              	%T_A or Teq [K]  Measured As<T_A<Af or Ms<Teq<Af
mparams(13) = 274.0;            	%s_hi [MPa] Loading Transf. Point at T_hi>T_lo
mparams(14) = 293.0;             	%T_hi [K]  Temperature for s_hi
mparams(15) = 200.58;           	%s_lo [MPa] (scale measured tensile data to shear)
mparams(16) = 243.0;             	%T_lo [K]  Temperature for s_lo
mparams(17) = 52^2; 	            %Gap Variance [MPa^2] {52 default)
mparams(18) = 55^2;  	            %Effective Stress Variance [MPa^2] {45 default}
mparams(19) = 0.1*mparams(6);        %Specific Heat Difference (dc) [J/(kg*K)]
%  End of matparam.m
"""
# Constants used in Boltzman integrals
KB = 1.3806503e-23        # Boltzman constant [J/K]
VOL =  1e-0*KB            # Lattice layer volume [m^3]
MASS = 3e3*VOL**(1/3)      # Mass of a layer [kg]


class Model_Parameters(object):

    def __init__(self):
        """
        Initialize all parameters as a dictionary. 
        Default values from MATLAB code. Will want to use nitinol-specific values.
        """
        self.val                        = {} # initialize dictionary
        self.val["EA"]                  = 37.2e3 # EA [MPa]
        self.val["EM"]                  = 25.5e3 # EM [MPa]
        self.val["eT"]                  = 0.0300 # eT [] ???????
        self.val["dsig"]                = 100.00 # Mean Gap Stress [MPa]
        self.val["CHT"]                 = 1.0 # Convective Heat Transfer Coeff [W/m^2/K]
        self.val["HTC"]                 = 2e-3*(1e3/5e2+1/7)*self.val["CHT"] # Geometry-dependent HTC [MW/m^3/K].   Formula: surf_area/volume * h.
        self.val["rho"]                 = 6.45e3 # Volume density [kg/m^3]
        self.val["c"]                   = 2.30e2 # Specific heat (constant volume) [J/(kg*K)]
        self.val["pcv"]                 = self.val["rho"]*self.val["c"]*1.0e-6 # pcv [MJ/m^3/K]
        self.val["LA"]                  = 11.0 # lamA [1e-6/K]  Austenite Thermal Expansion
        self.val["LM"]                  = 6.6 # lamM [1e-6/K]  Martensite Thermal Expansion
        self.val["RA"]                  = 100.0 # rhoA [1e-8 Ohm*m]  Austenite Resistivity
        self.val["RM"]                  = 80.0 # rhoM [1e-8 Ohm*m]  Martensite Resistivity
        self.val["TempM"]               = 248.56 # T_M [K]  Measured Mf<T_M<Ms
        self.val["TempA"]               = 260.0 # T_A or Teq [K]  Measured As<T_A<Af or Ms<Teq<Af
        self.val["s_hi"]                = 274.0 # s_hi [MPa] Loading Transf. Point at T_hi>T_lo
        self.val["T_hi"]                = 293.0 # T_hi [K]  Temperature for s_hi
        self.val["s_lo"]                = 200.59 # s_lo [MPa] (scale measured tensile data to shear)
        self.val["T_lo"]                = 243.0 # T_lo [K]  Temperature for s_lo
        self.val["gap_var"]             = 52.0**2# Gap Variance [MPa^2] {52 default)
        self.val["stress_var_eff"]      = 55.0**2# Effective Stress Variance [MPa^2] {45 default}
        self.val["deltac"]              = 0.1*self.val["pcv"] # Specific Heat Difference (dc) [J/(kg*K)]

    def compute_parameters(self):
        """
        These parameters are computed using the tabulated values above.

        """
        
        self.cval               = {}
        self.cval["Er"]         = self.val["EM"]/self.val["EA"] # relative E
        self.cval[""] = 
        self.cval[""] = 
        self.cval[""] = 
        self.cval[""] = 
        self.cval["T_A"] = 
        self.cval["T_A_sup"] = np.# eT highest temp beyond superelasticity 
        self.cval["qA"] = sqrt((VOL/KB)/self.val["EA"]/2)*1e3;               		# qA [sqrt(K)/MPa]
        self.cval["qM"] = sqrt((VOL/KB)/self.val["EM"]/2)*1e3;               		# qM [sqrt(K)/MPa]
        self.cval["relax_tau"] =  math.pi*math.sqrt(MASS / self.val["EA"] / VOL**(1/3))*1e-3 # relaxation tau [seconds]
        self.cval["eleven"] = 0
        self.cval["twelve"] = 1


    def G(self, T,dU,dS,T_R,E):
        """
        Returns the Gibbs free energy, accounting for chemical free energy.

        """
        return dU - T*dS + self.val["dc"]*(T - T_R - T*math.log(T/T_R)) - E