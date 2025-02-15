import math
import sympy as sp
from sympy import *
from sympy.solvers import solve

class Brayton:
    def __init__(self, in_vel, in_temp, in_pres, in_area_ratio, comp_pres_ratio = 1, heat_addition = 1, turb_pres_ratio = 1, dif_area_ratio = 1, gamma = 1.4):
        """
        """
        
        #Inlet Params
        self.in_vel = in_vel
        self.in_temp = in_temp
        self.in_pres = in_pres
        self.in_area_ratio = in_area_ratio

        #Compressor Params

        self.comp_pres_ratio = comp_pres_ratio
        
        #Combuston Params
        self.heat_addition = heat_addition

        #Turbine Params
        self.turb_pres_ratio = turb_pres_ratio

        #Diffusor Params
        self.dif_area_ratio = dif_area_ratio

        #Constants
        self.gamma = gamma
        self.cp = 1004

        print("done")
    def inletOutput(self):
        
        print("hi")
        V2, P2, T2 = symbols("V2, P2, T2")
        massflowrate = sp.Eq(self.in_pres*self.in_area_ratio*self.in_vel/self.in_temp, P2*V2/T2)
        isentropic = sp.Eq(T2/self.in_temp, (P2/self.in_pres)**((self.gamma - 1)/self.gamma))
        energy = sp.Eq(0, self.cp * (T2 - self.in_temp)+(V2 ** 2 - self.in_vel**2)/2)
        sol = nsolve([massflowrate,isentropic, energy], [V2,P2,T2] , [self.in_vel*self.dif_area_ratio, self.in_pres, self.in_temp], prec = 10)
        

        self.in_vel2 = sol[0]
        self.in_pres2 = sol[1]
        self.in_temp2 = sol[2] 
        

A = Brayton(10, 273, 101325, 0.1)
A.inletOutput()