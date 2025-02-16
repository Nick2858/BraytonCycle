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
        self.R = 287

        print("done")
    def inletOutput(self, inlet = True):
        
        V2, P2, T2 = symbols("V2, P2, T2")

        h1 = 1/3000 * self.in_temp ** 3 + 287/10000 * self.in_temp ** 2 + 24394/25 * self.in_temp
        massflowrate = sp.Eq(self.in_pres*self.in_area_ratio*self.in_vel/self.in_temp, P2*V2/T2)
        isentropic = sp.Eq(T2/self.in_temp, (P2/self.in_pres)**(self.R/(0.0001*((T2 + self.in_temp)/2)**2 + 0.0574 * ((T2 + self.in_temp)/2) * 975.76  )))
        energy = sp.Eq(0, 1/3000 * T2 ** 3 + 287/10000 * T2 ** 2 + 24394/25 * T2 - h1 +(V2 ** 2 - self.in_vel**2)/2)
        sol = nsolve([massflowrate,isentropic, energy], [V2,P2,T2] , [self.in_vel*self.dif_area_ratio, self.in_pres, self.in_temp], prec = 10)
        
        self.in_vel2 = float(sol[0])
        self.in_pres2 = float(sol[1])
        self.in_temp2 = float(sol[2])
       

        print(self.in_vel2, self.in_pres2,self.in_temp2)

    def comp(self):
        
        print("hi")
        self.comp_pres = self.comp_pres_ratio * self.in_pres2

        T2 = (self.comp_pres_ratio)**(self.R/(0.0001*((T2 + self.in_temp2)/2)**2 + 0.0574 * ((T2 + self.in_temp2)/2) * 975.76  )) * self.in_temp2
        w = (0.0001*((T2 + self.in_temp2)/2)**2 + 0.0574 * ((T2 + self.in_temp2)/2) * 975.76 )

        print(w, T2)
        


A = Brayton(30, 273, 101325, 0.1, 10)

A.inletOutput()
A.comp()