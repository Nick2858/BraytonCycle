import math
import sympy as sp
from sympy import *
from sympy.solvers import solve


class Brayton:
    def __init__(self, in_vel, in_temp, in_pres, in_area = 0.5, in_area_ratio = 0.1, comp_pres_ratio = 1, fuel_flow= 1, turb_pres_ratio = 0.1, noz_area_ratio = 10):
        """
        """
        
        #Inlet Params
        self.in_vel = in_vel
        self.in_temp = in_temp
        self.in_pres = in_pres
        self.in_area = in_area
        self.in_area_ratio = in_area_ratio

        #Compressor Params

        self.comp_pres_ratio = comp_pres_ratio
        
        #Combuston Params
        

        #Turbine Params
        self.turb_pres_ratio = turb_pres_ratio

        #nozfusor Params
        self.noz_area_ratio = noz_area_ratio

        #Constants
        self.R = 287

        #Flow Rate

        self.mass_flow = self.in_pres*self.in_area*self.in_vel/ (self.R * self.in_temp)
        self.fuel_flow = self.mass_flow/fuel_flow 


    def inletOutput(self, inlet = True):
        
        V2, P2, T2 = symbols("V2, P2, T2")

        h1 = 1/3000 * self.in_temp ** 3 + 287/10000 * self.in_temp ** 2 + 24394/25 * self.in_temp
        massflowrate = sp.Eq(self.in_pres*self.in_area_ratio*self.in_vel/self.in_temp, P2*V2/T2)
        isentropic = sp.Eq(T2/self.in_temp, (P2/self.in_pres)**(self.R/(0.0001*((T2 + self.in_temp)/2)**2 + 0.0574 * ((T2 + self.in_temp)/2) + 975.76  )))
        energy = sp.Eq(0, 1/3000 * T2 ** 3 + 287/10000 * T2 ** 2 + 24394/25 * T2 - h1 + (V2 ** 2 - self.in_vel**2)/2)
        sol = nsolve([massflowrate,isentropic, energy], [V2,P2,T2] , [self.in_vel*self.in_area_ratio, self.in_pres, self.in_temp], prec = 5)
        
        self.in_vel2 = float(sol[0])
        self.in_pres2 = float(sol[1])
        self.in_temp2 = float(sol[2])
       

    def compression(self):
        
        self.comp_pres = self.comp_pres_ratio * self.in_pres2

        w, T2 = symbols("w, T2")
        h1 = 1/3000 * self.in_temp2 ** 3 + 287/10000 * self.in_temp2 ** 2 + 24394/25 * self.in_temp2
        energy = sp.Eq(w, 1/3000 * T2 ** 3 + 287/10000 * T2 ** 2 + 24394/25 * T2 - h1)
        isentropic = sp.Eq(T2/self.in_temp2, (self.comp_pres_ratio)**(self.R/(0.0001*((T2 + self.in_temp2)/2)**2 + 0.0574 * ((T2 + self.in_temp2)/2) + 975.76  )))
        sol = nsolve([energy, isentropic], [w,T2], [1, 1], prec = 10)

        self.comp_work = sol[0] * self.mass_flow
        self.comp_temp = sol[1]
        
    
    def combustion(self):
        """
        Combustion of Kerosene (assume it's dodecane) 
        2 C12H16 (l) + 37 O2 (g) --> 24 CO2 (g) +26 H2O (g) Enthalpy of Combustion = -7901.74 kJ/mol 
        """ 
        self.outlet_flow = self.fuel_flow + self.mass_flow
        self.emissions = self.fuel_flow * 0.012/ (12 * 0.17034) 

        h1 = 1/3000 * self.comp_temp ** 3 + 287/10000 * self.comp_temp ** 2 + 24394/25 * self.comp_temp

        T2 = symbols("T2")
        Q = sp.Eq(self.outlet_flow * (1/3000 * T2 ** 3 + 287/10000 * T2 ** 2 + 24394/25 * T2 - h1), 7901740 * self.fuel_flow / 0.17034 )
        sol = nsolve([Q], [T2], [self.comp_temp + 1000], prec = 10)
        
        self.comb_temp = sol[0]
        self.comb_heat = 7901740 * self.fuel_flow / 0.17034 

        print(self.comp_temp, self.comb_temp)
    

    def turbine(self):
        self.turb_pres2 = self.turb_pres_ratio * self.comp_pres

        w, T2 = symbols("w, T2")
        h1 = 1/3000 * self.comb_temp ** 3 + 287/10000 * self.comb_temp ** 2 + 24394/25 * self.comb_temp
        energy = sp.Eq(w, 1/3000 * T2 ** 3 + 287/10000 * T2 ** 2 + 24394/25 * T2 - h1)
        isentropic = sp.Eq(T2/self.comb_temp, (self.turb_pres_ratio)**(self.R/(0.0001*((T2 + self.comb_temp)/2)**2 + 0.0574 * ((T2 + self.comb_temp)/2) + 975.76  )))
        sol = nsolve([energy, isentropic], [w,T2], [1, 1], prec = 10)

        self.turb_work = sol[0] * self.outlet_flow
        self.turb_temp = sol[1]
        

    def nozzle(self):
        
        V2, P2, T2 = symbols("V2, P2, T2")

        h1 = 1/3000 * self.turb_temp ** 3 + 287/10000 * self.turb_temp ** 2 + 24394/25 * self.turb_temp
        massflowrate = sp.Eq(self.turb_pres2*self.noz_area_ratio*self.in_vel2/self.turb_temp, P2*V2/T2)
        isentropic = sp.Eq(T2/self.turb_temp, (P2/self.turb_pres2)**(self.R/(0.0001*((T2 + self.turb_temp)/2)**2 + 0.0574 * ((T2 + self.turb_temp)/2) + 975.76  )))
        energy = sp.Eq(0, 1/3000 * T2 ** 3 + 287/10000 * T2 ** 2 + 24394/25 * T2 - h1 +(V2 ** 2 - self.in_vel2**2)/2)
        sol = nsolve([massflowrate,isentropic, energy], [V2,P2,T2] , [self.in_vel2*self.noz_area_ratio, self.turb_pres2, self.turb_temp], prec = 10)
        
        self.noz_vel2 = float(sol[0])
        self.noz_pres2 = float(sol[1])
        self.noz_temp2 = float(sol[2])
        self.noz_exit_area = self.outlet_flow*self.R*self.noz_temp2/(self.noz_pres2*self.noz_vel2)


    def work(self):
        self.net_work = self.turb_work + self.comp_work
        self.eff = -self.net_work/self.comb_heat
        print(self.turb_temp)
        print(self.mass_flow, self.eff, self.net_work, self.comp_work, self.turb_work, self.comb_heat)


A = Brayton(in_vel=10, in_temp = 280, in_pres = 101325, in_area= 0.5, in_area_ratio= 1, comp_pres_ratio= 40, fuel_flow= 90, turb_pres_ratio= 1/100, noz_area_ratio=1/1)

A.inletOutput()
A.compression()
A.combustion()
A.turbine()
A.nozzle()
A.work()