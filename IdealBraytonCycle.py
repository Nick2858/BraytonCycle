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

       

    def inletOutput(self):
        
        V2, P2, T2 = symbols("V2, P2, T2")

        h1 = 1/3000 * self.in_temp ** 3 + 287/10000 * self.in_temp ** 2 + 24394/25 * self.in_temp
        massflowrate = sp.Eq(self.in_pres*self.in_vel/self.in_temp, P2*V2*self.in_area_ratio/T2)
        isentropic = sp.Eq(T2/self.in_temp, (P2/self.in_pres)**(self.R/(0.0001*((T2 + self.in_temp)/2)**2 + 0.0574 * ((T2 + self.in_temp)/2) + 975.76  )))
        energy = sp.Eq(0, 1/3000 * T2 ** 3 + 287/10000 * T2 ** 2 + 24394/25 * T2 - h1 + (V2 ** 2 - self.in_vel**2)/2)
        sol = nsolve([massflowrate,isentropic, energy], [V2,P2,T2] , [self.in_vel*self.in_area_ratio, self.in_pres, self.in_temp], prec = 50)
        
        self.in_vel2 = float(sol[0])
        self.in_pres2 = float(sol[1])
        self.in_temp2 = float(sol[2])
       
        
    def compression(self):
       
        self.comp_pres = self.comp_pres_ratio * self.in_pres2

        w, T2= symbols("w, T2")
        h1 = 1/3000 * self.in_temp2 ** 3 + 287/10000 * self.in_temp2 ** 2 + 24394/25 * self.in_temp2
        energy = sp.Eq(w, 1/3000 * T2 ** 3 + 287/10000 * T2 ** 2 + 24394/25 * T2 - h1)
        isentropic = sp.Eq(T2/self.in_temp2, (self.comp_pres_ratio)**(self.R/(0.0001*((T2 + self.in_temp2)/2)**2 + 0.0574 * ((T2 + self.in_temp2)/2) + 975.76  )))
        sol = nsolve([energy, isentropic], [w,T2], [1, 1], prec = 50)

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

        T2,V2 = symbols("T2, V2")
        Q = sp.Eq(self.outlet_flow * (1/3000 * T2 ** 3 + 287/10000 * T2 ** 2 + 24394/25 * T2 - h1 + (V2**2-self.in_vel2**2)/2), 7901740 * self.fuel_flow / 0.17034 )
        massflow = sp.Eq(V2/T2, self.in_vel2/self.comp_temp )
        sol = nsolve([Q,massflow], [T2, V2], [self.comp_temp + 1000, self.in_vel2], prec = 50)
        
        self.comb_temp = sol[0]
        self.comb_heat = 7901740 * self.fuel_flow / 0.17034 
        self.comb_exit_vel = sol[1]
       

    def turbine(self):

       
        self.turb_pres2 = self.turb_pres_ratio * self.comp_pres

        w, T2 = symbols("w, T2")
        h1 = 1/3000 * self.comb_temp ** 3 + 287/10000 * self.comb_temp ** 2 + 24394/25 * self.comb_temp
        energy = sp.Eq(w, 1/3000 * T2 ** 3 + 287/10000 * T2 ** 2 + 24394/25 * T2 - h1)
        isentropic = sp.Eq(T2/self.comb_temp, (self.turb_pres_ratio)**(self.R/(0.0001*((T2 + self.comb_temp)/2)**2 + 0.0574 * ((T2 + self.comb_temp)/2) + 975.76  )))
        sol = nsolve([energy, isentropic], [w,T2], [1, 1], prec = 50)

        self.turb_work = sol[0] * self.outlet_flow
        self.turb_temp = sol[1]


        #find max turbine pressure ratio
        w, P2 = symbols("w, P2")
        h1 = 1/3000 * self.comb_temp ** 3 + 287/10000 * self.comb_temp ** 2 + 24394/25 * self.comb_temp
        energy = sp.Eq(w, 1/3000 * self.in_temp ** 3 + 287/10000 * self.in_temp ** 2 + 24394/25 * self.in_temp - h1)
        isentropic = sp.Eq(self.in_temp/self.comb_temp, (P2/ self.comp_pres)**(self.R/(0.0001*((self.in_temp + self.comb_temp)/2)**2 + 0.0574 * ((self.in_temp + self.comb_temp)/2) + 975.76  )))
        sol = nsolve([energy, isentropic], [w,P2], [1, 1], prec = 50)
        

    def nozzle(self):
        
        V2, T2 = symbols("V2, T2")

        h1 = 1/3000 * self.turb_temp ** 3 + 287/10000 * self.turb_temp ** 2 + 24394/25 * self.turb_temp
        massflowrate = sp.Eq(self.turb_pres2*self.noz_area_ratio*self.in_vel2/self.turb_temp, self.in_pres*V2/T2)
        energy = sp.Eq(0, 1/3000 * T2 ** 3 + 287/10000 * T2 ** 2 + 24394/25 * T2 - h1 + (V2 ** 2 - self.comb_exit_vel**2)/2)
        sol = nsolve([massflowrate, energy], [V2,T2] , [self.comb_exit_vel/self.noz_area_ratio, self.in_temp], prec = 50)
        
        self.noz_vel2 = float(sol[0])
        self.noz_temp2 = float(sol[1])
        self.noz_exit_area = self.outlet_flow*self.R*self.noz_temp2/(self.in_pres*self.noz_vel2)
        self.thrust = self.outlet_flow*(self.noz_vel2 - self.in_vel)

    def table(self):
        self.dictionary = {
            "----------------Environment Properties----------------":{
                "Air Velocity [m/s]": self.in_vel,
                "Ambient Temperature [°C]": round(self.in_temp -273),
                "Ambient Pressure [kPa]": round(self.in_pres/1000,2),
                "Air Gas Constant [kJ/kgK]": 0.287,

            },
            "----------------Intake/Diffuser----------------":{
                "Diameter [m]": round(math.sqrt(self.in_area*4/math.pi),2),
                "Area Ratio (1:x)": self.in_area_ratio,
                "Exit Velocity [m/s]": round(self.in_vel2,2),
                "Exit Pressure [kPa]": round(self.in_pres2/1000,2),
                "Exit Temperature [°C]": round(self.in_temp2 - 273),
                "Mass Flow Rate [kg/s]": "{:.2e}".format(self.mass_flow)
            },
            "----------------Compressor----------------":{
                "Pressure Ratio (1:x)": round(self.comp_pres_ratio),
                "Exit Pressure [kPa]": round(self.comp_pres/1000, 2),
                "Exit Temperature [°C]": round(self.comp_temp - 273),
                "Power Required [W]":  "{:.2e}".format(self.comp_work)
            },
            "----------------Combustor----------------":{
                "Air - Fuel Ratio (x:1)": round(self.mass_flow/self.fuel_flow), 
                "Fuel Mass Flow Rate [kg/s]": "{:.2e}".format(self.fuel_flow),
                "Dodecane ΔHcomb [kJ/kg]": round(7901.740 / 0.17034),
                "Heat Added [W]": "{:.2e}".format(self.comb_heat),
                "Exit Temperature [°C]": round(self.comb_temp - 273),
                "Exit Velocity [m/s]": round(self.comb_exit_vel,2),
                "Exhaust Mass Flow Rate [kg/s]": "{:.2e}".format(self.outlet_flow),
                "CO2 Emissions [kg/s]": "{:.2e}".format(self.emissions),
            },
            "----------------Turbine----------------":{
                "Pressure Ratio (1:x)": round(1/self.turb_pres_ratio),
                "Exit Pressure [kPa]": round(self.turb_pres2/1000, 2),
                "Exit Temperature [°C]": round(self.turb_temp - 273),
                "Power Generated [W]": "{:.2e}".format(-self.turb_work)
            },
            "----------------Outlet/Nozzle----------------":{
                "Area Ratio (x:1)": round(self.noz_area_ratio),
                "Outlet Area [m²]": "{:.2e}".format(self.noz_exit_area),
                "Exit Pressure [kPa]": round(self.in_pres/1000,2),
                "Exit Temperature [°C]": round(self.noz_temp2 - 273),
                "Exit Velocity [m/s]": round(self.noz_vel2,2),
                "Thrust [N]": "{:.2e}".format(self.thrust),
            },
            "----------------Power/Efficiency----------------":{
                "Net Power [W]": "{:.2e}".format(- self.net_work),
                "Thermal Efficiency [%]":  round(self.eff * 100)
            }
        }


    def work(self):
        self.inletOutput()
        self.compression()
        self.combustion()
        self.turbine()
        self.nozzle()

        self.net_work = self.turb_work + self.comp_work
        self.eff = -self.net_work/self.comb_heat

        self.table()

if __name__ == "__main__":
    A = Brayton(in_vel=10.0, in_temp = 293.0, in_pres = 101325.0, in_area= 0.7853981633974483, in_area_ratio= 1.0, comp_pres_ratio= 30.0, fuel_flow= 100, turb_pres_ratio= 1/60, noz_area_ratio=1/1)
    A.work()