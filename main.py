# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import math
import sympy

#SLS format
M_sls = [88,49] # in KNm
N_sls = [35,29] # in KN
V_sls = [21.5,7] # in KN
SLS = {'M':M_sls,'N':N_sls,'V':V_sls}
M_uls = [120,18.5,76] # in KNm
N_uls = [30,-5,32] # in KN
V_uls = [28,2,12] # in KN
ULS = {'M':M_uls,'N':N_uls,'V':V_uls}
Dimensions = {'D':0.5,'W':1} # in metres
allowable_end_pressure = 150 # in KPa
density_concrete = 25 # in KN/m3
overburden = 0.2 * 18 # KPa from soil above (depth x unit weight)
Friction_angle = 24 # in degrees
cohesion = 5 # in KPa

def get_int(string):
    while True:
        try:
            integer = int(input(string))
            break
        except ValueError:
            print("Invalid input: please input integer")
            pass
    return integer
def overturning():
    #TODO: check if element passes overturning check
    return

def service_bearing():
    #TODO: check if element passes service bearing check
    return

def sliding():
    #TODO: check if element passes sliding check
    return

def main():
    #TODO take user inputs
    allowable = get_int("Input allowable end bearing pressure in kPa: ")
    width = get_int("Input set width in metres or type N: ")

    print(allowable)


    #TODO Add user inputs into database

    #TODO Check pad footing
    return


def pad_footing(SLS,ULS,Dimensions,allowable_end_pressure,density_concrete,overburden,Friction_angle,cohesion):
    #Steps are:
    #1. Check servicibility bearing pressures
    #2. Check overturning moment
    L2 = [0]* len(ULS['M'])
    for i in range(len(SLS['M'])):
        e = SLS['M'][i]/SLS['N'][i]
        L1 = 0
        L =  [0]
        selfweight = density_concrete*Dimensions['D']
        while True:
            if e < L[0]/6:
                x = sympy.Symbol('x')
                L = sympy.solve(sympy.Eq((SLS['N'][i] + selfweight * L[0] * Dimensions['W']+ overburden * L[0] * Dimensions['W'])/(Dimensions['D']*x) + (6*SLS['M'][i])/(Dimensions['W']*x**2) , allowable_end_pressure)) #Round L to nearest 0.1 afterwards
                L[0] = L[1]

                e = SLS['M'][i]/(SLS['N'][i]+ selfweight * L[0] * Dimensions['W']+ overburden * L[0] * Dimensions['W'])
                a = L[0]
            elif e >= L[0]/6:
                x = sympy.Symbol('x')
                L = sympy.solve(sympy.Eq(2*(SLS['N'][i] + selfweight * L[0] * Dimensions['W']+ overburden * L[0] * Dimensions['W'])/(Dimensions['W']*(3*(x/2 - e))),allowable_end_pressure))

                e = SLS['M'][i] / (SLS['N'][i] + selfweight * L[0] * Dimensions['W']+ overburden * L[0] * Dimensions['W'])
                a = 3*(L[0]/2 - e)
            if L1 == round(L[0],3):
                print('L > '+ str( round(L[0],2)))
                L2[i] = L1
                break
            elif L1 != round(L[0],3):
                print(L[0])
                L1 = round(L[0],3)
    L1 = max(L2)
    Ast = [0]*len(ULS['M'])
    for i in range(len(ULS['M'])):
        while True:
            #Check overturning moments
            e = ULS['M'][i] / (
                    ULS['N'][i] + selfweight * L1 * Dimensions['W'] + overburden * L1 * Dimensions['W'])
            selfweight = density_concrete * Dimensions['D']
            if e < L1 / 6:
                    x = sympy.Symbol('x')
                    Qmax = sympy.solve(sympy.Eq((ULS['N'][i]+selfweight * L1 * Dimensions['W'] + overburden * L1 * Dimensions['W'])/(Dimensions['W']*L1) + 6*ULS['M'][i]/(Dimensions['W']*L1**2),x))  # Round L to nearest 0.1 afterwards
                    a = L1
            elif e >= L1 / 6:
                    x = sympy.Symbol('x')
                    Qmax = sympy.solve(sympy.Eq(
                        2 * (ULS['N'][i] + selfweight * L1 * Dimensions['W'] + overburden * L1 * Dimensions['W']) / (
                                    Dimensions['W'] * (3 * (L1 / 2 - e))), x))
                    a = 3 * (L1 / 2 - e)
            Moment_check = ULS['M'][i] - (
                        ULS['N'][i] + selfweight * L1 * Dimensions['W'] + overburden * L1 * Dimensions['W']) * L1 / 2
            if Qmax[0] > 2*allowable_end_pressure:
                L1 += 0.1

            elif Moment_check > 0:
                print('Overturning Fail', Moment_check)
                L1 += 0.1
            elif Moment_check < 0 and Qmax[0] < 2*allowable_end_pressure:
                print('Passes overturning check',Moment_check)
                L2[i] = L1
                break
        #Check sliding failure
        F = (ULS['N'][i] + selfweight * L1 * Dimensions['W']+ overburden * L1 * Dimensions['W'])*math.tan(Friction_angle*0.75*math.pi/180) + cohesion*a
        if F > ULS['V'][i]:
            print('Passes Shear')
        elif F < ULS['V'][i]:
            print('Fails Shear')
        #Check required reinforcement for toe
        Qmax = Qmax[0]
        if a > L1/2:
            max_moment = ((0.5*a*Qmax) - (a - L1/2)*Qmax*0.5) * (0.5*(Qmax - 0.5*L1/a * Qmax)*(0.5*L1)**2*2/3 + (0.5*L1)**2/(2*a)*Qmax)/((0.5*L1)**2/a*Qmax + 0.5*(Qmax - 0.5*L1/a*Qmax)*0.5*L1)
            max_shear = ((0.5*a*Qmax) - (a - L1/2)*Qmax*0.5)
        elif a <= L1/2:
            max_moment = 0.5*a*Qmax * (L1/2 - 1/3*a)
            max_shear = 0.5*a*Qmax
        #Determine required Area of steel
        fc = 32
        alpha = max(0.85 - 0.0015*fc,0.67)
        Ast[i] = sympy.solve(sympy.Eq(x*500*(Dimensions['D']- 0.075 - ((x*500)/(alpha*Dimensions['W']*32*10**6))/2),max_moment/0.65*1000))
        print(Ast,'Area of steel required')
        #Check required reinforcement for heel
        max_moment_heel = (selfweight * Dimensions['W'] + overburden * Dimensions['W']) * (L1/2)**2/2
        max1 = sympy.solve(sympy.Eq(x*500*(Dimensions['D']- 0.075 - ((x*500)/(alpha*Dimensions['W']*32*10**6))/2),max_moment_heel))
        print(max1,'max heel')
        #Check Shear requirements

    #Check minimum requirements
    Pmin = 0.19*(Dimensions['D']/(Dimensions['D']-0.075))**2*(0.6*math.sqrt(fc)/500)*Dimensions['D']*Dimensions['W']*10**6
    print(Pmin,'Pmin')
    print(round(L1,2))

#pad_footing(SLS,ULS,Dimensions,allowable_end_pressure,density_concrete,overburden,Friction_angle,cohesion)

if __name__ == "__main__":
    main()