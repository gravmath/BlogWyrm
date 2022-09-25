#Compressible Fluid Flow Utility Functions
import numpy as np

#Basic thermodynamic functions anchored to the stagnation resevoir
#see http://underthehood.blogwyrm.com/?p=1402 for formulae
def P0_P(M,gamma):
    body = T0_T(M,gamma)
    return body**(gamma/(gamma-1.0))

def P_P0(M,gamma):
    return 1.0/P0_P(M,gamma)

def P_Pstar(M,gamma):
    return P_P0(M,gamma)*P0_P(1.0,gamma)

def rho0_rho(M,gamma):
    body = T0_T(M,gamma)
    return body**(1.0/(gamma-1.0))

def rho_rho0(M,gamma):
    return 1.0/rho0_rho(M,gamma)

def rho_rhostar(M,gamma):
    return rho_rho0(M,gamma)*rho0_rho(1.0,gamma)

def T0_T(M,gamma):
    return 1.0 + (gamma - 1.0)/2.0*M**2

def T_T0(M,gamma):
    return 1.0/T0_T(M,gamma)

def T_Tstar(M,gamma):
    return T_T0(M,gamma)*T0_T(1.0,gamma)

def AMR(M,gamma):
    body  = T0_T(M,gamma)
    term1 = (1.0 + (gamma-1.0)/2.0)**(0.5*(gamma+1.0)/(1.0-gamma))
    term2 = body**(0.5*(gamma+1.0)/(gamma-1.0))
    return term1*term2/M

def inv_AMR(A_Astar,gamma,range):

    #setup bisection
    counter = 0
    tol     = 1e-6
    if range == 'subsonic':
        M_a   = 0.0001
        M_b   = 0.9999
    elif range == 'supersonic':
        M_a = 1.0001
        M_b = 1000.0
    AAS_a  = AMR(M_a,gamma)
    AAS_b  = AMR(M_b,gamma)
    if AAS_a - A_Astar > 0.0:
        sign_a = 1
    else:
        sign_a = -1
    
    while(counter < 40):
        M_t    = (M_a+M_b)/2.0
        AAS_t  = AMR(M_t,gamma)
        if AAS_t - A_Astar > 0:
            sign_t = 1
        else:
            sign_t = -1
        if abs(A_Astar-AAS_t) < tol:
            break
        elif sign_t == sign_a:
            M_a  = M_t
        else:
            M_b = M_t
        counter += 1
    return M_t, counter

#Normal shock relations with u = upstream and d = downstream
def Mach_jump(Mu,gamma):
    coeff = (gamma - 1.0)/2.0
    num = 1.0 + coeff*Mu**2
    den = gamma*Mu**2 - coeff 
    return np.sqrt(num/den)

def Pd_Pu(Mu,gamma):
    return 1.0 + (2.0*gamma)/(gamma+1.0)*(Mu**2 - 1.0)

def rhod_rhou(Mu,gamma):
    num = (gamma + 1.0)*Mu**2
    den = 2.0 + (gamma - 1.0)*Mu**2
    return num/den

def Td_Tu(Mu,gamma):
    coeff = Pd_Pu(Mu,gamma)
    return coeff/rhod_rhou (Mu,gamma)

#cases checked on 9/21/2022 against the VT compressible fluid flow calculator
if __name__ == '__main__':
    print('**************************')
    gamma = 1.4
    M     = 0.3
    print('Case 1: gamma = 1.4 & M = 0.3')
    print('P/P0:',P_P0(M,gamma),   'rho/rho0:',rho_rho0(M,gamma),   'T/T0:',T_T0(M,gamma))
    print('P/P*:',P_Pstar(M,gamma),'rho/rho*:',rho_rhostar(M,gamma),'T/T*:',T_Tstar(M,gamma),'A/A*:',AMR(M,gamma))
    print('Confirmed against the VT Compressible Flow Calculator')
    print('**************************')    
    gamma = 1.4
    M     = 2.3
    print('Case 2: gamma = 1.4 & M = 2.3')
    print('P/P0:',P_P0(M,gamma),   'rho/rho0:',rho_rho0(M,gamma),   'T/T0:',T_T0(M,gamma))
    print('P/P*:',P_Pstar(M,gamma),'rho/rho*:',rho_rhostar(M,gamma),'T/T*:',T_Tstar(M,gamma),'A/A*:',AMR(M,gamma))
    print('Confirmed against the VT Compressible Flow Calculator')    
    print('**************************')    
    gamma = 1.4
    M     = 1.1
    print('Case 3: gamma = 1.4 & M = 1.1; Normal Shock Conditions')
    print('M1:',M,'M2',Mach_jump(M,gamma))
    print('P2/P1:',Pd_Pu(M,gamma),'rho2/rho1:',rhod_rhou(M,gamma),'T2/T1:',Td_Tu(M,gamma))
    print('Confirmed against the VT Compressible Flow Calculator')    
    print('**************************')    

    gamma   = 1.4
    A_Astar = 2.03506526
    print(inv_AMR(A_Astar,gamma,'subsonic'))
    A_Astar = 2.19313081
    print(inv_AMR(A_Astar,gamma,'supersonic'))

    import numpy as np
    import matplotlib.pyplot as plt

    M    = np.arange(0.01,7,0.01)
    A_As = AMR(M,1.4)

    plt.plot(A_As,M)
    plt.show()
