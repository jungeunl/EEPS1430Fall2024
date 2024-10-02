import numpy as np

def satvap (T):

    '''
    Calculates saturation vapor pressure
    Input is temperature in K
    Returns saturation vapor pressure in Pa
    '''

    T0=273.15
    Tc=T-T0  # Note that the equation below uses temperature in Celcius
    es=6.11*np.exp(17.502*Tc/(Tc+240.97))*100. #Pa
    return es

def lapse_dry(T):

    '''
    Calculates dry adiabatic lapse rate 
    Input is temperature in K
    Returns lapse rate K/m
    '''

    Cp=1004.    # specific heat when pressure is constant J/K/kg
    g=9.81      # gravitation accelration

    gammad=g/Cp

    return gammad

def lapse_moist(T):

    '''
    Calculates moist adiabatic lapse rate at 1000 hPa
    Input is temperature in K
    Returns lapse rate K/m
    '''

    es=satvap(T) # Pa

    L=2.5e6   # Latent heat of vaporization J/kg
    Rd=287.     # gas constant for dry air
    Rv=461.5   # gas constant for water vapor

    eps= Rd/Rv

    ws=eps*es/(100000.-es)   # 100000. is 1000 hPa

    Cp=1004.    # specific heat when pressure is constant J/K/kg
    g=9.81      # gravitation accelration

    num = 1.+ L*ws/(Rd*T)
    denum = 1. + L**2*eps*ws/Rd/Cp/T**2 

    gammas=lapse_dry(T)*num/denum # moist adiabatic lapse rate
    return gammas

