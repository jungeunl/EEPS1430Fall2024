#Library of functions for computing solar flux distribution for
#circular and elliptical orbits
#
#ToDo:  
#       Phi is sometimes used for latitude, in some argument lists.
#
#       The function names could use some improvement
#
#       Check: The annual mean flux seems independent of
#       precession, even for high eccentricity. Is this right?

from ClimateUtilities import *
#----------------------------------------------------------------
#Basic functions to get latitude of the Sun, zenith angle, and
#so forth. Input angles are all in radians

from math import *
#Latitude of the Sun
def sunlat( season, gamma):
	return asin( cos(season)*sin(gamma))
#Hour angle giving the terminator position
def h( phi, delta):
	temp = -sin(phi)*sin(delta)/(cos(phi)*cos(delta));
	if temp < -1. :
		return pi
	elif temp > 1.:
		 return 0.
	else:
		 return acos(temp)

#Flux factor as a function of latitude and lat of the Sun
def flux( phi, delta):
	return (2.*sin(h(phi,delta))*cos(phi)*cos(delta) +
	     2.*h(phi,delta)*sin(phi)*sin(delta))/(2.*3.14159)

#  A convenience function, with a calling sequence simpler to use
#  but otherwise, it does the same as the function flux(phi,delta)
def solar(phi,season,obliquity):
	return flux(phi,sunlat(season,obliquity))

# The following function computes the annual
# mean flux, but it is only valid for a circular orbit
#ToDo: Replace integration with Romberg integrator from ClimateUtilities
def AnnMeanFluxCirc( phi, obliquity):
	avg = 0.
	n = 180 # n should be even for symmetry
	ds = 2.*pi/n
	for season in [ds/2. + i*ds for i in range(n)]:
		avg += flux(phi,sunlat(season,obliquity))
	return avg/180.

#----------------------------------------------------------------
#
#Functions for computing fluxes for an elliptical orbit
#
#Angular velocity function. e is the eccentricity,
#J is the orbital angular momentum and a is the
#semi-major axis
def kappadot(t,kappa,params):
    e = params.e
    J = params.J
    a = params.a
    return (J/a**2)*(1+e*cos(kappa))**2/(1-e**2)**2

#Distance as function of semi-major axis and angle
def r(a,e,kappa):
    return a*(1.-e*e)/(1.+e*cos(kappa))

#Compute kappa(t)
#    Return Curve object with t,kappa(t),r(t)
#    Distance is normalized so that semi-major
#    axis = 1, and time such that angular momentum = 1. 
#yday=number of days per year
def orbit(e,yday,kappaStart=0.):
    #Set up the parameter object
    params = Dummy()
    params.e = e
    params.J = 1.
    params.a = 1.
    #Set up the integrator and install parameter object
    dt = 2.*pi/yday
    m = integrator(kappadot,0.,kappaStart,dt)
    m.setParams(params)
    tList = []
    rList = []
    kappaList = []
    kappa = kappaStart
    while kappa < 2.*pi +kappaStart:
        t,kappa = m.next()
        tList.append(t)
        kappaList.append(kappa)
        rList.append(r(params.a,params.e,kappa))
    #
    #Normalize time here?
    #Stuff it into a Curve() and return
    c= Curve()
    c.addCurve(tList,'t','Time')
    c.addCurve(kappaList,'kappa','Season angle')
    c.addCurve(rList,'r','Distance from Sun')
    return c

#Return an array giving the flux factor as a function
#of time, for the specified orbit (returned by the orbit() function)
def orbitFlux(lat,obl,precess,orb):
    L = 1./orb['r']**2 #Uses array arithmetic
    fluxes =[
        flux(lat,sunlat(kappa-precess,obl)) for kappa in orb['kappa']]
    return Numeric.array(fluxes)*L

#Return the annual mean flux for a given orbit
def AnnMeanFlux(lat,obl,precess,orb):
    flux = orbitFlux(lat,obl,precess,orb)
    #Integrate to get annual average (trapezoidal rule)
    sum = 0.
    for i in range(len(flux)-1):
        sum += .5*(flux[i+1]+flux[i])*(orb['t'][i+1]-orb['t'][i])
    sum = sum/(orb['t'][-1]-orb['t'][0])
    return sum
    
#Function to make a map of the flux vs latitude and time
#Input angles are in degrees
#yday is number of days per year
def fluxmap(e,obl,precess,lats,yday):
    degrad = pi/180.
    fluxes = []
    orb = orbit(e,yday,precess*degrad)#Start at the solstice, to make the graphs nicer
    fluxes = [orbitFlux(lat*degrad,obl*degrad,precess*degrad,orb) for lat in lats]
    return numpy.array(fluxes)
