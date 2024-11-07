def orbital_parameter(kyear):
  import numpy as np
# Calculate orbital parameter
# Required file: orbital_parameter_data.txt

# Load orbital parameters (given each kyr for 0-5Myr)===
  orbitParams = np.loadtxt('orbital_parameter_data.txt')
# print orbitParams.shape

  time0=-orbitParams[:,0]
  e0=orbitParams[:,1]
  pre0=orbitParams[:,2]
  obl0=orbitParams[:,3]
  e=np.interp(kyear,time0,e0,left=None,right=None)
  obl=np.interp(kyear,time0,obl0,left=None,right=None)
  pre=np.interp(kyear,time0,pre0,left=None,right=None)
  return obl,e,pre


