# Notes
# - CHECK substrate absorbtion. Comment and value don't agree.
# - Which of the 4 effects contribute to HR-RoC change and thermal lensing?
# - Is the expression correct given how the input powers are given?
# - Is it the path length or physical distance that should be used in the end to
#   model the new effective RoCs?

##################################################
# To Valeria
##################################################
# - Check and comment the function hellovinet() down until the line where 

import numpy as np
import pykat
import matplotlib.pyplot as plt
import scipy.special as sp
import scipy.optimize as so
import pylab
from scipy.optimize import minimize
import pykat.exceptions as pkex
    

def hellovinet(P_coat, P_sub_in, P_sub_out, HR_RoC = None, AR_RoC = None, HR_zOff=None, AR_zOff=None, **kwargs):
    """
    Computes the effective change in optical path length in a mirror due to thermal lensing
    and suface deformation effects from the power and absorbtion coefficients.
    
    If the initial HR_RoC and/or AR_RoC are specified as inputs, spherical surfaces are
    fitted to the new deformed surfaces, and the resulting change in RoC is extraced. If also
    the initial z-offsets are specified, both the RoCs and offsets are fitted.
    
    Based on Hello & Vinet, J. Phys. France 51 (1990) 1267-1282. Function written by Valeria
    Sequino and Daniel Toyra, translated from Matlab function by ?.


    Power definitions:
    ------------------
       
           Mirror
            
    AR               HR 
     |  P_sub_in      |      
     |  ---------->   |     P_coat
     |                |   <----------  
     |  P_sub_out     |         
     |  <----------   |    
     |                |      


    RoC and z-offset signs and definitions:
    ---------------------------------------
    If r is the transverse coordinate, and if z (optical axis) increases in the direction AR --> HR,
    then the z-axis passes through the mirror surface (AR or HR, they are fitted independently)
    at the point (r, z) = (0, zOffset), and the surface's center of curvature is located at
    (r, Z) = (0, RoC + zOffset). E.g., for Advanced Virgo IMs, both HR and AR surfaces have positive RoC.
    
    Inputs:
    -------
    P_coat            - Power on HR coating from vacuum side. See figure above.
    P_sub_in          - Power propagaint from AR-surface to HR-surface. See figure above.
    P_sub_out         - Power propagating from HR-surface to AR-surface. See figure above. 
    HR_RoC            - Initial RoC of the HR surface [m]. If specified, the new RoC will be
                        fitted to the deformed surface. For sign, see above. 
    AR_RoC            - Same as for HR_RoC, but for the AR-surface.
    HR_zOff           - Initial guess for the new HR-surface's offset along the optical axis [m].
                        If given a numeric value, the offset will be fitted along with the RoC,
                        which affects the RoC fitting. Set this to None (no fitting) or 0 (the
                        offset should always be small).
    AR_zOff           - Same as for HR_zOff above, but for the AR-surface.

    mirror_properties - Dictionary where all the below properties can be specified directly.
    
    thickness         - Mirror thickness [m] (default 0.2 m)
    aCoat             - Coating power absorption (default 2.5 ppm)
    aSub              - Substrate power absorptio [1/m] default (30 ppm/m)
    n                 - Mirror index of refraction (default SiO2, 1.452).
    a                 - Mirror radius [m] (default 0.175 m)
    w                 - Gaussian spot size [m] (default 0.049 m)
    K                 - Mirror thermal conductivity (default 1.380)
    T0                - Surrounding temperature [K] (default 295.0 K)
    emiss             - Mirror emissivity (default 0.89)
    alpha             - Mirror thermal expansion coefficient (default 0.54e-6)
    sigma             - Mirror Poisson's ratio (default 0.164)
    dndT              - Mirror index of refraction change with temperature (default 8.7 ppm)
    N                 - Number of data points along the mirror radius a (default 176).
    nScale            - If set to true, the optical path length data is scaled to physical distance.

    Returns:
    --------

    r              - Array with distances from the optical axis [m]. 0 <= r <= a
    fitData        - List that contains the RoC and z-offsets for the new fitted
                     spherical surfaces, if the initial values for these parameters
                     have been given as inputs.
                     fitData = [ np.array(HR RoC, HR z-offset),
                                 np.array(AR RoC, AR z-offset) ]
    oplData        - List containing the four contributions to extra optical path length [m]
                     as functions of r.
                     oplData = [ coating thermo-optic,   coating thermo-elastic,
                                 substrate thermo-optic, substrate thermo-elastic ]
    """ 


    #############################################
    # Default values
    #############################################    
    sigmab = 5.670367e-8   # Stephan-Boltzmann constant [W m**(-2) K**(-4)]

    h=0.2          # test mass thickness [m]
    aCoat=2.5e-6   # coating absorption 
    aSub=3.0e-5    # substrate absorption  [1/m] [TDR, table 2.6]
    n=1.452        # SiO2 refraction index
    a=0.175        # test mass radius [m]
    w=0.049        # Spot size at Virgo input mirrors
    K=1.380        # SiO2 thermal conductivity
    T0=295.0       # room temperature
    emiss=0.89     # SiO2 emissivity
    alpha=0.54e-6  # SiO2 thermal expansion coefficient
    sigma=0.164    # SiO2 Poisson's ratio
    dndT=8.7e-6    # Mirror index of refraction change with temperature (default 8.7 ppm)
    N = 176        # Number of data points along the mirror radius
    scale = 1.0    # Scales path-length data before fitting RoCs
    nScale = False # If true, scale is set to 1/n

    #Pin=125;       # IFO input power
    #RecG=37.5;     # recycling gain
    #Finesse=443;   # F-P cavity finesse not used



    # Updating values specified as a dictionary. 
    if 'mirror_properties' in kwargs:
        for k,v in kwargs['mirror_properties'].items():
            if k == 'aCoat':
                aCoat = v
            elif k == 'thickness':
                h = v
            elif k == 'aSub':
                aSub = v
            elif k == 'n':
                n = v
            elif k == 'a':
                a = v
            elif k == 'w':
                w = v
            elif k == 'K':
                K = v
            elif k == 'T0':
                T0 = v
            elif k == 'emiss':
                emiss = v
            elif k == 'alpha':
                alpha = v
            elif k == 'sigma':
                sigma = v
            elif k == 'dndT':
                dndT = v
            elif k == 'N':
                N = v
            elif k == 'nScale':
                nScale = v

    # Updating values specified as arguments. These have precedence
    # over the parameters specified in the dictionary. 
    for k, v in kwargs.items():
        if k == 'aCoat':
            aCoat = v
        elif k == 'thickness':
            h = v
        elif k == 'aSub':
            aSub = v
        elif k == 'n':
            n = v
        elif k == 'a':
            a = v
        elif k == 'w':
            w = v
        elif k == 'K':
            K = v
        elif k == 'T0':
            T0 = v
        elif k == 'emiss':
            emiss = v
        elif k == 'alpha':
            alpha = v
        elif k == 'sigma':
            sigma = v
        elif k == 'dndT':
            dndT = v
        elif k == 'N':
            N = v
        elif k == 'nScale':
            nScale = v

    

    # Scales the optical path length data to physical distances
    if nScale:
        scale = 1.0/n

    # Total power going into the coating
    Pc = P_coat + P_sub_in

    # Total power going through the substrate
    Ps = P_sub_in + P_sub_out

    ######################################################################################
    # To Valeria:
    # ----------
    # If it's possible for you, it would be great if you could check, correct and comment
    # what's happening in the various steps below.
    ######################################################################################

    # Radiative heat losses
    chi=4*emiss*sigmab*T0**3*a/K 

    # What is this function? General solution to differential equation I think... Reference?
    def g(x):
        return x*sp.jv(1,x) - chi*sp.jv(0,x)

    # Initial root guesses of g(x)
    i = np.linspace(1,51,51)
    x0s = (i-1.0+1.0/4.0)*np.pi
        
    # Finding roots of g(x)
    zetha = np.zeros(len(x0s))
    for k,x0 in enumerate(x0s):
        out = so.fsolve(g, x0, xtol = 1e-8, full_output = True, maxfev = 1000)
        # Printing message if no root found.
        if out[2] != 1:
            print(out[3])
            print(x0)
        zetha[k] = float(out[0])

    # What are these 5 lines? Comments? Reference?
    gamma = h/(2.0*a)*zetha
    A = 1.0/(2.0*(zetha*np.sinh(gamma)+chi*np.cosh(gamma)))
    B = 1.0/(2.0*(zetha*np.cosh(gamma)+chi*np.sinh(gamma)))
    beta = 1.0/8.0*w**2/a**2*zetha**2
    p = 1.0/(np.pi*a**2)*zetha**2/((zetha**2+chi**2)*(sp.jv(0,zetha))**2)*np.exp(-beta)

    # Array with distances from the optical axis
    r = np.linspace(0, a, N)
    # Unused? Remove?
    z = np.linspace(-h/2, h/2, int(np.round(h/0.001))+1)

    # What are these? Comments? Reference?
    oos = np.zeros([len(r), len(zetha)])
    ooc = np.zeros([len(r), len(zetha)])
    for i in range(len(r)):  
        for k in range(len(zetha)):
                oos[i,k] = p[k]/zetha[k]**2*(1-2*chi*A[k]/gamma[k]*np.sinh(gamma[k]))*sp.jv(0,zetha[k]*r[i]/a)
                ooc[i,k] = p[k]/zetha[k]*2*A[k]*np.sinh(gamma[k])*sp.jv(0,zetha[k]*r[i]/a)

    # Thermo-optic effect due to coating absorption
    OPLc=Pc*aCoat*a**2/K*dndT*np.sum(ooc,1)
    # Thermo-optic effect due to substrate absorption
    OPLs=Ps*aSub*h*a**2/K*dndT*np.sum(oos,1)
    # Thermo-elastic effect due to coating absorption
    OPLtec=alpha*(sigma+1)*(n-1)*Pc*aCoat*a**2/K*np.sum(ooc,1)
    # Thermo-elastic effect due to substrate absorption
    OPLtes=alpha*(sigma+1)*(n-1)*Ps*aSub*h*a**2/K*np.sum(oos,1)
    # Total thermal effect
    OPLTM=OPLc+OPLs+OPLtec+OPLtes    

    ######################################################################################
    # To Valeria:
    # ----------
    # I take responsibility for everything that below this line, so you don't
    # need to check anything below here. /Daniel
    ######################################################################################
    
    oplData = [OPLc, OPLtec, OPLs, OPLtes]

    # In Finesse we want to separate this into two effects:
    # 1. Effective HR-surface deformation seen by the intra cavity field.
    # 2. Effective thermal lens seen by beams passing through the substrate.
    # Question: How do we achieve this?
    # Guess: I think only the thermo-elastic effect due to coating absorption
    #        should affect 1. The other three contributions gets added as a
    #        thermal lens.

    # Setting HR-deformation
    OPL_HR = OPLtec*scale
    # Setting AR-thermal
    OPL_AR = (OPLs + OPLtes + OPLc)*scale

    fitData = []
    # Computing new RoC for HR-surface, if initial HR_RoC was given. 
    if not HR_RoC is None:
        # Creating initial mirror surface
        Z_HR0 = createSurface(r, HR_RoC, HR_zOff)
        # Adding the distortion
        Z_HR1 = Z_HR0 + OPL_HR
        # Fitting
        out = fit_circle(r, Z_HR1, Rc0 = HR_RoC, zOff0 = HR_zOff, w = w)
        fitData.append(out)
        #if isinstance(out, list):
        #    HR_RoC1 = out[0]
        #    HR_zOff1 = out[1]
        #else:
        #    HR_RoC1 = out

    # Computing new effective RoC for AR-surface, if initial AR_RoC was given.
    if not AR_RoC is None:
        ######################
        # AR-surface
        ######################
        # Creating initial mirror surface
        Z_AR0 = createSurface(r, AR_RoC, AR_zOff)
        # Adding the distortion
        Z_AR1 = Z_AR0 - OPL_AR
        # Fitting
        out = fit_circle(r, Z_AR1, Rc0 = AR_RoC, zOff0 = AR_zOff, w = w)
        fitData.append(out)

    return r, fitData, oplData



def createSurface(r,Rc,zOffset=None):
    '''
    Creating circular surface with radius of curvature Rc. The z-axis (optical axis)
    passes through the surface at the point (r, Z) = (0, zOffset), and the circles
    center of curvature is located at (r, Z) = (0, Rc+zOffset). 

    Inputs
    ------
    Rc      - Radius of curvature, and center of sphere on z-axis in case zOffset=0 [m].
    r       - Array with distanecs from the optical axis [m].
    zOffset - Surface center offset [m].

    Returns:
    -------
    Z       - Array of surface positions [m], Z(X). 
    '''

    if Rc == np.inf:
        Rc = None
    elif Rc == 0:
        raise pkex.BasePyKatException("Rc cannot be 0")
    elif np.abs(Rc) < np.abs(r).max():
        raise pkex.BasePyKatException("abs(Rc) must be >= max(abs(X))")

    # Adding offset
    if zOffset is None:
        Z = 0
    else:
        Z = zOffset

    # Adding spherical shape.
    if not Rc is None:
        Z = Z + Rc - np.sign(Rc)*np.sqrt(Rc**2 - r**2)
    else:
        Z = Z*np.ones(len(r), dtype=float)

    return Z
    
def fit_circle(r, z, Rc0=None, zOff0=None, w=None, maxfev=2000):

    '''
    Fits circle segment to data. By default, only the radius of curvature is fitted,
    the circle-segment crosses the z-axis (optical axis) at z=0, and the center of
    curvature is at (r, z) = (0, Rc). 

    If zOff0 is set, the offset along the z-axis is fitted as well, the circle-segment
    crosses the z-axis (optical axis) at z = zOff, and the center of curvature is at
    (r, z) = (0, Rc + zOff).

    If w is set, gaussian weights are used.

    Inputs:
    -------
    r        - Array with distances away from the z-axis [m].
    z        - Array with data [m]
    Rc0      - Initial guess of the radius of curvature [m]. 
    w        - Gaussian weighting parameter. The distance from the z-axis where the
               weigths have decreased by a factor of exp(-2) compared to on-axis.
               Should normally be equal to the beam radius at the mirror [m].
    zOff0    - Initial guess of the z-offset [m]. Fits the z-offset if set. 
    
    Returns:
    --------
    Rc       - Radius of curvature of the fitted circle-segment [m].
    zOff     - z-offset of the circle segment [m]. Only returned if zOff0 was set. 
    '''
    
    # Initial guesses of radius of curvature and z-offset.
    p = []
    if Rc0 is None:
        p.append(0)
    else:
        p.append(Rc0)
    if not zOff0 is None:
        p.append(zOff0)

    
    # Cost-function to minimize.
    def costFunc(p, zOff=0):
        Rc = p[0]
        if len(p) == 2:
            zOff = p[1]
        Z = createSurface(r,Rc,zOff)
        if w is None:
            # Mean squared difference between map and the created sphere.
            res = np.sqrt( ((Z - z)**2).sum() )/(len(z)*np.mean(z))
        else:
            # Weights centered around the center of the mirror xy-plane.
            weight = np.exp(-2*r**2/w**2)
            # Weighted mean squared difference between map and the created sphere.
            res = np.sqrt( ( weight*( (Z - z)**2 )).sum() )/np.abs(weight*z).sum()
            return res

    opts = {'xtol': 1.0e-5, 'ftol': 1.0e-9, 'maxiter': 10000, 'maxfev': maxfev, 'disp': False}

    if len(p) < 2:
        # Removing z-offset from z-array if z-offset is not being fitted
        zOff0 = z[np.abs(r).argmin()]
        z = z - zOff0
        
    out = minimize(costFunc, p, method='Nelder-Mead', options=opts)

    if not out['success']:
        msg = '  Warning: ' + out['message'].split('.')[0] + ' (nfev={:d}).'.format(out['nfev'])

    return out['x']


