"""
---------------------------------------------------------
Example that shows how to:
1. Read the file formats below.
- metroPro (Binary files ending with .dat. Used for the
  current LIGO mirrors.)
- zygo (both zygo files ending with .asc and .xyz)
- ligo (files ending with _asc.dat)
- finesse (map files written by simtools or pykat, usually
  ends with _finesse.txt)
2. Use the standard map preparation procedure.
3. Create a procedure of your own by using existing
   surfacemap methods.
---------------------------------------------------------
"""

from pykat.optics.maps import *
import pylab

# --------------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------------
def read_all():
    # Reading metroPro map
    smap = read_map('ETM08_S1_-power160.dat', mapFormat='metroPro')
    smap.plot()
    # Reading Zygo xyz-map
    smap1 = read_map('ITM04_XYZ.xyz', mapFormat='zygo')
    smap1.plot()
    # Reading Zygo ascii-map
    smap2 = read_map('CVI1-S1-BC.asc', mapFormat='zygo')
    smap2.plot()
    # Reading ligo ascii-map
    smap3 = read_map('ETM02-S2_asc.dat', mapFormat='ligo')
    smap3.plot()
    # Writing map to file
    filename = smap.name + '_finesse.txt'
    smap.write_map(filename)
    # Reading finesse map (that was just written to file)
    smap4 = read_map(filename,mapFormat='finesse')
    smap4.plot()

def autoProcess():
    # Reading map
    smap = read_map('ETM08_S1_-power160.dat', mapFormat='metroPro')
    # Uncomment next line to crop the outermost 0.5 cm of the mirror. Run
    # without cropping first. To replicate the mirror map figueres in the
    # LIGO figure measurement reports, use: smap.crop(0.0802)
    '''
    smap.crop(0.155)
    '''
    # Showing unprocessed map
    smap.plot()
    # The standard automatised preparation procedure is about to take place.
    # The map is cropped, and the curvature, offset, and tilts are removed,
    # in that order.
    # There are three parameters that can be set here though. First, and
    # most important, the weighting parameter w. If not set (or set to None),
    # the whole map is convolved with Zernike-polynomials. If set, surfaces
    # are fitted to the map instead by gaussian weighting with the radius [m]
    # given in w. Splitting up the map into two so both cases can be
    # demonstrated.
    smap2 = deepcopy(smap)
    # xyOffset should be set if the beam spot is off the mirror centre. Doesn't
    # affect the processing, just sets the origo off the mirror centre
    # afterwards. Verbose should be set if we want extra information about the
    # processing, which we want here.
    amap = smap.preparePhaseMap(w=None, xyOffset=None, verbose=True)
    smap.plot()
    amap2 = smap2.preparePhaseMap(w=0.062, xyOffset=None, verbose=True)
    smap2.plot()
    # Plotting the aperture map created in the preparePhaseMap procedure.
    amap.plot()

def manualProcess():
    # Reading map and plotting
    smap = read_map('ETM08_S1_-power160.dat', mapFormat='metroPro')
    smap.plot()
    # If we want to get rid of the outher edge, use the crop-method.
    # Doing this here just because the figure gains resolution. To
    # replicate the mirror maps in the official figure measurement
    # reports, use: smap.crop(0.0802)
    
    # smap.crop(0.155)
    # smap.crop(0.0802)
    smap.crop(0.062)
    
    smap.plot()
    # Recentering is useful. To show the effect of xyOffset first:
    smap.xyOffset = (0.02,0.05)
    # In the figure the mirror centre should be at (-0.02, 0.05).
    smap.plot()
    # And now we recentering. The origo is again at the mirror centre.
    smap.recenter()
    smap.plot()# Splitting into two versions: One where we process by convolving the mirror surface
    # with Zernike-polynomials, and one where we process by fitting surfaces to the
    # mirror. In the latter case Gaussian weighting is used to make the centre of the
    # mirror more important.

    # To fit sphere, use this line:
    w = 0.062
    # To use Zernike-polynomials use:
    '''
    w = None
    '''
    if w is None:
        '''
        Using Zernike polynomials to process the map.
        '''
        # Removing curvature by using the parabolic Z(2,0) mode. Rc is the spherical
        # radius of curvature (the parabola is approximately spherical at the centre).
        # znm is dictionary containg the zernike polynomials and amplitudes that have
        # been removed.

        Rc, znm = smap.remove_curvature(method='zernike', zModes = 'defocus')
        
        # In case we want to remove astigmatism instead, use:
        '''
        Rc, znm = smap.remove_curvature(method='zernike', zModes = 'astigmatism')
        '''
        # In case we want to remove both defocus and astigmatism, use:
        '''
        Rc, znm = smap.remove_curvature(method='zernike', zModes = 'all')
        '''

        # Printing stuff
        # --------------------------
        print('Curvatures removed')
        if len(smap.zernikeRemoved)==1:
            print(' Rc = {:.1f} m,  Zernike(m,n,amp) = ({:d}, {:d}, {:.2f} nm)'.
                  format(Rc,znm['02'][0],znm['02'][1],znm['02'][2]))
        elif len(smap.zernikeRemoved)==2:
             print(' ROC(min,max) = ({:.1f},{:.1f}) m'.format(Rc[0],Rc[1]))
             print(' Zernike(m,n,amp) = ({:d}, {:d}, {:.2f} nm)'.format(znm['-22'][0],znm['-22'][1],znm['-22'][2]))
             print(' Zernike(m,n,amp) = ({:d}, {:d}, {:.2f} nm)'.format(znm['22'][0],znm['22'][1],znm['22'][2]))
        elif len(smap.zernikeRemoved)==3:
             print(' ROC(min,max) = ({:.1f},{:.1f}) m'.format(Rc[0],Rc[1]))
             print(' Zernike(m,n,amp) = ({:d}, {:d}, {:.2f} nm)'.format(znm['02'][0],znm['02'][1],znm['02'][2]))
             print(' Zernike(m,n,amp) = ({:d}, {:d}, {:.2f} nm)'.format(znm['-22'][0],znm['-22'][1],znm['-22'][2]))
             print(' Zernike(m,n,amp) = ({:d}, {:d}, {:.2f} nm)'.format(znm['22'][0],znm['22'][1],znm['22'][2]))
        # ----------------------------
        smap.plot()
        
        # Removing offset by convolving with Z(0,0). zOff is the amplitude removed.
        zOff = smap.removeOffset(None)
        # Printing stuff
        # --------------------------
        print('z-offset removed')
        print(' Z(m,n,amp) = ({:d}, {:d}, {:.2f} nm)'.
              format(smap.zernikeRemoved['00'][0], smap.zernikeRemoved['00'][1],
                     smap.zernikeRemoved['00'][2]) )
        # --------------------------
        smap.plot()
        
        # Removing tilts. A1 is a list with the amplitudes of the two Zernike
        # polynomials involved, Z(-1,1) and Z(1,1).
        A1,xbeta,ybeta = smap.rmTilt(method='zernike')
        # Printing stuff
        # --------------------------
        print('Tilts removed')
        print(' Z(m,n,amp) = ({:d}, {:d}, {:.2f} nm)'.
              format(smap.zernikeRemoved['-11'][0], smap.zernikeRemoved['-11'][1],
                     smap.zernikeRemoved['-11'][2]))
        print(' Z(m,n,amp) = ({:d}, {:d}, {:.2f} nm)'.
              format(smap.zernikeRemoved['11'][0], smap.zernikeRemoved['11'][1],
                     smap.zernikeRemoved['11'][2]))
        print(' or')
        print(' xbeta = {:.2e} rad'.format(xbeta))
        print(' ybeta = {:.2e} rad'.format(ybeta))
        # --------------------------
        smap.plot()

       
    else:
        '''
        Fitting surfaces to process the map.
        '''
        
        # Fitting sphere to the mirror surface. Also converts to the
        # equivalent Zernike(2,0) mode amplitude. 
        Rc, zOff, A20 = smap.remove_curvature(method='sphere', w=w)
        smap.plot()
        
        print('Curvature removed')
        print(' Removed Rc = {0:.2f} m'.format(Rc) )
        print(' or')
        print(' Z(n=2,m=0) amplitude A20 = {0:.2f} nm'.format(A20))
        
        # Fitting flat surface to the mirror.
        zOff = smap.removeOffset(w)
        smap.plot()
        print('Offset removed')
        print(' Removed z-offset (A00) = {0:.3f} nm'.format(zOff))
        
        # Fitting tilted surface to the map. Also calculating
        # equivalent Zernike-amplitude.
        A1,xbeta,ybeta,zOff = smap.rmTilt(method='fitSurf', w=w)
        smap.plot()
        print('Tilted surface removed:')
        print(' xbeta    = {:.2e} rad'.format(xbeta))
        print(' ybeta    = {:.2e} rad'.format(ybeta))
        print(' z-offset = {:.2e} nm'.format(zOff))
        print('Equivalent Zernike amplitudes:')
        print(' A(1,-1) = {:.2f} nm'.format(A1[0]))
        print(' A(1, 1) = {:.2f} nm'.format(A1[1]))
        
    # Creating aperture/absorption map
    amap = aperturemap(smap.name, smap.size, smap.step_size,
                       smap.find_radius(method='min',unit='meters'), smap.center)
    amap.plot()
    print('Aperture map created')
    filename = smap.name + '_finesse.txt'
    smap.write_map(filename)
    print(' Phase map written to file {:s}'.format(filename))
    filename = amap.name + '_aperture.txt'
    amap.write_map(filename)
    print(' Aperture map written to file {:s}'.format(filename))
    
# --------------------------------------------------------------------------------
# End of functions
# --------------------------------------------------------------------------------


pylab.close('all')
# --------------------------------------------------------------------------------
# 1. Read all different map formats. Set isReadAll to True, and
#    see code in the read_all function above.
# --------------------------------------------------------------------------------
isReadAll = False
if isReadAll:
    read_all()
# --------------------------------------------------------------------------------
# 2. Automatically prepare phase map for finesse (in two different ways).
#    Set isAutoProcessing to True and see and play with code in the
#    autoProcess function above.
# --------------------------------------------------------------------------------
isAutoProcessing = True
if isAutoProcessing:
    autoProcess()
# --------------------------------------------------------------------------------
# 3. Manually prepare phase map for finesse. Set isManualProcessing to
#    True, and see and play with code in the manualProcess function
#    above.
# --------------------------------------------------------------------------------
isManualProcessing = False
if isManualProcessing:
    manualProcess()


