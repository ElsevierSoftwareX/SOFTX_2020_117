# This all needs to be reworked to fit in with new plotting
# methods. Daniel 1/12/2015


# # ------------------------------------------------------
# # This file is very much in flux at the moment and will
# # undergo many significant changes without notice.
# # Don't use this for anything but playing/testing yet.
# # Andreas 06.03.2014
# # ------------------------------------------------------
#
# import numpy as np
# import matplotlib
# import matplotlib.backends.backend_pdf
# from matplotlib import cm
# BACKEND = 'Qt4Agg'
# matplotlib.use(BACKEND)
# from matplotlib import rc
# import matplotlib.pyplot as plt
#
# formatter = matplotlib.ticker.EngFormatter(unit='', places=0)
# formatter.ENG_PREFIXES[-6] = 'u'
#
#
# def argmaxM(data):
#     return(np.unravel_index(np.argmax(data),data.shape))
#
#
# def printPDF(fig, filename):
#     if filename !='' and filename != None:
#         if not filename.lower().endswith(('.pdf')):
#             filename=filename+'.pdf'
#         pdfp = matplotlib.backends.backend_pdf.PdfPages(filename)
#         pdfp.savefig(fig, dpi=300,bbox_inches='tight')
#         pdfp.close()
#
# def plot_field(field):
#     ax,fig=plot_setup()
#     im = ax.imshow(np.abs(field),origin='lower', aspect='auto')
#     cb = fig.colorbar(im, format="%.4g")
#     #cb.set_clim(-1.0*zlimit, zlimit)
#     ax.autoscale(False)
#     plt.show(block=0)
#
# def plot_propagation(field, grid, Lambda, axis, normalise, Lrange, nr):
#     # axis = 0 for x and =1 for y crosssection
#     # if normalise = 1, each slice is normalises to peak intensity
#     [n,m]=field.shape
#     slices = np.zeros((n,np.size(Lrange)))
#     from pykat.optics.fft import FFT_propagate
#     for idx, L in enumerate(Lrange):
#         field2 = FFT_propagate(field,grid, Lambda, L, nr)
#         if axis==0:
#             slices[:,idx]=np.abs(field2[:,int(m/2)])**2
#         else:
#             slices[:,idx]=np.abs(field2[int(n/2),:])**2
#         if normalise==1:
#             peak_intensity= np.abs(field2[int(n/2),int(m/2)])**2
#             slices[:,idx]=slices[:,idx]/peak_intensity
#     ax,fig=plot_setup()
#     im = ax.imshow(slices, aspect='auto')
#     cb = fig.colorbar(im, format="%.4g")
#     #cb.set_clim(-1.0*zlimit, zlimit)
#     ax.autoscale(False)
#     plt.show(block=0)
#
#
# def plot_power_contour(x,y,z,xlabel, ylabel, clabel, title='', filename=''):
#     ax,fig=plot_setup()
#     xm,ym=np.meshgrid(x,y)
#     extent = [x[0], x[-1], y[0], y[-1]]
#     mycm='YlOrRd_r'
#     #mycm='hot'
#     #im = ax.imshow(z,origin='lower',extent=extent,cmap=mycm, aspect='auto', interpolation='nearest')
#     im = ax.imshow(z,origin='lower',extent=extent,cmap=mycm, aspect='auto')
#     #im = ax.imshow(z,origin='lower',cmap=mycm, aspect='auto')
#     cb = fig.colorbar(im, format="%.4g")
#     #cb.set_clim(-1.0*zlimit, zlimit)
#     ax.autoscale(False)
#     #ct = ax.contour(xm,ym,z, zdir='z', levels = [0])
#     #ax.clabel(ct,inline=1,fontsize=10)
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel(ylabel)
#     cb.set_label(clabel)
#     ax.xaxis.set_major_formatter(formatter)
#     ax.yaxis.set_major_formatter(formatter)
#     fig.canvas.manager.set_window_title(title)
#     plt.show(block=0)
#     printPDF(fig, filename)
#
#
#
# def plot_error_contour(x,y,z,xlabel, ylabel, clabel, title='', filename=''):
#     global fig, ax, cb, ct, data
#     rc('font',**pp.font)
#     rc('xtick',labelsize=pp.TICK_SIZE)
#     rc('ytick',labelsize=pp.TICK_SIZE)
#     rc('text', usetex=pp.USETEX)
#     rc('axes', labelsize = pp.LABEL_SIZE)
#     fig, ax =plt.subplots()
#     fig.set_size_inches(pp.fig_size)
#     fig.set_dpi(pp.FIG_DPI)
#     xm,ym=np.meshgrid(x,y)
#     RGB1 = 255*np.array([0.23, 0.299, 0.754])
#     RGB2 = 255*np.array([0.706,0.016, 0.15])
#     cm1 = CM(RGB1, RGB2)
#     cm_values = cm1.getMap()
#     #cm1.showMap()
#     mycm = matplotlib.colors.ListedColormap(cm_values)
#     extent = [x[0], x[-1], y[0], y[-1] ]
#     data=z
#     # make symmetric color display
#     zlimit=np.max([abs(data.max()),abs(data.min())])
#     im = ax.imshow(data,origin='lower',extent=extent,cmap=mycm, aspect='auto', interpolation='nearest')
#     #im = ax.imshow(data,origin='lower',extent=extent,cmap=mycm, aspect='auto')
#     cb = fig.colorbar(im, format="%.4g")
#     cb.set_clim(-1.0*zlimit, zlimit)
#     ax.autoscale(False)
#     ct = ax.contour(xm,ym,data, zdir='z', levels = [0])
#     ax.clabel(ct,inline=1,fontsize=10)
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel(ylabel)
#     cb.set_label(clabel)
#     ax.xaxis.set_major_formatter(formatter)
#     fig.canvas.manager.set_window_title(title)
#     plt.show()
#     printPDF(fig, filename)
#
# def plot_setup():
#     rc('font',**pp.font)
#     rc('xtick',labelsize=pp.TICK_SIZE)
#     rc('ytick',labelsize=pp.TICK_SIZE)
#     rc('text', usetex=pp.USETEX)
#     rc('axes', labelsize = pp.LABEL_SIZE)
#     fig, ax =plt.subplots()
#     fig.set_size_inches(pp.fig_size)
#     fig.set_dpi(pp.FIG_DPI)
#     return ax,fig
#
# class pp():
#     # set some gobal settings first
#     BACKEND = 'Qt4Agg' # matplotlib backend
#     FIG_DPI=90 # DPI of on sceen plot
#     # Some help in calculating good figure size for Latex
#     # documents. Starting with plot size in pt,
#     # get this from LaTeX using \showthe\columnwidth
#     fig_width_pt = 484.0
#     inches_per_pt = 1.0/72.27  # Convert TeX pt to inches
#     golden_mean = (np.sqrt(5)-1.0)/2.0   # Aesthetic ratio
#     fig_width = fig_width_pt*inches_per_pt  # width in inches
#     fig_height = fig_width*golden_mean      # height in inches
#     fig_size = [fig_width,fig_height]
#     # some plot options:
#     LINEWIDTH = 1 # linewidths of traces in plot
#     AA = True # antialiasing of traces
#     USETEX = False # use Latex encoding in text
#     SHADOW = False # shadow of legend box
#     GRID = True # grid on or off
#     # font sizes for normal text, tick labels and legend
#     FONT_SIZE = 10 # size of normal text
#     TICK_SIZE = 10 # size of tick labels
#     LABEL_SIZE = 10 # size of axes labels
#     LEGEND_SIZE = 10 # size of legend
#     # font family and type
#     font = {'family':'sans-serif','sans-serif':['Helvetica'],'size':FONT_SIZE}
#     DPI=300 # DPI for saving via savefig
#     # print options given to savefig command:
#     print_options = {'dpi':DPI, 'transparent':True, 'bbox_inches':'tight', 'pad_inches':0.1}
#     # for Palatino and other serif fonts use:
#     #font = {'family':'serif','serif':['Palatino']}
#     SCREEN_TITLE = True # show title on screen?
#     PRINT_TITLE = False # show title in saved file?
#
# # =============================================================================
# # ====================== The Class ColorMapCreator ============================
# # =============================================================================
#
# class CM:
#     """
#     Class ColorMapCreator:
#     Create diverging colormaps from RGB1 to RGB2 using the method of Moreland
#     or a simple CIELAB-interpolation. numColors controls the number of color
#     values to output (odd number) and divide gives the possibility to output
#     RGB-values from 0.0-1.0 instead of 0-255. If a filename different than
#     "" is given, the colormap will be saved to this file, otherwise a simple
#     output using print will be given.
#     """
#
#     # ======================== Global Variables ===============================
#
#     # Reference white-point D65
#     Xn, Yn, Zn = [95.047, 100.0, 108.883] # from Adobe Cookbook
#
#     # Transfer-matrix for the conversion of RGB to XYZ color space
#     transM = np.array([[0.4124564, 0.2126729, 0.0193339],
#                         [0.3575761, 0.7151522, 0.1191920],
#                         [0.1804375, 0.0721750, 0.9503041]])
#
#
#     # ============================= Functions =================================
#
#
#     def __init__(self, RGB1, RGB2, numColors = 257., divide = 255.,
#                   method = "moreland", filename = ""):
#
#         # create a class variable for the number of colors
#         self.numColors = numColors
#
#         # assert an odd number of points
#         assert np.mod(numColors,2) == 1, \
#             "For diverging colormaps odd numbers of colors are desireable!"
#
#         # assert a known method was specified
#         knownMethods = ["moreland", "lab"]
#         assert method in knownMethods, "Unknown method was specified!"
#
#         if method == knownMethods[0]:
#             #generate the Msh diverging colormap
#             self.colorMap = self.generateColorMap(RGB1, RGB2, divide)
#         elif method == knownMethods[1]:
#             # generate the Lab diverging colormap
#             self.colorMap = self.generateColorMapLab(RGB1, RGB2, divide)
#
#     def getMap(self):
#         return self.colorMap
#
#     def showMap(self):
#         #rc('text', usetex=False)
#         a=np.outer(np.arange(0,1,0.01),np.ones(10))
#         fig=plt.figure(99,figsize=(10,2))
#         plt.axis("off")
#         cm = matplotlib.colors.ListedColormap(self.colorMap)
#         pm=plt.imshow(a,aspect='auto',cmap=cm,origin="lower")
#         plt.clim(0,1)
#         fig.colorbar(pm)
#         plt.draw()
#
#
#     def rgblinear(self, RGB):
#         """
#         Conversion from the sRGB components to RGB components with physically
#         linear properties.
#         """
#
#         # initialize the linear RGB array
#         RGBlinear = np.zeros((3,))
#
#         #  calculate the linear RGB values
#         for i,value in enumerate(RGB):
#             value = float(value) / 255.
#             if value > 0.04045 :
#                 value = ( ( value + 0.055 ) / 1.055 ) ** 2.4
#             else :
#                 value = value / 12.92
#             RGBlinear[i] = value * 100.
#         return RGBlinear
#     #-
#
#     def sRGB(self, RGBlinear):
#         """
#         Back conversion from linear RGB to sRGB.
#         """
#
#         # initialize the sRGB array
#         RGB = np.zeros((3,))
#
#         #  calculate the sRGB values
#         for i,value in enumerate(RGBlinear):
#             value = float(value) / 100.
#
#             if value > 0.00313080495356037152:
#                 value = (1.055 * np.power(value,1./2.4) ) - 0.055
#             else :
#                 value = value * 12.92
#
#             RGB[i] = round(value * 255.)
#         return RGB
#     #-
#
#     def rgb2xyz(self, RGB):
#         """
#         Conversion of RGB to XYZ using the transfer-matrix
#         """
#         return np.dot(self.rgblinear(RGB), self.transM)
#     #-
#
#     def xyz2rgb(self, XYZ):
#         """
#         Conversion of RGB to XYZ using the transfer-matrix
#         """
#         #return np.round(np.dot(XYZ, np.array(np.matrix(transM).I)))
#         return self.sRGB(np.dot(XYZ, np.array(np.matrix(self.transM).I)))
#     #-
#
#     def rgb2Lab(self, RGB):
#         """
#         Conversion of RGB to CIELAB
#         """
#
#         # convert RGB to XYZ
#         X, Y, Z = (self.rgb2xyz(RGB)).tolist()
#
#         # helper function
#         def f(x):
#             limit = 0.008856
#             if x> limit:
#                 return np.power(x, 1./3.)
#             else:
#                 return 7.787*x + 16./116.
#
#         # calculation of L, a and b
#         L = 116. * ( f(Y/self.Yn) - (16./116.) )
#         a = 500. * ( f(X/self.Xn) - f(Y/self.Yn) )
#         b = 200. * ( f(Y/self.Yn) - f(Z/self.Zn) )
#         return np.array([L, a, b])
#     #-
#
#     def Lab2rgb(self, Lab):
#         """
#         Conversion of CIELAB to RGB
#         """
#
#         # unpack the Lab-array
#         L, a, b = Lab.tolist()
#
#         # helper function
#         def finverse(x):
#             xlim = 0.008856
#             a = 7.787
#             b = 16./116.
#             ylim = a*xlim+b
#             if x > ylim:
#                 return np.power(x, 3)
#             else:
#                 return ( x - b ) / a
#
#         # calculation of X, Y and Z
#         X = self.Xn * finverse( (a/500.) + (L+16.)/116. )
#         Y = self.Yn * finverse( (L+16.)/116. )
#         Z = self.Zn * finverse( (L+16.)/116. - (b/200.) )
#
#         # conversion of XYZ to RGB
#         return self.xyz2rgb(np.array([X,Y,Z]))
#     #-
#
#     def Lab2Msh(self, Lab):
#         """
#         Conversion of CIELAB to Msh
#         """
#
#         # unpack the Lab-array
#         L, a, b = Lab.tolist()
#
#         # calculation of M, s and h
#         M = np.sqrt(np.sum(np.power(Lab, 2)))
#         s = np.arccos(L/M)
#         h = np.arctan2(b,a)
#         return np.array([M,s,h])
#     #-
#
#     def Msh2Lab(self, Msh):
#         """
#         Conversion of Msh to CIELAB
#         """
#
#         # unpack the Msh-array
#         M, s, h = Msh.tolist()
#
#         # calculation of L, a and b
#         L = M*np.cos(s)
#         a = M*np.sin(s)*np.cos(h)
#         b = M*np.sin(s)*np.sin(h)
#         return np.array([L,a,b])
#     #-
#
#     def rgb2Msh(self, RGB):
#         """ Direct conversion of RGB to Msh. """
#         return self.Lab2Msh(self.rgb2Lab(RGB))
#     #-
#
#     def Msh2rgb(self, Msh):
#         """ Direct conversion of Msh to RGB. """
#         return self.Lab2rgb(self.Msh2Lab(Msh))
#     #-
#
#     def adjustHue(self, MshSat, Munsat):
#         """
#         Function to provide an adjusted hue when interpolating to an
#         unsaturated color in Msh space.
#         """
#
#         # unpack the saturated Msh-array
#         Msat, ssat, hsat = MshSat.tolist()
#
#         if Msat >= Munsat:
#             return hsat
#         else:
#             hSpin = ssat * np.sqrt(Munsat**2 - Msat**2) / \
#                     (Msat * np.sin(ssat))
#             if hsat > -np.pi/3:
#                 return hsat + hSpin
#             else:
#                 return hsat - hSpin
#     #-
#
#     def interpolateColor(self, RGB1, RGB2, interp):
#         """
#         Interpolation algorithm to automatically create continuous diverging
#         color maps.
#         """
#
#         # convert RGB to Msh and unpack
#         Msh1 = self.rgb2Msh(RGB1)
#         M1, s1, h1 = Msh1.tolist()
#         Msh2 = self.rgb2Msh(RGB2)
#         M2, s2, h2 = Msh2.tolist()
#
#         # If points saturated and distinct, place white in middle
#         if (s1>0.05) and (s2>0.05) and ( np.abs(h1-h2) > np.pi/3. ):
#             Mmid = max([M1, M2, 88.])
#             if interp < 0.5:
#                 M2 = Mmid
#                 s2 = 0.
#                 h2 = 0.
#                 interp = 2*interp
#             else:
#                 M1 = Mmid
#                 s1 = 0.
#                 h1 = 0.
#                 interp = 2*interp-1.
#
#         # Adjust hue of unsaturated colors
#         if (s1 < 0.05) and (s2 > 0.05):
#             h1 = self.adjustHue(np.array([M2,s2,h2]), M1)
#         elif (s2 < 0.05) and (s1 > 0.05):
#             h2 = self.adjustHue(np.array([M1,s1,h1]), M2)
#
#         # Linear interpolation on adjusted control points
#         MshMid = (1-interp)*np.array([M1,s1,h1]) + \
#                  interp*np.array([M2,s2,h2])
#
#         return self.Msh2rgb(MshMid)
#     #-
#
#     def generateColorMap(self, RGB1, RGB2, divide):
#         """
#         Generate the complete diverging color map using the Moreland-technique
#         from RGB1 to RGB2, placing "white" in the middle. The number of points
#         given by "numPoints" controls the resolution of the colormap. The
#         optional parameter "divide" gives the possibility to scale the whole
#         colormap, for example to have float values from 0 to 1.
#         """
#
#         # calculate
#         scalars = np.linspace(0., 1., self.numColors)
#         RGBs = np.zeros((self.numColors, 3))
#         for i,s in enumerate(scalars):
#             RGBs[i,:] = self.interpolateColor(RGB1, RGB2, s)
#         return RGBs/divide
#     #-
#
#     def generateColorMapLab(self, RGB1, RGB2, divide):
#         """
#         Generate the complete diverging color map using a transition from
#         Lab1 to Lab2, transitioning true RGB-white. The number of points
#         given by "numPoints" controls the resolution of the colormap. The
#         optional parameter "divide" gives the possibility to scale the whole
#         colormap, for example to have float values from 0 to 1.
#         """
#
#         # convert to Lab-space
#         Lab1 = self.rgb2Lab(RGB1)
#         Lab2 = self.rgb2Lab(RGB2)
#         LabWhite = np.array([100., 0., 0.])
#
#         # initialize the resulting arrays
#         Lab = np.zeros((self.numColors ,3))
#         RGBs = np.zeros((self.numColors ,3))
#         N2 = np.floor(self.numColors/2.)
#
#         # calculate
#         for i in range(3):
#             Lab[0:N2+1, i] = np.linspace(Lab1[i], LabWhite[i], N2+1)
#             Lab[N2:, i] = np.linspace(LabWhite[i], Lab2[i], N2+1)
#         for i,l in enumerate(Lab):
#             RGBs[i,:] = self.Lab2rgb(l)
#         return RGBs/divide
#     #-
#
