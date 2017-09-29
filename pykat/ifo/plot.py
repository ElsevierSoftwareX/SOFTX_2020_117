
def cavity_mode_comparison(kat, node, ax=None, show=True, reference=None, scaling=None, ppm=None):
    import pykat
    import matplotlib.pyplot as plt
    import numpy as np

    from . import mismatch_scan_RoC
    
    mmx, mmy, qs = pykat.ifo.mismatch_cavities(kat, node)
    
    if ax is None:
        plt.figure()
        ax = plt.subplot(111)
    
    avgz = []
    avgzr = []
    
    if reference is None:
        _zx = 0
        _zrx = 0
        _zy = 0
        _zry = 0
    else:
        found = False
        for c, qx, qy in qs:
            if reference in c:
                _zx = qx.z
                _zrx = qx.zr
                _zy = qy.z
                _zry = qy.zr
                
                found = True
                break
                
        if not found:
            raise Exception("Couldn't find reference cavity named '%s'"%reference)
    
    if scaling is None:
        fzx = fzrx = fzy = fzry = lambda x: x
        
        ax.set_xlabel('$z$')
        ax.set_ylabel('$z_r$')
        
    elif scaling == "difference":
        if reference is None:
            raise Exception("reference must be set for scaling")
            
        fzx  = lambda z: (z-_zx)
        fzrx = lambda zr: (zr-_zrx)
        fzy  = lambda z: (z-_zy)
        fzry = lambda zr: (zr-_zry)
        
        ax.set_xlabel("$z - z'$")
        ax.set_ylabel("$z_r - z_r'$")
        
    elif scaling == "relative":
        if reference is None:
            raise Exception("reference must be set for scaling")
            
        fzx  = lambda z: (z-_zx)/_zrx
        fzrx = lambda zr: (zr-_zrx)/_zrx
        fzy  = lambda z: (z-_zy)/_zry
        fzry = lambda zr: (zr-_zry)/_zry
        
        ax.set_xlabel("$\\frac{z-z'}{z_r'}$")
        ax.set_ylabel("$\\frac{z_r-z_r'}{z_r'}$")
        
    for c, qx, qy in qs:
        c = c.strip("cav")
        p1 = ax.scatter(fzx(qx.z), fzrx(qx.zr), marker='x')
        p2 = ax.scatter(fzy(qy.z), fzry(qy.zr),c=p1.get_facecolor()[0], marker='o', label=c)
    
    if ppm is not None:
        if reference is None:
            raise Exception("reference must be set for ppm plotting")
        
        for __z, __zr, ls in ((_zx, _zrx, '-'), (_zy, _zry, '--')):
            z = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100)
            zr = np.linspace(ax.get_ylim()[0], ax.get_ylim()[1], 100)
            
            if scaling is None:
                Z,ZR = np.meshgrid(z,zr)
            elif scaling == "difference":
                Z,ZR = np.meshgrid(z+__z, zr+__zr)
            elif scaling == "relative":
                Z,ZR = np.meshgrid(__zr*z + __z, __zr*zr + __zr)

            Q1 = Z + 1j*ZR

            Z_ = np.ones_like(Z) * __z
            ZR_ = np.ones_like(ZR) * __zr

            Q2 = Z_ + 1j*ZR_

            # 2 times 1-overlap roughly gives ppm loss for ppm << 1
            Y = 2*(1-abs(4*Q1.imag * Q2.imag)/abs(Q1.conjugate()-Q2)**2)/1e-6

            CS = ax.contour(z, zr, Y, levels=ppm, linestyles=ls)
            ax.clabel(CS, CS.levels[::2], inline=True, fmt="%g", fontsize=10)
        
        
    l1 = ax.legend([p1, p2], ["X", "Y"], fontsize=10, bbox_to_anchor=(1.04,0), loc="lower left")
    l2 = ax.legend(fontsize=10, bbox_to_anchor=(1.04,1), loc="upper left")
    
    l1.get_frame().set_alpha(1)
    l2.get_frame().set_alpha(1)
    
    ax.add_artist(l1)
    
    plt.tight_layout()
    
    if show: plt.show()
    
def roc_vs_cavity_overlap(base, cav1, cav2, node, mirror, lower, upper, steps, plot=True, getData=False):
    """
    This function should be used to find the optimum cavity overlap between two cavities
    by varying the radius of curvture of a mirror. From this you can visually see
    if a minimum is found. The function returns the radii of curvature that give the minimum
    overlap in both the X and Y directions.
    
    Function will alter the current curvature between -lower and +upper for number of steps
    provided.
    
    base: kat object
    cav1: cavity for comparison (str or cavity object)
    cav2: cavity for comparison (str or cavity object)
    node: node to trace beam parameter to (str or node object)
    mirror: mirror or beamsplitter to vary curvature of (name str or object)
    lower: lower change in radii of curvature
    upper: upper change in radii of curvature
    steps: number of points to compute
    plot: Set whether to plot or not
    
    return: (x_min_overlap_RoC, y_min_overlap_RoC)
    """
    import pykat
    import matplotlib.pyplot as plt
    import numpy as np

    from . import mismatch_scan_RoC
    
    _kat = base.deepcopy()
    
    # User could pass pykat object or string names in, here we convert to string to get
    # name of object to use to select new deepcopied version
    cav1 = _kat.commands[str(cav1)]
    cav2 = _kat.commands[str(cav2)]
    node = str(node)
    mirror = str(mirror)
    
    qxs = []
    qys = []
    cavs = []
    
    for cav in [cav1, cav2]:
        for _ in _kat.getAll(pykat.commands.cavity):
            _.enabled = False
            
        cav.enabled = True
        
        qx, qy = mismatch_scan_RoC(_kat, node, mirror, lower, upper, steps)
        
        cavs.append(cav.name)
        qxs.append(qx)
        qys.append(qy)
    
    x = np.linspace(lower, upper, steps+1)
    y1 = 1 - pykat.BeamParam.overlap(qxs[0], qxs[1])
    
    if plot:
        plt.figure()
        plt.semilogy(x, y1, label='X')
    
    qx_min = x[y1.argmin()]

    y2 = 1 - pykat.BeamParam.overlap(qys[0], qys[1])
    
    if plot: plt.semilogy(x, y2, label='Y')
    
    qy_min = x[y2.argmin()]
    
    if plot:
        plt.title(mirror + ": "+ " vs ".join(cavs))
        plt.legend(loc=0)
        plt.xlabel(mirror + ' dRoC [m]')
        plt.ylabel('1-Overlap')
        plt.tight_layout()
        plt.show()
    
    if getData:
        return x, y1, y2
    else:
        return qx_min, qy_min
        
def roc_vs_cavity_overlap_2D(base, cav1, cav2, node, mirror1, lower1, upper1, steps1,
                             mirror2, lower2, upper2, steps2, plot=True, getData=False):
    """
    This function should be used to find the optimum cavity overlap between two cavities
    by varying the radius of curvture of two mirrors. From this you can visually see
    if a minimum is found. The function returns the radii of curvature that give the minimum
    overlap in both the X and Y directions.
    
    Function will alter the current curvature between -lower and +upper for number of steps
    provided.
    
    base: kat object
    cav1: cavity for comparison (str or cavity object)
    cav2: cavity for comparison (str or cavity object)
    node: node to trace beam parameter to (str or node object)
    mirror[1/2]: mirror or beamsplitter to vary curvature of (name str or object)
    lower[1/2]: lower change in radii of curvature
    upper[1/2]: upper change in radii of curvature
    steps[1/2]: number of points to compute
    plot: Set whether to plot or not
    getData: If true returns the data rather than min mismatch RoCs
                             
    return: ((x_min_RoC1, x_min_RoC2), (y_min_RoC1, y_min_RoC2)
    
    if getData=True
    
    return: ((x_min_RoC1, x_min_RoC2), (y_min_RoC1, y_min_RoC2), dRoC1, dRoc2, overlap_x, overlap_y
    """
    import pykat
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import LogNorm
    
    from pykat.ifo import mismatch_scan_RoC_2D
    
    _kat = base.deepcopy()
    
    # User could pass pykat object or string names in, here we convert to string to get
    # name of object to use to select new deepcopied version
    cav1 = _kat.commands[str(cav1)]
    cav2 = _kat.commands[str(cav2)]
    node = str(node)
    mirror1 = str(mirror1)
    mirror2 = str(mirror2)
    
    qxs = []
    qys = []
    cavs = []
    
    for cav in [cav1, cav2]:
        for _ in _kat.getAll(pykat.commands.cavity):
            _.enabled = False
            
        cav.enabled = True
        
        qx, qy = mismatch_scan_RoC_2D(_kat, node, mirror1, lower1, upper1, steps1, mirror2, lower2, upper2, steps2)
        
        cavs.append(cav.name)
        qxs.append(qx)
        qys.append(qy)
    
    x = np.linspace(lower1, upper1, steps1+1)
    y = np.linspace(lower2, upper2, steps2+1)
    
    zx = 1 - pykat.BeamParam.overlap(qxs[0], qxs[1])
    zy = 1 - pykat.BeamParam.overlap(qys[0], qys[1])
    
    if plot:
        norm = LogNorm(vmin=min((zx.min(), zy.min())), vmax=max((zx.max(), zy.max())))
        plt.figure()
        plt.subplot(121)
        p = plt.pcolormesh(x, y, zx, norm=norm)
        p.set_rasterized(True)
        plt.xlabel("dRoC %s [m]" % mirror1, fontsize=10)
        plt.ylabel("dRoC %s [m]" % mirror2, fontsize=10)
        plt.title("X", fontsize=10)
        cb = plt.colorbar()
        cb.solids.set_rasterized(True) 
        
        plt.subplot(122)
        p = plt.pcolormesh(x, y, zy, norm=norm)
        p.set_rasterized(True)
        plt.xlabel("dRoC %s [m]" % mirror1, fontsize=10)
        plt.ylabel("dRoC %s [m]" % mirror2, fontsize=10)
        plt.title("Y", fontsize=10)
        cb = plt.colorbar()
        cb.solids.set_rasterized(True) 
        
        cb.set_label("1-Overlap(%s,%s)"% (cav1.name, cav2.name))
        plt.tight_layout()
        plt.show()
    
    xb,xa = np.unravel_index(zx.argmin(), zx.shape)
    yb,ya = np.unravel_index(zy.argmin(), zy.shape)
    
    if getData:
        return ((x[xa],y[xb]), (x[ya],y[yb])), x,y,zx,zy
    else:
        return ((x[xa],y[xb]), (x[ya],y[yb]))