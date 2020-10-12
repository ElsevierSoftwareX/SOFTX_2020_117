import pykat
import pykat.exceptions as pkex
import copy
import numpy as np
from warnings import warn

def plot_beam_trace(_kat, from_node, to_node,
                    fig=None,size_axes=True,gouy_axes=True,return_data=False):
    """Plots the beam radius between two nodes.

    Args:
        _kat (pykat.finesse.kat): The kat object containing the parsed katcode
        from_node (str): Node at which to start plotting
        to_node (str): Node to end plotting
        fig (matplotlib.figure): Figure to plot onto. None creates a new axis.
        size_axes (matplotlib.axes): Beam shape axes. True creates a axes, False disables.
        gouy_axes (matplotlib.axes): Gouy phase axes. True creates a axes, False disables.
        return_data (bool): Return the figure object and data.
    """
    import matplotlib.pyplot as plt
    
    if _kat == None:
        raise pkex.BasePyKatException('kat object in None')

    kat = copy.deepcopy(_kat)
        
    components = kat.nodes.getComponentsBetween(from_node, to_node)

    if len(components) == 0:
        raise pkex.BasePyKatException('Could not find a direct path between nodes %s and %s'%(from_node, to_node))
        
    nodes = []

    for cn in range(len(components)):
        for n in components[cn].nodes:
    
            cs = np.array(n.components)

            if cn < len(components)-1 and components[cn] in cs and components[cn+1] in cs:
                nodes.append(n)
                break
            elif cn == len(components)-1 and components[cn-1] not in cs and components[cn] in cs:
                nodes.append(n)
                break
        
    v = [c for c in zip(components, nodes) if isinstance(c[0], pykat.components.space)]
    spaces = [d[0] for d in v]
    nodes2  = [d[1] for d in v]
    
    # If at least one axes is enabled and no figure object is availale
    # create a new figure object
    # (Is not False needs to be used because we need to use *is* not *==*)
    plot_enabled = (size_axes is not False or gouy_axes is not False)
    if (fig is None) and plot_enabled:
        fig = plt.figure()
    
    if len(spaces) > 0:
        print("Plotting beam along spaces", ", ".join([s.name for s in spaces]))
    
        for d in kat.detectors.values():
            d.remove()

        if hasattr(kat, "xaxis"):
            kat.xaxis.remove()

        if hasattr(kat, "x2axis"):
            kat.x2axis.remove()

        if hasattr(kat, "x3axis"):
            kat.x3axis.remove()

        kat.noxaxis = False
        kat.maxtem = 0
    
        currL = 0
        
        L = []
        wx = []
        gx = []
        wy = []
        gy = []
            
        for n in range(len(spaces)):
            s = spaces[n]
            Lmax = s.L
            N = 100
            node = None

            cmds = """
            gouy gx1 x {spaces}
            bp bpx x w {node}*
            gouy gy1 y {spaces}
            bp bpy y w {node}*
            xaxis {space} L lin 0 {Lmax} {N}
            """.format(space=s.name, Lmax=Lmax, N=int(N), spaces=" ".join([s.name for s in spaces[0:n+1]]), node=nodes2[n])
        
            k = copy.deepcopy(kat)
            k.parse(cmds)
            k.verbose = False
            out = k.run()
        
            L.extend(currL + out.x)
            wx.extend(out['bpx'])
            wy.extend(out['bpy'])
            gx.extend(out['gx1'])
            gy.extend(out['gy1'])
            currL += Lmax
        
        wx = np.array(wx)/1e-2
        wy = np.array(wy)/1e-2
        
        if size_axes is True and gouy_axes is True:
            nplot=2
        else:
            nplot=1
        
        if size_axes is not False:
            if size_axes is True:
                size_axes = fig.add_subplot(nplot,1,1)
            size_axes.plot(L, wx, 'b', L, wy, 'r', L, -wx, 'b', L, -wy, 'r')
            size_axes.set_title("Beam size between %s and %s" % (from_node, to_node))
            size_axes.set_ylabel("Beam size [cm]")
            size_axes.set_xlabel("Distance from %s to %s [m]" % (from_node, to_node))
            size_axes.set_xlim(min(L), max(L))
    
        if gouy_axes is not False:
            if gouy_axes is True:
                gouy_axes = fig.add_subplot(nplot,1,nplot)
            gouy_axes.plot(L, gx, 'b', L, gy, 'r')
            gouy_axes.set_title("Gouy phase accumulation between %s and %s" % (from_node, to_node))
            gouy_axes.set_ylabel("Phase accumulation [deg]")
            gouy_axes.set_xlabel("Distance [m]")
            gouy_axes.set_xlim(min(L), max(L))
            
        if plot_enabled:
            fig.legend(['x','y'], loc='upper right')
            fig.tight_layout()
        
        if return_data:
            return fig,{'wx':wx,'wy':wy,'L':L,'gx':gx,'gy':gy}
        elif plot_enabled:
            fig.show()
        else:
            warn('Data and plotting disabled, no output will be returned.')
        
