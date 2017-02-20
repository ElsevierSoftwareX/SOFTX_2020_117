import pykat
import pykat.exceptions as pkex
import copy
import numpy as np

def plot_beam_trace(_kat, from_node, to_node):
    import pylab
    
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
    
    pylab.figure()
    
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
            k.parseCommands(cmds)
            k.verbose = False
            out = k.run()
            pylab.subplot(2,1,1)
        
            L.extend(currL + out.x)
            wx.extend(out['bpx'])
            wy.extend(out['bpy'])
            gx.extend(out['gx1'])
            gy.extend(out['gy1'])
            currL += Lmax
            
        pylab.subplot(2,1,1)
        wx = np.array(wx)/1e-2
        wy = np.array(wy)/1e-2
        pylab.plot(L, wx, 'b', L, wy, 'r', L, -wx, 'b', L, -wy, 'r')
        pylab.title("Beam size between %s and %s" % (from_node, to_node))
        pylab.ylabel("Beam size [cm]")
        pylab.xlabel("Distance from %s to %s [m]" % (from_node, to_node))
        pylab.xlim(min(L), max(L))
    
        pylab.subplot(2,1,2)
        pylab.plot(L, gx, 'b', L, gy, 'r')
        pylab.title("Gouy phase accumulation between %s and %s" % (from_node, to_node))
        pylab.ylabel("Phase accumulation [deg]")
        pylab.xlabel("Distance [m]")
        pylab.xlim(min(L), max(L))
        pylab.legend(['x','y'], loc=0)
        pylab.tight_layout()
        pylab.show()
        