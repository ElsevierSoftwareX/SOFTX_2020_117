from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import numpy as np
import math
import copy
import warnings
import cmath
import inspect
import six 

from itertools import chain

from pykat import finesse
from pykat.finesse import BlockedKatFile
from pykat.ifo import *

import pykat.components
import pykat.exceptions as pkex
import pykat.external.peakdetect as peak
import pkg_resources

from scipy.constants import c as clight
from scipy.optimize import fmin

class VOYAGER_IFO(IFO):   
    """
    This contains Voyager specific methods for computing interferometer
    variables.
    
    Functions that operate on the kat/IFO objects, manipulate their
    structure and return a new one, should not be included here. They should
    be separate functions that are called by the user. 
    
    The functions here should be those that update the kat with information
    from the IFO object or vice-versa.
    """

    def __init__(self, kat, tuning_keys_list, tunings_components_list):
        IFO.__init__(self, kat, tuning_keys_list, tunings_components_list)
        self._f1 = np.nan
        self._f2 = np.nan
        self._f3 = np.nan
        self._f36M = np.nan

    @property
    def f1(self):
        return self._f1
    @f1.setter
    def f1(self, value):
        self._f1 = float(value)
        self.f36M = self.f2 - self.f1
        # Updating ports
        if hasattr(self, 'LSC_DOFs'):
            for a in self.LSC_DOFs:
                if a.port.name[-2:] == 'f1':
                    a.port.f = self.f1
                
    @property
    def f2(self):
        return self._f2
    @f2.setter
    def f2(self, value):
        self._f2 = float(value)
        self.f36M = self.f2 - self.f1
        # Updating ports
        if hasattr(self, 'LSC_DOFs'):
            for a in self.LSC_DOFs:
                if a.port.name[-2:] == 'f2':
                    a.port.f = self.f2
        
    @property
    def f3(self):
        return self._f3
    @f3.setter
    def f3(self, value):
        self._f3 = float(value)
        # Updating ports
        if hasattr(self, 'LSC_DOFs'):
            for a in self.LSC_DOFs:
                if a.port.name[-2:] == 'f3':
                    a.port.f = self.f3
        
    @property
    def f36M(self):
        return self._f36M
    @f36M.setter
    def f36M(self, value):
        self._f36M = float(value)
        # Updating ports
        if hasattr(self, 'LSC_DOFs'):
            for a in self.LSC_DOFs:
                if a.port.name[-4:] == 'f36M':
                    a.port.f = self.f36M
        
    
    def compute_derived_resonances(self):
        self.fsrX = 0.5 * clight / float(self.kat.LX.L)
        self.fsrY = 0.5 * clight / float(self.kat.LY.L)
        self.fsrPRC = 0.5 * clight / self.lPRC
        self.fsrSRC = 0.5 * clight / self.lSRC
        self.f1_PRC = 3.5 * self.fsrPRC
    
    def compute_derived_lengths(self, verbose=False):
        """
        Compute derived length from individual space components.
        Design values are currently:
        lPRC = 57.656, lSRC = 56.008, lSchnupp = 0.08
        and the individual lengths:
        PRC: Lp1 16.6107, Lp2 16.1647, Lp3 19.5381 
        SRC: Ls1 15.7586, Ls2 15.4435, Ls3 19.3661
        """
        # distances between HR surfaces:
        self.lpr = self.kat.lp1.L + self.kat.lp2.L + self.kat.lp3.L
        self.lx = self.kat.lx1.L + self.kat.BSsub1.L * self.kat.BSsub1.n + self.kat.ITMXsub.L * self.kat.ITMXsub.n
        self.ly = self.kat.ly1.L + self.kat.ITMYsub.L * self.kat.ITMYsub.n
        self.lsr = self.kat.ls1.L + self.kat.ls2.L + self.kat.ls3.L + self.kat.BSsub2.L * self.kat.BSsub2.n
        # resulting combined distances (single, not roundtrip)
        self.lMI =  0.5 * (self.lx + self.ly)
        self.lPRC = self.lpr + self.lMI
        self.lSRC = self.lsr + self.lMI
        self.lSchnupp = self.lx - self.ly
    
        self.compute_derived_resonances()
    
    def suspend_mirrors_z(self):
        """
        Suspends the main mirrors in an aLIGO model in z in the
        supplied kat object.
    
        Returns the commands used for reference.
        """
    
        code = """
        attr ETMY mass 203
        attr ETMX mass 203
        attr ITMY mass 203
        attr ITMX mass 203
        """
    
        self.kat.parse(code)
    
        return code

    
    def fix_mirrors(self, z=True, pitch=True, yaw=True):
        """
        This function will iterate through the main mirrors
        and remove any suspension settings on them. This can be
        done individuallly or for z, pitch, and yaw.
        """
    
        for mirror in ["ETMY","ETMX","ITMY","ITMX","PRM","PR2","PR3","SRM","SR2","SR3","BS"]:
            mirror = self.kat.components[mirror]
        
            if z:
                mirror.mass = None
                mirror.zmech = None
            
            if pitch:
                mirror.Iy = None
                mirror.rymech = None
        
            if yaw:
                mirror.Ix = None
                mirror.rxmech = None
       
    def cavity_status(self):
        kat = self.kat.deepcopy()
        kat.verbose = False
        kat.noxaxis = True
        kat.IFO.fix_mirrors()
        kat.IFO.remove_modulators()
        kat.removeBlock('locks',False)
        kat.maxtem = 0

        kat.parse('yaxis abs')

        cavs = kat.getAll(pykat.commands.cavity)   
        optics = [x for x in kat.getAll(pykat.components.mirror) + kat.getAll(pykat.components.beamSplitter) if 'AR' not in x.name]

        cp_cmds = ''
        for cav in cavs:
            cp_cmds += '''cp {0}_x_g {0} x stability
                        cp {0}_y_g {0} y stability
                        cp {0}_x_B {0} x B
                        cp {0}_y_B {0} y B
                        cp {0}_x_fsr {0} x fsr
                        cp {0}_y_fsr {0} y fsr
                        '''.format(cav)
    
        bp_cmds = ''
        for optic in optics:
            bp_cmds += '''bp {0}_x_w x w {1}
                        bp {0}_y_w y w {1}
                        '''.format(optic,optic.nodes[0].name)

        kat.parse(cp_cmds)
        kat.parse(bp_cmds)
    
        out = kat.run()

        print(" .-------------------------------------------------.")
        print(" | optic    | w_x[cm]    | w_y[cm]                 |")
        print(" +----------|------------|-------------------------+")
        for optic in optics:
            print(' |  {:8}| {:10.4g} | {:10.4g}              |'\
                  .format(optic.name,100*float(out[optic.name+'_x_w']),100*float(out[optic.name+'_x_w'])))
        print(" `-------------------------------------------------'")
        print()
        print(" .----------------------------------------------------------.")
        print(" | cavity   | g_x[0->1]| psi_RT[deg]| delta_f[kHz]| FSR[kHz]|")
        print(" +----------|----------|------------|-------------|---------+")
        for cav in cavs:
            # shift g from -1:1 to 0:1
            g_x = (out[cav.name+'_x_g']+1)/2
            B_x = out[cav.name+'_x_B']
            psi_rt_x = 2*np.arccos(np.sign(B_x)*np.sqrt(g_x))
            FSR = out[cav.name+'_x_fsr']
    
            if np.abs(psi_rt_x/(2*np.pi)) > 1/2:
                delta_f_x = psi_rt_x/(2*np.pi) * FSR - FSR*np.sign(psi_rt_x)
            else:
                delta_f_x = psi_rt_x/(2*np.pi) * FSR
            print(' |  {:8}| {:8.4g} | {:10.4g} | {:9.4g}   | {:6.4g}  |'\
                  .format(cav.name,g_x,180/np.pi*float(psi_rt_x),1e-3*float(delta_f_x),1e-3*float(FSR)))
        print(" `--------------------------------------------------------'")
        print()
        print(" .----------------------------------------------------------.")
        print(" | cavity   | g_y[0->1]| psi_RT[deg]| delta_f[kHz]| FSR[kHz]|")
        print(" +----------|----------|------------|-------------|---------+")
        for cav in cavs:
            # shift g from -1:1 to 0:1
            g_y = (out[cav.name+'_y_g']+1)/2
            B_y = out[cav.name+'_y_B']
            psi_rt_y = 2*np.arccos(np.sign(B_y)*np.sqrt(g_y))
            FSR = out[cav.name+'_y_fsr']       
            if np.abs(psi_rt_y/(2*np.pi)) > 1/2:
                delta_f_y = psi_rt_y/(2*np.pi) * FSR - FSR*np.sign(psi_rt_y)
            else:
                delta_f_y = psi_rt_y/(2*np.pi) * FSR
            print(' |  {:8}| {:8.4g} | {:10.4g} | {:9.4g}   | {:6.4g}  |'\
                  .format(cav.name,g_y,180/np.pi*float(psi_rt_y),1e-3*float(delta_f_y),1e-3*float(FSR)))
        print(" `-----------------------------------------------------------'")
        
     
    def lengths_status(self):
        self.compute_derived_lengths()
        
        print(" .--------------------------------------------------.")
        print("| - arm length:                                     |")
        print("| Lx   = {:11.7}m, Ly   = {:11.7}m          |".format(float(self.kat.LX.L), float(self.kat.LY.L)))
        print("| - small MI and recycling lengths:                 | ")
        print("| lx   = {:11.7}m, ly   = {:11.7}m          |".format(self.lx, self.ly))
        print("| lpr  = {:11.7}m, lsr  = {:11.7}m          |".format(self.lpr, self.lsr))
        print("| lMI  = {:11.7}m, lSchnupp = {:11.5}m      |".format(self.lMI, self.lSchnupp))
        print("| lPRC = {:11.7}m, lSRC = {:11.7}m          |".format(self.lPRC, self.lSRC))
        print("+---------------------------------------------------+")
        print("| - associated cavity frequencies [Hz]:             |")
        print("| fsrx   = {:11.8}, fsry   = {:11.8}        |".format(self.fsrX, self.fsrY))
        print("| fsrPRC = {:11.8}, fsrSRC = {:11.8}        |".format(self.fsrPRC, self.fsrSRC))
        print("| f1_PRC = {:11.8}                              |".format(self.f1_PRC))
        print("| f1     = {:11.8}, f2     = {:11.9}        |".format(self.f1, self.f2))
        print(" `--------------------------------------------------'")
    
    def remove_modulators(self):
        """
        Removes the input modulators and reconnects the input laser to the PRC reflection node.
        
        This function alters the kat object directly.
        """
        self.kat.remove("lmod1", "mod1", "lmod2", "mod2", "lmod3")    # Remove modulators
        # Set output node of laser block to be on the laser
        self.kat.nodes.replaceNode(self.kat.L0, 'n0', 'nLaserOut')
        
    def remove_IMC_HAM2(self, removeIMC, removeHAM2):
        """
        For use with files that have the IMC and HAM2 blocks.
        
        Removes the IMC and HAM2 blocks if not required in the model. Reconnects
        spaces between the laser and HAM2 and PRC. Assumes spaces exists
        with name and node:
            sHAM2in and node nIMCout
            sPRCin  and node nHAM2out
        
        
        This function alters the kat object directly.
        """
        
        if removeHAM2 and not removeIMC:
            raise pkex.BasePyKatException("Must remove IMC if removing HAM2 block")
        
        if removeIMC:
            self.kat.removeBlock("IMC")
            self.kat.cavIMC.remove()
            self.kat.nodes.replaceNode(self.kat.sHAM2in, 'nIMCout', 'nLaserOut')
        
        if removeHAM2:
            self.kat.removeBlock("HAM2")
            self.kat.nodes.replaceNode(self.kat.sPRCin, 'nHAM2out', 'nLaserOut')


    def remove_FI_OMC(self, removeFI=True, removeOMC=True):
        raise NotImplementedError("Disabled for Voyager modelling as we need both LMC and OMC for BHD")

    def adjust_PRC_length(self, verbose=False):
        """
        Adjust PRC length so that it fulfils the requirement
        lPRC = (N+1/2) * c/(2*f1), see [1] equation C.1
        In the current design N=3.
    
        This function directly alters the lengths of the associated kat object.
        """
        kat = self.kat
        
        vprint(verbose, "-- adjusting PRC length")
        ltmp = 0.5 * clight / kat.IFO.f1
        delta_l = 3.5 * ltmp - kat.IFO.lPRC
        vprint(verbose, "   adusting kat.lp1.L by {:.4g}m".format(delta_l))
        kat.lp1.L += delta_l
    
        kat.IFO.compute_derived_lengths(kat)

    def apply_lock_feedback(self, out, idx=None):
        """
        This function will apply the lock values that have been calculated
        in a previous kat run. This should bring the kat object closer to an
        initial lock point so that the lock commands do not need to be run
        on startup.
    
        out: kat run object containing data on lock outputs
        idx: the step in the output array to use
    
        This function directly alters the tunings of the associated kat object.
        """
    
        tuning = self.kat.IFO.get_tunings()

        if "ETMX_lock" in out.ylabels:
            if idx is None:
                tuning["ETMX"] += float(out["ETMX_lock"])
            else:
                tuning["ETMX"] += float(out["ETMX_lock"][idx])
        else:
            pkex.printWarning("could not find ETMX lock")
    
        if "ETMY_lock" in out.ylabels:
            if idx is None:
                tuning["ETMY"] += float(out["ETMY_lock"])
            else:
                tuning["ETMY"] += float(out["ETMY_lock"][idx])
        else:
            pkex.printWarning("could not find ETMY lock")
    
        if "PRCL_lock" in out.ylabels:
            if idx is None:
                tuning["PRM"]  += float(out["PRM_lock"])
            else:
                tuning["PRM"]  += float(out["PRM_lock"][idx])
        else:
            pkex.printWarning("could not find PRCL lock")
    
        if ("MICH_lock" in out.ylabels) and ("ITMY_lock" in out.ylabels):
            if idx is None:
                tuning["ITMX"] += float(out["ITMX_lock"])
                tuning["ITMY"] += float(out["ITMY_lock"])
            else:
                tuning["ITMX"] += float(out["ITMX_lock"][idx])
                tuning["ITMY"] += float(out["ITMY_lock"][idx])
        else:
            pkex.printWarning("could not find MICH (ITMY) lock")
    
        if "SRCL_lock" in out.ylabels:
            if idx is None:
                tuning["SRM"]  += float(out["SRM_lock"])
            else:
                tuning["SRM"]  += float(out["SRM_lock"][idx])
        else:
            pkex.printWarning("could not find SRCL lock")
            
        self.kat.IFO.apply_tunings(tuning)
    
    def add_errsigs_block(self, noplot=True):
        """
        Creates and adds the 'errsigs' block to the kat object based on the
        DARM, CARM, PRCL, MICH and SRCL DOF objects
        
        Removes exisiting errsigs block if present.
        
        Returns the commands added for reference.
        """
        kat = self.kat
        
        sigDARM = kat.IFO.DARM.signal()
        sigCARM = kat.IFO.CARM.signal()
        sigPRCL = kat.IFO.PRCL.signal()
        sigMICH = kat.IFO.MICH.signal()
        sigSRCL = kat.IFO.SRCL.signal()
    
        code2 = ""
        for _ in [sigDARM, sigCARM, sigPRCL, sigMICH, sigSRCL]:
            code2 += "\n".join(_) + "\n"

        code3= ""
    
        if noplot:
            nameDARM = kat.IFO.DARM.signal_name()
            nameCARM = kat.IFO.CARM.signal_name()
            namePRCL = kat.IFO.PRCL.signal_name()
            nameMICH = kat.IFO.MICH.signal_name()
            nameSRCL = kat.IFO.SRCL.signal_name()
                    
        cmds = "".join([code2, code3])
        kat.removeBlock("errsigs", False)
        kat.parse(cmds, addToBlock="errsigs")
        
        return cmds
        
    def add_locks_block(self, lock_data, verbose=False):
        """
        Accepts a dictionary describing the lock gains and accuracies, e.g.:
            data = {
                "DARM": {"accuracy":1, "gain":1},
                "CARM": {"accuracy":1, "gain":1},
                "PRCL": {"accuracy":1, "gain":1},
                "MICH": {"accuracy":1, "gain":1},
                "SRCL": {"accuracy":1, "gain":1},
            }
        
        This then generates the lock block and adds it to the kat object in the 'locks' block.
        
        Removes exisiting locks block if present.
        
        Returns the commands added for reference.
        """
        
        DOFs = ["DARM", "CARM", "PRCL", "MICH", "SRCL"]
        
        names = [getattr(self, _).signal_name() for _ in DOFs]
        accuracies = [lock_data[_]['accuracy'] for _ in DOFs]
        gains = [lock_data[_]['gain'] for _ in DOFs]
        
        code1 = ("set DARM_err {} re\n"
                 "set CARM_err {} re\n"
                 "set PRCL_err {} re\n"
                 "set MICH_err {} re\n"
                 "set SRCL_err {} re\n").format(*names)

        code2 = ("lock DARM_lock $DARM_err {:8.2g} {:8.2g}\n"
                 "lock CARM_lock $CARM_err {:8.2g} {:8.2g}\n"
                 "lock PRCL_lock $PRCL_err {:8.2g} {:8.2g}\n"
                 "lock MICH_lock $MICH_err {:8.2g} {:8.2g}\n"
                 "lock SRCL_lock $SRCL_err {:8.2g} {:8.2g}\n"
                 ).format(*chain.from_iterable(zip(gains, accuracies)))

        # TODO: Use DOF optics and factors to define this. 
        code3 = ("func ETMX_lock = (-1.0) * $CARM_lock - 0.5 * $MICH_lock - $DARM_lock\n"
                 "func ETMY_lock = (-1.0) * $CARM_lock + 0.5 * $MICH_lock + $DARM_lock\n"
                 "func ITMX_lock = (-0.5) * $MICH_lock\n"
                 "func ITMY_lock = 0.5 * $MICH_lock\n"
                 "func PRM_lock = 1.0 * $PRCL_lock\n"
                 "func SRM_lock = (-1.0) * $SRCL_lock\n"

                 "put* PRM     phi     $PRM_lock\n"
                 "put* ITMX    phi     $ITMX_lock\n"
                 "put* ITMY    phi     $ITMY_lock\n"
                 "put* ETMX    phi     $ETMX_lock\n"
                 "put* ETMY    phi     $ETMY_lock\n"
                 "put* SRM     phi     $SRM_lock\n"
                 "put* PRM     phi     $PRM_lock\n"

                 "noplot PRCL_lock\n"
                 "noplot SRCL_lock\n"
                 "noplot MICH_lock\n"
                 "noplot DARM_lock\n"
                 "noplot CARM_lock\n"
                 "noplot ETMX_lock\n"
                 "noplot ETMY_lock\n"
                 "noplot ITMX_lock\n"
                 "noplot ITMY_lock\n"
                 "noplot PRM_lock\n"
                 "noplot SRM_lock\n"
                 )

        if verbose:
            print(" .--------------------------------------------------.")
            print(" | Lock commands used:                              |")
            print(" +--------------------------------------------------+")
            for l in code2.splitlines():
                print (" | {:49}|".format(l))
            print(" `--------------------------------------------------'")

        cmds = "".join([code1, code2, code3])
        
        self.kat.removeBlock("locks", False) # Remove existing block if exists
        self.kat.parse(cmds, addToBlock="locks")
        
        return cmds
    
    def add_REFL_gouy_telescope(self, loss=0, gouy_REFL_BS=0, gouy_A=0, gouy_B=90):
        """
        Adds in the gouy phase telescope for WFS detectors and the IFO port objects.
        Commands added into block "REFL_gouy_tele". This attaches to the
        nREFL node which should be from an isolator on the input path.
        
        Also adds the relevant IFO port objects for generating detectors:
            * ASC_REFL9A, ASC_REFL9B
            * ASC_REFL45A, ASC_REFL45B
            * ASC_REFL36A, ASC_REFL36B
        
        These ports are associated with the block "REFL_gouy_tele".
        
        loss: Total loss accumulated along telescope up to the WFS BS [0 -> 1]
        gouy_REFL_BS:  Gouy phase along path from isolator to WFS BS [deg]
        gouy_A: Gouy phase along A path from BS to WFS [deg]
        gouy_B: Gouy phase along B path from BS to WFS [deg]
        """
        
        self.kat.removeBlock("REFL_gouy_tele", False) # Remove old one
        
        self.kat.parse("""
        s  sFI_REFL_WFS_LOSS 0 nREFL nREFL_loss1
        m2 mREFL_WFS_loss 0 {} 0 nREFL_loss1 nREFL_loss2
        s  sFI_REFL_WFS 0 nREFL_loss2 nREFL_WFS_BS1
        bs WFS_REFL_BS 0.5 0.5 0 0 nREFL_WFS_BS1 nREFL_WFS_BS2 nREFL_WFS_BS3 dump
        s  sWFS_REFL_A  0 nREFL_WFS_BS3 nREFL_WFS_A
        s  sWFS_REFL_B  0 nREFL_WFS_BS2 nREFL_WFS_B
        """.format(loss), addToBlock="REFL_gouy_tele", exceptionOnReplace=True)
        
        self.set_REFL_gouy_telescope_phase(gouy_REFL_BS, gouy_A, gouy_B)
        
        self.kat.IFO.ASC_REFL9A   = Output(self.kat.IFO, "ASC_REFL9A",  "nREFL_WFS_A", "f1", block="REFL_gouy_tele")
        self.kat.IFO.ASC_REFL9B   = Output(self.kat.IFO, "ASC_REFL9B",  "nREFL_WFS_B", "f1", block="REFL_gouy_tele")

        self.kat.IFO.ASC_REFL45A  = Output(self.kat.IFO, "ASC_REFL45A",  "nREFL_WFS_A", "f2", block="REFL_gouy_tele")
        self.kat.IFO.ASC_REFL45B  = Output(self.kat.IFO, "ASC_REFL45B",  "nREFL_WFS_B", "f2", block="REFL_gouy_tele")
        
        self.kat.IFO.ASC_REFL36A  = Output(self.kat.IFO, "ASC_REFL36A",  "nREFL_WFS_A", "36M", block="REFL_gouy_tele")
        self.kat.IFO.ASC_REFL36B  = Output(self.kat.IFO, "ASC_REFL36B",  "nREFL_WFS_B", "36M", block="REFL_gouy_tele")
        
        self.update()
        
    def set_REFL_gouy_telescope_phase(self, gouy_REFL_BS, gouy_A, gouy_B):
        """
        Sets the gouy phase from the the FI to the REFL WFS BS, and then
        the gouy on each path to the A and B detectors. Units all in degrees.
        """
        
        if "REFL_gouy_tele" in self.kat.getBlocks():
            self.kat.sFI_REFL_WFS.gouy = gouy_REFL_BS
            self.kat.sWFS_REFL_A.gouy = gouy_A
            self.kat.sWFS_REFL_B.gouy = gouy_B
        else:
            raise pkex.BasePyKatException("\033[91mREFL Gouy phase telescope isn't in the kat object, see kat.IFO.add_REFL_gouy_telescope()\033[0m")
        
    def scan_REFL_gouy_telescope_gouy_cmds(self, start, end, steps=20, xaxis=1, AB_gouy_diff=None, relative=False):
        """
        This will return commands to scan the REFL gouy telescope gouy phase of the A and B paths.
        """
        if "REFL_gouy_tele" not in self.kat.getBlocks():
            raise pkex.BasePyKatException("\033[91mREFL Gouy phase telescope isn't in the kat object, see kat.IFO.add_REFL_gouy_telescope()\033[0m")
        
        if xaxis not in [1, 2]:
            raise pkex.BasePyKatException("xaxis value must be 1 or 2")
        elif xaxis == 1:
            xaxis_cmd = "xaxis"
        elif xaxis == 2:
            xaxis_cmd = "x2axis"
            
        if AB_gouy_diff is None:
            AB_gouy_diff = self.kat.sWFS_REFL_B.gouy - self.kat.sWFS_REFL_A.gouy
            
        if relative:
            put = "put*"
        else:
            put = "put"
            
        cmds = ("var REFL_GOUY_SCAN 0\n"
        "{xaxis} REFL_GOUY_SCAN re lin {start} {end} {steps}\n"
        "{put} sWFS_REFL_A gx $x{axis}\n"
        "{put} sWFS_REFL_A gy $x{axis}\n"
        "func REFL_SCAN_B = $x{axis} + {AB_gouy_diff}\n"
        "{put} sWFS_REFL_B gx $REFL_SCAN_B\n"
        "{put} sWFS_REFL_B gy $REFL_SCAN_B\n").format(xaxis=xaxis_cmd, axis=xaxis, start=start, end=end, steps=steps, AB_gouy_diff=AB_gouy_diff, put=put)
        
        return cmds
    
    def update(self):
        """
        Iterates through the IFO and updates the DOFs and Outputs dictionaries with the latest ports and DOFs that have
        been added to the interferometer object.
        """
        self.DOFs = {}
    
        for _ in inspect.getmembers(self, lambda x: isinstance(x, DOF)):
            self.DOFs[_[0]] = _[1]
        
        self.Outputs = {}
    
        for _ in inspect.getmembers(self, lambda x: isinstance(x, Output)):
            self.Outputs[_[0]] = _[1]
            
def assert_voyager_ifo_kat(kat):
    if not isinstance(kat.IFO, VOYAGER_IFO):
        raise pkex.BasePyKatException("\033[91mkat file is not an VOYAGER_IFO compatiable kat\033[0m")
              
def make_kat(name="voyager_BSAR_LO", katfile=None, verbose = False, debug=False, keepComments=False, preserveConstants=False):
    """
    Returns a kat object and fills in the kat.IFO property for storing
    the associated interferometer data.
    
    name: Model to load
        - "voyager" base voyager model with IFO, OMC, LMC
    
    keepComments: If true it will keep the original comments from the file
    preserveComments: If true it will keep the const commands in the kat
    """
    names = ['voyager_BSAR_LO', 'voyager_POP_LO']
    
    if debug:
        kat = finesse.kat(tempdir=".",tempname="test")
    else:
        kat = finesse.kat()
    
    kat.verbose=verbose
    
    # Create empty object to just store whatever DOFs, port, variables in
    # that will be used by processing functions
    kat.IFO = VOYAGER_IFO(kat,
                        # Define which keys are used for a tuning description
                        ["maxtem", "phase"],
                        # Define which mirrors create the tuning description
                        ["PRM", "ITMX", "ETMX", "ITMY", "ETMY", "BS", "SRM"])
    
    kat.IFO._data_path=pkg_resources.resource_filename('pykat.ifo', os.path.join('voyager','files'))

    kat.IFO.rawBlocks = BlockedKatFile()
    
    if katfile:
        kat.load(katfile, keepComments=keepComments, preserveConstants=preserveConstants)
        kat.IFO.rawBlocks.read(katfile)
    else:
        if name not in names:
            pkex.printWarning("aLIGO name `{}' not recognised, options are {}, using default 'design'".format(name, names))
        
        katkile = os.path.join(kat.IFO._data_path, name+".kat")
        
        kat.load(katkile, keepComments=keepComments, preserveConstants=preserveConstants)
        kat.IFO.rawBlocks.read(katkile)
    
    # ----------------------------------------------------------------------
    # get and derive parameters from the kat file
    
    # get main frequencies
    if "f1" in kat.constants.keys():
        kat.IFO.f1 = float(kat.constants["f1"].value)
    else:
        kat.IFO.f1 = 9099471.0
        
    if "f2" in kat.constants.keys():
        kat.IFO.f2 = float(kat.constants["f2"].value)
    else:
        kat.IFO.f2 = 5.0 * kat.IFO.f1
        
    if "f3" in kat.constants.keys():
        kat.IFO.f3 = float(kat.constants["f3"].value)
    
    kat.IFO.f36M = kat.IFO.f2 - kat.IFO.f1
        
    # TODO add else here!
    # check modultion frequencies
    if (5 * kat.IFO.f1 != kat.IFO.f2):
        print(" ** Warning: modulation frequencies do not match: 5*f1!=f2")
    
    # defining a dicotionary for the main mirror positions (tunings),
    # keys should include maxtem, phase and all main optics names
    #kat.IFO.tunings = get_tunings(dict.fromkeys(["maxtem", "phase", "PRM", "ITMX", "ETMX", "ITMY", "ETMY", "BS", "SRM"]))
    kat.IFO.compute_derived_lengths()
        
    # ----------------------------------------------------------------------
    # define ports and signals 
    
    # useful ports
    kat.IFO.POP_f1  = Output(kat.IFO, "POP_f1",  "nPOP",  "f1", phase=101)
    kat.IFO.POP_f2  = Output(kat.IFO, "POP_f2",  "nPOP",  "f2", phase=13)
    kat.IFO.REFL_f1 = Output(kat.IFO, "REFL_f1", "nREFL", "f1", phase=101)
    kat.IFO.REFL_f2 = Output(kat.IFO, "REFL_f2", "nREFL", "f2", phase=14)
    
        
    kat.IFO.AS_f1  = Output(kat.IFO, "AS_f1",  "nSRM2", "f1", phase=101)
    kat.IFO.AS_f2  = Output(kat.IFO, "AS_f2",  "nSRM2", "f2", phase=14)
    kat.IFO.AS_f36 = Output(kat.IFO, "AS_f36", "nSRM2", "f36M", phase=14)
    
    # AS_DC refers to what is coming out of the OMC now
    kat.IFO.AS_DC   = Output(kat.IFO, "AS_DC", "nOMC_OCc")
    kat.IFO.POW_BS  = Output(kat.IFO, "PowBS", "nPRBS*")
    kat.IFO.POW_X   = Output(kat.IFO, "PowX",  "nITMX2")
    kat.IFO.POW_Y   = Output(kat.IFO, "PowY",  "nITMY2")
    kat.IFO.POW_SRC = Output(kat.IFO, "PowSRC",  "nSRM1")
    kat.IFO.POW_PRC = Output(kat.IFO, "PowPRC",  "nPRM2")

    # pretune LSC DOF
    kat.IFO.preARMX =  DOF(kat.IFO, "ARMX", kat.IFO.POW_X,   "", "ETMX", 1, 1.0, sigtype="z")
    kat.IFO.preARMY =  DOF(kat.IFO, "ARMY", kat.IFO.POW_Y,   "", "ETMY", 1, 1.0, sigtype="z")
    kat.IFO.preMICH =  DOF(kat.IFO, "AS"  , kat.IFO.AS_DC,   "", ["ITMX", "ETMX", "ITMY", "ETMY"], [1,1,-1,-1], 6.0, sigtype="z")
    kat.IFO.prePRCL =  DOF(kat.IFO, "PRCL", kat.IFO.POW_BS,  "", "PRM",  1, 10.0, sigtype="z")
    kat.IFO.preSRCL =  DOF(kat.IFO, "SRCL", kat.IFO.AS_DC,   "", "SRM",  1, 10.0, sigtype="z")
    
    # control scheme as in [1] Table C.1. Due to Finesse
    # conventions, the overall factor for all but PRCL are multiplied by -1
    # compared to the LIGO defintion, to match the same defintion. 
    kat.IFO.PRCL =  DOF(kat.IFO, "PRCL", kat.IFO.POP_f1,  "I", "PRM", 1, 100.0, sigtype="z")
    kat.IFO.MICH =  DOF(kat.IFO, "MICH", kat.IFO.POP_f2,  "Q", ["ITMX", "ETMX", "ITMY", "ETMY"], [-0.5,-0.5,0.5,0.5], 100.0, sigtype="z")
    kat.IFO.CARM =  DOF(kat.IFO, "CARM", kat.IFO.REFL_f1, "I", ["ETMX", "ETMY"], [-1, -1], 1.5, sigtype="z")
    kat.IFO.DARM =  DOF(kat.IFO, "DARM", kat.IFO.AS_f2,   "Q", ["ETMX", "ETMY"], [-1,1], 1.0, sigtype="z")
    kat.IFO.SRCL =  DOF(kat.IFO, "SRCL", kat.IFO.REFL_f2, "I", "SRM", -1, 1e2, sigtype="z")

    kat.IFO.DARM_h =  DOF(kat.IFO, "DARM_h", None, "", ["LY", "LX"], [-1,1], 1.0, sigtype="phase")

    kat.IFO.LSC_DOFs = (kat.IFO.PRCL, kat.IFO.MICH, kat.IFO.CARM, kat.IFO.DARM, kat.IFO.SRCL)
    kat.IFO.CAV_POWs = (kat.IFO.POW_X, kat.IFO.POW_Y, kat.IFO.POW_BS)
    
    # Pitch DOfs
    # There is a difference in the way LIGO and Finesse define positive and negative
    # rotations of the cavity mirrors. For LIGO the rotational DOFs assume ITM + rotation
    # is clockwise and ETM + rotation is anticlockwise.
    # I'll be explict here for future reference.
    cav_mirrors = ["ETMX", "ETMXAR", "ETMY", "ETMYAR", "ITMX", "ITMXAR", "ITMY", "ITMYAR"]

    # LIGO definitions
    # Based on figure 7 in T0900511-v4
    CHARD_factors   = np.array([ 1, 1, 1, 1,-1,-1,-1,-1])
    DHARD_factors   = np.array([ 1, 1,-1,-1,-1,-1, 1, 1])
    CSOFT_factors   = np.array([-1,-1,-1,-1,-1,-1,-1,-1])
    # DSOFT_factors   = np.array([-1,-1, 1, 1, 1, 1,-1,-1])   # Wrong!
    DSOFT_factors   = np.array([-1,-1, 1, 1,-1,-1, 1, 1])
    
    # Finesse definitions
    # negative for ITM rotations
    ITMS = np.in1d(cav_mirrors, np.array(["ITMX", "ITMXAR", "ITMY", "ITMYAR"]))
    CHARD_factors[ITMS] *= -1
    DHARD_factors[ITMS] *= -1
    CSOFT_factors[ITMS] *= -1
    DSOFT_factors[ITMS] *= -1

    kat.IFO.CHARD_P = DOF(kat.IFO, "CHARD_P", None , None, cav_mirrors, CHARD_factors, 1, sigtype="pitch")
    kat.IFO.DHARD_P = DOF(kat.IFO, "DHARD_P", None , None, cav_mirrors, DHARD_factors, 1, sigtype="pitch")
    kat.IFO.CSOFT_P = DOF(kat.IFO, "CSOFT_P", None , None, cav_mirrors, CSOFT_factors, 1, sigtype="pitch")
    kat.IFO.DSOFT_P = DOF(kat.IFO, "DSOFT_P", None , None, cav_mirrors, DSOFT_factors, 1, sigtype="pitch")
    kat.IFO.PRM_P   = DOF(kat.IFO, "PRM_P"  , None , None, ["PRM", "PRMAR"], [1,1], 1, sigtype="pitch")
    kat.IFO.PRC2_P  = DOF(kat.IFO, "PRC2_P" , None , None, ["PR2"], [1], 1, sigtype="pitch")
    kat.IFO.PRC3_P  = DOF(kat.IFO, "PRC3_P" , None , None, ["PR3"], [1], 1, sigtype="pitch")
    kat.IFO.SRM_P   = DOF(kat.IFO, "SRM_P"  , None , None, ["SRM", "SRMAR"], [1,1], 1, sigtype="pitch")
    kat.IFO.SRC2_P  = DOF(kat.IFO, "SRC2_P" , None , None, ["SR2"], [1], 1, sigtype="pitch")
    kat.IFO.SRC3_P  = DOF(kat.IFO, "SRC3_P" , None , None, ["SR3"], [1], 1, sigtype="pitch")
    kat.IFO.MICH_P  = DOF(kat.IFO, "MICH_P" , None , None, ["BS", "BSAR1", "BSAR2"], [1,1,1], 1, sigtype="pitch")
    
    kat.IFO.ASC_P_DOFs = (kat.IFO.CHARD_P, kat.IFO.DHARD_P,
                          kat.IFO.CSOFT_P, kat.IFO.DSOFT_P,
                          kat.IFO.PRM_P, kat.IFO.PRC2_P,
                          kat.IFO.PRC3_P, kat.IFO.SRM_P,
                          kat.IFO.SRC2_P, kat.IFO.SRC3_P,
                          kat.IFO.MICH_P)
    
    kat.IFO.update()

    kat.IFO.lockNames = None
    
    return kat
    

    
def scan_to_precision(kat, DOF, pretune_precision, minmax="max", phi=0.0, precision=60.0):
    assert_voyager_ifo_kat(kat)
    
    while precision > pretune_precision * DOF.scale:
        out = scan_DOF(kat, DOF, xlimits = [phi-1.5*precision, phi+1.5*precision])
        phi, precision = find_peak(out, DOF.port.name, minmax=minmax)
        
    return phi, precision
    
    
def pretune(_kat, pretune_precision=1.0e-4, verbose=False):
    assert_voyager_ifo_kat(_kat)
    
    # This function needs to apply a bunch of pretunings to the original
    # kat and associated IFO object passed in
    IFO = _kat.IFO
    
    print("-- pretuning interferometer to precision {0:2g} deg = {1:2g} m".format(pretune_precision, pretune_precision*_kat.lambda0/360.0))
    
    kat = _kat.deepcopy()
    kat.removeBlock("locks", False)
    
    vprint(verbose, "   scanning X arm (maximising power)")
    
    make_transparent(kat, ["PRM", "SRM"])
    make_transparent(kat, ["ITMY", "ETMY"])
    
    kat.BS.setRTL(0.0, 1.0, 0.0) # set BS refl. for X arm
    
    phi, precision = scan_to_precision(kat, IFO.preARMX, pretune_precision)
    phi = round(phi/pretune_precision)*pretune_precision
    phi = round_to_n(phi,5)
    
    vprint(verbose, "   found max/min at: {} (precision = {:2g})".format(phi, precision))
    
    IFO.preARMX.apply_tuning(phi)

    vprint(verbose, "   scanning Y arm (maximising power)")
    kat = _kat.deepcopy()
    kat.removeBlock("locks", False)
    
    make_transparent(kat,["PRM","SRM"])
    make_transparent(kat,["ITMX", "ETMX"])
    kat.BS.setRTL(1.0,0.0,0.0) # set BS refl. for Y arm
    phi, precision = scan_to_precision(kat, IFO.preARMY, pretune_precision)
    phi=round(phi/pretune_precision)*pretune_precision
    phi=round_to_n(phi,5)
    vprint(verbose, "   found max/min at: {} (precision = {:2g})".format(phi, precision))
    IFO.preARMY.apply_tuning(phi)

    vprint(verbose, "   scanning MICH (minimising power)")
    kat = _kat.deepcopy()
    kat.removeBlock("locks", False)
    
    make_transparent(kat,["PRM","SRM"])
    phi, precision = scan_to_precision(kat, IFO.preMICH, pretune_precision, minmax="min", precision=30.0)
    phi=round(phi/pretune_precision)*pretune_precision
    phi=round_to_n(phi,5)
    vprint(verbose, "   found max/min at: {} (precision = {:2g})".format(phi, precision))
    IFO.preMICH.apply_tuning(phi, add=True)

    vprint(verbose, "   scanning PRCL (maximising power)")
    kat = _kat.deepcopy()
    kat.removeBlock("locks", False)
    make_transparent(kat,["SRM"])
    phi, precision = scan_to_precision(kat, IFO.prePRCL, pretune_precision)
    phi=round(phi/pretune_precision)*pretune_precision
    phi=round_to_n(phi,5)
    vprint(verbose, "   found max/min at: {} (precision = {:2g})".format(phi, precision))
    IFO.prePRCL.apply_tuning(phi)

    vprint(verbose, "   scanning SRCL (maximising carrier power, then adding 90 deg)")
    kat = _kat.deepcopy()
    kat.removeBlock("locks", False)
    
    phi, precision = scan_to_precision(kat, IFO.preSRCL, pretune_precision, phi=0, precision = 10)
    phi=round(phi/pretune_precision)*pretune_precision
    phi=round_to_n(phi,4)-90.0
    
    vprint(verbose, "   found max/min at: {} (precision = {:2g})".format(phi, precision))
    IFO.preSRCL.apply_tuning(phi)
    
    print("   ... done")
    


def pretune_status(_kat):
    assert_voyager_ifo_kat(_kat)
    
    kat = _kat.deepcopy()
    kat.verbose = False
    kat.noxaxis = True
    kat.yaxis = "abs"
    
    pretune_DOFs = [kat.IFO.preARMX, kat.IFO.preARMY, kat.IFO.prePRCL, kat.IFO.preMICH, kat.IFO.preSRCL]
    
    _detStr=""
    
    for dof in pretune_DOFs:
        dof.add_signal()
        
    out = kat.run()
    Pin = float(kat.L0.P)

    tunings = kat.IFO.get_tunings()
    
    if tunings['keys']["maxtem"] == -1:
        _maxtemStr="off"
    else:
        _maxtemStr = "{:3}".format(tunings['keys']["maxtem"])
        
    print(" .--------------------------------------------------.")
    print(" | pretuned for maxtem = {}, phase = {:2}            |".format(_maxtemStr, int(kat.phase)))
    
    keys_t = list(tunings.keys())
    keys_t.remove("keys")
    
    print(" .--------------------------------------------------.")
    print(" | port   power[W] pow. ratio | optics   tunings    |")
    print(" +----------------------------|---------------------+")
    
    idx_p = 0
    idx_t = 0
    
    while (idx_p < len(pretune_DOFs) or idx_t < len(keys_t)):
        if idx_p < len(pretune_DOFs):
            p = pretune_DOFs[idx_p]
            print(" | {:5}: {:9.4g} {:9.4g} |".format(p.name, float(out[p.port.name]), float(out[p.port.name])/Pin),end="")
            idx_p +=1
        else:
            print(" |                            |", end="")
            
        if idx_t < len(keys_t):
            t=keys_t[idx_t]
            print(" {:5}: {:9.3g}    |".format(t, float(tunings[t])))
            idx_t +=1
        else:
            print("                     |")
            
    print(" `--------------------------------------------------'")

# probably extra and can be removed
def power_ratios(_kat):
    assert_voyager_ifo_kat(_kat)
    
    kat = _kat.deepcopy()
    kat.verbose = False
    kat.noxaxis = True

    ports = [kat.IFO.POW_X, kat.IFO.POW_Y, kat.IFO.AS_DC, kat.IFO.POW_BS]
    _detStr = ""
    
    for p in ports:
        _sigStr = p.signal(kat)
        _detStr = "\n".join([_detStr, _sigStr])
    
    kat.parse(_detStr)
    
    out = kat.run()
    
    Pin = float(kat.L0.P)

    print("-- power ratios (Pin = {0:.3g} W)".format(Pin))
    
    for p in ports:
        print(" {0:6} = {1:8.3g} W ({0:6}/Pin = {2:8.2g})" .format(p.name, float(out[p.name]), float(out[p.name])/Pin))


def generate_locks(kat, gainsAdjustment = [0.5, 0.005, 1.0, 0.5, 0.025],
                    gains=None, accuracies=None,
                    rms=[1e-13, 1e-13, 1e-12, 1e-11, 50e-11], verbose=True,
                    useDiff = True):
    """
    gainsAdjustment: factors to apply to loop gains computed from optical gains
    gains:           override loop gain [W per deg]
    accuracies:      overwrite error signal threshold [W]
    useDiff:         use diff command instead of fsig to compute optical gains
                    
    rms: loop accuracies in meters (manually tuned for the loops to work
         with the default file)
         to compute accuracies from rms, we convert
         rms to radians as rms_rad = rms * 2 pi/lambda
         and then multiply by the optical gain.
                    
    NOTE: gainsAdjustment, gains, accuracies and rms are specified in the order of DARM, CARM, PRCL, MICH, SRCL.
    """
    assert_voyager_ifo_kat(kat)
        
    # optical gains in W/rad
    
    ogDARM = optical_gain(kat.IFO.DARM, kat.IFO.DARM, useDiff=useDiff)
    ogCARM = optical_gain(kat.IFO.CARM, kat.IFO.CARM, useDiff=useDiff)
    ogPRCL = optical_gain(kat.IFO.PRCL, kat.IFO.PRCL, useDiff=useDiff)
    ogMICH = optical_gain(kat.IFO.MICH, kat.IFO.MICH, useDiff=useDiff)
    ogSRCL = optical_gain(kat.IFO.SRCL, kat.IFO.SRCL, useDiff=useDiff)

    if gains is None:            
        # manually tuning relative gains
        factor = -1.0 * 180 / math.pi # convert from rad/W to -1 * deg/W
        
        gainDARM = round_to_n(gainsAdjustment[0] * factor / ogDARM, 2) # manually tuned
        gainCARM = round_to_n(gainsAdjustment[1] * factor / ogCARM, 2) # factor 0.005 for better gain hirarchy with DARM
        gainPRCL = round_to_n(gainsAdjustment[2] * factor / ogPRCL, 2) # manually tuned
        gainMICH = round_to_n(gainsAdjustment[3] * factor / ogMICH, 2) # manually tuned
        gainSRCL = round_to_n(gainsAdjustment[4] * factor / ogSRCL, 2) # gain hirarchy with MICH
        
        gains = [ gainDARM, gainCARM, gainPRCL, gainMICH, gainSRCL]
    
    if accuracies is None:
        factor = 2.0 * math.pi / kat.lambda0 # convert from m to radians
        
        accDARM = round_to_n(np.abs(factor * rms[0] * ogDARM), 2) 
        accCARM = round_to_n(np.abs(factor * rms[1] * ogCARM), 2)
        accPRCL = round_to_n(np.abs(factor * rms[2] * ogPRCL), 2)
        accMICH = round_to_n(np.abs(factor * rms[3] * ogMICH), 2)
        accSRCL = round_to_n(np.abs(factor * rms[4] * ogSRCL), 2)

        accuracies = [accDARM, accCARM, accPRCL, accMICH, accSRCL]
        
    factor1 = 2.0 * math.pi / 360.0 
    factor2 = 2.0 * math.pi / kat.lambda0 
    factor3 = 360.0  / kat.lambda0
    factor4 = -1.0 * 180 / math.pi 

    if verbose:
        print(" .--------------------------------------------------.")
        print(" | Parameters for locks:                            |")
        print(" +--------------------------------------------------+")
        print(" | -- optical gains [W/rad], [W/deg] and [W/m]:     |")
        print(" | DARM: {:12.5}, {:12.5}, {:12.5}   |".format(ogDARM, ogDARM*factor1, ogDARM*factor2))
        print(" | CARM: {:12.5}, {:12.5}, {:12.5}   |".format(ogCARM, ogCARM*factor1, ogCARM*factor2))
        print(" | PRCL: {:12.5}, {:12.5}, {:12.5}   |".format(ogPRCL, ogPRCL*factor1, ogPRCL*factor2))
        print(" | MICH: {:12.5}, {:12.5}, {:12.5}   |".format(ogMICH, ogMICH*factor1, ogMICH*factor2))
        print(" | SRCL: {:12.5}, {:12.5}, {:12.5}   |".format(ogSRCL, ogSRCL*factor1, ogSRCL*factor2))
        print(" +--------------------------------------------------+")
        print(" | -- defult loop accuracies [deg], [m] and [W]:    |")
        print(" | DARM: {:12.6}, {:12.6}, {:12.6}   |".format(factor3*rms[0], rms[0], np.abs(rms[0]*ogDARM*factor2)))
        print(" | CARM: {:12.6}, {:12.6}, {:12.6}   |".format(factor3*rms[1], rms[1], np.abs(rms[1]*ogCARM*factor2)))
        print(" | PRCL: {:12.6}, {:12.6}, {:12.6}   |".format(factor3*rms[2], rms[2], np.abs(rms[2]*ogPRCL*factor2)))
        print(" | MICH: {:12.6}, {:12.6}, {:12.6}   |".format(factor3*rms[3], rms[3], np.abs(rms[3]*ogMICH*factor2)))
        print(" | SRCL: {:12.6}, {:12.6}, {:12.6}   |".format(factor3*rms[4], rms[4], np.abs(rms[4]*ogSRCL*factor2)))
        print(" +--------------------------------------------------+")
        print(" | -- extra gain factors (factor * 1/optical_gain): |")
        print(" | DARM: {:5.4} * {:12.6} = {:12.6}        |".format(gainsAdjustment[0],factor4/ogDARM, gainsAdjustment[0]*factor4/ogDARM))
        print(" | CARM: {:5.4} * {:12.6} = {:12.6}        |".format(gainsAdjustment[1],factor4/ogCARM, gainsAdjustment[1]*factor4/ogCARM))
        print(" | PRCL: {:5.4} * {:12.6} = {:12.6}        |".format(gainsAdjustment[2],factor4/ogPRCL, gainsAdjustment[2]*factor4/ogPRCL))
        print(" | MICH: {:5.4} * {:12.6} = {:12.6}        |".format(gainsAdjustment[3],factor4/ogMICH, gainsAdjustment[3]*factor4/ogMICH))
        print(" | SRCL: {:5.4} * {:12.6} = {:12.6}        |".format(gainsAdjustment[4],factor4/ogSRCL, gainsAdjustment[4]*factor4/ogSRCL))
        print(" `--------------------------------------------------'")
        
    data = {
        "DARM": {"accuracy": accuracies[0], "gain": gains[0]},
        "CARM": {"accuracy": accuracies[1], "gain": gains[1]},
        "PRCL": {"accuracy": accuracies[2], "gain": gains[2]},
        "MICH": {"accuracy": accuracies[3], "gain": gains[3]},
        "SRCL": {"accuracy": accuracies[4], "gain": gains[4]},
    }
    
    return data
