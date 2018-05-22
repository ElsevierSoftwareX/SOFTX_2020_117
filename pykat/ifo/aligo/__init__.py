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

class ALIGO_IFO(IFO):   
    """
    This contains aLIGO specific methods for computing interferometer
    variables.
    
    Functions that operate on the kat/IFO objects, manipulate their
    structure and return a new one, should not be included here. They should
    be separate functions that are called by the user. 
    
    The functions here should be those that update the kat with information
    from the IFO object or vice-versa.
    """

    def __init__(self, kat, tuning_keys_list, tunings_components_list):
        IFO.__init__(self, kat, tuning_keys_list, tunings_components_list)
        self._f36M = np.nan
    
    @property
    def DCoffset(self):
        if 'DCoffset' not in self.kat.data:
            return 0
        else:
            return float(self.kat.data['DCoffset'])

    @DCoffset.setter
    def DCoffset(self, value):
        self.kat.data['DCoffset'] = float(value)
    
    @property
    def DCoffsetW(self):
        if 'DCoffsetW' not in self.kat.data:
            return 0
        else:
            return float(self.kat.data['DCoffsetW'])

    @DCoffsetW.setter
    def DCoffsetW(self, value):
        self.kat.data['DCoffsetW'] = float(value)

    @property
    def f1(self):
        return self.kat.mod1.f.value
        
    @f1.setter
    def f1(self, value):
        self.kat.mod1.f.value = value
                
    @property
    def f2(self):
        return self.kat.mod2.f.value
        
    @f2.setter
    def f2(self, value):
        self.kat.mod2.f.value = value
        
    @property
    def f36M(self):
        return self.f2 - self.f1
        
    def createPorts(self):
        # useful ports
        self.POP_f1  = Output(self, "POP_f1",  "nPOP",  "f1", phase=101)
        self.POP_f2  = Output(self, "POP_f2",  "nPOP",  "f2", phase=13)
        self.REFL_f1 = Output(self, "REFL_f1", "nREFL", "f1", phase=101)
        self.REFL_f2 = Output(self, "REFL_f2", "nREFL", "f2", phase=14)
        self.AS_DC   = Output(self, "AS_DC", "nAS")
        self.POW_BS  = Output(self, "PowBS", "nPRBS*")
        self.POW_X   = Output(self, "PowX",  "nITMX2")
        self.POW_Y   = Output(self, "PowY",  "nITMY2")
    
    def compute_derived_resonances(self):
        clight
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
        attr ETMY mass 40
        attr ETMX mass 40
        attr ITMY mass 40
        attr ITMX mass 40
    
        attr PRM  mass 2.9
        attr PR2  mass 2.9
        attr PR3  mass 12
    
        attr SRM  mass 2.9
        attr SR2  mass 2.9
        attr SR3  mass 12
    
        attr BS   mass 14
        """
    
        self.kat.parse(code)
    
        return code


    def suspend_mirrors_pitch(self):
        """
        Suspends the main mirrors in an aLIGO model in pitch in the
        supplied kat object.
    
        Returns the commands used for reference.
        
        TODO: Assumes all suspensions are QUADS currently.
        """
    
        code = """
        tf2 QUAD 2.38663 0.0 {-0.0050+86.8639i,-0.0050+61.0536i,-0.0050+32.0042i,-0.0050+21.3735i,-0.0050+20.6567i,-0.0050+19.0823i,-0.0050+22.3646i,-0.0050+17.2518i,-0.0050+16.5670i,-0.0050+15.0288i,-0.0050+12.4591i,-0.0050+13.1589i,-0.0050+10.0625i,-0.0050+8.4105i,-0.0050+8.4829i,-0.0050+6.2308i,-0.0050+6.5431i,-0.0050+5.5092i,-0.0050+2.7083i,-0.0050+3.2843i,-0.0050+2.8957i,-0.0050+3.7645i,-0.0050+14.0137i,-0.0050+3.4691i} {-0.0050+86.8639i,-0.0050+61.0536i,-0.0050+32.0042i,-0.0050+21.3735i,-0.0050+20.6566i,-0.0050+19.0823i,-0.0050+17.2493i,-0.0050+16.5665i,-0.0050+22.3646i,-0.0050+15.0288i,-0.0050+12.4591i,-0.0050+13.1589i,-0.0050+9.4995i,-0.0050+8.4829i,-0.0050+5.5072i,-0.0050+6.2177i,-0.0050+6.7464i,-0.0050+6.5428i,-0.0050+2.7591i,-0.0050+2.8957i,-0.0050+3.7645i,-0.0050+14.0137i,-0.0050+3.4691i}
    
        attr ETMY iy 1 rymech QUAD
        attr ETMX iy 1 rymech QUAD
        attr ITMY iy 1 rymech QUAD
        attr ITMX iy 1 rymech QUAD
                              
        attr PRM  iy 1 rymech QUAD
        attr PR2  iy 1 rymech QUAD
        attr PR3  iy 1 rymech QUAD
                              
        attr SRM  iy 1 rymech QUAD
        attr SR2  iy 1 rymech QUAD
        attr SR3  iy 1 rymech QUAD
                              
        attr BS   iy 1 rymech QUAD
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
        
    def lengths_status(self):
        self.compute_derived_lengths()
        
        print("  .-------------------------------------------------.")
        print(" | - arm length:                                    |")
        print(" | Lx   = {:11.7}m, Ly   = {:11.7}m         |".format(float(self.kat.LX.L), float(self.kat.LY.L)))
        print(" | - small MI and recycling lengths:                | ")
        print(" | lx   = {:11.7}m, ly   = {:11.7}m         |".format(self.lx, self.ly))
        print(" | lpr  = {:11.7}m, lsr  = {:11.7}m         |".format(self.lpr, self.lsr))
        print(" | lMI  = {:11.7}m, lSchnupp = {:11.5}m     |".format(self.lMI, self.lSchnupp))
        print(" | lPRC = {:11.7}m, lSRC = {:11.7}m         |".format(self.lPRC, self.lSRC))
        print(" +--------------------------------------------------+")
        print(" | - associated cavity frequencies [Hz]:            |")
        print(" | fsrx   = {:11.8}, fsry   = {:11.8}       |".format(self.fsrX, self.fsrY))
        print(" | fsrPRC = {:11.8}, fsrSRC = {:11.8}       |".format(self.fsrPRC, self.fsrSRC))
        print(" | f1_PRC = {:11.8}                             |".format(self.f1_PRC))
        print(" | f1     = {:11.8}, f2     = {:11.9}       |".format(self.f1, self.f2))
        print(" `-------------------------------------------------'")
    
    def remove_modulators(self):
        """
        Removes the input modulators and reconnects the input laser to the PRC reflection node.
        
        This function alters the kat object directly.
        """
        comps, nodes = self.kat.nodes.getComponentsBetween(self.kat.L0.nodes[0].name, self.kat.mod2.nodes[1].name, getNodes=True)
        
        for i,c in enumerate(comps):
            if c.name == "lmod1":
                i -= 1
                break
        
        self.kat.remove("lmod1", "mod1", "lmod2", "mod2")    # Remove modulators
                
        # Set output node of laser block to be on the laser
        self.kat.nodes.replaceNode(comps[i], nodes[i][-1].name, 'nLaserOut')
        
            
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
        """
        Method for removing the OMC and the FI blocks in kat-objects having these
        included. The FI block contains an ideal Faraday isolator as well as the
        path from the isolator to the OMC, which is used to mode match the OMC to
        the interferometer. The node nAS is re-set such that it always corresponds
        to the "last" node of the output path (dark port, asymmetric port, etc). 

        Parameters
        ----------
        removeFI  : Boolean
                    If True, the Faraday isolator is removed along with the path
                    to the OMC.
        removeOMC : Boolean
                    If True, the OMC is removed. Must be True if removeFI = True.
        """
        
        if removeFI and not removeOMC:
            raise pkex.BasePyKatException("Must remove OMC if removing FI")
        if removeFI:
            self.kat.nodes.replaceNode(self.kat.sSRM_FI, 'nFI2a', 'nAS')
            self.kat.removeBlock('FI')
            self.kat.removeBlock('OMC')
            self.kat.cavOMC.remove()
        elif removeOMC:
            self.kat.nodes.replaceNode(self.kat.sOM3_OMC, 'nOMC_ICa', 'nAS')
            self.kat.removeBlock('OMC')
            self.kat.cavOMC.remove()

    def adjust_PRC_length(self, verbose=False):
        """
        Adjust PRC length so that it fulfils the requirement
        lPRC = (N+1/2) * c/(2*f1), see [1] equation C.1
        In the current design N=3.
    
        This function directly alters the lengths of the associated kat object.
        """
        kat = self.kat
        self.compute_derived_lengths()
        
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
        last = np.size(out.x)
        if idx==None and last>1:
            #print("switching to last value of out")
            idx = last-1

        if "ETMX_lock" in out:
            if idx is None:
                tuning["ETMX"] += float(out["ETMX_lock"].real)
            else:
                tuning["ETMX"] += float(out["ETMX_lock"][idx].real)
        else:
            pkex.printWarning("could not find ETMX lock")
    
        if "ETMY_lock" in out:
            if idx is None:
                tuning["ETMY"] += float(out["ETMY_lock"].real)
            else:
                tuning["ETMY"] += float(out["ETMY_lock"][idx].real)
        else:
            pkex.printWarning("could not find ETMY lock")
    
        if "PRM_lock" in out:
            if idx is None:
                tuning["PRM"]  += float(out["PRM_lock"].real)
            else:
                tuning["PRM"]  += float(out["PRM_lock"][idx].real)
        else:
            pkex.printWarning("could not find PRCL lock")
    
        if ("ITMX_lock" in out) and ("ITMY_lock" in out):
            if idx is None:
                tuning["ITMX"] += float(out["ITMX_lock"].real)
                tuning["ITMY"] += float(out["ITMY_lock"].real)
            else:
                tuning["ITMX"] += float(out["ITMX_lock"][idx].real)
                tuning["ITMY"] += float(out["ITMY_lock"][idx].real)
        else:
            pkex.printWarning("could not find MICH (ITMX, ITMY) lock")
    
        if "SRM_lock" in out:
            if idx is None:
                tuning["SRM"]  += float(out["SRM_lock"].real)
            else:
                tuning["SRM"]  += float(out["SRM_lock"][idx].real)
        else:
            pkex.printWarning("could not find SRCL lock")
            
        self.kat.IFO.apply_tunings(tuning)
    
    def set_DC_offset(self, DCoffset=None, verbose=False):
        """
        Sets the DC offset for this inteferometer.
        
        This function directly alters the tunings of the associated kat object.
        
        DCoffset - degrees of additional DC offset to apply current tunings
        """
        _kat = self.kat
        
        if DCoffset is not None:
            self.DCoffset = DCoffset
        
            vprint(verbose, "-- applying user-defined DC offset:")
        
            tunings = self.get_tunings()
        
            tunings["ETMY"] += self.DCoffset
            tunings["ETMX"] -= self.DCoffset
        
            self.apply_tunings(tunings)        
        
            # Compute the DC offset powers
            kat = _kat.deepcopy()
        
            signame = kat.IFO.AS_DC.add_signal()
            kat.parse("ad adp00 0 0 0 nAS")
        
            kat.noxaxis=True
        
            out = kat.run(cmd_args=["-cr=on"])
            
            TEM00_DC = abs(out['adp00'])**2
            
            _kat.IFO.DCoffsetW = float(out[signame])
        else:
            # Finding light power in AS port (mostly due to RF sidebands now)
            kat = _kat.deepcopy()
        
            signame = kat.IFO.AS_DC.add_signal()
        
            kat.noxaxis=True
        
            out = kat.run()
        
            print("-- adjusting DCoffset based on light in dark port:")
        
            waste_light = round(float(out[signame]),1)
        
            print("   waste light in AS port of {:2} W".format(waste_light))
        
            #kat_lock = _kat.deepcopy()
        
            self.find_DC_offset(_kat, 2*waste_light)
        
            _kat.parse("ad adp00 0 0 0 nAS")
            _kat.noxaxis=True
        
            out = _kat.run(cmd_args=["-cr=on"])
            
            TEM00_DC = abs(out['adp00'])**2
            
        vprint(verbose, "   DCoffset = {:6.4} deg ({:6.4}m)".format(self.DCoffset, self.DCoffset / 360.0 * _kat.lambda0 ))
        vprint(verbose, "   at dark port power: {:6.4}W".format(self.DCoffsetW))
        vprint(verbose, "   at dark port power (TEM00 0Hz): {:6.4}W".format(TEM00_DC))

    def find_DC_offset(self, AS_power, precision=1e-4, verbose=False):
        """
        Returns the DC offset of DARM that corresponds to the specified power in the AS power.
        
        This function directly alters the tunings of the associated kat object.
        """
        vprint(verbose, "   finding DC offset for AS power of {:3g} W".format(float(AS_power)))
    
        _kat = self.kat
        
        kat = _kat.deepcopy()
        kat.verbose = False
        kat.noxaxis = True
        
        kat.removeBlock("locks", False)
        kat.removeBlock("errsigs", False)
        
        kat.IFO.AS_DC.add_signal()
    
        Xphi = float(kat.ETMX.phi)
        Yphi = float(kat.ETMY.phi)

        def powerDiff(phi):
            kat.ETMY.phi = Yphi + phi
            kat.ETMX.phi = Xphi - phi
        
            out = kat.run()
            print("   ! ", out[self.AS_DC.get_signal_name()], phi)
            
            return np.abs(out[self.AS_DC.get_signal_name()] - AS_power)

        vprint(verbose, "   starting peak search...")
        out = fmin(powerDiff, 0, xtol=precision, ftol=1e-3, disp=verbose)
    
        vprint(verbose, "   ... done")
        vprint(verbose, "   DC offset for AS_DC={} W is: {}".format(AS_power, out[0]))
    
        self.DCoffset = round(out[0], 6)
        self.DCoffsetW = AS_power
    
        tunings = self.get_tunings()
        tunings["ETMY"] += self.DCoffset
        tunings["ETMX"] -= self.DCoffset
    
        self.apply_tunings(tunings)
        
        return self.DCoffset

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
        
            # code3 = """
            #         noplot {}
            #         noplot {}
            #         noplot {}
            #         noplot {}
            #         noplot {}""".format(nameDARM, nameCARM, namePRCL, nameMICH, nameSRCL).replace("  ","")
                    
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

        if "DC" in self.kat.IFO.DARM.port.name:
            code2 = ("lock DARM_lock $DARM_err {:8.2} {:8.2} {DC}\n"
                     "lock CARM_lock $CARM_err {:8.2g} {:8.2g}\n"
                     "lock PRCL_lock $PRCL_err {:8.2g} {:8.2g}\n"
                     "lock MICH_lock $MICH_err {:8.2g} {:8.2g}\n"
                     "lock SRCL_lock $SRCL_err {:8.2g} {:8.2g}\n"
                     ).format(*chain.from_iterable(zip(gains, accuracies)),
                              DC=-self.kat.IFO.DCoffsetW)
        else:
            code2 = ("lock DARM_lock $DARM_err {:8.2} {:8.2}\n"
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
        
        self.kat.IFO.ASC_REFL9A   = Output(self.kat.IFO, "ASC_REFL9A",  "nREFL_WFS_A",  "f1", block="REFL_gouy_tele")
        self.kat.IFO.ASC_REFL9B   = Output(self.kat.IFO, "ASC_REFL9B",  "nREFL_WFS_B",  "f1", block="REFL_gouy_tele")

        self.kat.IFO.ASC_REFL45A  = Output(self.kat.IFO, "ASC_REFL45A",  "nREFL_WFS_A",  "f2", block="REFL_gouy_tele")
        self.kat.IFO.ASC_REFL45B  = Output(self.kat.IFO, "ASC_REFL45B",  "nREFL_WFS_B",  "f2", block="REFL_gouy_tele")
        
        self.kat.IFO.ASC_REFL36A  = Output(self.kat.IFO, "ASC_REFL36A",  "nREFL_WFS_A",  "f36M", block="REFL_gouy_tele")
        self.kat.IFO.ASC_REFL36B  = Output(self.kat.IFO, "ASC_REFL36B",  "nREFL_WFS_B",  "f36M", block="REFL_gouy_tele")
        
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
            
def assert_aligo_ifo_kat(kat):
    if not isinstance(kat.IFO, ALIGO_IFO):
        raise pkex.BasePyKatException("\033[91mkat file is not an ALIGO_IFO compatiable kat\033[0m")
              
def make_kat(name="design", katfile=None, verbose = False, debug=False, use_RF_DARM_lock=False,
             keepComments=False, preserveConstants=False):
    """
    Returns a kat object and fills in the kat.IFO property for storing
    the associated interferometer data.
    
    The `name` argument selects from default aLIGO files included in Pykat:
    
        - design: A file based on the design parameters for the final aLIGO setup.
          125W input, T_SRM = 20%.
    
        - design_low_power: A file based on the design parameters for the final aLIGO setup.
          20W input, T_SRM = 35%. The higher SRM transmission mirror is used for low power
          operation. 20W input power from O1 observation.
        
        - design_with_IMC_HAM2: A file based on `design` but has the IMC and HAM2 blocks
          which contain design parameter input optics
    
        - design_with_IMC_HAM2_FI_OMC: A file with the OMC and IMC, most complete file
    
    keepComments: If true it will keep the original comments from the file
    preserveComments: If true it will keep the const commands in the kat
    """
    names = ['design', 'design_low_power', 'design_with_IMC_HAM2', 'design_with_IMC_HAM2_FI_OMC']
    
    if debug:
        kat = finesse.kat(tempdir=".",tempname="test")
    else:
        kat = finesse.kat()
    
    kat.verbose=verbose
    
    # Create empty object to just store whatever DOFs, port, variables in
    # that will be used by processing functions
    kat.IFO = ALIGO_IFO(kat,
                        # Define which keys are used for a tuning description
                        ["maxtem", "phase"],
                        # Define which mirrors create the tuning description
                        ["PRM", "ITMX", "ETMX", "ITMY", "ETMY", "BS", "SRM"])
    
    kat.IFO._data_path=pkg_resources.resource_filename('pykat.ifo', os.path.join('aligo','files'))

    kat.IFO.rawBlocks = BlockedKatFile()
    
    if katfile:
        kat.load(katfile, keepComments=keepComments, preserveConstants=preserveConstants)
        kat.IFO.rawBlocks.read(katfile)
    else:
        if name not in names:
            pkex.printWarning("aLIGO name `{}' not recognised, options are {}, "
                              "using default 'design'".format(name, names))
        
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
    
    # If we don't have an OMC then we need to attach
    # directly to the AS node. Otherwise use OMC refl
    if "OMC" in kat.getBlocks():
        nAS_RF = ["nOMC_ICb","nAS"] # Output() class accepts list of node names and match to the first one it finds
    else:
        nAS_RF = "nAS"
        
    kat.IFO.AS_f1  = Output(kat.IFO, "AS_f1",  nAS_RF,  "f1", phase=101)
    kat.IFO.AS_f2  = Output(kat.IFO, "AS_f2",  nAS_RF,  "f2", phase=14)
    kat.IFO.AS_f36 = Output(kat.IFO, "AS_f36", nAS_RF, "f36M", phase=14)

    kat.IFO.AS_DC   = Output(kat.IFO, "AS_DC", "nAS")
    kat.IFO.POW_BS  = Output(kat.IFO, "PowBS", "nPRBS*")
    kat.IFO.POW_X   = Output(kat.IFO, "PowX",  "nITMX2")
    kat.IFO.POW_Y   = Output(kat.IFO, "PowY",  "nITMY2")
    kat.IFO.TRX     = Output(kat.IFO, "TRX",   "nETMX2")
    kat.IFO.TRY     = Output(kat.IFO, "TRY",  " nETMY2")

    # pretune LSC DOF
    kat.IFO.preARMX  =  DOF(kat.IFO, "ARMX", kat.IFO.POW_X,   "", "ETMX", 1, 1.0, sigtype="z")
    kat.IFO.preARMY  =  DOF(kat.IFO, "ARMY", kat.IFO.POW_Y,   "", "ETMY", 1, 1.0, sigtype="z")
    kat.IFO.preMICH  =  DOF(kat.IFO, "AS"  , kat.IFO.AS_DC,   "", ["ITMX", "ETMX", "ITMY", "ETMY"], [1,1,-1,-1], 6.0, sigtype="z")
    kat.IFO.prePRCL  =  DOF(kat.IFO, "PRCL", kat.IFO.POW_BS,  "", "PRM",  1, 10.0, sigtype="z")
    kat.IFO.preSRCL  =  DOF(kat.IFO, "SRCL", kat.IFO.AS_DC,   "", "SRM",  1, 10.0, sigtype="z")
     
    # Used by new pretuning scripts DOFs - based on lock aquisition stuff in martynov thesis
    kat.IFO._preMICH =  DOF(kat.IFO, "preMICH", kat.IFO.AS_f2,   "Q", ["ITMX", "ETMX", "ITMY", "ETMY"], [1,1,-1,-1], 1.0, sigtype="z")
    kat.IFO._preSRCL =  DOF(kat.IFO, "preSRCL", kat.IFO.AS_f2,   "I", "SRM", 1, 1.0, sigtype="z")
    kat.IFO._prePRCL =  DOF(kat.IFO, "prePRCL", kat.IFO.REFL_f1, "I", "PRM", 1, 1.0, sigtype="z")
    
    kat.IFO._preALSX =  DOF(kat.IFO, "ALSX", kat.IFO.POW_X,   "", "ETMX", 1, 1.0, sigtype="z")
    kat.IFO._preALSY =  DOF(kat.IFO, "ALSY", kat.IFO.POW_Y,   "", "ETMY", 1, 1.0, sigtype="z")
    
    # control scheme as in [1] Table C.1. Due to Finesse
    # conventions, the overall factor for all but PRCL are multiplied by -1
    # compared to the LIGO defintion, to match the same defintion. 
    kat.IFO.PRCL =  DOF(kat.IFO, "PRCL", kat.IFO.POP_f1,  "I", "PRM", 1, 100.0, sigtype="z")
    kat.IFO.MICH =  DOF(kat.IFO, "MICH", kat.IFO.POP_f2,  "Q", ["ITMX", "ETMX", "ITMY", "ETMY"], [-0.5,-0.5,0.5,0.5], 100.0, sigtype="z")
    kat.IFO.CARM =  DOF(kat.IFO, "CARM", kat.IFO.REFL_f1, "I", ["ETMX", "ETMY"], [-1, -1], 1.5, sigtype="z")
    
    if use_RF_DARM_lock:
        kat.IFO.DARM =  DOF(kat.IFO, "DARM", kat.IFO.AS_f2, "Q", ["ETMX", "ETMY"], [-1,1], 1.0, sigtype="z")
    else:
        kat.IFO.DARM =  DOF(kat.IFO, "DARM", kat.IFO.AS_DC, "",  ["ETMX", "ETMY"], [-1,1], 1.0, sigtype="z")
                            
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
        
    
def pretune(_kat, pretune_precision=1.0e-4, verbose=False, debug={}):
    assert_aligo_ifo_kat(_kat)
    
    # This function needs to apply a bunch of pretunings to the original
    # kat and associated IFO object passed in
    IFO = _kat.IFO
    
    vprint(verbose, "-- pretuning interferometer to precision {0:2g} deg = {1:2g} m".format(pretune_precision, pretune_precision*_kat.lambda0/360.0))
    
    kat = _kat.deepcopy()
    kat.removeBlock("locks", False)
    
    vprint(verbose, "   scanning X arm (maximising power)")
    
    make_transparent(kat, ["PRM", "SRM"])
    make_transparent(kat, ["ITMY", "ETMY"])
    
    kat.BS.setRTL(0.0, 1.0, 0.0) # set BS refl. for X arm
    
    phi, precision = scan_to_precision(kat.IFO.preARMX, pretune_precision)
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
    phi, precision = scan_to_precision(kat.IFO.preARMY, pretune_precision)
    phi=round(phi/pretune_precision)*pretune_precision
    phi=round_to_n(phi,5)
    vprint(verbose, "   found max/min at: {} (precision = {:2g})".format(phi, precision))
    IFO.preARMY.apply_tuning(phi)

    vprint(verbose, "   scanning MICH (minimising power)")
    kat = _kat.deepcopy()
    kat.removeBlock("locks", False)
    
    make_transparent(kat,["PRM","SRM"])
    phi, precision = scan_to_precision(kat.IFO.preMICH, pretune_precision, minmax="min", precision=30.0)
    phi=round(phi/pretune_precision)*pretune_precision
    phi=round_to_n(phi,5)
    vprint(verbose, "   found max/min at: {} (precision = {:2g})".format(phi, precision))
    IFO.preMICH.apply_tuning(phi, add=True)

    vprint(verbose, "   scanning PRCL (maximising power)")
    kat = _kat.deepcopy()
    kat.removeBlock("locks", False)
    make_transparent(kat,["SRM"])
    phi, precision = scan_to_precision(kat.IFO.prePRCL, pretune_precision)
    phi=round(phi/pretune_precision)*pretune_precision
    phi=round_to_n(phi, 5)
    vprint(verbose, "   found max/min at: {} (precision = {:2g})".format(phi, precision))
    IFO.prePRCL.apply_tuning(phi)

    vprint(verbose, "   scanning SRCL (maximising carrier power, then adding 90 deg)")
    kat = _kat.deepcopy()
    kat.removeBlock("locks", False)
    
    phi, precision = scan_to_precision(kat.IFO.preSRCL, pretune_precision, phi=0, precision=90.0, debug=("SRCL" in debug))
    phi=round(phi/pretune_precision)*pretune_precision
    phi=round_to_n(phi, 5) - 90.0
    
    vprint(verbose, "   found max/min at: {} (precision = {:2g})".format(phi, precision))
    IFO.preSRCL.apply_tuning(phi)
    
    vprint(verbose,"   ... done")
    

def pretune_2(_kat, pretune_precision=1.0e-4, verbose=False, debug={}):
    """
    This will pretune the arms and PRC cavities. This should be used in conjunction with the
    pretune_SRCL function after the DC offset has been set.
    """
    assert_aligo_ifo_kat(_kat)
    
    # This function needs to apply a bunch of pretunings to the original
    # kat and associated IFO object passed in
    IFO = _kat.IFO
    
    vprint(verbose, "-- pretuning interferometer to precision {0:2g} deg = {1:2g} m".format(pretune_precision, pretune_precision*_kat.lambda0/360.0))
    
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
    phi=round_to_n(phi, 5)
    vprint(verbose, "   found max/min at: {} (precision = {:2g})".format(phi, precision))
    IFO.prePRCL.apply_tuning(phi)
    
    vprint(verbose, "   ... done")
    
    
def pretune_SRCL(_kat, verbose=False, debug=False):
    """
    This pretunes the SRC cavity to find a reasonable operating point.
    It works slightly better with distorted interferometer models compared
    to just using pretune().
    """
    assert_aligo_ifo_kat(_kat)
    
    IFO = _kat.IFO
    
    vprint(verbose, "   scanning SRCL (finding TEM00 antiresonance)")
    base1 = _kat.deepcopy()

    # add some detectors to make a metric
    base1.removeBlock("locks", False)
    base1.parse("""
    pd Psrm nSRM1
    """)
    base1.parse(base1.IFO.SRCL.signal())

    # Set some detectors to only see the TEM00 mode
    for n in range(_kat.maxtem+1):
        for m in range(_kat.maxtem+1):
            base1.detectors[base1.IFO.SRCL.signal_name()].mask(n, m, int(n==m==0))
            base1.Psrm.mask(n, m, int(n==m==0))

    out = pykat.ifo.scan_DOF(base1, base1.IFO.SRCL, xlimits = [-10, 180], steps=1000)
    # find resonance tuning
    mx_deg_Psrm = out.x[out['Psrm'].argmax()]

    c = ((mx_deg_Psrm + 90) % 360) - 180

    out = pykat.ifo.scan_DOF(base1, base1.IFO.SRCL, xlimits = [c-180, c+180], steps=1000)
    
    # A weighting to favour the minimum value between the two SRC resonances
    # highlighting the anti-resonance
    W_Psrm     = 1+np.cos(np.deg2rad((out.x-mx_deg_Psrm)*2))

    # Error signal should be small as well as small amount of power
    SRCL_metric = abs(out["REFL_f2_I"]*out["Psrm"])
    SRCL_metric /= SRCL_metric.max()
    # favour points further away from resonance
    SRCL_metric += W_Psrm

    if debug:
        plt.plot(out.x, SRCL_metric)
        plt.plot(out.x, out["REFL_f2_I"])
    
    deg = (out.x[SRCL_metric.argmin()] % 360) - 180
    IFO.preSRCL.apply_tuning(-deg)
    
    
def setup(base, old=True, DC_offset_pm=20, verbose=False):
    """
    Runs a preparation routine to produce a LIGO model at a resonable operating point.
    This uses the pretune2 and pretune_SRCL methods which allow you 
    """
    assert_aligo_ifo_kat(base)
    
    base = base.deepcopy()
    base.verbose = False

    base.removeBlock('locks', False)
    base.removeBlock('errsigs', False)
    base.removeBlock('powers', False)

    base.phase = 2
    base.IFO.fix_mirrors()

    kat = base.deepcopy()
    kat.IFO.remove_modulators()

    if old:
        pretune(kat,pretune_precision=1e-3, verbose=verbose)
    else:
        pretune_2(kat, pretune_precision=1e-3, verbose=verbose)

    # Apply the tunings to our base kat file
    base.IFO.apply_tunings(kat.IFO.get_tunings())
    base.IFO.adjust_PRC_length(verbose=verbose)

    if verbose:
        pretune_status(base)
        base.IFO.lengths_status()

    # Set DC offset and plot
    DCoffset = DC_offset_pm*1e-12 / base.lambda0 * 180.0
    base.IFO.set_DC_offset(DCoffset=DCoffset, verbose=verbose)

    if not old:
        pretune_SRCL(base, verbose=verbose)

    errsigs_cmds = base.IFO.add_errsigs_block()

    # Generates a dictionary of the lock values to use
    locks = generate_locks(base, verbose=verbose)

    # Takes these values and then generates the commands and adds them to
    # the lock block of the kat
    lock_cmds = base.IFO.add_locks_block(locks, verbose=verbose)

    base.SRCL_lock.accuracy /= 10
    
    return base






def pretune_status(_kat):
    assert_aligo_ifo_kat(_kat)
    
    kat = _kat.deepcopy()
    kat.verbose = False
    kat.noxaxis = True
    
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
    assert_aligo_ifo_kat(_kat)
    
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
    assert_aligo_ifo_kat(kat)
        
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

def setup(base, DC_offset_pm=20, verbose=False, debug=False):
    """
    Runs a preparation routine to produce a LIGO model at a resonable operating point.
    
    Returns a copy of the base file provided with all the locks and error signals
    added.
    
    base - Base aLIGO model to tune and find an operating for
    """
    
    # Will change later when this works with the 
    old = True
    
    assert_aligo_ifo_kat(base)
    
    base = base.deepcopy()
    base.verbose = False

    base.removeBlock('locks',   False)
    base.removeBlock('errsigs', False)
    base.removeBlock('powers',  False)

    base.phase = 2
    base.IFO.fix_mirrors()

    kat = base.deepcopy()
    kat.IFO.remove_modulators()

    if old:
        pretune(kat, pretune_precision=1e-3, verbose=verbose)
    else:
        _pretune_ARM(kat, pretune_precision=1e-3, verbose=verbose)
        _pretune_MICH(kat, pretune_precision=1e-3, verbose=verbose)
        
    # Apply the tunings to our base kat file
    base.IFO.apply_tunings(kat.IFO.get_tunings())
    base.IFO.adjust_PRC_length(verbose=verbose)
    
    if not old:
        _pretune_PRCL(base, verbose=verbose, debug=debug)
    
    DCoffset = DC_offset_pm*1e-12 / base.lambda0 * 180.0
    base.IFO.set_DC_offset(DCoffset=DCoffset, verbose=verbose)
    
    if not old:
        _pretune_SRCL(base, verbose=verbose, debug=debug)
        
    if verbose:
        pretune_status(base)
        base.IFO.lengths_status()
        
    errsigs_cmds = base.IFO.add_errsigs_block()

    #Generates a dictionary of the lock values to use
    locks = generate_locks(base, verbose=verbose)

    #Takes these values and then generates the commands and adds them to
    #the lock block of the kat
    lock_cmds = base.IFO.add_locks_block(locks, verbose=verbose)
    
    return base