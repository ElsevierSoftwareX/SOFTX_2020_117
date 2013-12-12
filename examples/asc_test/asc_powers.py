import copy
from pykat import finesse

def run(tmpkat):

    kat = copy.deepcopy(tmpkat)
    
    code1 = """
    ad EOM_up 9M nEOM1
    ad EOM_low -9M nEOM1
    pd cav_pow nITM2
    ad cav_c 0 nITM2
    ad WFS1_u  9M nWFS1
    ad WFS1_l -9M nWFS1
    ad WFS1_c  0  nWFS1
    ad WFS2_u  9M nWFS2
    ad WFS2_l -9M nWFS2
    ad WFS2_c   0 nWFS2
    noxaxis
    """

    kat.parseKatCode(code1)

    out = kat.run(printout=0,printerr=0)

    code1 = code1.split("\n")
    for i in range(len(out.y)):
        print " %8s: %.4e" % (out.ylabels[i], out.y[i])
 
