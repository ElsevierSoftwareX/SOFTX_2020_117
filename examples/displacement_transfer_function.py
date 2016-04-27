# -*- coding: utf-8 -*-
"""
Transfer function of mirror rotation to photodiode signal for a Fabry-Perot
cavity.

      ]------=----/--(=========================)
                  |
                  Y

    laser   EOM      ITM      10m cavity      ETM

              Photodiode

This simulation sets up an optical environment in FINESSE and then creates
a signal to rotate the ETM. This signal is then attached to the x-axis and
added as a demodulation frequency to the photodiode to allow for a transfer
function to be obtained.

The magnitude and phase are then plotted, and the figure is saved as a PDF.

Sean Leavey
s.leavey.1@research.gla.ac.uk

February 2014
"""

import pykat
import pylab

# instantiate PyKat object
kat = pykat.finesse.kat()

# create a FINESSE environment with an impedence matched Fabry-Perot cavity of length 10 m
kat.parseCommands("""
l laser 30.0 0.0 0.0 n_LASER_OUT
s space_laser_to_eom 1.0 n_LASER_OUT n_EOM_IN
mod eom 10M 0.3 2 pm 0.0 n_EOM_IN n_EOM_OUT
s space_EOM_to_BS 0.01 n_EOM_OUT n_BS_A
bs beamsplitter 0.1 0.9 0.0 0.0 n_BS_A dump n_BS_C n_BS_D
s space_BS_to_ITM 0.1 n_BS_C n_ITM_IN
m M_ITM_AR 0.001 0.999 0.0 n_ITM_IN n_ITM_BULK_IN
attr M_ITM_AR Rcx -1.7763
attr M_ITM_AR Rcy -1.7763
attr M_ITM_AR r_ap 0.04645
s M_ITM_BULK 0.027 1.45 n_ITM_BULK_IN n_ITM_BULK_OUT
m M_ITM_HR 0.995 0.005 0.0 n_ITM_BULK_OUT n_ITM_OUT
attr M_ITM_HR Rcx -5.7
attr M_ITM_HR Rcy -5.7
attr M_ITM_HR r_ap 0.04645
s space_cavity 10 n_ITM_OUT n_ETM_IN
m M_ETM_HR 0.995 0.005 0.0 n_ETM_IN n_ETM_BULK_IN
attr M_ETM_HR Rcx 5.7
attr M_ETM_HR Rcy 5.7
attr M_ETM_HR r_ap 0.04645
s M_ETM_BULK 0.027 1.45 n_ETM_BULK_IN n_ETM_BULK_OUT
m M_ETM_AR 0.001 0.999 0.0 n_ETM_BULK_OUT dump
attr M_ETM_AR Rcx 1.7763
attr M_ETM_AR Rcy 1.7763
attr M_ETM_AR r_ap 0.04645
pd2 PD_PDH_i 10M max 1 n_BS_D
cav cavity M_ITM_HR n_ITM_OUT M_ETM_HR n_ETM_IN
""")

##############################
# define what we want to see #
##############################

# set up a mirror rotation signal
kat.signals.f = 1
kat.signals.apply(kat.M_ETM_HR.phi, 1, 180) # amplitude = 1, phase = 180 degrees

# manually add yaxis command, as this is not implemented yet as of PyKat 0.3.1
# this plots both magnitude and phase
kat.parseKatCode("""
yaxis lin abs:deg
""")

# set the x-axis to plot frequency between 1 Hz and 1 MHz
kat.add(pykat.commands.xaxis('log', [1, 1e6], kat.signals.f, 1000))

# now, we need to demodulate the photodiode at the frequency of the signal at each point on the x-axis
kat.PD_PDH_i.f2.put(kat.xaxis.x)

# turn off higher order TEM modes, we don't need these for this transfer function
kat.maxtem = "off"

##################
# run simulation #
##################

result = kat.run()

###########
# plot it #
###########

# create figure with 2 subplots, sharing the x-axis
figure, axes = pylab.subplots(2, sharex = True)

# plot amplitude on first subplot
axes[0].loglog(result.x, result.y[:, 0])
axes[0].set_title('Amplitude')
axes[0].set_ylabel(result.ylabels[0])
axes[0].grid(True)

# plot phase on second subplot
axes[1].semilogx(result.x, result.y[:, 1])
axes[1].set_title('Phase')
axes[1].set_xlabel(result.xlabel)
axes[1].set_ylabel(result.ylabels[1])
axes[1].set_ylim([10, -100])
axes[1].grid(True)

# space everything properly
figure.tight_layout()

# save figure as a PDF
pylab.savefig('displacement_transfer_function', format = 'PDF')

# show
pylab.show()
