# Example on running multiple kat objects at once using the parakat feature.
#
# Firstly you need to start an ipython cluster on your computer. To do this open
# a new terminal and type the command:
#
#   ipcluster start --n=4
#
# (or something similar)
# This will start a cluster with 4 workers. You should set this number to how many
# cores you have.
#
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pykat
from pykat.parallel import parakat

# Create a connection to the server
pk = parakat()

# Create a bunch of kat objects to run
kat1 = pykat.finesse.kat()
kat2 = pykat.finesse.kat()
kat3 = pykat.finesse.kat()

kat1.parseCommands("""
l l1 1 0 n0
s s1 1 n0 n1
pd P n1
noxaxis
""")

kat2.parseCommands("""
l l1 2 0 n0
s s1 1 n0 n1
pd P n1
noxaxis
""")

kat3.parseCommands("""
l l1 3 0 n0
s s1 1 n0 n1
pd P n1
noxaxis
""")

# Each kat object that is started straightaway when called like this.
# You can add as many as you want, the number that will
# run at once depends on the number of workers in the cluster
pk.run(kat1, cmd_args=["-cr=on"])
pk.run(kat2, cmd_args=["-cr=on"])
pk.run(kat3, cmd_args=["-cr=on"])

# Now you can get the output objects for each of the files run
# They will be returned in order. So you can get each of the outputs using:
out1, out2, out3 = pk.getResults()

print(out1["P"])
print(out2["P"])
print(out3["P"])

# or you could also get the outputs in a list:
outs = pk.getResults()
# so you'd then access the ouputs like: outs[0], outs[1], etc.
# This option is useful if you want to iterate over each output
# result in a for loop for example:

for out in outs:
    print(out["P"])

# Call `clear` to remove the last output results
# After this you could run more in parallel and get the outputs
# for those.
pk.clear()

# When finally done you should close the connection to the cluster
pk.close()
