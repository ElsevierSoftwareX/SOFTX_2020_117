# Submitted by Sebastian.
# Checks that trace data output returns the same beamsize when a different default wavelength is used

from pykat import finesse

kat_code = """
l laser 1 0 0 n1
s s1 10 1 n1 n2

m2 mITM 0.9 0.0 0 n2 n3
attr mITM Rc -2.0
s scav 1 n3 n4
m2 mETM 0.9 0.0 0 n4 n5
attr mETM Rc 2.0

cav myCav mITM n3 mETM n4

s s3 1 n5 n6

pd pd_trans n6

noxaxis
yaxis abs:deg                                 # move detector frequency with xaxis
"""

kat = finesse.kat()
kat.verbose = False
kat.lambda0 = 1550e-9
kat.parse(kat_code)
kat.maxtem = 0
kat.trace = 2

out, T = kat.run(getTraceData=True)

print (out.stdout)

bp = T[0]['n4'][0]

# this is not overwritten to 1550nm as above
print('beam_param.wavelength: {:.0f}nm'.format(bp.wavelength*1e9))
# therefore this is wrong
print('w0 from beam_param.w0: {:.2f}um'.format(bp.w0*1e6))
# and this does not really work as the wavelength cancels out
# for z=0 and therefore the waist does not change
print('w0 from beam_param.beamsize(): {:.2f}um'.format(bp.beamsize(z=-bp.z)*1e6))

assert(bp.w0 == bp.beamsize(z=-bp.z))