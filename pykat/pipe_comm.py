import time
import sys
import struct
import pykat.exceptions as pkex
import os

from struct import calcsize as cs

def open_pipes(pipe_name, start_kat, duration):
    fifo_r = None
    fifo_w = None
    
    while fifo_r is None or fifo_w is None:
        try:
            if time.time() < start_kat + duration:
                if fifo_w is None and os.path.exists(pipe_name+".pf"):
                    fifo_w = open(pipe_name+".pf", "wb", 0)
                
                if fifo_r is None and os.path.exists(pipe_name+".fp"):
                    fifo_r = open(pipe_name+".fp", "rb", 0)
            else:
                raise pkex.BasePyKatException("Could not connect to pykat pipe in {0} seconds. Ensure you are using Finesse >= v2.1 and Pykat >= v1.0.0. Or set usePipe=False when making kat object.".format(duration))
                
        except FileNotFoundError as ex:
            sys.stdout.flush()
            
    return fifo_r, fifo_w
                    
def send_test(fifo_w, a, b, c):
    if fifo_w is not None:
        data = [a,b,c]
        bdata = struct.pack('ddi', *data)
        
        fifo_w.write(("<test>").encode())
        fifo_w.write(struct.pack('b', len(bdata)))
        fifo_w.write(bdata)
        fifo_w.flush()
        print("Sent",len(bdata), "bytes")
        
def send_test_str(fifo_w, a):
    if fifo_w is not None:
        bdata = struct.pack('20s', a.encode('ascii'))
        
        fifo_w.write(("<test_str>").encode())
        fifo_w.write(struct.pack('b', len(bdata)))
        fifo_w.write(bdata)
        fifo_w.flush()
        print("Sent",len(bdata), "bytes")
        
def send_finished(fifo_w):
    if fifo_w is not None:
        fifo_w.write(("<finished>").encode())
        fifo_w.write(struct.pack('b', 0))
        fifo_w.flush()
        
def send_do_step(fifo_w):
    if fifo_w is not None:
        fifo_w.write(("<do_step>").encode())
        fifo_w.write(struct.pack('b', 0))
        fifo_w.flush()
        
def send_do_axis(fifo_w, idx, _from, to, N):
    if fifo_w is not None:
        fifo_w.write(("<do_axis>").encode())
        fifo_w.write(struct.pack('b', cs('i') + cs('d') + cs('d') + cs('i')))
        fifo_w.write(struct.pack('i', idx))
        fifo_w.write(struct.pack('d', _from))
        fifo_w.write(struct.pack('d', to))
        fifo_w.write(struct.pack('I', N))
        fifo_w.flush()
        
def send_update(fifo_w, idx, value):
    if fifo_w is not None:
        bdata = struct.pack('id', idx, value)
        
        fifo_w.write(("<update>").encode())
        fifo_w.write(struct.pack('b', struct.calcsize('i') + struct.calcsize('d')))
        fifo_w.write(struct.pack('i', idx))
        fifo_w.write(struct.pack('d', value))
        fifo_w.flush()