#!/usr/bin/env python2
'''
Just plot first particle only trajectory.

Usage:
  ./animate.py [options] <infile>
  ./animate.py [options]

Options:
  -e END                Read to index END.
  -s STEP               Read with STEP.
  -o OUTPUT             Output to OUTPUT, if not chosen, just display.
  --unitless -u         Plot for the unitless case.
  --dt DT -t DT         Set dt factor for unitless. [default: 1e4].
'''
import math as m;
import matplotlib.pyplot as plt;
import matplotlib.animation as anim;
import numpy as np;
import pickle;
import sys;
from docopt import docopt;

#arguments
opts=docopt(__doc__,help=True);
end  = int(opts['-e']) if opts['-e'] else None;
step = int(opts['-s']) if opts['-s'] else None;
output = opts['-o'];

#IMPORTANT: we use microns and seconds here.
lm = 8e-1;
c = 2.99792458e14;


if opts['<infile>']:
    with open(opts['<infile>'],'rb') as f:
        t,d=pickle.load(f);

else:
    t,d=pickle.loads(sys.stdin.read());

#first particle
d = d[:,:,0]

x=d[0,:]; y=d[1,:]; z=d[2,:]
del d;
    
x = x[:end:step];
y = y[:end:step];
t = t[:end:step];

if opts['--unitless']:
    units = lm/(2*np.pi);
    x*=units;
    y*=units;
    t*=lm/(2*np.pi*c);
else:
    x*=1e4;
    y*=1e4;

xmin = m.floor((x.min())/(lm/4))*lm/4;
xmax =  m.ceil(x.max()/(lm/4))*lm/4;
ymin = y.min();
ymax = y.max();

fig = plt.figure();
ax = fig.add_subplot(111);
p, = ax.plot(y,x,marker='o',label='r',c='r');

if output:
    plt.savefig(output);
else:
    plt.show();
