#!/usr/bin/env python2
'''
Scatter plot in 3d.

Usage:
  ./plot3d.py [options] <infile>
  ./plot3d.py [options]

Options:
  -e END         Read to index END.
  -s STEP        Read with STEP. [default: 100]
  -o OUTPUT      Output to OUTPUT, if not chosen, just display.
'''
import math as m;
import matplotlib.pyplot as plt;
from mpl_toolkits.mplot3d import Axes3D
import numpy as np;
import cPickle as pickle;
import sys;
from docopt import docopt;

#arguments
opts=docopt(__doc__,help=True);
end  = int(opts['-e']) if opts['-e'] else None;
step = int(opts['-s']);
output = opts['-o'];

if opts['<infile>']:
    with open(opts['<infile>'],'rb') as f:
        t,d,_=pickle.load(f);
else:
    t,d,_=pickle.loads(sys.stdin.read());

#first particle
d = d[:,:,0]

x=d[0,:]; y=d[1,:]; z=d[2,:]
del d;

x = x[:end:step];
y = y[:end:step];
z = z[:end:step];
t = t[:end:step];

fig = plt.figure(1);
ax = fig.add_subplot(111,projection='3d');
ax.scatter(x,y,z,marker='o',label='r',c='b');
plt.show();

if output:
    a.save(output,fps=15);
else:
    plt.show();
