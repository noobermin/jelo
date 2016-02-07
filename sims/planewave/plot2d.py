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
  --X -x                Plot x component.
  --Y -y                Plot y component.
  --Z -z                Plot z component.
  --velocity -v         Plot velocity instead of position.
  --momentum -m         Plot momentum instead of position.          
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
        _,d=pickle.load(f);
else:
    _,d=pickle.loads(sys.stdin.read());

#first particle
d = d[:,:,0];

select = []
if opts['--X']: select.append(0);
if opts['--Y']: select.append(1);
if opts['--Z']: select.append(2);

if len(select) == 0:
    select = [0,1];
elif len(select) != 2:
    print("select two dimensions or none (default: xy)");
    quit(-1);

#only handle unitless sims right now.    
if opts['--velocity']:
    labels = [r'$\beta_'+l for l in ['x$','y$','z$']];
    x = d[3+select[0],:];
    y = d[3+select[1],:];
    xmin,xmax = x.min(),x.max();
    ymin,ymax = y.min(),y.max();
elif opts['--momentum']:
    labels = [r'$\gamma\beta_'+l for l in ['x$','y$','z$']];
    x = d[3+select[0],:];
    y = d[3+select[1],:];
    gamma = 1/np.sqrt(1-d[3,:]**2-d[4,:]**2-d[5,:]**2);
    print(min(gamma));
    x*=gamma;
    y*=gamma;
    xmin,xmax = x.min(),x.max();
    ymin,ymax = y.min(),y.max();
else:
    labels = ['$x$','$y$','$z$'];
    x = d[select[0],:];
    y = d[select[1],:];
    units = lm/(2*np.pi);
    x*=units;
    y*=units;
    xmin = x.min(); #m.floor((x.min())/(lm/4))*lm/4;
    xmax = x.max(); #m.ceil(x.max()/(lm/4))*lm/4;
    ymin = y.min();
    ymax = y.max();
del d;
xlabel = labels[select[0]];
ylabel = labels[select[1]];
x = x[:end:step];
y = y[:end:step];


fig = plt.figure();
ax = fig.add_subplot(111);
p, = ax.plot(y,x,marker='o',label='r',c='r');
plt.xlabel(xlabel);
plt.ylabel(ylabel);
plt.xlim(ymin,ymax);
plt.ylim(xmin,xmax);

if output:
    plt.savefig(output);
else:
    plt.show();
