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
  --unitless -u         Do not scale.
'''
import math as m;
import matplotlib.pyplot as plt;
import numpy as np;
import pickle;
import sys;
import re;
from docopt import docopt;

#arguments
opts=docopt(__doc__,help=True);
end  = int(opts['-e']) if opts['-e'] else None;
step = int(opts['-s']) if opts['-s'] else None;
output = opts['-o'];


lm = 8e-7;
c = 2.99792458e8;
e0 = 8.854e-12
me_si = 9.10938291e-31;
e  = 1.602e-19


if opts['<infile>']:
    with open(opts['<infile>'],'rb') as f:
        _,d,cl=pickle.load(f);
else:
    _,d,cl=pickle.loads(sys.stdin.read());

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
    xmin = x.min(); 
    xmax = x.max(); 
    ymin = y.min();
    ymax = y.max();
del d;
xlabel = labels[select[0]];
ylabel = labels[select[1]];
x = x[:end:step];
y = y[:end:step];

#generate analytic solution plot
#obtain b0 from the comment lines (HACKING)
b0x_s = [s for s in cl if re.match(r'^# *b0',s)][0];
b0x,b0y,b0z = eval(re.sub(r'^# *b0 *= *(.*)$',r'[\1]',b0x_s));
b0 = np.sqrt(b0x**2+b0y**2+b0z**2);

#get quiver velocity
a0_s = [s for s in cl if re.match(r'^# *a0',s)][0];
a0 = eval(re.sub(r'^# *a0 *= *(.*)$',r'\1',a0_s));
#get T
a0_s = [s for s in cl if re.match(r'^# *T',s)][0];
T = eval(re.sub(r'^# *T *= *(.*)$',r'\1',a0_s));

#redshift factor
rs = np.sqrt(1-b0**2)/(1-b0x);

h = np.linspace(0,T*2*np.pi, len(x));
kx = h*b0x/(1-b0x) + a0**2/4*rs**2*(h - np.sin(2*h)/2) + a0*(np.cos(h)-1)*rs*b0y/(1-b0x)
ky = h*b0y/(1-b0x) + a0*rs*(np.cos(h)-1);
kz = h*b0z/(1-b0x);

ax = kx;
ay = ky;

fig = plt.figure();
plt.plot(y, x, marker='o',label='simulation',c='b');
plt.plot(ky, kx, marker='+',label='analytic',c='r');
plt.legend(loc="upper left");
plt.xlabel(ylabel);
plt.ylabel(xlabel);

if output:
    plt.savefig(output);
else:
    plt.show();
