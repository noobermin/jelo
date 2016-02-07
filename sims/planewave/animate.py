#!/usr/bin/env python2
'''
Animate for plane wave output, first particle only.

Usage:
  ./animate.py [options] <infile>
  ./animate.py [options]

Options:
  -e END                Read to index END.
  -s STEP               Read with STEP.
  -o OUTPUT             Output to OUTPUT, if not chosen, just display.
  --intensity=I -I I    Set the intensity to I in W/cm^2. [default: 1e18]
  --Ephi=Ephi           Run with a phase of Ephi*pi for the electric field. [default: 1.5].
  --Bphi=Bphi           Run with a phase of Bphi*pi for the magnetic field. [default: 1.5].
  --no-B                Don't plot the B-field.
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
Ephi = float(opts['--Ephi']);
Bphi = float(opts['--Ephi']);

I     = float(opts['--intensity'])
E_0   = np.log10(2*np.sqrt(4*np.pi*I*1e11/c));

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

Y,X = np.mgrid[ ymin : ymax : 16j,
                xmin : xmax : 16j];
zeros = np.zeros(X.shape);
Ey = np.cos(2*np.pi*X/lm+Ephi*np.pi)*E_0;

#initial plotting
fig = plt.figure();
ax = fig.add_subplot(111);
p, = ax.plot(y[0:1],x[0:1],marker='o',label='r',c='r');
if not opts['--no-B']:
    Yh,Xh = np.mgrid[ ymin : ymax : 64j,
                      xmin : xmax : 64j];
    Bz = np.cos(2*np.pi/lm*(c*0-Xh) + Bphi*np.pi)*E_0
    Bz=Bz[::-1,::1];
    #Bz=Bz[:-1,:-1];
    pc = ax.pcolormesh(Yh,Xh,Bz,vmin=-E_0,vmax=E_0);

q = ax.quiver(Y,X,Ey,zeros,scale=E_0*20);
ax.set_ylim(xmin,xmax);
ax.set_xlim(ymin,ymax);
ax.set_ylabel('x ($\mu$m)');
ax.set_xlabel('y ($\mu$m)');

def animate(ii):
    i,t=ii;
    Ey = np.cos(2*np.pi/lm*(c*t-X) + Ephi*np.pi)*E_0;
    p.set_data(y[0:i],x[0:i]);
    q.set_UVC(Ey,zeros);
    if not opts['--no-B']:
        Bz = np.cos(2*np.pi/lm*(c*t-Xh) + Bphi*np.pi)*E_0
        Bz=Bz[:1:-1,:-1:1];
        pc.set_array(Bz.ravel());
    return p,q;

ts = list(enumerate(t));
a = anim.FuncAnimation(fig,animate, ts, interval=0.01);
if output:
    a.save(output,fps=15);
else:
    plt.show();
