#!/usr/bin/env python2
'''
Animate for standing wave output, first particle only.

Usage:
  ./animate.py [options] <infile>
  ./animate.py [options]

Options:
  -e END                Read to index END.
  -s STEP               Read with STEP.
  -o OUTPUT             Output to OUTPUT, if not chosen, just display.
  --intensity=I -I I    Set the intensity to I in W/cm^2. [default: 1e18]
  --ET=ET               Run with a time phase of ET*pi for the electric field. [default: 0.0].
  --BT=BT               Run with a time phase of BT*pi for the magnetic field. [default: 0.0].
  --ES=ES               Run with a space phase of ES*pi for the electric field. [default: 1.5].
  --BS=BS               Run with a space phase of BS*pi for the magnetic field. [default: 1.5].
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
Et = float(opts['--ET']);
Es = float(opts['--ES']);
Bt = float(opts['--BT']);
Bs = float(opts['--BS']);
I     = float(opts['--intensity']) if opts['--intensity'] else 1e18;
E_0   = 2*np.sqrt(4*np.pi*I*1e11/c);

if opts['<infile>']:
    with open(opts['<infile>'],'rb') as f:
        t,d=pickle.load(f);

else:
    t,d=pickle.loads(sys.stdin.read());

#first particle
d = d[:,:,0]

x=d[0,:]; y=d[1,:]; z=d[2,:]
del d;

    
y = y[:end:step];
z = z[:end:step];
t = t[:end:step];

if opts['--unitless']:
    units = lm/(2*np.pi);
    y*=units;
    z*=units;
    t*=lm/(2*np.pi*c)#*float(opts['--dt']);
else:
    y*=1e4;
    z*=1e4;

zmin = m.floor(z.min()/(lm/2))*lm/2;
zmax = m.ceil(z.max()/(lm/2))*lm/2;
ymin = y.min();
ymax = y.max();

Y,Z = np.mgrid[ ymin : ymax : 16j,
                zmin : zmax : 16j];
zeros = np.zeros(Z.shape);
Ez = np.cos(2*np.pi*Z/lm+Es*np.pi)*E_0;

#initial plotting
fig = plt.figure();
ax = fig.add_subplot(111);
p, = ax.plot(y[0:1],z[0:1],marker='o',label='r',c='r');
if not opts['--no-B']:
    Yh,Zh = np.mgrid[ ymin : ymax : 64j,
                      zmin : zmax : 64j];
    By = -np.cos(2*np.pi*Zh/lm+Bs*np.pi)*E_0*0.0;
    By=By[:-1,:-1];
    pc = ax.pcolormesh(Yh,Zh,By,vmin=-E_0,vmax=E_0);
q = ax.quiver(Y,Z,Ez,zeros,scale=max(E_0, 5e10));
ax.set_ylim(zmin,zmax);
ax.set_xlim(y.min(),y.max());
ax.set_ylabel('z ($\mu$m)');
ax.set_xlabel('y ($\mu$m)');

def animate(ii):
    i,t=ii;
    Ez = np.cos(2*np.pi*Z/lm+Es*np.pi)*np.cos(2*np.pi/lm*c*t+Et*np.pi)*E_0;
    p.set_data(y[0:i],z[0:i]);
    q.set_UVC(Ez,zeros);
    if not opts['--no-B']:
        By = -np.cos(2*np.pi*Zh/lm+Bs*np.pi)*np.cos(2*np.pi/lm*c*t+Bt*np.pi)*E_0;
        By=By[:-1,:-1];
        pc.set_array(By.ravel());
    return p,q;

ts = list(enumerate(t));
a = anim.FuncAnimation(fig,animate, ts, interval=0.01);
if output:
    a.save(output,fps=5);
else:
    plt.show();
