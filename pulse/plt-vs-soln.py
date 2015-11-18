#!/usr/bin/env python2
# vim: set fileencoding=utf-8 :
'''
Plot momentum vs analytic

Usage:
  ./animate.py [options] <infile>
  ./animate.py [options]

Options:
  -e END                 Read to index END.
  -s STEP                Read with STEP.
  -o OUTPUT              Output to OUTPUT, if not chosen, just display.
  --title=TITLE -t TITLE Set the title.
'''
import math as m;
import matplotlib.pyplot as plt;
import matplotlib.animation as anim;
import numpy as np;
import pickle;
import sys;
from docopt import docopt;
import pulsesoln as soln;
import re;

#arguments
opts=docopt(__doc__,help=True);
end  = int(opts['-e']) if opts['-e'] else None;
step = int(opts['-s']) if opts['-s'] else None;
output = opts['-o'];

if opts['<infile>']:
    with open(opts['<infile>'],'rb') as f:
        t,d,c=pickle.load(f);
else:
    t,d,c=pickle.loads(sys.stdin.read());

#first particle
d = d[:,:,0];
q = lambda i:d[i,:][:end:step];

x  = q(0);
y  = q(1);
bx = q(3);
by = q(4);

gamma = 1/np.sqrt(1-d[3,:]**2-d[4,:]**2-d[5,:]**2);
px=bx*gamma;
py=by*gamma;

px = px[:end:step];
py = py[:end:step];

fig = plt.figure();

#getting sim info from the comment lines
def get_comment(quantity,list=False):
    match_r = "^# *"+quantity;
    get_r   = match_r+" *= *(.*)$";
    l = [s for s in c if re.match(match_r,s)][0];
    return eval(
        re.sub(get_r,
               r'[\1]' if list else r'\1',
               l));

b0 = get_comment("b0",True);
b0 = np.array(b0);
b0mod = np.sqrt(b0[0]**2+b0[1]**2+b0[2]**2);
gm0 = 1.0/np.sqrt(1-np.dot(b0,b0));
al = get_comment("α");
a0 = get_comment("a0");
t0 = get_comment("t0");
n  = get_comment("n");
#now, for the solution.
if n < 1: n =1;
h = np.linspace(t0, 2*np.pi*al*n, 500);
rho = b0*gm0/abs(a0)
twt = soln.wt(h,a0,rho) - t0;
tpx = soln.px(h,a0,rho);
tpy = soln.py(h,a0,rho);
tkx = soln.kx(h,a0,rho);
tky = soln.ky(h,a0,rho);
plt.subplot(411);
if opts["--title"]:
    plt.title(opts['--title']);
plt.plot(t,px,marker='o',label='simulation',c='b');
plt.plot(twt,tpx,marker='x',label='analytic',c='r');
plt.ylabel(r"$\gamma\beta_x$",size=22);
plt.legend();
plt.subplot(412);
plt.plot(  t,  py,marker='o',label='r',c='b')
plt.plot(twt, tpy,marker='x',label='r',c='r');
plt.ylabel(r"$\gamma\beta_y$",size=22);

plt.subplot(413);
plt.plot(t,x,marker='o',label='r',c='b');
plt.plot(twt,tkx, marker='x',label='r',c='r');
plt.ylabel(r"$kx$",size=22);

plt.subplot(414);
plt.plot(t,y,marker='o',label='r',c='b');
plt.plot(twt,tky, marker='x',label='r',c='r')
plt.ylabel(r"$ky$",size=22);
plt.xlabel(r"$\omega t$",size=22);
if output:
    plt.savefig(output);
else:
    plt.show();
