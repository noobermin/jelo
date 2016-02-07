#!/usr/bin/env python2
# vim: set fileencoding=utf-8 :
'''
Plot momentum vs analytic

Usage:
  ./check.py [options] <infile>
  ./check.py [options]

Options:
  -b BEGIN               Read to index BEGIN.
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
import re;

#arguments
opts=docopt(__doc__,help=True)
begin = int(opts['-b']) if opts['-b'] else None;
end   = int(opts['-e']) if opts['-e'] else None;
step  = int(opts['-s']) if opts['-s'] else None;
output = opts['-o'];

if opts['<infile>']:
    with open(opts['<infile>'],'rb') as f:
        t,d,c=pickle.load(f);
else:
    t,d,c=pickle.loads(sys.stdin.read());
#first particle
d = d[:,:,0];
q = lambda i:d[i,:][begin:end:step];
t = t[begin:end:step];
x  = q(0);
y  = q(1);
bx = q(3);
by = q(4);
bz = q(5);
gamma = 1/np.sqrt(1-bx**2-by**2-bz**2);
px=bx*gamma;
py=by*gamma;

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

plt.subplot(411);
if opts["--title"]:
    plt.title(opts['--title']);
plt.plot(t,px,marker='o',label='$p_x$',c='b');
plt.plot(t,py,marker='o',label='$p_y$',c='r');
plt.legend();
plt.ylabel(r"$\gamma\beta_{(x,y)}$",size=22);
plt.subplot(412);
plt.plot(t,y,marker='o',label='r',c='b');
plt.ylabel(r"$ky$",size=22);
plt.subplot(413);
plt.plot(t,x,marker='o',label='r',c='b');
plt.ylabel(r"$kx$",size=22);
plt.xlabel(r"$\omega t$",size=22);
plt.subplot(414);
plt.plot(y,x,marker='o',label='r',c='b');
plt.ylabel(r"$kx$",size=22);
plt.xlabel(r"$ky$",size=22);
if output:
    plt.savefig(output);
else:
    plt.show();
