#!/usr/bin/env python2
'''
Check the if the particles do not exceed c.

Usage:
  ./txt2pickle [options] <infile>
  ./txt2pickle [options]

Options:
  --natural-units -n        Use natural units instead of cgs.
'''
from docopt import docopt;
import numpy as np;
import cPickle;
import fileinput;
import itertools as it;
import sys;

opts=docopt(__doc__,help=True);
c = 1.0 if opts['--natural-units'] else 2.99792458e10;

if opts['<infile>']:
    with open(opts['<infile>'],'rb') as f:
        t,d=cPickle.load(f);
else:
    t,d = cPickle.loads(sys.stdin.read());

v = np.sqrt(d[3,:,:]**2 + d[4,:,:]**2 +d[5,:,:]**2).T #makes the first index be particles
del d;
errmsg = "exceeded speed of light at time {} s of particle {} (speed is {}, v-c={})"
#iterate over particles
for i,cv in enumerate(v):
    import matplotlib.pyplot as plt;
    plt.hist(cv,bins=100);
    plt.show();
    bad = cv > c;
    bad_times = t[bad];
    bad_speeds = cv[bad];
    for bv,bt in it.izip(bad_speeds,bad_times):
        print(errmsg.format(bt,i,bv,bv-c));
    pass;
pass;
