#!/usr/bin/env python2
'''
Convert generic output to a pickle of data.

Usage:
  ./txt2p.py <infile> <outfile>
  ./txt2p.py
'''
from docopt import docopt;
import numpy as np;
import re;
import cPickle as pickle;
from string import strip;
import fileinput;

opts=docopt(__doc__,help=True);
inname = opts['<infile>'];

#get data and split
if inname:
    with open(inname,'r') as f:
        lines = map(strip,f.readlines());
else:
    lines = list(fileinput.input());
# skip lines starting with #
comment = re.compile(r"^ *#");
datalines = [l for l in lines if not comment.match(l)];
commentlines = [l for l in lines if comment.match(l)];
#count number of particles
N=len(datalines[0].split(';'));
r = re.compile(r'[:,\,,;]');
datalines[:] = map(lambda i: r.sub(' ',i), datalines);

lines = np.array(map(lambda i: np.fromstring(i,sep=' '), datalines));
time = lines[:,0];
 #shape now is [quantity, time, particle];
lines = lines[:,1:].T.reshape(6,-1,N);
if inname:
    with open(opts['<outfile>'],'wb') as f:
        pickle.dump((time,lines,commentlines),f,2);
else:
    print(pickle.dumps((time,lines,commentlines)));
pass;

