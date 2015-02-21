#!/usr/bin/env python2
'''
Test plot.

Usage:
  testplot.py [options]

Options:
  --es=ES           Set E space phase. [default: 0.0].
  --bs=BS           Set B space phase. [default: 0.0].
  --et=ET           Set E space phase. [default: 1.5].
  --bt=ES           Set E space phase. [default: 1.5].
  --I=I             Set intensity in W/cm^2 [default: 1e16].
  --L=l             Set wavelength in cm [default: 800e-7].
'''
from docopt import docopt;
from born import born;
import numpy as np;
import matplotlib.pyplot as plt;
opts = docopt(__doc__);
es = float(opts['--es'])*np.pi;
et = float(opts['--et'])*np.pi;
bs = float(opts['--bs'])*np.pi;
bt = float(opts['--bt'])*np.pi;

I = float(opts['--I'])*1e4;

e_si   = 1.602176565e-19;
c_si   = 2.99792458e8;
m_e_si = 9.10938291e-31;
e0     = 8.854187817e-12;
w      = 2*np.pi*c_si/(float(opts['--L'])*1e-2);

E_0 = np.sqrt(2*I/(c_si*e0));
a0  = e_si*E_0/(m_e_si*c_si*w);
print(a0);
#fucking matlab
t = np.linspace(0,2*np.pi*10,10000);
x  = np.array([0.0, 0.0, 0.125*2*np.pi]);
b0 = np.array([0.0, 0.0, 0.0]);
beta, v = born(x,b0,
               a0,es,et,bs,bt);


_,T = np.meshgrid(x,t);
db = beta(T);

plt.plot(t,db[:,1]);
plt.show();
