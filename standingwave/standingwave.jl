#!/usr/bin/env julia
doc="
Run the standingwave simulation.

Usage:
  ./standingwave.jl [options]

Options:
  --help -h            Print this help.
  --periods=T -n T     Run for T periods. [default: 1].
  --intensity=I -I I   Use an intensity of I in W/cm^2. [default: 1e18].
  --starting=Z  -Z X   Use the starting fraction of wavelength. [default: 0.125].
  --wavelength=L -l L  Set the wavelength in cm. [default: 800e-7].
  --ET=ET              Run with a time phase of ET*pi for the electric field. [default: 0.0].
  --BT=BT              Run with a time phase of BT*pi for the magnetic field. [default: 0.0].
  --ES=ES              Run with a space phase of ES*pi for the electric field. [default: 1.5].
  --BS=BS              Run with a space phase of BS*pi for the magnetic field. [default: 1.5].
"
using DocOpt
opts=docopt(doc);
using jelo

const c_si = 2.99792458e8;
const c_cgs = 2.99792458e10;
#const e_cgs = 4.80320425e-10;
const m_e = 9.10938291e-28;
const m_e_si = 9.10938291e-31;
const e_si = 1.602176565e-19;
const e0 = 8.854187817e-12;
const dt_factor = 10000;

I  = float(opts["--intensity"]);
l  = float(opts["--wavelength"]);
Z  = float(opts["--starting"]);
T  = float(opts["--periods"]);
Et = float(opts["--ET"]);
Bt = float(opts["--BT"]);
Es = float(opts["--ES"]);
Bs = float(opts["--BS"]);

# Calculate a_0

# E_0 = sqrt(2I/(c_mks*e0));
# qmr = -(E_0*l)/(pi*.511e6);
# dt  = abs(1.0/qmr/dt_factor);
# zi  = 2pi*Z;
# N   = int(round(2pi*T/dt));
# E(v,u)=[1.0,0.0,0.0]*cos(u+pi*Et)*cos(v[3]+pi*Es);
# B(v,u)=[0.0,1.0,0.0]*cos(u+pi*Bt)*cos(v[3]+pi*Bs);

E_0 = sqrt(2I*1e4/(c_si*e0));
period = l/c_cgs;
w = 2pi/period;
# the natural timescale
tau = m_e_si*c_si/(2e_si*E_0);
# and a_0
a_0 = 1/(2tau*w)
dt = 1/a_0/dt_factor;
zi  = 2pi*Z;
N   = int(round(2pi*T/dt));
# E and B fields

E(v,u) = cos(v[3]+Es)*cos(u+Et)*[0.0,1.0,0.0];
B(v,u) = cos(v[3]+Bs)*cos(u+Bt)*[1.0,0.0,0.0];
# function E(x,t)
#     [1.0,0.0,0.0]*2.0*E_0*cos(2*pi/l*c*t+pi*Et)*cos(2*pi/l*x[3]+pi*Es);
# end
# function B(x,t)
#     [0.0,1.0,0.0]*2.0*E_0*cos(2*pi/l*c*t+pi*Bt)*cos(2*pi/l*x[3]+pi*Bs);
# end

#initializing
j=Jelo(E, B, dt)
add(j, particle([0,0,zi],[0.0,0.0,0.0],2*a_0))

#running
for x = 1:N
    step(j)
    println(output(j))
end
