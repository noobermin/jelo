#!/usr/bin/env julia
doc="
Run the planewave simulation.

Usage:
  ./planewave.jl [options]

Options:
  --help -h            Print this help.
  --periods=T -n T     Run for T periods. [default: 1].
  --intensity=I -I I   Use an intensity of I in W/cm^2. [default: 1e18].
  --starting=X -Z X    Use the starting fraction of wavelength. [default: 0.125].
  --wavelength=L -l L  Set the wavelength in cm. [default: 800e-7].
  --Ephi=ET            Run with a time phase of ET*pi for the electric field. [default: 1.5].
  --Bphi=BT            Run with a time phase of BT*pi for the magnetic field. [default: 1.5].
  --b0=B0              Set the initial velocity as a comma separated list.
"
using DocOpt
opts=docopt(doc);
using jelo

const c_si = 2.99792458e8;
const c_cgs = 2.99792458e10;
const m_e = 9.10938291e-28;
const m_e_si = 9.10938291e-31;
const e_si = 1.602176565e-19;
const e0 = 8.854187817e-12;
const dt_factor = 10000;

I    = float(opts["--intensity"]);
l    = float(opts["--wavelength"]);
Z    = float(opts["--starting"]);
T    = float(opts["--periods"]);
Ephi = float(opts["--Ephi"]);
Bphi = float(opts["--Bphi"]);


# calculate a_0

E_0 = sqrt(2I*1e4/(c_si*e0));
period = l/c_cgs;
w = 2pi/period;
# the natural timescale
tau = m_e_si*c_si/(e_si*E_0);
# and a_0
a_0 = 1/(tau*w)
dt = 1/a_0/dt_factor;
zi  = 2pi*Z;
N   = int(round(2pi*T/dt));

# parsing b0
if typeof(opts["--b0"]) == Nothing
    b0 = [-a_0^2/(a_0^2+4), 0.0, 0.0];
else
    b0 = map(float,split(opts["--b0"],","));
end
# E and B fields
E(v,u) = cos(u-v[3]+Ephi)*[0.0,1.0,0.0];
B(v,u) = cos(u-v[3]+Bphi)*[0.0,0.0,1.0];

#initializing
j=Jelo(E, B, dt)

add(j, particle([0,0,zi], b0, -a_0))

#running
for x = 1:N
    step(j)
    println(output(j))
end
