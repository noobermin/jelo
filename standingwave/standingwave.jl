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
  --ET=ET              Run with a time phase of ET*pi for the electric field. [default: 0.0].
  --BT=BT              Run with a time phase of BT*pi for the magnetic field. [default: 0.0].
  --ES=ES              Run with a space phase of ES*pi for the electric field. [default: 1.5].
  --BS=BS              Run with a space phase of BS*pi for the magnetic field. [default: 1.5].
"
using DocOpt
opts=docopt(doc);
using jelo

const l = 800e-9;
const c_mks = 2.99792458e8;
const e_cgs = 4.80320425e-10;
const m_cgs = 9.10938291e-28;
const e0 = 8.854187817e-12;
const dt_factor = 100;

I  = float(opts["--intensity"]);
Z  = float(opts["--starting"]);
T  = float(opts["--periods"]);
Et = float(opts["--ET"]);
Bt = float(opts["--BT"]);
Es = float(opts["--ES"]);
Bs = float(opts["--BS"]);


# E_0 = sqrt(2I/(c_mks*e0));
# qmr = -(E_0*l)/(pi*.511e6);
# dt  = abs(1.0/qmr/dt_factor);
# zi  = 2pi*Z;
# N   = int(round(2pi*T/dt));
# E(v,u)=[1.0,0.0,0.0]*cos(u+pi*Et)*cos(v[3]+pi*Es);
# B(v,u)=[0.0,1.0,0.0]*cos(u+pi*Bt)*cos(v[3]+pi*Bs);

E_0 = sqrt(4pi*I*1e7/c);
qmr = -e/m_cgs;
period = l*1e2/c;
dt  = min(abs(c/(2qmr*E_0)/dt_factor),period/10000);
zi  = l*Z*100;
N   = int(round(period*T/dt));

function E(x,t)
    [1.0,0.0,0.0]*2.0*E_0*cos(2*pi/l*c*t+pi*Et)*cos(2*pi/l*x[3]+pi*Es);
end
function B(x,t)
    [0.0,1.0,0.0]*2.0*E_0*cos(2*pi/l*c*t+pi*Bt)*cos(2*pi/l*x[3]+pi*Bs);
end


#initializing
j=Jelo(E, B, dt)
add(j, particle([0,0,zi],[0.0,0.0,0.0],qmr))

#running
for x = 1:N
    step(j)
    println(output(j))
end
