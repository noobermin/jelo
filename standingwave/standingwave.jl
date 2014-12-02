#!/usr/bin/env julia
doc="
Run the standingwave simulation.

Usage:
  ./standingwave.jl [options]

Options:
  --help -h            Print this help.
  --periods=T -n T     Run for T periods. [default: 2].
  --intensity=I -I I   Use an intensity of I in W/cm^2. [default: 1e18].
  --starting=Z  -Z X   Use the starting fraction of wavelength. [default: 0.125].
  --ET=ET              Run with a time phase of ET*pi for the electric field.
  --BT=BT              Run with a time phase of BT*pi for the magnetic field.
  --ES=ES              Run with a space phase of ES*pi for the electric field.
  --BS=BS              Run with a space phase of BS*pi for the magnetic field.
  --unitless           Run the unitless simulation.
"
using DocOpt
opts=docopt(doc);
using jelo

const l = 800e-9;
const c_mks = 2.99792458e8;
const c_cgs = 2.99792458e10;
const e_cgs = 4.80320425e-10;
const m_cgs = 9.10938291e-27;
const e0 = 8.854187817e-12;
const dt_factor = 100;

I=float(opts["--intensity"]);
Z=float(opts["--starting"]);
T=float(opts["--periods"]);
if opts["--unitless"]
    E_0 = sqrt(2I/(c_mks*e0));
    qmr = -(E_0*l)/(pi*.511e6);
    dt  = abs(1.0/qmr/dt_factor);
    zi  = 2pi*Z;
    N   = 2pi*T;
    function E(v,u)
        [1.0,0.0,0.0]*cos(u+pi*Et)*cos(v[3]+pi*Es);
    end
    function B(v,u)
        [0.0,1.0,0.0]*cos(u+pi*Bt)*cos(v[3]+pi*Bs);
    end
else
    E_0 = sqrt(4pi*I*1e-7/c_cgs);
    qmr = -e_cgs/m_cgs;
    dt  = abs(c_cgs/(2qmr*E_0)/dt_factor);
    zi  = l*Z*100;
    N   = l*1e2/c_cgs*T;
    function E(x,t)
        [1.0,0.0,0.0]*2.0*E_0*cos(2*pi/l*c*t+pi*Et)*cos(2*pi/l*x[3]+pi*Es);
    end
    function B(x,t)
        [0.0,1.0,0.0]*2.0*E_0*cos(2*pi/l*c*t+pi*Bt)*cos(2*pi/l*x[3]+pi*Bs);
    end
end

#initializing
j=Jelo(E, B, dt)
add(j, particle([0,0,zi],[0.0,0.0,0.0],qmr))

#running
for x = 1:N
    step(j)
    println(output(j))
end
