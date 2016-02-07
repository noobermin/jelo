#!/usr/bin/env julia
doc="
Run the planewave simulation.

Usage:
  ./planewave.jl [options]

Options:
  --help -h            Print this help.
  --periods=T -n T     Run for wt = 2pi*T. [default: 1].
  --intensity=I -I I   Use an intensity of I in W/cm^2. [default: 1e18].
  --starting=X -X X    Use the starting fraction of wavelength. [default: 0.0].
  --wavelength=L -l L  Set the wavelength in cm. [default: 800e-7].
  --Ephi=ET            Run with a time phase of ET*pi for the electric field. [default: 0.0].
  --Bphi=BT            Run with a time phase of BT*pi for the magnetic field. [default: 0.0].
  --b0=B0 -b BO        Set the initial velocity as a comma separated list.
  --dt=DT -t DT        Manually enter a timestep.
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
const dt_factor = 2000;

I    = float(opts["--intensity"]);
l    = float(opts["--wavelength"]);
X    = float(opts["--starting"]);
T    = float(opts["--periods"]);
Ephi = float(opts["--Ephi"]);
Bphi = float(opts["--Bphi"]);


# calculate a0

E_0 = sqrt(2I*1e4/(c_si*e0));
period = l/c_cgs;
w = 2pi/period;
# the natural timescale
tau = m_e_si*c_si/(e_si*E_0);
# and a_0
a0 = 1/(tau*w);
if typeof(opts["--dt"]) == Nothing
    dt = 1/(a0^2*dt_factor);
else
    dt = float(opts["--dt"]);
end

xi  = 2pi*X;

println("# a0 = $(a0)");
println("# dt = $(dt)");
println("# T = $(T)");
# parsing b0
if typeof(opts["--b0"]) == Nothing
    b0 = [-a0^2/(a0^2+4), 0.0, 0.0];
else
    b0 = map(float,split(opts["--b0"],","));
end
println("#b0 = $(b0[1]),$(b0[2]),$(b0[3])");

N   = int(round(2pi*T/dt));
println("# N = $(N)");
# parsing dt

# E and B fields
E(v,u) = cos(u-v[1]+Ephi*pi)*[0.0,1.0,0.0];
B(v,u) = cos(u-v[1]+Bphi*pi)*[0.0,0.0,1.0];

#initializing
j=Jelo(E, B, dt)

add(j, particle([xi,0,0], b0, -a0))

#running
for x = 1:N
    step(j)
    println(output(j))
end
