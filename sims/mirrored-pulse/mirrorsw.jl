#!/usr/bin/env julia
doc="
Run the pulse simulation.

Usage:
  ./planewave.jl [options]

Options:
  --help -h              Print this help.
  --pulses=N -n N        Run for wt = 2pi*α. [default: 4].
  --intensity=I -I I     Use an intensity of I in W/cm^2. [default: 1e18].
  --start-time=T -T T    Start at the time T*alpha*2pi [default: 0.0].
  --wavelength=L -l L    Set the wavelength in cm. [default: 800e-7].
  --b0=B0 -b BO          Set the initial velocity as a comma separated list.
  --dt=DT -t DT          Manually enter a timestep.
  --pulse-period=P -p P  Provide the pulse period in seconds. [default: 60e-15]
  --X=X -x X             Set the starting place as x=-X*2pi*α [default: 0.5]
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
const dt_factor = 2500;

I    = float(opts["--intensity"]);
l    = float(opts["--wavelength"]);
period = l/c_cgs;
n    = float(opts["--pulses"]);
P    = float(opts["--pulse-period"]);
α = P/period;
t0   = float(opts["--start-time"])*2pi*α;
xi   =-float(opts["--X"])*2pi*α;

# calculate a0...
E_0 = sqrt(2I*1e4/(c_si*e0));
w = 2pi/period;
# the natural timescale
tau = m_e_si*c_si/(e_si*E_0);
# and a_0
a0 = -1/(tau*w);
if typeof(opts["--dt"]) == Nothing
    dt = 1/(4*a0^2*dt_factor);
else
    dt = float(opts["--dt"]);
end

println("# a0 = $(a0)");
println("# dt = $(dt)");
println("#  n = $(n)");
println("# t0 = $(t0)");
# parsing b0
if typeof(opts["--b0"]) == Nothing
    b0 = [0.0, 0.0, 0.0];
else
    b0 = map(float,split(opts["--b0"],","));
end
println("#b0 = $(b0[1]),$(b0[2]),$(b0[3])");
println("#α = $(α)");
N   = int(round(2pi*α*n/dt));
println("# N = $(N)");
# parsing dt
unitstep(x) = x > 0.0 ? 1.0 : 0.0;
ηi(x,t) = t-x[1]-α*2pi;
ηr(x,t) = t+x[1]-α*2pi;

Pi(x,t) = sin(ηi(x,t)/2α)*unitstep(ηi(x,t))*unitstep(α*2pi-ηi(x,t))*unitstep(-x[1])
Pr(x,t) = sin(ηr(x,t)/2α)*unitstep(ηr(x,t))*unitstep(α*2pi-ηr(x,t))*unitstep(-x[1])

fi(x,t) = Pi(x,t)*sin(ηi(x,t))*a0
fr(x,t) = -Pr(x,t)*sin(ηr(x,t))*a0
# E and B fields
E(x,t) = (fi(x,t)+fr(x,t))*[0.0,1.0,0.0];
B(x,t) = (fi(x,t)-fr(x,t))*[0.0,0.0,1.0];

#initializing
j=Jelo(E, B, dt)
add(j, particle([xi,0,0], b0, 1.0))

#running
for x = 1:N
    step(j)
    println(output(j))
end
