#!/usr/bin/env julia
using jelo
using HDF5
using Base.Test
using DocOpt
doc="
Run the magnetic field test Btest.

Usage:
    ./Btest.jl [options]

Options:
    --help -h     Print this help.
    --B=B -B B    Use this magnetic field strength in Gauss. [default: 1.0]
    --qmr=Q -q Q  Use this charge-to-mass ratio. [default: 1.0]
    --b=B -b      Set the initial velocity as a comma separated list.
    --output=O    Output to this file.
    --print       Print instead of running a test.
"
opts=docopt(doc);
chkopt(s)=opts[s] != Void()
truopt(s)=chkopt(s) && opts[s]
b0 = if !chkopt("--b")
    [0.5, 0.0, 0.5]
else
    map(float, split(opts["--b"],","))
end

if dot(b0,b0) > 1.0
    error("Unphysical speed suggested.")
end

Bc = float(opts["--B"])
qmr= float(opts["--qmr"])

#creating fields
E(x,t)=[0,0,0]
B(x,t)=[0,0,Bc]

γ = 1.0/sqrt(1-dot(b0,b0));
ωg = Bc*qmr/γ;

dt = 1/(200ωg)
#run for 10 periods
totalt = round(Int,(5/(ωg/(2pi)))/dt)
j=Jelo(E, B, dt,
       Dict{Any,Any}("imethod"=>"rk4"))
jelo.add(j, particle([0,0,0], b0, qmr))

# See Btest.pdf, this is a freshman undergraduate
# problem. This is the anayltic solution.

function b_sol(t)
    [ b0[2]*sin(ωg*t) + b0[1]*cos(ωg*t),
      b0[2]*cos(ωg*t) - b0[1]*sin(ωg*t),
      b0[3] ];
end
function x_sol(t)
    [ b0[2]*(1-cos(ωg*t)) + b0[1]*sin(ωg*t),
      b0[2]*sin(ωg*t)     + b0[1]*(cos(ωg*t)-1),
      b0[3]*ωg*t ] / ωg;
end

getv(j::Jelo)=j.particles[1].v
getx(j::Jelo)=j.particles[1].x

out=Array(Float64,totalt,7)
#running
for i = 1:totalt
    jelo.step(j)
    out[i,1]   = j.t;
    out[i,2:4] = abs((b_sol(j.t) - getv(j))./b_sol(j.t))
    out[i,5:7] = abs((x_sol(j.t) - getx(j))./x_sol(j.t))
end
#outputting
if chkopt("--output")
    h5open(opts["--output"], "w") do f
        f["/Btest"]=out
    end
elseif truopt("--print")
    for i = 1:totalt
        println("$(out[i,0]):"*
                jelo.vect2s(getx(j))*" "jelo.vect2s(getv(j)))
    end
else
    #gathering max
    maxes = mapslices(maximum, out, 1)[:,2:end]
    maxes = reshape(maxes,6);
    println("errors:")
    println(maxes)
    const tol=5e-5
    println("tolerance: $(tol)");
    @test minimum(maxes .< tol*ones(6))
end
