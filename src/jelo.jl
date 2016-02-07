# Copyright (C) 2014,2015  Gregory K. Ngirmang
#
# This file is part of Jelo.
#
# Jelo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Jelo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Jelo.  If not, see <http://www.gnu.org/licenses/>.
#
include("integrate.jl");



#unfortunately, I don't know any other way around making this
#"unitful" or not.
module lorentz_force

export c_cgs, lorentz, relativistic_lorentz
export lorentz_cgs, relativistic_lorentz_cgs
const c_cgs = 2.99792458e10;
const c_si  = 2.99792458e8;

function lorentz(x::Vector,b::Vector,t::Real,qmr::Real,E,B)
    qmr*(E(x,t)+cross(b,B(x,t)))
end
function relativistic_lorentz(x::Vector,b::Vector,t::Real,qmr::Real,E,B)
    M = eye(3)-[b[1]*b[1] b[1]*b[2] b[1]*b[3]
                b[1]*b[2] b[2]*b[2] b[2]*b[3]
                b[1]*b[3] b[3]*b[2] b[3]*b[3]]
    M*lorentz(x,b,t,qmr,E,B) * sqrt(1-dot(b,b))
end

#cgs
function lorentz_cgs(x::Vector,v::Vector,t::Real,qmr::Real,E,B)
    lorentz(x,v/c_cgs,t,qmr,E,B)
end
function relativistic_lorentz_cgs(x::Vector,v::Vector,t::Real,qmr::Real,E,B)
    relativistic_lorentz(x,v/c_cgs,t,qmr,E,B)
end

#si
function lorentz_si(x::Vector,v::Vector,t::Real,qmr::Real,E,B)
    lorentz(x,v,t,qmr,E,B)
end
function relativistic_lorentz_si(x::Vector,v::Vector,t::Real,qmr::Real,E,B)
    b=v/c_si;
    M = eye(3)-[b[1]*b[1] b[1]*b[2] b[1]*b[3]
                b[1]*b[2] b[2]*b[2] b[2]*b[3]
                b[1]*b[3] b[3]*b[2] b[3]*b[3]]
    M*lorentz_si(x,v,t,qmr,E,B) * sqrt(1-dot(b,b))
end



end

module jelo;
export particle
export Jelo, add, step, output_str;

using integrate;
using lorentz_force;
export c_cgs, c_si;

type particle
    x::Vector
    v::Vector
    qmr::Real
end

type Jelo
    particles::Vector{particle};
    #E and B fields
    E::Function
    B::Function
    force::Function
    imethod::Function
    t::Real
    dt::Real
    function Jelo(E::Function,B::Function,dt::Real,
                  options::Dict{Any,Any})
        userunits = if haskey(options, "units")
            options["units"] else "unitless" end
        usermethod= if haskey(options, "imethod")
            options["imethod"] else "leapfrog"  end
        force = if userunits == "cgs"
            relativistic_lorentz_cgs
        elseif userunits == "si"
            relativistic_lorentz_si    
        elseif userunits == "unitless"
            relativistic_lorentz
        else
            error("units must be either \"cgs\",\"si\", or \"unitless\"")
        end
        method = if usermethod == "rk4"
            integrate.rk4
        elseif usermethod == "leapfrog"
            integrate.leapfrog
        else usermethod end

        new(particle[],E,B,force,method,0,dt);
    end
            
    function Jelo(E::Function,B::Function,dt::Real,units::AbstractString)
        Jelo(E,B,dt,Dict{AbstractString,Any}("units"=>units))
    end
    Jelo(E::Function,B::Function,dt::Real)=Jelo(E,B,dt,"unitless")
end
#adding particles to the simulation
add(j::Jelo,p::particle) = push!(j.particles,p)
#timestep

function step(j::Jelo)
    function step_one(p::particle)
        a(x,v,t::Real) = j.force(x,v,t::Real,
                                 p.qmr,j.E,j.B)
        p.x,p.v = j.imethod(p.x,p.v,a,j.t,j.dt);
        p
    end
    map!(step_one,j.particles);
    j.t+=j.dt
end


vect2s(v::Vector)="$(x[1]),$(x[2]),$(x[3])"

function output_str(j::Jelo)
    out_one(p::particle,last=false)=
        vect2s(p.x)*" "*vect2s(p.v)*(last?"":";");
    if length(j.particles)>1
        out_first(p::particle)=out_one(p);
        en = endof(j.particles);
        strs = map(out_first,j.particles[1:en-1])
        push!(strs, out_one(j.particles[en],true));
        "$(j.t):"*reduce(*,strs);
    else
        "$(j.t):"*out_one(j.particles[1],true);
    end
end

end

