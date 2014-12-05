module integrate
export leapfrog


function leapfrog(x,v,a::Function,t::Real,dt::Real)
    a1 = a(x,v,t)
    x += (v + a1*dt/2.0)*dt;
    a2 = a(x,v,t+dt)
    v += (a1+a2)/2.0*dt
    x,v
end

function rk4(x,v,a::Function,t::Real,dt::Real)
    sample(x,v,t)= v, a(x,v,t);
    xk1,vk1 = sample(x,            v,            t);
    xk2,vk2 = sample(x+xk1*dt/2.0, v+vk1*dt/2.0, t+dt/2.0);
    xk3,vk3 = sample(x+xk2*dt/2.0, v+vk2*dt/2.0, t+dt/2.0);
    xk4,vk4 = sample(x+xk3*dt,     v+vk3*dt,     t+dt);
    x+(xk1/2.0+xk2+xk3+xk4/2.0)*dt/3.0,v+(vk1/2.0+vk2+vk3+vk4/2.0)*dt/3.0
end

end



module jelo

export lorentz,relativistic_lorentz
export particle
export Jelo, add, step, output;
export e,c;

using integrate
#constants
global c,e;
const c = 2.99792458e10;
const e = 4.80320425e-10;
# const c = 1.0;
# const e = 0.302822;

function lorentz(x::Vector,v::Vector,t::Real,qmr::Real,E,B)
    qmr*(E(x,t)+cross(v/c,B(x,t)))
end

function relativistic_lorentz(x::Vector,v::Vector,t::Real,qmr::Real,E,B)
    b = v/c;
    M = eye(3)-[b[1]*b[1] b[1]*b[2] b[1]*b[3]
                b[1]*b[2] b[2]*b[2] b[2]*b[3]
                b[1]*b[3] b[3]*b[2] b[3]*b[3]]
    M*lorentz(x,v,t,qmr,E,B) * sqrt(1-dot(b,b))
end

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
    t::Real
    dt::Real
    Jelo(E::Function,B::Function,dt::Real)=new(particle[],E,B,0,dt);
end
#adding particles to the simulation
add(j::Jelo,p::particle) = push!(j.particles,p)
#timestep

function step(j::Jelo)
    function step_one(p::particle)
        a(x,v,t::Real) = relativistic_lorentz(x,v,t::Real,
                                              p.qmr,j.E,j.B)
        p.x,p.v = integrate.leapfrog(p.x,p.v,a,j.t,j.dt);
        p
    end
    map!(step_one,j.particles);
    j.t+=j.dt
end

function output(j::Jelo)
    out_one(p::particle,last=false)=
        "$(p.x[1]),$(p.x[2]),$(p.x[3]) $(p.v[1]),$(p.v[2]),$(p.v[3])"*(last?"":";");
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
