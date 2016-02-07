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
module integrate
export leapfrog, rk4

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
