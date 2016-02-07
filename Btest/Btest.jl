#!/usr/bin/env julia
using jelo

Bc=1.0

E(x,t)=[0,0,0]
B(x,t)=[0,0,Bc]

#initializing
b0 = [0.5, 0, 0];
qmr= 1.0;

j=Jelo(E,B,0.01)
jelo.add(j, particle([0,0,0],[0.5,0,0],qmr))

#see Btest.pdf, this is a freshman undergraduate problem.
#anayltic solution.
γ = 1.0/sqrt(1-dot(b0,b0));
ωg = Bc*qmr/γ;
b_sol(t) = b0+[sin(ωg*t + p), cos(ωg*t + p), 1];
x_sol(t) = [[sin(ωg*t + p), sin(ωg*t + p), 1], b0[3]*t]
          

#running
for x = 1:100000
    jelo.step(j)
    println(output_str(j))
end
