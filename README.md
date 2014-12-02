jelo
====

A Julia implementation of the EM solver electro, that is <b>J</b>ulia
<b>El</b>ectr<b>o</b> or `jelo`. Mind you, this was written within about
a week of studying the [Julia docs](http://julia.readthedocs.org), which
is a testament to how shit easy it is to learn/use Julia. The code itself was
written over, I don't know, two days or less?

Usage / Configuration
---------------------
To use, it import jelo.

```
using jelo
```

Define electric and magnetic field functions, then create the `jelo` object.
```
j=Jelo(myE,myB, dt)
```
Add particles to the simulation using the `add` method,
```
#               position  velocity    charge-mass ratio
add(j, particle([0,0,0],  [0.5c,0,0],  e/.511e6))
```
Then, loop, using the `step` method to step the solver,
```
for i=1:10000
    step(j)
```
and the `output` method to output the result of the step
```
    println(output(j))
end
```

Technical
---------
As with [Electro](http://github.com/noobermin/electro), the integration method is the "leapfrog" algorithm, although rk4 is included (as it is in `electro`). It uses the relativistic stepper described in this [file](http://github.com/noobermin/electro/blob/master/docs/equationsofmotion.pdf).
