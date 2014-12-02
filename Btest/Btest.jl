using jelo

Bc=1.0

E(x,t)=[0,0,0]
B(x,t)=[Bc,0,0]

#initializing
j=Jelo(E,B,0.0001)
add(j, particle([0,0,0],[0,0,0.5c],1.0))

#running
for x = 1:100000
    step(j)
    println(output(j))
end
