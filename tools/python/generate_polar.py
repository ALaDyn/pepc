import random as rnd

number_particles = 1000

print number_particles+2

print -5.0, 0.0, 0.0,  1e4, 1e3
print  5.0, 0.0, 0.0, -1e4, 1e3

for p in range(number_particles/2):
    print -0.5+1*rnd.random(), -0.5+1*rnd.random(), -0.5+1*rnd.random(),  1.0, 1.0
    print -0.5+1*rnd.random(), -0.5+1*rnd.random(), -0.5+1*rnd.random(), -1.0, 1.0
