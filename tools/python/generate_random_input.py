import random as rnd

number_particles = 1000000

print number_particles
for p in range(number_particles):
    print rnd.random(), rnd.random(), rnd.random(), 2*rnd.random()-1, 0.9+0.2*rnd.random()
