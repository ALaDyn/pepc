#!/usr/bin/python

from math import *


CRYSTAL_TYPE_MADELUNG         = 1


# setup type
type = CRYSTAL_TYPE_MADELUNG
# lattice cells per edge of cubic simulation box
boxes_per_edge = 32
# cubic simulation box edge length
edge_length = 1.0
# particle charge and mass
qelectron =     1.
qion      =    -1.*qelectron
melectron =     1.
mion      =  1836.*melectron


def crystal_madelung(i, j, k):
    x = (i+0.5) * delta
    y = (j+0.5) * delta
    z = (k+0.5) * delta
    if ((i+j+k)%2) == 0:
        q = qelectron
        m = melectron
    else:
        q = qion
        m = mion

    return x, y, z, q, m


get_crystal = { CRYSTAL_TYPE_MADELUNG    : crystal_madelung }

number_particles_per_edge = 2*boxes_per_edge
number_particles = number_particles_per_edge**3
delta            = edge_length/(1.*number_particles_per_edge)


print number_particles

for i in range(0,number_particles_per_edge):
    for j in range(0,number_particles_per_edge):
        for k in range(0,number_particles_per_edge):
            x, y, z, q, m =  get_crystal[type](i,j,k)
            print x, y, z, q, m
