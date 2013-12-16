#!/usr/bin/python

import random as rnd
from math import *
from numpy import *

number_particles = 1000000

add_grid = False

SPHERE_TYPE_HOLLOW             = 1
SPHERE_TYPE_HOMOGENEOUS        = 2
SPHERE_TYPE_DEGRADING_DENSITY  = 3

# setup type
type = SPHERE_TYPE_HOMOGENEOUS
# particle charge and mass
q = 1.
m = 1.
# sphere radius
r_sphere   = 1.
neutral    = False
grid_count = 400
grid_max   = 10.


def r_uniform():
    myr = r_sphere + 1.
    while myr > r_sphere:
        x = rnd.uniform(-r_sphere, r_sphere)
        y = rnd.uniform(-r_sphere, r_sphere)
        z = rnd.uniform(-r_sphere, r_sphere)
        myr = sqrt(x*x+y*y+z*z)

    return x, y, z

def hollow_sphere(radius):
    # see http://mathworld.wolfram.com/SpherePointPicking.html
    theta = rnd.uniform(0, 2*pi)
    u     = rnd.uniform(-1, 1)
    return radius*sqrt(1-u*u)*cos(theta), radius*sqrt(1-u*u)*sin(theta), radius*u

def r_hollow():
    return hollow_sphere(r_sphere)

def rho(S,mu):
    # profile according to DOI: 10.1103/PhysRevLett.91.143401, eq. (9)
    return (S*S) * (3./4.*pi)/(1+S**(3*mu))**(1/mu+1)

def r_degrading():
    # use profile according to DOI: 10.1103/PhysRevLett.91.143401, eq. (9)
    mu   =  1.
    maxS = 10.
    # first choose S
    myS    = 0.
    maxrho = max(rho(arange(0.,maxS,0.05), mu))
    myrho  = maxrho + 1.
    while myrho > rho(myS, mu):
        myS   = rnd.uniform(0., maxS)
        myrho = rnd.uniform(0., maxrho)

    r = r_sphere * myS
    # pick other coordinates as for the hollow sphere
    return hollow_sphere(r)
    
    
get_r = { SPHERE_TYPE_HOLLOW             : r_hollow,
          SPHERE_TYPE_HOMOGENEOUS        : r_uniform,
          SPHERE_TYPE_DEGRADING_DENSITY  : r_degrading }


if add_grid:
	print number_particles + grid_count
else:
	print number_particles

if add_grid:
    for x in linspace(-grid_max, grid_max, grid_count):
        if (x<0):
            print x, 0.0, 0.0, 0.0, -2.0
        else:
            print x, 0.0, 0.0, 0.0, -1.0

for p in range(number_particles):
    x,y,z = get_r[type]()
    
    if neutral:
        print x, y, z, (2*(p%2)-1)*q, m
    else:
        print x, y, z, q, m

