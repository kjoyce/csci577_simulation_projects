from Container import Container
from Force import SledForces
from Integrator import VerletIntegrator
from DistanceMatrix import DistanceMatrix
from numpy import linspace,sqrt,mod,inf,sin,pi
from pylab import rand,randint
from Container import Container
class ParticleInitialize(object):
  def __init__(self,case):
######   Params   ################
    xlim = (10,60)
    ylim = (-1,10)
    n_floor = 100
    n_sled = 13
    start_sled = 20
    pot_energy_lim = (-10,30)
    kin_energy_lim = (40,70)
    tot_energy_lim = (40,70)
    pressure_lim = (40,70)
    r = 2**(1./6.)
    dims = 2
    dt = .01
##################################
    initialization = case
    dy = r*sqrt(3)
    dx = r
    xinit = (2*start_sled+1)*.5*dx + n_sled*dx
    self.distance_matrix = DistanceMatrix()
    self.force = SledForces(dims,self.distance_matrix,xinit)
    self.integrate = VerletIntegrator(dt,self.force)
    self.c = Container(self.integrate)
    c = self.c
    for i in range(n_floor):
      c.addParticle(float(i)*dx,0.,0.,0.,1.)

    for i in range(n_sled):
      x = (2*start_sled+1)*.5*dx + i*dx
      y = dy*(i%2) + r
      c.addParticle(x,y,0.,0.,1.)
      print "{}, {}".format(x,y)

    self.xlim            = xlim
    self.ylim            = ylim
    self.pot_energy_lim  = pot_energy_lim
    self.kin_energy_lim  = kin_energy_lim
    self.tot_energy_lim  = tot_energy_lim
    self.pressure_lim    = pressure_lim
  def __call__(self):
    return self.c,self.distance_matrix,self.force,self.integrate,self.xlim,self.ylim,self.pot_energy_lim,self.kin_energy_lim,self.tot_energy_lim,self.pressure_lim


