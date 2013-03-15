from Container import Container
from Force import SledForces
from Integrator import VerletIntegrator
from DistanceMatrix import DistanceMatrix
from numpy import linspace,sqrt,mod,inf,sin,pi
from pylab import rand,randint
from Container import Container
class ParticleInitialize(object):
  def __init__(self,case):
    dims = 2
    dt = .01
    self.distance_matrix = DistanceMatrix()
    self.force		 = SledForces(dims,self.distance_matrix)
    self.integrate	 = VerletIntegrator(dt,self.force)
    self.c		 = Container(self.integrate)
    c = self.c
    initialization = case
    xlim = (0,18)
    ylim = (-1,5)
    pot_energy_lim = (-10,30)
    kin_energy_lim = (40,70)
    tot_energy_lim = (40,70)
    pressure_lim = (40,70)
    
    dx = 2**(1./3.)
    dy = sin(pi/3.)*dx
    for i in range(100):
      c.addParticle(float(i)*dx,0.,0.,0.,1.)

    for i in range(100,113):
      print "{}, {}".format((i-99)*dx,dy*i%2 + dx)
      c.addParticle((i-99)*dx,dy*(i%2) + dx,0.,0.,1.)

    self.xlim            = xlim
    self.ylim            = ylim
    self.pot_energy_lim  = pot_energy_lim
    self.kin_energy_lim  = kin_energy_lim
    self.tot_energy_lim  = tot_energy_lim
    self.pressure_lim    = pressure_lim
  def __call__(self):
    return self.c,self.distance_matrix,self.force,self.integrate,self.xlim,self.ylim,self.pot_energy_lim,self.kin_energy_lim,self.tot_energy_lim,self.pressure_lim


