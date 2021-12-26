# -*- coding: utf-8 -*-
"""
Created on Sun Dec 26 14:16:23 2021

@author: grego
"""
import numpy as np

#Constants
G = 6.67 * 1e-20

class Star():
    """
    Class desribes star. 
    Placed in (0,0) and has solar mass, if doesn't mentioned another
    """
    def __init__(self, mass: float = 1e30, xy = (0,0)):
        self.mass = mass
        self.vec_p = np.array(xy)
    
    def getMass(self):
        return self.mass
    
    def getPosition(self):
        return self.vec_p

class CosmicBody():
    """
    Class desribes any CosmicBody. 
    Placed in (1000,1000), has 100 kg mass and zero velosity, if doesn't mentioned another
    Mass [kg]
    P[km, km]
    V[km/s]
    """
    def __init__(self, mass: float = 1e4, vec_p = (0, 1e5), vec_v = (0,0)):
        self.mass = mass
        self.vec_p = np.array(vec_p)
        self.vec_v = np.array(vec_v)
        
    def destroy(self):
        self.mass = 0
        self.vec_p = 0
        self.vec_v = 0
    
    def move(self, dt):
        self.vec_p = self.vec_p + self.vec_v * dt
        
    def grav(self, Cosmic1):
        self.vec_v = self.vec_v - G * Cosmic1.getMass * vec_p / ( np.sum(vec_p * vec_p) )^(3/2)
        
    def getMass(self):
        return self.mass
    
    def getPosition(self):
        return self.vec_p
    
    def getVelosity(self):
        return self.vec_v
    
    
if __name__ == '__main__':    

    def test_Star():
        Sun = Star(35, (3, 9))
        assert Sun.getMass() == 35
        # print (Sun.getPosition())
        # assert Sun.getPosition() == (3,9)
        print("test_Star is Ok")        
    
    def test_Body():
        Tesla_Roadster = CosmicBody(1814, (0, 6400), (7.91, 0))
        assert Tesla_Roadster.getMass() == 1814
        # print (Tesla_Roadster.getPosition())
        # print (Tesla_Roadster.getVelosity())
        print("Test Tesla Motorspots is Ok")
    
    def straight_Motion():
        """
        Discrete motion. During time dt moves with constant velosity. Then momentum changes according to the Gravitational Law.
        Time [sec]
        Mass [kg]
        P[km, km]
        V[km/sec]
        """
        Sun = Star()
        Aster = CosmicBody()
        dt = 60
        n = 4000
        t = 0
        while t < dt*n:
            
        
        
    test_Star()
    test_Body()
