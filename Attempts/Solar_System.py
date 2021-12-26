# -*- coding: utf-8 -*-
"""
Created on Sun Dec 26 14:16:23 2021

@author: grego
"""
import numpy as np


class Star():
    """
    Class desribes star. 
    Situated in (0,0) and has solar mass, if doesn't mentioned another
    """
    def __init__(self, mass: float = 1e30, xy = (0,0)):
        self.mass = mass
        self.xy = np.array(xy)
    
    def getMass(self):
        return self.mass
    
    def getPosition(self):
        return self.xy

class CosmicBody():
    """
    Class desribes any CosmicBody. 
    Situated in (1000,1000), has 100 kg mass and zero velosity, if doesn't mentioned another
    """
    def __init__(self, mass: float = 100, xy = (1000, 1000), vec_v = (0,0)):
        self.mass = mass
        self.xy = np.array(xy)
        self.vec_v = np.array(vec_v)
    
    def getMass(self):
        return self.mass
    
    def getPosition(self):
        return self.xy
    
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
        Tesla_Roadster = CosmicBody(1814, (5400, 6400), (100, 200))
        assert Tesla_Roadster.getMass() == 1814
        # print (Tesla_Roadster.getPosition())
        # print (Tesla_Roadster.getVelosity())
        print("Test Tesla Motorspots is Ok")
        
    test_Star()
    test_Body()
