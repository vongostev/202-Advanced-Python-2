# -*- coding: utf-8 -*-
"""
Created on Sun Dec 26 14:16:23 2021

@author: grego
"""

class Star(mass: float):
    """
    Class desribes star. 
    Situated in (0,0) and has solar mass, if doesn't mentioned another
    """
    def __init__(self, mass = 10^30, xy = (0, 0)):
        self.mass = mass
        self.xyz = xy
    
    def getMass(self):
        return self.mass
    
    def getPosition(self):
        return self.xy


    
