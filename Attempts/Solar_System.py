# -*- coding: utf-8 -*-
"""
Created on Sun Dec 26 14:16:23 2021

@author: grego
"""

class Star():
    """
    Class desribes star. 
    Situated in (0,0) and has solar mass, if doesn't mentioned another
    """
    def __init__(self, mass: float = 1e30, xy = (0, 0)):
        self.mass = mass
        self.xy = xy
    
    def getMass(self):
        return self.mass
    
    def getPosition(self):
        return self.xy


if __name__ == '__main__':    

    def test_Star():
        Sun = Star(35, (3, 9))
        assert Sun.getMass() == 35
        assert Sun.getPosition() == (3,9)
        print("test_Star is Ok")        
        
    test_Star()
