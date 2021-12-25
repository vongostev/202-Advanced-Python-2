# -*- coding: utf-8 -*-

from random import random
import numpy


SimulationSettings = {}
SimulationSettings.update({     "SimulationRadius"          : 30.0E+12    })
SimulationSettings.update({     "ObjectMaxInitialSpeed"     : 15.0E+10    })
SimulationSettings.update({     "ObjectMaxInitialAngVel"    : 62.0000     })
SimulationSettings.update({     "ObjectMaxInitialMass"      : 2.00E+27    })
SimulationSettings.update({     "ObjectMaxInitialRadius"    : 8.00E+7     })


class ObjectBox:
    
    echo = False
    
    def __init__(self):
        self.pos = numpy.array([0.0,0.0,0.0])
        self.spd = numpy.array([0.0,0.0,0.0])
    
        self.axs = numpy.array([0.0,1.0,0.0])
    
        self.mas = 1
        self.rad = 1
        
        
        
        VecLength   = random() * SimulationSettings.get("SimulationRadius")
        VecVerAngle = random() * 2.0 * numpy.pi
        VecAziAngle = random() * 2.0 * numpy.pi
        
        self.pos[0] = VecLength * numpy.cos(VecVerAngle) * numpy.cos(VecAziAngle)
        self.pos[1] = VecLength * numpy.cos(VecVerAngle) * numpy.sin(VecAziAngle)
        self.pos[2] = VecLength * numpy.sin(VecVerAngle)
        
        
        
        VecLength   = random() * SimulationSettings.get("ObjectMaxInitialSpeed")
        VecVerAngle = random() * 2.0 * numpy.pi
        VecAziAngle = random() * 2.0 * numpy.pi
        
        self.spd[0] = VecLength * numpy.cos(VecVerAngle) * numpy.cos(VecAziAngle)
        self.spd[1] = VecLength * numpy.cos(VecVerAngle) * numpy.sin(VecAziAngle)
        self.spd[2] = VecLength * numpy.sin(VecVerAngle)
        
        
        
        VecLength   = random() * SimulationSettings.get("ObjectMaxInitialAngVel")
        VecVerAngle = random() * 2.0 * numpy.pi
        VecAziAngle = random() * 2.0 * numpy.pi
        
        self.axs[0] = VecLength * numpy.cos(VecVerAngle) * numpy.cos(VecAziAngle)
        self.axs[1] = VecLength * numpy.cos(VecVerAngle) * numpy.sin(VecAziAngle)
        self.axs[2] = VecLength * numpy.sin(VecVerAngle)
        
        
        
        self.mas = random() * SimulationSettings.get("ObjectMaxInitialMass")
        
        self.rad = random() * SimulationSettings.get("ObjectMaxInitialRadius")
        
        
        if(self.echo):
            print("ObjectBox ", self, " - has been sucessfully created with:")
            print("\tPosition: \t{", self.pos[0], ",", self.pos[1], ",", self.pos[2], "}")
            print("\tVelocity: \t{", self.spd[0], ",", self.spd[1], ",", self.spd[2], "}")
            print("\tSelfAxis: \t{", self.axs[0], ",", self.axs[1], ",", self.axs[2], "}")
            print("\tMass: \t\t\t", self.mas)
            print("\tSelfRadius: \t", self.rad)
    
        
    
    def CalcInteractionForce(self, other):
        VecForce = numpy.array([0,0,0])
        
        VecDistance = self.pos - other.pos
        
        Force = (6.67E-11 * self.mas * other.mas) / (VecDistance * VecDistance)
                 
        VecForce = -( VecDistance / numpy.linalg.norm(VecDistance) ) * Force
        
        return VecForce
        


a = ObjectBox()
b = ObjectBox()

print("Gravity a-b:", a.CalcInteractionForce(b));
print("Gravity b-a:", b.CalcInteractionForce(a));