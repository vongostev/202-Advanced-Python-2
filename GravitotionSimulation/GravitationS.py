# -*- coding: utf-8 -*-

from random import random
import numpy


SimulationSettings = {}
SimulationSettings.update({     "SimulationRadius"          : 30.0E+12    })
SimulationSettings.update({     "ObjectMaxInitialSpeed"     : 15.0E+10    })
SimulationSettings.update({     "ObjectMaxInitialAngVel"    : 62.0000     })
SimulationSettings.update({     "ObjectMaxInitialMass"      : 2.00E+27    })
SimulationSettings.update({     "ObjectMaxInitialRadius"    : 8.00E+7     })


class Object:
    
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
        Distance    = numpy.linalg.norm(VecDistance)
        
        Force = (6.67E-11 * self.mas * other.mas) / (Distance * Distance)
                 
        VecForce = -( VecDistance / Distance ) * Force
        
        return VecForce
        

a = Object()
b = Object()

#print("Gravity a-b:", a.CalcInteractionForce(b));
#print("Gravity b-a:", b.CalcInteractionForce(a));




class ObjectSystem:
    
    def __init__(self):
        self.objects = []
        
        self.energy = 0
    
        
    def add_object(self, some_object):
        self.objects.append(some_object)
        
    def add_object_by_parameters(self, pos, spd, axs, mas, rad):
        temp = Object()
        
        temp.pos = numpy.array(pos)
        temp.spd = numpy.array(spd)
        temp.axs = numpy.array(axs)
        temp.mas = mas
        temp.rad = rad
        
        self.add_object(temp)
        
    def add_object_by_dict(self, dic):
        temp = Object();
        
        if(dic.get("pos") != None):
            temp.pos = numpy.array(dic.get("pos"))
        if(dic.get("spd") != None):
            temp.spd = numpy.array(dic.get("spd"))
        if(dic.get("axs") != None):
            temp.axs = numpy.array(dic.get("axs"))
        if(dic.get("mas") != None):
            temp.mas = dic.get("mas")  
        if(dic.get("rad") != None):
            temp.rad = dic.get("rad")
            
        self.add_object(temp)
        
        
    def Calc_Energy(self):
        energy = 0;
        for i in range(len(self.objects)):
            energy = energy + 0.5 * self.objects[i].mas * (numpy.linalg.norm(self.objects[i].spd)**2)
            for q in range(i + 1, len(self.objects)):
                energy = energy - 6.67E-11 * self.objects[i].mas * self.objects[q].mas / numpy.linalg.norm(self.objects[i].pos - self.objects[q].pos)
        return energy
    
    def Calc_Forces(self):
        forces = [numpy.array([0.0,0.0,0.0])] * len(self.objects)
        
        for i in range(len(self.objects)):
            for q in range(i + 1, len(self.objects)):
                #print(i, " - ", q)
                IForce = self.objects[i].CalcInteractionForce(self.objects[q])
                forces[i] = forces[i] + IForce
                forces[q] = forces[q] - IForce
                
        return forces
    
    
    def Calc_Iteration(self):
        BegEnergy = self.Calc_Energy()
        
        dt = 1E-4
        InterForces = self.Calc_Forces()
        
        for i in range(len(self.objects)):
            self.objects[i].spd = self.objects[i].spd + (InterForces[i] / self.objects[i].mas) * dt
            self.objects[i].pos = self.objects[i].pos + self.objects[i].spd * dt
            
        EndEnergy = self.Calc_Energy()
        
        #print("Energy factor:", BegEnergy - EndEnergy)
        
        
        
        
    
    
    
    
System = ObjectSystem();

System.add_object_by_parameters([0,0,0], [0,0,0], [0,0,1E-5], 1E+30, 1E+9)
System.add_object_by_parameters([0,1E+10,0], [1E+5,0,0], [0,1E-5,1E-5], 1E+20, 1E+5)
System.add_object_by_parameters([0,-1E+10,0], [-1E+5,0,0], [0,1E-5,1E-5], 1E+20, 1E+5)

#print("All Interforces:", System.Calc_Forces())
#print("System Energy:  ", System.Calc_Energy())

for i in range(10):
    System.Calc_Iteration()
    print(System.objects[0].pos, System.objects[1].pos, System.objects[2].pos)
