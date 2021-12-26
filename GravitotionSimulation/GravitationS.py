# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from celluloid import Camera

fig = plt.figure(figsize=(10, 10), dpi=100)
ax = fig.add_subplot(111, projection='3d')
camera = Camera(fig)


from random import random
import numpy


import math

def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
  return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))



SimulationSettings = {}
SimulationSettings.update({     "SimulationRadius"          : 30.0E+10    })
SimulationSettings.update({     "ObjectMaxInitialSpeed"     : 15.0E+6    })
SimulationSettings.update({     "ObjectMaxInitialAngVel"    : 62.0000     })
SimulationSettings.update({     "ObjectMaxInitialMass"      : 2.00E+25    })
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
        
        
        
        self.trajectoryX = [];
        self.trajectoryY = [];
        self.trajectoryZ = [];
        
        
        if(self.echo):
            print("ObjectBox ", self, " - has been sucessfully created with:")
            print("\tPosition: \t{", self.pos[0], ",", self.pos[1], ",", self.pos[2], "}")
            print("\tVelocity: \t{", self.spd[0], ",", self.spd[1], ",", self.spd[2], "}")
            print("\tSelfAxis: \t{", self.axs[0], ",", self.axs[1], ",", self.axs[2], "}")
            print("\tMass: \t\t\t", self.mas)
            print("\tSelfRadius: \t", self.rad)
     
        
    def update(self):
        self.trajectoryX.append(self.pos[0])
        self.trajectoryY.append(self.pos[1])
        self.trajectoryZ.append(self.pos[2])
        
    
    def CalcInteractionForce(self, other):
        
        VecForce = numpy.array([0,0,0])
        
        VecDistance = self.pos - other.pos
        Distance    = numpy.linalg.norm(VecDistance)
        
        Force = (6.67E-11 * self.mas * other.mas) / (Distance * Distance)
                 
        VecForce = -( VecDistance / Distance ) * Force
        
        return VecForce
        
    
    def data(self):
        print("ObjectBox ", self, " - has been sucessfully created with:")
        print("\tPosition: \t{", self.pos[0], ",", self.pos[1], ",", self.pos[2], "} :", numpy.linalg.norm(self.pos))
        print("\tVelocity: \t{", self.spd[0], ",", self.spd[1], ",", self.spd[2], "} :", numpy.linalg.norm(self.spd))
        print("\tSelfAxis: \t{", self.axs[0], ",", self.axs[1], ",", self.axs[2], "} :", numpy.linalg.norm(self.axs))
        print("\tMass: \t\t\t", self.mas)
        print("\tSelfRadius: \t", self.rad)

class ObjectSystem:
    
    def __init__(self):
        self.objects = []
    
        
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
        self.add_object(temp)
        
    def add_object_period(self, period, RN, mas):
        temp = Object()
        
        M1 = 0;
        M2 = numpy.array([0, 0, 0]);
        for obj in self.objects:
            M1 = M1 + obj.mas
            M2 = M2 + obj.mas * obj.pos
        
        mass_pos = M2 / M1
        
        PosRad = numpy.abs( 0.5*(RN**2)*6.67E-11*M1*((period/numpy.pi)**2)**(1.0/3.0) )
        VelRad = RN * numpy.sqrt(6.67E-11*2*M1/PosRad)
        
        AngleVer    = random() * 2.0 * numpy.pi
        AngleAzi    = random() * 2.0 * numpy.pi
        RandomVec   = numpy.array([numpy.cos(AngleVer) * numpy.cos(AngleAzi), numpy.cos(AngleVer) * numpy.sin(AngleAzi), numpy.sin(AngleVer)])
        
        temp.pos = PosRad * RandomVec + mass_pos
        
        AngleVer    = random() * 2.0 * numpy.pi
        AngleAzi    = random() * 2.0 * numpy.pi
        RandomVec   = numpy.array([numpy.cos(AngleVer) * numpy.cos(AngleAzi), numpy.cos(AngleVer) * numpy.sin(AngleAzi), numpy.sin(AngleVer)])
        
        VelocityVec = numpy.cross(temp.pos, RandomVec)
        VelocityVec = VelocityVec / numpy.linalg.norm(VelocityVec)
        VelocityVec = VelocityVec * VelRad
        
        temp.spd = VelocityVec
        temp.mas = mas
        
        self.add_object(temp)
        
    def add_object_period2(self, period, Beta, mas):
        temp = Object()
        
        M1 = 0;
        M2 = numpy.array([0, 0, 0]);
        for obj in self.objects:
            M1 = M1 + obj.mas
            M2 = M2 + obj.mas * obj.pos
        
        mass_pos = M2 / M1
        
        PosRad = (Beta*Beta*6.67E-11*M1*((0.5*period/numpy.pi)**2))**(1.0/3.0)
        VelRad = 2.0 * numpy.pi * PosRad / period
        
        AngleVer    = random() * 2.0 * numpy.pi
        AngleAzi    = random() * 2.0 * numpy.pi
        RandomVec   = numpy.array([numpy.cos(AngleVer) * numpy.cos(AngleAzi), numpy.cos(AngleVer) * numpy.sin(AngleAzi), numpy.sin(AngleVer)])
        
        temp.pos = PosRad * RandomVec + mass_pos
        
        AngleVer    = random() * 2.0 * numpy.pi
        AngleAzi    = random() * 2.0 * numpy.pi
        RandomVec   = numpy.array([numpy.cos(AngleVer) * numpy.cos(AngleAzi), numpy.cos(AngleVer) * numpy.sin(AngleAzi), numpy.sin(AngleVer)])
        
        VelocityVec1 = numpy.cross(temp.pos, RandomVec)
        VelocityVec1 = VelocityVec1 / numpy.linalg.norm(VelocityVec1)
        
        VelocityVec2 = (temp.pos - mass_pos) / numpy.linalg.norm(temp.pos - mass_pos)
        
        Rn1 = random()
        Rn2 = random()
        
        VelocityVec1 = VelocityVec1 * Rn1
        VelocityVec2 = VelocityVec2 * Rn2
        
        VelocityVec = VelocityVec1 + VelocityVec2
        VelocityVec = VelocityVec / numpy.linalg.norm(VelocityVec)
        
        temp.spd = VelocityVec * VelRad
        
        temp.mas = mas
        
        self.add_object(temp)
        
        
        
        
    def object_list(self):
        for obj in self.objects:
            obj.data()
            print("\n\n")
        
        
    def Calc_Energy(self):
        energy = 0;
        for i in range(len(self.objects)):
            energy = energy + 0.5 * self.objects[i].mas * (numpy.linalg.norm(self.objects[i].spd)**2)
            for q in range(i + 1, len(self.objects)):
                energy = energy - 6.67E-11 * self.objects[i].mas * self.objects[q].mas / numpy.linalg.norm(self.objects[i].pos - self.objects[q].pos)
        return energy
    
    def Calc_Impulse(self):
        impulse = numpy.array([0,0,0])
        for obj in self.objects:
            impulse = impulse + obj.mas * obj.spd
        return impulse;
            
    
    def Calc_Forces(self):
        forces = [numpy.array([0.0,0.0,0.0])] * len(self.objects)
        
        for i in range(len(self.objects)):
            for q in range(i + 1, len(self.objects)):
                IForce = self.objects[i].CalcInteractionForce(self.objects[q])
                forces[i] = forces[i] + IForce
                forces[q] = forces[q] - IForce
                
        return forces
    
    
    def Calc_Iteration(self, dt):
        BegEnergy = self.Calc_Energy()
        
        InterForces = self.Calc_Forces()
        
        for i in range(len(self.objects)):
            self.objects[i].spd = self.objects[i].spd + (InterForces[i] / self.objects[i].mas) * dt
            self.objects[i].pos = self.objects[i].pos + self.objects[i].spd * dt + 0.5 * (InterForces[i] / self.objects[i].mas) * dt * dt
            self.objects[i].update()
        
        EndEnergy = self.Calc_Energy()
        
        EnergyError = 2.0 * (BegEnergy - EndEnergy) / (BegEnergy + EndEnergy)
        
        assert EnergyError <= 1E+2
        
    def Calc_time(self, Time, dt):
        t = 0;
        
        while t <= Time:
            t = t + dt;
            self.Calc_Iteration(dt)
    
    
    def Create_Animation(self, Step, End, Begin):
        
        colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:olive", "tab:gray", "tab:cyan"]
        
        for a in range(Begin, End, Step):
            for i in range(len(self.objects)):
                X = self.objects[i].trajectoryX[::-1][:a]
                Y = self.objects[i].trajectoryY[::-1][:a]
                Z = self.objects[i].trajectoryZ[::-1][:a]
                ax.plot3D(X, Y, Z,"--", color=colors[i%10], linewidth=1.5)
                ax.scatter(X[::-1][:2], Y[::-1][:2], Z[::-1][:2], s=40, color=colors[i%10])
            
            ax.view_init(elev=25, azim=45+a)
            
            camera.snap()
        
        animation = camera.animate(blit=False,interval=100)
        animation.save("filename.gif")
    
    def Create_NoTrack_Animation(self, Step, End, Begin):
        
        colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:olive", "tab:gray", "tab:cyan"]
        
        X = [0, 0]
        Y = [0, 0]
        Z = [0, 0]
        
        for a in range(Begin, End, Step):
            for i in range(len(self.objects)):
                X[1] = X[0]
                X[0] = self.objects[i].trajectoryX[a]
                
                Y[1] = Y[0]
                Y[0] = self.objects[i].trajectoryY[a]
                
                Z[1] = Z[0]
                Z[0] = self.objects[i].trajectoryZ[a]
                
                ax.scatter(X, Y, Z, s=40, color=colors[i%10])
            
            ax.view_init(elev=25, azim=45+a)
            
            camera.snap()
        
        animation = camera.animate(blit=False,interval=100)
        animation.save("filename.gif")
    
    def draw(self, a = 25, b = -45):
        ax.clear()
        
        ax.view_init(elev=a, azim=b)
        
        colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:olive", "tab:gray", "tab:cyan"]
        
        for i in range(len(self.objects)):
            ax.plot3D(self.objects[i].trajectoryX, self.objects[i].trajectoryY, self.objects[i].trajectoryZ,"--", color=colors[i%10], linewidth=1.5)
            ax.scatter(self.objects[i].trajectoryX[::-1][:2], self.objects[i].trajectoryY[::-1][:2], self.objects[i].trajectoryZ[::-1][:2], s=40, color=colors[i%10])
            
        fig.canvas.draw()
    
    
    

    
System = ObjectSystem()

System.add_object_by_parameters([0,0,0], [0,0,0], [0,1,0], 2E+10,     1E+9)

System.add_object_period(30,  0.707, 7E+4)
System.add_object_period(20,  0.9, 7E+4)
System.add_object_period(30,  0.7, 7E+8)

# System.add_object_by_parameters([+0.5,0,0], [0,+1.08,0], [0,1,0], 2E+10,     1E+9)
# System.add_object_by_parameters([-0.5,0,0], [0,-1.08,0], [0,1,0], 2E+10,     1E+9)

# System.add_object_period(5,   0.707, 1E+3)
# System.add_object_period(10,  0.707, 5E+3)
# System.add_object_period(20,  0.707, 4E+4)
# System.add_object_period(30,  0.707, 7E+4)
# System.add_object_period(60,  0.707, 6E+8)


System.object_list()
print("System Energy:  ", System.Calc_Energy())
print("System Impulse: ", System.Calc_Impulse())


print("Angle:", angle(System.objects[1].pos, System.objects[1].spd))

System.Calc_time(120, 0.002 )

System.draw()
   
plt.show()

#System.Create_Animation(250, 30000, 0)
System.Create_NoTrack_Animation(300, 60001, 0)

print("System Energy:  ", System.Calc_Energy())
print("System Impulse: ", System.Calc_Impulse())

print("Finish", len(System.objects[0].trajectoryX))


