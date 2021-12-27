# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt

from celluloid import Camera
from IPython.display import Image

import random
import numpy as np

import math

import colorsys

from statistics import mean

import time as tm



GravityConst = 6.67E-11



class Object:
    
    #Праметры, общие для всех экземпляров класса Object, которые можно трогать пользователю
    Echo = False
    MaxDevAng = (np.pi / 180.0) * 1.0
    
    #Праметры, общие для всех экземпляров класса Object, которые вычисляются автоматически и недоступны для редактирования
    CosMaxDevAng = math.cos(MaxDevAng)
    
    
    def __init__(self, position = [0.0, 0.0, 0.0], velocity = [0.0, 0.0, 0.0], mass = 1.0, radius = 1.0):
        self.Position = np.array(position)
        self.Velocity = np.array(velocity)
        
        self.Mass = mass
        
        self.Radius = radius
        
        self.Trajectory = [[],[],[],[]]
        
        if self.Echo:
            print("Добавлен новый объект с параметрами:")
            print("\tПоложение в пространстве:", self.Position, "\tРасстояние до центра -", np.linalg.norm(self.Position))
            print("\tВектор скорости:         ", self.Velocity, "\tМодуль вектора -      ", np.linalg.norm(self.Velocity))
            print("\tРадиус: ", self.Radius)
            print("\tМасса:  ", self.Mass)
            print("")
            
    def CalcForce(self, other):
        ForceVec = np.array([0.0, 0.0, 0.0])
        
        DistanceVec = (self.Position - other.Position)
        DistanceMod = np.linalg.norm(DistanceVec)
        
        ForceVec = -GravityConst * (self.Mass * other.Mass * DistanceVec) / (DistanceMod * DistanceMod * DistanceMod)
        
        return ForceVec
    
    def CalcEnergy(self, other):
        Energy = 0.0
        
        DistanceMod = np.linalg.norm(self.Position - other.Position)
        
        Energy = -GravityConst * self.Mass * other.Mass / DistanceMod
        
        return Energy
    
    def CalcCollision(self, other):
        DistanceMod = np.linalg.norm(self.Position - other.Position)
        DistanceMax = self.Radius + other.Radius
        
        if DistanceMod <= DistanceMax:
            return True
        return False
    
    def UpdateTrajectory(self, time):
        LPN = len(self.Trajectory[0]) - 1
        
        if LPN <= 0:
            self.Trajectory[0].append(time)
            self.Trajectory[1].append(self.Position[0])
            self.Trajectory[2].append(self.Position[1])
            self.Trajectory[3].append(self.Position[2])
            return
        
        Point0 = np.array([self.Trajectory[1][LPN - 1], self.Trajectory[2][LPN - 1], self.Trajectory[3][LPN - 1]])
        Point1 = np.array([self.Trajectory[1][LPN - 0], self.Trajectory[2][LPN - 0], self.Trajectory[3][LPN - 0]])
        Point2 = self.Position
        
        Vector0 = Point1 - Point0
        Vector1 = Point2 - Point1
        
        Vector0 = Vector0 / np.linalg.norm(Vector0)
        Vector1 = Vector1 / np.linalg.norm(Vector1)
        
        CosDevAng = np.dot(Vector0, Vector1)
        
        if CosDevAng < self.CosMaxDevAng:
            self.Trajectory[0].append(time)
            self.Trajectory[1].append(self.Position[0])
            self.Trajectory[2].append(self.Position[1])
            self.Trajectory[3].append(self.Position[2])
         
    def ForceUpdateTrajectory(self, time):
        self.Trajectory[0].append(time)
        self.Trajectory[1].append(self.Position[0])
        self.Trajectory[2].append(self.Position[1])
        self.Trajectory[3].append(self.Position[2])
        
    def print(self):
        print("Объект:")
        print("\tПоложение в пространстве:", self.Position, "\tРасстояние до центра -", np.linalg.norm(self.Position))
        print("\tВектор скорости:         ", self.Velocity, "\tМодуль вектора -      ", np.linalg.norm(self.Velocity))
        print("\tРадиус: ", self.Radius)
        print("\tМасса:  ", self.Mass)
        print("")