# -*- coding: utf-8 -*-


from ObjectSystem import ObjectSystem


import random
import time

random.seed(5484)

EarthSunSystem = ObjectSystem()
EarthSunSystem.Normal_Accuracy()

#Добавим солнце в нуль координат
EarthSunSystem.Add_Object_parameters([0.,0.,0.], [0.,0.,0.], 2E+30, 1.)

#Скажем системе самой создать планету на расстояние 1.5E+11 от центра масс (солнца) 
EarthSunSystem.Add_Planet_orbit(1.5E+11, 1.0, 0.78)
EarthSunSystem.Add_Planet_orbit(1.6E+11, 1.0, 0.76)
EarthSunSystem.Add_Planet_orbit(1.7E+11, 1.0, 0.74)
#EarthSunSystem.Add_Planet_orbit(1.8E+11, 1.0, 0.72)
#EarthSunSystem.Add_Planet_orbit(1.9E+11, 1.0, 0.70)

#Запустим симуляцию на два года
EarthSunSystem.IterateTime(63.2E+6 / 10.0)

EarthSunSystem.draw()


T1 = time.perf_counter()

for i in range(1000):
    EarthSunSystem.CalcEnergy()

T2 = time.perf_counter()


print("\nOne oper time:", ((T2 - T1) / 1000.0)*1000000.0, "ns")