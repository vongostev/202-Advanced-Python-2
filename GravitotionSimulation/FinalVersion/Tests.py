# -*- coding: utf-8 -*-


from ObjectSystem import ObjectSystem


import random
import time as time

random.seed(time.asctime())

EarthSunSystem = ObjectSystem()
EarthSunSystem.Normal_Accuracy()

#Добавим солнце в нуль координат
EarthSunSystem.Add_Object_parameters([0.,0.,0.], [0.,0.,0.], 2E+30, 1.)

#Скажем системе самой создать планету на расстояние 1.5E+11 от центра масс (солнца) 
EarthSunSystem.Add_Planet_orbit(1.5E+11, 1.0, 0.707)

#Запустим симуляцию на два года
EarthSunSystem.IterateTime(63.2E+6)

EarthSunSystem.draw()