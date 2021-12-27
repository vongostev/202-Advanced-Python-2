# -*- coding: utf-8 -*-


from ObjectSystem import ObjectSystem


import random
import time

random.seed(256)

Sys = ObjectSystem()
Sys.Normal_Accuracy()

#Добавим солнце в нуль координат
Sys.Add_Object_parameters([0.,0.,0.], [0.,0.,0.], 2E+30, 1.)

for i in range(10):
    Sys.Add_Planet_orbit((random.random() + 1.0)*1.0E+11, 1.0, 0.707 + (random.random() - 0.5) / 5.0)

#Sys.Add_Planet_orbit((random.random() + 1.0)*1.0E+11, 1.0, 0.6)

#Sys.Add_Planet_orbit((random.random() + 1.0)*1.0E+11, 1.0, 0.707)
#Sys.objects[len(Sys.objects)-1].Velocity *= 1.5

#EarthSunSystem.IterateTime(63.2E+6)

print("\n")
N = 1000
AvTimeStep = 0;
for n in range(N):
    Time1 = time.perf_counter()
    for i in range(int(1.0E+3)):
        Sys.Iteration()
    Time2 = time.perf_counter()
    if AvTimeStep > 0:
        AvTimeStep = (4.0 * AvTimeStep + 1.0 * (Time2 - Time1)) / 5.0
    else:
        AvTimeStep = (Time2 - Time1)
    print("Escaped Time:", int(10.0 * AvTimeStep * (N - n - 1)) / 10.0, "seconds")
    
    
    

Sys.draw()

Sys.AnalizeOrbits()

Sys.animate(0, Sys.T, 2000, 30, 1E+7)