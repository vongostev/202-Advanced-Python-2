# -*- coding: utf-8 -*-


from ObjectSystem import ObjectSystem


import random
import time

random.seed(256)

EarthSunSystem = ObjectSystem()
EarthSunSystem.Normal_Accuracy()

#Добавим солнце в нуль координат
EarthSunSystem.Add_Object_parameters([0.,0.,0.], [0.,0.,0.], 2E+30, 1.)

for i in range(8):
    EarthSunSystem.Add_Planet_orbit((random.random() + 1.0)*1.0E+11, 1.0, 0.707 + (random.random() - 0.5) / 5.0)


#EarthSunSystem.IterateTime(63.2E+6)
#EarthSunSystem.IterateTime(800.0E+6)

print("\n")
N = 400
AvTimeStep = 0;
for n in range(N):
    Time1 = time.perf_counter()
    for i in range(int(1.0E+3)):
        EarthSunSystem.Iteration()
    Time2 = time.perf_counter()
    if AvTimeStep > 0:
        AvTimeStep = (4.0 * AvTimeStep + 1.0 * (Time2 - Time1)) / 5.0
    else:
        AvTimeStep = (Time2 - Time1)
    print("Escaped Time:", int(10.0 * AvTimeStep * (N - n - 1)) / 10.0, "seconds")
    
    
    
    

EarthSunSystem.draw()






# print("Calculating time...")
# T1 = time.perf_counter()

# for i in range(1000):
#     EarthSunSystem.Iteration()

# T2 = time.perf_counter()


# print("\nOne oper time:", ((T2 - T1) / 1000.0)*1000000.0, "ns")