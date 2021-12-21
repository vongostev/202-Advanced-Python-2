import unittest

import numpy as np
from Tests import Tests
from cosmicbody import CosmicBody2D, CosmicBody3D
from Star import Star



t_end = 360*24*60
G = 6.674e-11
M = 1.98847e30
m = 5.9722e24
R = 696e6
r = 1496e8
r_p = 147098291e3
v = np.sqrt(2 * G * M / r) ##вторая космическая
u = np.sqrt(M*G*(1+0.017)/r_p) ##скорость Земли отн Солнца в перигелии

"Модель: Земля - Солнце"
star = Star(M,R)
a1 = CosmicBody2D(m, np.array([0, u]), np.array([r_p, 0]))
print(a1.orbit_type(star))
CosmicBody2D.graph(star, t_end, [a1])



"Проверка известных траекторий тел, вылетающих с орбиты земли с разными космическими скоростями "

t_end = 3000
earth = Star(m)
v_2 = np.sqrt(G * m / 6371e3)

body1 = CosmicBody3D(1e3, np.array([7.91e3 / np.sqrt(2), 7.91e3 / np.sqrt(2), 0]),
                     np.array([0, 0, 6371e3]))
body2 = CosmicBody3D(1e3, np.array([v_2, v_2, 0]),
                     np.array([0, 0, 6371e3]))
body3 = CosmicBody3D(1e3, np.array([16.7e3 / np.sqrt(2), 16.7e3 / np.sqrt(2), 0]),
                     np.array([0, 0, 6371e3]))
CosmicBody3D.graph(earth, t_end, [body1,body2,body3])



body1 = CosmicBody2D(1e22, np.array([7.91e3, 0]),
                     np.array([0, 6371e3]))
body2 = CosmicBody2D(1e22, np.array([np.sqrt(2 * G * m / 6371e3), 0]),
                     np.array([0, 6371e3]))
body3 = CosmicBody2D(1e22, np.array([16.7e3, 0]),
                     np.array([0, 6371e3]))
print(body1.orbit_type(earth))
print(body2.orbit_type(earth))
print(body3.orbit_type(earth))
CosmicBody2D.graph(earth, t_end, [body1,body2,body3])

"Рандом 2D"
t_end_random = 100
star1 = Star(M)
a_1 = CosmicBody2D.create_random_cosmic_body()
CosmicBody2D.graph(star1, t_end_random, [a_1])

"Рандом 3D"

star2 = Star(M)
a_2 = CosmicBody3D.create_random_cosmic_body()
CosmicBody3D.graph(star2, t_end_random, [a_2])

if __name__ == '__main__':
    unittest.main()