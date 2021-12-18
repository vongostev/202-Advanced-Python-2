import numpy as np
import matplotlib.pyplot as plt
import random
from mpl_toolkits import mplot3d

G = 0.1
dt = 0.01

class Star:
    def __init__(self, mass:float):
        self.mass = mass
    "Здесь описывается класс звезды"

class CosmicBody:
    def __init__(self, mass: float, vec_v: np.ndarray, radius: np.ndarray):
        self.mass = mass
        self.vec_v = vec_v
        self.radius = radius
    "класс тела космического"

    def gravitate(self, Star : Star):
        return - G * self.mass * Star.mass * self.radius / np.linalg.norm(self.radius) ** 3
    
    def destroy(self):
        #if self.radius[0] < 0.01 and self.radius[1] < 0.01:
        if abs(self.radius[0]) < 0.01 and abs(self.radius[1]) < 0.01:
            self.mass = 0
            return 1
        "При радиусе = 0 происходит столкновение"
        
    def track(self, Star : Star):
        if (self.destroy() != 1):
            delta_v = self.gravitate(Star) * dt / self.mass
            self.radius = self.radius + self.vec_v * dt + delta_v * dt * dt / 2
            self.vec_v = delta_v + self.vec_v
        
        "по хорошему лучше здесь фунцию gravitate переписать и не пересчитвать массу"
        
    def trajectory_2D(self, Star : Star, time_beggining, time_stop):
        position_x = list()
        position_y = list()
        time = list()
        while(time_beggining < time_stop):
            position_x.append(self.radius[0])
            position_y.append(self.radius[1])
            time.append(time_beggining)
            #print(self.radius)
            self.track(Star)
            time_beggining += dt #вроде как глобальная time_beggining должна меняться, но это не проблема
            E = self.mass * np.linalg.norm(self.vec_v) * np.linalg.norm(self.vec_v) / 2 - G * self.mass * Star.mass / np.linalg.norm(self.radius)
        return [position_x, position_y, time]
    
    def graf_trajectory_2D(self, Star : Star, time_beggining, time_stop):
        trajectory = self.trajectory_2D(Star, time_beggining, time_stop)  
        plt.plot(trajectory[0], trajectory[1], 'ro')
        plt.show()
    
    def graf_trajectory_2D_more(self, Star : Star, time_beggining, time_stop, list_Cosmos_Badys = []):
        k = 0
        color_list = ['mediumorchid', 'blueviolet', 'navy', 'royalblue', 'darkslategrey', 'limegreen', 'darkgreen']
        trajectory = self.trajectory_2D(Star, time_beggining, time_stop)
        plt.plot(trajectory[0], trajectory[1], 'ro', color = color_list[k])
        k += 1
        for i in list_Cosmos_Badys:
            trajectory_i = i.trajectory_2D(Star, float(input()), time_stop)
            plt.plot(trajectory_i[0], trajectory_i[1], 'ro', color = color_list[k])
            k += 1
        plt.show()
            
        
        
    def what_trajectory(self, Star : Star):
        E = self.mass * np.linalg.norm(self.vec_v) * np.linalg.norm(self.vec_v) / 2 - G * self.mass * Star.mass / np.linalg.norm(self.radius)
        if (E > 1e-5):
            return "hyperbole"
        if (E < 1e-5):
            return "ellipse"
        return "parabol"
    
    def trajectory_3D(self, Star : Star, time_beggining, time_stop):
        position_x = list()
        position_y = list()
        position_z = list()
        time = list()
        while(time_beggining < time_stop):
            position_x.append(self.radius[0])
            position_y.append(self.radius[1])
            position_z.append(self.radius[2])
            time.append(time_beggining)
            #print(self.radius)
            self.track(Star)
            time_beggining += dt #вроде как глобальная time_beggining должна меняться, но это не проблема
            E = self.mass * np.linalg.norm(self.vec_v) * np.linalg.norm(self.vec_v) / 2 - G * self.mass * Star.mass / np.linalg.norm(self.radius)
        return [position_x, position_y, position_z, time]
    
    def graf_trajectory_3D_more(self, Star : Star, time_beggining, time_stop, list_Cosmos_Badys = []):
        k = 0
        color_list = ['mediumorchid', 'blueviolet', 'navy', 'royalblue', 'darkslategrey', 'limegreen', 'darkgreen']
        trajectory = self.trajectory_3D(Star, time_beggining, time_stop)
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.plot3D(trajectory[0], trajectory[1], trajectory[2], 'red')
        k += 1
        for i in list_Cosmos_Badys:
            trajectory_i = i.trajectory_3D(Star, float(input()), time_stop)
            ax.plot3D(trajectory_i[0], trajectory_i[1], trajectory_i[2], color = color_list[k])
            k += 1    
        plt.show()
        
        
"""
a = CosmicBody(2, np.array([1, 2]), np.array([3, 4]))
star = Star(5)
print(a.gravitate(star))

e = np.array([1, 2])
d = np.array([6, 10])
print(d - e)
g = np.array([3, 4])
print(g / np.linalg.norm(g) ** (3 / 2))
print(2 * [1, 2])

"""

"тесты"

"падение на звезду"
G = 0.1
dt = 0.001
t_beg = 0
t_end = 20
a1 = CosmicBody(2, np.array([0, 0]), np.array([3, 4]))
star = Star(5)
#print(a1.what_trajectory(star))
a1.graf_trajectory_2D(star, t_beg, t_end)

"""
a = CosmicBody(2, np.array([0.03, 0.1]), np.array([3, 4]))
print(a.what_trajectory(star))
a.trajectory_2D(star, t_beg, t_end)
"""


"""
"Земля и Солнце"
G = 6.6743e-11
dt = 60 * 60 * 24 / 100
t_beg = 0
t_end = 365 * 60 * 60 * 24
Sun = Star(1.98892e30)
Earth = CosmicBody(5.9722e24, np.array([0, 30270]) ,np.array([147098290000, 0]))
print(Earth.what_trajectory(Sun))
Earth.trajectory_2D(Sun, t_beg, t_end)
"""




"3D"
"v = 0"
G = 0.1
dt = 0.001
t_beg = 0
t_end = 40
a1 = CosmicBody(2, np.array([0, 0, 0]), np.array([3, 4, 5]))
star = Star(5)
#print(a1.what_trajectory(star))
a1.graf_trajectory_3D_more(star, t_beg, t_end)


"Эллипс"
G = 0.1
dt = 0.001
t_beg = 0
t_end = 220
a2 = CosmicBody(11, np.array([0, -0.3, -0.1]), np.array([-5, -2, 3]))
a2.graf_trajectory_3D_more(star, t_beg, t_end)
print(a2.what_trajectory(star))

"Гипербола"
G = 0.1
dt = 0.001
t_beg = 0
t_end = 220
a2 = CosmicBody(11, np.array([0, -0.3, -0.4]), np.array([-5, -2, 3]))
a2.graf_trajectory_3D_more(star, t_beg, t_end)
print(a2.what_trajectory(star))

"несколько тел"
G = 0.1
dt = 0.001
t_beg = 0
t_end = 220
a2 = CosmicBody(11, np.array([0, -0.3, -0.1]), np.array([-5, -2, 3]))
a3 = CosmicBody(11, np.array([0, -0.3, -0.4]), np.array([-5, -2, 3]))
a2.graf_trajectory_3D_more(star, t_beg, t_end, [a3])

"2D"
"эллипс"
G = 0.1
dt = 0.001
t_beg = 0
t_end = 110
a2 = CosmicBody(11, np.array([0, -0.3]), np.array([-5, -2]))
a2.graf_trajectory_2D(star, t_beg, t_end)
print(a2.what_trajectory(star))

"гипербола"
G = 0.1
dt = 0.001
t_beg = 0
t_end = 20
a3 = CosmicBody(11, np.array([-1, 1/3]), np.array([10, -2]))
a3.graf_trajectory_2D(star, t_beg, t_end)
print(a3.what_trajectory(star))


"Рандомим"
G = random.random() * 10
dt = 0.001
t_beg = 0
t_end = 10
star = Star(random.random() * 10)
a4 = CosmicBody(random.random() * 10, np.array([random.random(), random.random()]), np.array([random.random() * 2, random.random() * 2]))
a4.graf_trajectory_2D(star, t_beg, t_end)
print(a4.what_trajectory(star))

"Несколько тел"
G = 0.1
dt = 0.001
star = Star(5)
t_beg = 0
t_end = 110
a2 = CosmicBody(11, np.array([0, -0.3]), np.array([-5, -2]))
a3 = CosmicBody(11, np.array([-1, 1/3]), np.array([10, -2]))
a2.graf_trajectory_2D_more(star, t_beg, t_end, [a3])



    
        
