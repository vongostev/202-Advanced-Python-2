import numpy as np
import matplotlib.pyplot as plt
import random
from mpl_toolkits import mplot3d
from numba import njit, prange
import unittest
#import numba
G = 0.1
dt = 0.001

class Star:
    def __init__(self, mass:float):
        self.mass = mass
    "Здесь описывается класс звезды"

class CosmicBody:
    def __init__(self, mass: float, vec_v: np.ndarray, vec_r: np.ndarray):
        self.mass = mass
        self.vec_v = vec_v
        self.vec_r = vec_r
    "класс тела космического"
    
    #@nb.njit(fastmath=True)
    def gravitate(self, Star : Star):
        return - G * self.mass * Star.mass * self.vec_r / np.linalg.norm(self.vec_r) ** 3
    
    #@njit(fastmath=True, cache=True,parallel=True)
    def destroy(self):
        if abs(self.vec_r[0]) < 0.01 and abs(self.vec_r[1]) < 0.01:
            self.mass = 0
            return 1
        "При радиусе = 0 происходит столкновение"
        
    #@njit(fastmath=True, cache=True,parallel=True)
    def track(self, Star : Star):
        if (self.destroy() != 1):
            delta_v = self.gravitate(Star) * dt / self.mass
            self.vec_r = self.vec_r + self.vec_v * dt + delta_v * dt * dt / 2
            self.vec_v = delta_v + self.vec_v
        "функция описывает изменнеие параметров космического тела за время dt"
        "по хорошему лучше здесь фунцию gravitate переписать и не пересчитвать массу"
    
    #@njit(fastmath=True, cache=True,parallel=True)
    def trajectory_2D(self, Star : Star, time_beggining, time_stop):
        position_x = list()
        position_y = list()
        #while(time_beggining < time_stop):
            #position_x.append(self.vec_r[0])
            #position_y.append(self.vec_r[1])
            #time.append(time_beggining)
            #self.track(Star)
            #time_beggining += dt #вроде как глобальная time_beggining должна меняться, но это не проблема
        for i in np.arange(time_beggining, time_stop, dt):
            position_x.append(self.vec_r[0])
            position_y.append(self.vec_r[1])
            self.track(Star)
        #было заменено с while на arange, так быстрее получилось
        return [position_x, position_y] 
    "функция описывает траекторию 2D тела за [t1, t2] секунд"
    
    #@njit(fastmath=True, cache=True,parallel=True)
    def graf_trajectory_2D(self, Star : Star, time_beggining, time_stop):
        trajectory = self.trajectory_2D(Star, time_beggining, time_stop)  
        plt.plot(trajectory[0], trajectory[1], 'ro')
        plt.show()
    "строит график для одной кометы в 2D"
    
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
    "строит график для нескольких тел"
    "можно было бы и удалить функцию graf_trajectory_2D, тк эта функция включает ее сама по себе, но я начинал именно с graf_trajectory_2D"
        
    #@njit(fastmath=True, cache=True,parallel=True)
    def what_trajectory(self, Star : Star):
        E = self.mass * np.linalg.norm(self.vec_v) * np.linalg.norm(self.vec_v) / 2 - G * self.mass * Star.mass / np.linalg.norm(self.vec_r)
        if (E > 1e-5):
            return "hyperbole"
        if (E < 1e-5):
            return "ellipse"
        return "parabol"
    "определяет тип траектрии"
    
    #@njit(fastmath=True, cache=True,parallel=True)
    def trajectory_3D(self, Star : Star, time_beggining, time_stop):
        position_x = list()
        position_y = list()
        position_z = list()
        #while(time_beggining < time_stop):
            #position_x.append(self.vec_r[0])
            #position_y.append(self.vec_r[1])
            #position_z.append(self.vec_r[2])
            #time.append(time_beggining)
            #self.track(Star)
            #time_beggining += dt #вроде как глобальная time_beggining должна меняться, но это не проблема
        for i in np.arange(time_beggining, time_stop, dt):
            position_x.append(self.vec_r[0])
            position_y.append(self.vec_r[1])
            position_z.append(self.vec_r[2])
            self.track(Star) 
        #было заменено с while на arange, так быстрее получилось
        return [position_x, position_y, position_z]
    "функция описывает траекторию 3D тела за [t1, t2] секунд"
    
    def graf_trajectory_3D_more(self, Star : Star, time_beggining, time_stop, list_Cosmos_Badys = []):
        k = 0
        color_list = ['mediumorchid', 'blueviolet', 'navy', 'royalblue', 'darkslategrey', 'limegreen', 'darkgreen']
        trajectory = self.trajectory_3D(Star, time_beggining, time_stop)
        #fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.plot3D(trajectory[0], trajectory[1], trajectory[2], 'red')
        k += 1
        for i in list_Cosmos_Badys:
            trajectory_i = i.trajectory_3D(Star, float(input()), time_stop)
            ax.plot3D(trajectory_i[0], trajectory_i[1], trajectory_i[2], color = color_list[k])
            k += 1    
        plt.show()
    "строит график для нескольких тел в 3D"
    
    
    
"В этой части кода происходить будет ускорение кода, видимо в классах это нельзя реализовать"
"думал думал, только изменил while на np.range"


class TestGrafic(unittest.TestCase):
    def test_2D_destroy(self):
        G = 0.1
        dt = 0.001
        t_beg = 0
        t_end = 20
        a1 = CosmicBody(2, np.array([0, 0]), np.array([3, 4]))
        star = Star(5)
        a1.graf_trajectory_2D(star, t_beg, t_end)
        print("---------------")
    
    def test_2D_ellipse(self):
        G = 0.1
        dt = 0.001
        t_beg = 0
        t_end = 110
        star = Star(5)
        a2 = CosmicBody(11, np.array([0, -0.3]), np.array([-5, -2]))
        a2.graf_trajectory_2D(star, t_beg, t_end)
        print(a2.what_trajectory(star))
        print("---------------")
        
    def test_2D_hyperbole(self):
        star = Star(5)
        G = 0.1
        dt = 0.001
        t_beg = 0
        t_end = 20
        a3 = CosmicBody(11, np.array([-1, 1/3]), np.array([10, -2]))
        a3.graf_trajectory_2D(star, t_beg, t_end)
        print(a3.what_trajectory(star))
        print("---------------")
        
    def test_2D_random_CosmicBody(self):
        G = random.random() * 10
        dt = 0.001
        t_beg = 0
        t_end = 10
        star = Star(random.random() * 10)
        a4 = CosmicBody(random.random() * 10, np.array([random.random(), random.random()]), np.array([random.random() * 2, random.random() * 2]))
        a4.graf_trajectory_2D(star, t_beg, t_end)
        print(a4.what_trajectory(star))
        print("---------------")
        
    def test_2D_more_CosmicBody(self):
        G = 0.1
        dt = 0.001
        star = Star(5)
        t_beg = 0
        t_end = 110
        a2 = CosmicBody(11, np.array([0, -0.3]), np.array([-5, -2]))
        a3 = CosmicBody(11, np.array([-1, 1/3]), np.array([10, -2]))
        a2.graf_trajectory_2D_more(star, t_beg, t_end, [a3])
        "не забываем писать время, когда второе космическое тело появилось"
    
    def test_3D_destroy(self):
        G = 0.1
        dt = 0.001
        t_beg = 0
        t_end = 40
        a1 = CosmicBody(2, np.array([0, 0, 0]), np.array([3, 4, 5]))
        star = Star(5)
        a1.graf_trajectory_3D_more(star, t_beg, t_end)
        
    def test_3D_ellipse(self):
        star = Star(5)
        G = 0.1
        dt = 0.001
        t_beg = 0
        t_end = 220
        a2 = CosmicBody(11, np.array([0, -0.3, -0.1]), np.array([-5, -2, 3]))
        a2.graf_trajectory_3D_more(star, t_beg, t_end)
        print(a2.what_trajectory(star))
        print("---------------")
    
    def test_3D_hyperbole(self):
        star = Star(5)
        G = 0.1
        dt = 0.001
        t_beg = 0
        t_end = 220
        a2 = CosmicBody(11, np.array([0, -0.3, -0.4]), np.array([-5, -2, 3]))
        a2.graf_trajectory_3D_more(star, t_beg, t_end)
        print(a2.what_trajectory(star))
        print("---------------")
        
    def test_3D_more_CosmicBody(self):
        star = Star(5)
        G = 0.1
        dt = 0.001
        t_beg = 0
        t_end = 220
        a2 = CosmicBody(11, np.array([0, -0.3, -0.1]), np.array([-5, -2, 3]))
        a3 = CosmicBody(11, np.array([0, -0.3, -0.4]), np.array([-5, -2, 3]))
        a2.graf_trajectory_3D_more(star, t_beg, t_end, [a3])
        "не забываем писать время, когда второе космическое тело появилось"
        
        
        
if __name__ == "__main__":
    unittest.main()
    
    
    #сделал юнит тесты
"""
    "тесты"
    
    "2D"
    
    "падение на звезду"
    G = 0.1
    dt = 0.001
    t_beg = 0
    t_end = 20
    a1 = CosmicBody(2, np.array([0, 0]), np.array([3, 4]))
    star = Star(5)
    a1.graf_trajectory_2D(star, t_beg, t_end)
    
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
    
    
    "3D"
    "v = 0"
    G = 0.1
    dt = 0.001
    t_beg = 0
    t_end = 40
    a1 = CosmicBody(2, np.array([0, 0, 0]), np.array([3, 4, 5]))
    star = Star(5)
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
"""
    
        