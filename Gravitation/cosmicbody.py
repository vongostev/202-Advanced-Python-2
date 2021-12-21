import random
from abc import ABCMeta, abstractmethod
import numpy as np
import matplotlib.pyplot as plt
from numba import njit

from Star import Star

dt = 1000
G = 6.674e-11


@njit('float64(float64[:])', cache=True,
      nogil=False, fastmath=True)
def norm(vec: np.ndarray):
    return np.linalg.norm(vec)


class CosmicBody(metaclass=ABCMeta):
    @staticmethod
    @abstractmethod
    def create_random_cosmic_body():
        pass

    @abstractmethod
    def destroy(self, gravitate_star: Star):
        pass

    @staticmethod
    @abstractmethod
    def graph(gravitate_star: Star, time_stop, list_cosmic_bodies=None):
        if list_cosmic_bodies is None:
            list_cosmic_bodies = []

    @abstractmethod
    def locus(self, gravitate_star: Star, time_stop):
        pass

    def __init__(self, mass: float, velocity: np.ndarray, position: np.ndarray):
        self.mass = mass
        self.velocity = velocity
        self.position = position

    def gravitate(self, gravitate_star: Star):
        return - G * self.mass * gravitate_star.mass * self.position / np.power(np.linalg.norm(self.position), 3)

    def step(self, gravitate_star: Star):
        if self.destroy(gravitate_star) != 1:
            delta_v = self.gravitate(gravitate_star) * dt / self.mass
            self.position = self.position + self.velocity * dt + delta_v * np.float_power(dt, 2) / 2
            self.velocity = delta_v + self.velocity

    def orbit_type(self, gravitate_star: Star):
        e = np.float_power(np.linalg.norm(self.velocity),
                           2) / 2 - G * gravitate_star.mass / np.linalg.norm(self.position)
        print(e)
        if e > 0:
            return "hyperbola"
        if e < 0:
            return "ellipse"
        return "parabola"


class CosmicBody2D(CosmicBody):
    @staticmethod
    def create_random_cosmic_body():
        return CosmicBody2D(random.uniform(1e5, 1e7),
                            np.array([random.uniform(0, 100000),
                                      (random.uniform(0, 100000))]),
                            np.array([(random.uniform(-10000, 10000)),
                                      random.uniform(-10000, 10000)]))

    def destroy(self, gravitate_star: Star):
        if abs(self.position[0]) < gravitate_star.radius and abs(self.position[1]) < gravitate_star.radius:
            self.mass = 0
            return 1
        return 0

    def locus(self, gravitate_star: Star, time_stop):
        position_x = []
        position_y = []
        time = 0
        while time < time_stop:
            position_x.append(self.position[0])
            position_y.append(self.position[1])
            self.step(gravitate_star)
            if self.destroy(gravitate_star):
                break
            time += dt
        return [position_x, position_y]

    @staticmethod
    def graph(gravitate_star: Star, time_stop, list_cosmic_bodies=None):
        if list_cosmic_bodies is None:
            list_cosmic_bodies = []
        k = 0
        trajectory = []
        for i in list_cosmic_bodies:
            trajectory.append(i.locus(gravitate_star, time_stop))
            plt.scatter(trajectory[k][0][0], trajectory[k][1][0])
            plt.scatter(0, 0, color='darkgreen')
            plt.plot(trajectory[k][0], trajectory[k][1])
            k += 1
        plt.show()


class CosmicBody3D(CosmicBody):
    @staticmethod
    def create_random_cosmic_body():
        return CosmicBody3D(random.uniform(1e5, 1e7),
                            np.array([(random.uniform(0, 10000)),
                                      (random.uniform(0, 10000)), random.uniform(0, 10000)]),
                            np.array([(random.uniform(-10000, 10000)),
                                      (random.uniform(-10000, 10000)), random.uniform(-10000, 10000)]))

    def destroy(self, gravitate_star: Star):
        if abs(self.position[0]) < gravitate_star.radius and abs(self.position[1]) < gravitate_star.radius and abs(
                self.position[2]) < gravitate_star.radius:
            self.mass = 0
            return 1

    def locus(self, gravitate_star: Star, time_stop):
        position_x = []
        position_y = []
        position_z = []
        time = 0
        while time < time_stop:
            position_x.append(self.position[0])
            position_y.append(self.position[1])
            position_z.append(self.position[2])
            self.step(gravitate_star)
            time += dt
        return [position_x, position_y, position_z]

    @staticmethod
    def graph(gravitate_star: Star, time_stop, list_cosmic_bodies=None):
        if list_cosmic_bodies is None:
            list_cosmic_bodies = []
        k = 0
        plt.figure()
        ax = plt.axes(projection='3d')
        trajectory = []
        for i in list_cosmic_bodies:
            trajectory.append(i.locus(gravitate_star, time_stop))
            ax.plot3D(trajectory[k][0], trajectory[k][1], trajectory[k][2])
            k += 1
        plt.show()
