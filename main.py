import unittest

import numpy as np
import matplotlib.pyplot as plt


class Star:
    def __init__(self, mass: float, radius: float):
        self.mass = mass
        self.radius = radius


class CosmicBody:
    def __init__(self, mass: float, vec_v: np.ndarray, vec_p: np.ndarray):
        self.mass = mass
        self.vec_v = vec_v
        self.vec_p = vec_p

    def gravitate(self, Star: Star):
        return -1 * G * self.mass * Star.mass * self.vec_p / np.power(np.linalg.norm(self.vec_p), 3)

    def destroy(self, Star: Star):
        if np.linalg.norm(self.vec_p) < Star.radius:
            self.mass = 0
            return 1

    def move(self, Star: Star):
        if self.destroy(Star) != 1:
            delta_v = self.gravitate(Star) * dt / self.mass
            self.vec_p = self.vec_p + self.vec_v * dt + delta_v * dt / 2
            self.vec_v = delta_v + self.vec_v

    def kind_of_trajectory(self, Star: Star):
        energy = self.mass * np.power(np.linalg.norm(self.vec_v), 2) / 2 - G * self.mass * Star.mass / np.linalg.norm(
            self.vec_p)
        if (energy > 1e-6):
            return "hyperbola"
        if (energy < 1e-6):
            return "ellipse"
        return "parabola"

    def trajectory2D(self, Star: Star, start, end):
        x = []
        y = []
        for i in np.arange(start, end, dt):
            x.append(self.vec_p[0])
            y.append(self.vec_p[1])
            self.move(Star)
        return [x, y]

    def trajectory3D(self, Star: Star, start, end):
        x = []
        y = []
        z = []
        for i in np.arange(start, end, dt):
            x.append(self.vec_p[0])
            y.append(self.vec_p[1])
            z.append(self.vec_p[2])
            self.move(Star)
        return [x, y, z]

    def graph2D(self, Star: Star, start, end, bodies=None):
        track = self.trajectory2D(Star, start, end)
        plt.plot(track[0], track[1])
        if bodies != None:
            for body in bodies:
                track = body.trajectory2D(Star, start, end)
                plt.plot(track[0], track[1])
        plt.show()

    def graph3D(self, Star: Star, start, end, bodies=None):
        track = self.trajectory3D(Star, start, end)
        ax = plt.axes(projection='3d')
        ax.plot3D(track[0], track[1], track[2])
        if bodies != None:
            for body in bodies:
                track = body.trajectory3D(Star, start, end)
                ax.plot3D(track[0], track[1], track[2])
        plt.show()


class Tests(unittest.TestCase):
    def test_orbit_type(self):
        star = Star(1.9e30, 7e8)
        body1 = CosmicBody(3.33e23, np.array([29800. / np.sqrt(2), -29800. / np.sqrt(
            2)]), np.array([146e9 / np.sqrt(2), 146e9 / np.sqrt(2)]))
        body2 = CosmicBody(3.33e23, np.array([2980000. / np.sqrt(2), -29800. / np.sqrt(
            2)]), np.array([146e9 / np.sqrt(2), 146e9 / np.sqrt(2)]))
        self.assertEqual(body1.kind_of_trajectory(star), 'ellipse')
        self.assertEqual(body2.kind_of_trajectory(star), 'hyperbola')

    def test_destroy(self):
        star = Star(1.9e30, 7e8)
        body = CosmicBody(3.33e23, np.array([298. / np.sqrt(2), -29800. / np.sqrt(
            2)]), np.array([146e3 / np.sqrt(2), 146e8 / np.sqrt(2)]))
        body.graph2D(star, 0, 10000000)
        self.assertEqual(body.destroy(star), 1)

G = 6.67e-11
dt = 100
m = 5.97e24
r = 6.3e6
M = 1.989e30
R = 7e8

"Земля и Солнце"
star = Star(M, R)
earth = CosmicBody(m, np.array([0, 30.4e3]), np.array([1.47e11, 0]))
print(earth.kind_of_trajectory(star))
earth.graph2D(star, 0, 365 * 24 * 3600)

if __name__ == '__main__':
    unittest.main()






