import unittest

import numpy as np
import matplotlib.pyplot as plt

G = 6.67e-11
dt = 0.001

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
        star = Star(1.9e30, 7e4)
        body = CosmicBody(5.97e24, np.array([-5000, 0]), np.array([146e2, 0]))
        G = 6.674e-11
        body.graph2D(star, 0, 10)
        self.assertEqual(body.destroy(star), 1)



if __name__ == '__main__':
    unittest.main()
#G = 6.67e-11
#dt = 1000
#t_beg = 0
#t_end = 365.25 * 24 * 3600
#a1 = CosmicBody(3.33e23, np.array([29800. / np.sqrt(2), -29800. / np.sqrt(
#        2), 0]), np.array([146e9 / np.sqrt(2), 146e9 / np.sqrt(2), 0]))
#star = Star(1.9e30, 0)
#a1.graph2D(star, t_beg, t_end)
#a1.graph3D(star, t_beg, t_end)

#fish1 = CosmicBody(3.33e23, np.array([29800. / np.sqrt(2), -29800. / np.sqrt(
#        2)]), np.array([146e9 / np.sqrt(2), 146e9 / np.sqrt(2)]))
#fish2 = CosmicBody(3.33e23, np.array([29800. / np.sqrt(2), 29800. / np.sqrt(
#        2)]), np.array([146e9 / np.sqrt(2), -(146e9 / np.sqrt(2))]))
#t_beg = 0
#t_end = 100 * 24 * 3600
#fish1.graph2D(star, t_beg, t_end, [fish2])





