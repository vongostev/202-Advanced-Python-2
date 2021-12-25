import numpy as np
import matplotlib.pyplot as plt

dt = 1
G = 6.67e-11


class Star:
    def __init__(self, mass: float, radius: float):
        self.mass = mass
        self.radius = radius


class CosmicBody:
    def __init__(self, mass: float, vec_p: np.ndarray, vec_v: np.ndarray):
        self.mass = mass
        self.vec_v = vec_v
        self.vec_p = vec_p

    def gravitate(self, Star: Star):
        return G * self.mass * Star.mass * self.vec_p / np.power(np.linalg.norm(self.vec_p), 3)

    def destroy(self, Star: Star):
        if np.linalg.norm(self.vec_p) < Star.radius:
            self.mass = 0
            return True

    def move(self, Star: Star):
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
        return x, y

    def trajectory3D(self, Star: Star, start, end):
        x = []
        y = []
        z = []
        for i in np.arange(start, end, dt):
            x.append(self.vec_p[0])
            y.append(self.vec_p[1])
            z.append(self.vec_p[2])
            self.move(Star)
        return x, y, z











