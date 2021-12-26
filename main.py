import unittest

import numpy as np
import matplotlib.pyplot as plt


class Star:
    def __init__(self, mass: float, radius: float):
        """
        Создание звезды.

        Parameters
        ----------
        mass : float
            Масса звезды.
        radius : float
            Радиус звезды.
        """
        self.mass = mass
        self.radius = radius


class CosmicBody:
    def __init__(self, mass: float, vec_v: np.ndarray, vec_p: np.ndarray):
        """
        Создание звезды.

        Parameters
        ----------
        mass : float
            Масса звезды
        vec_v : np.ndarray
            Вектор скорости тела
        vec_p : np.ndarray
            Радиус-вектор положения тела
        """
        self.mass = mass
        self.vec_v = vec_v
        self.vec_p = vec_p

    def gravitate(self, Star):
        """
        Расчёт силы взаимодействия между звездой и телом.

        Parameters
        ----------
        Star : Star
            звезда.
        """
        return -1 * G * self.mass * Star.mass * self.vec_p / np.power(np.linalg.norm(self.vec_p), 3)

    def destroy(self, Star):
        """
        Уничтожение тела в результате падения на звезду.

        Parameters
        ----------
        Star : Star
            звезда.
        """
        if np.linalg.norm(self.vec_p) < Star.radius:
            self.mass = 0
            return 1

    def move(self, Star):
        """
        Обновление координат и скорости тела через время dt.

        Parameters
        ----------
        Star : Star
            звезда.
        """
        if self.destroy(Star) != 1:
            delta_v = self.gravitate(Star) * dt / self.mass
            self.vec_p = self.vec_p + self.vec_v * dt + delta_v * dt / 2
            self.vec_v = delta_v + self.vec_v

    def kind_of_trajectory(self, Star):
        """
        Определение вида траектории тела.

        Parameters
        ----------
        Star : Star
            звезда.
        """
        energy = self.mass * np.power(np.linalg.norm(self.vec_v), 2) / 2 - G * self.mass * Star.mass / np.linalg.norm(
            self.vec_p)
        if (energy > 1e-6):
            return "hyperbola"
        if (energy < 1e-6):
            return "ellipse"
        return "parabola"

    def trajectory2D(self, Star, start, end):
        """
        Получение 2D координат траектории тела.

        Parameters
        ----------
        Star : Star
            звезда.
        start : float
            начальное время.
        end : float
            конечное время.
        """
        x = []
        y = []
        for i in np.arange(start, end, dt):
            x.append(self.vec_p[0])
            y.append(self.vec_p[1])
            self.move(Star)
        return [x, y]

    def trajectory3D(self, Star, start, end):
        """
        Получение 3D координат траектории тела.

        Parameters
        ----------
        Star : Star
            звезда.
        start : float
            начальное время.
        end : float
            конечное время.
        """
        x = []
        y = []
        z = []
        for i in np.arange(start, end, dt):
            x.append(self.vec_p[0])
            y.append(self.vec_p[1])
            z.append(self.vec_p[2])
            self.move(Star)
        return [x, y, z]

    def graph2D(self, Star, start, end, bodies=None):
        """
        Построение 2D графика траектории тел.

        Parameters
        ----------
        Star : Star
            звезда.
        start : float
            начальное время.
        end : float
            конечное время.
        bodies : list
            список космических тел.
        """
        track = self.trajectory2D(Star, start, end)
        plt.plot(track[0], track[1])
        plt.scatter(0, 0, c='red')
        if bodies != None:
            for body in bodies:
                track = body.trajectory2D(Star, start, end)
                plt.plot(track[0], track[1])
        plt.show()

    def graph3D(self, Star, start, end, bodies=None):
        """
        Построение 3D графика траектории тел.

        Parameters
        ----------
        Star : Star
            звезда.
        start : float
            начальное время.
        end : float
            конечное время.
        bodies : list
            список космических тел.
        """
        track = self.trajectory3D(Star, start, end)
        ax = plt.axes(projection='3d')
        ax.scatter(0, 0, c='red')
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

"Первые 3 планеты солнечной системы"
star = Star(M, R)
mercury = CosmicBody(3.33e23, np.array([0, 49.4e3]), np.array([4.6e10, 0]))
venus = CosmicBody(4.87e24, np.array([0, 35e3]), np.array([1.07e11, 1]))
earth = CosmicBody(m, np.array([0, 30.4e3]), np.array([1.5e11, 0]))
CosmicBody.graph2D(mercury, star, 0, 365 * 24 * 900, [venus, earth])

"Разные тела, вылетающие с орбиты земли"
earth = Star(m, 0)
v_2 = np.sqrt(G * m / 6371e3)
body1 = CosmicBody(1e3, np.array([7.91e3 / np.sqrt(2), 7.91e3 / np.sqrt(2), 0]),
                   np.array([0, 0, 6371e3]))
body2 = CosmicBody(1e3, np.array([v_2, v_2, 0]),
                   np.array([0, 0, 6371e3]))
body3 = CosmicBody(1e3, np.array([16.7e3 / np.sqrt(2), 16.7e3 / np.sqrt(2), 0]),
                   np.array([0, 0, 6371e3]))
CosmicBody.graph3D(body1, earth, 0, 3000, [body2, body3])


if __name__ == '__main__':
    unittest.main()





