# -*- coding: utf-8 -*-
"""
Created on Sun Dec 26 14:16:23 2021

@author: grego
"""
import numpy as np
import matplotlib.pyplot as plt

"""
Constants and smth
"""
G = 6.67


def Norm(vec):
    """
    

    Parameters
    ----------
    vec : входной вектор.

    Returns
    -------
    Норму вектора.

    """
    return (np.sum(vec ** 2))**(1/2)


class Star():
    """
    Class desribes star. 
    Placed in (0,0) and has 10^3 kg mass, if doesn't mentioned another
    """

    def __init__(self, mass: float = 1e3, xy=(0, 0)):
        self.exist = 1
        self.mass = mass
        self.vec_p = np.array(xy)

    def getMass(self):
        return self.mass

    def getPosition(self):
        return self.vec_p


class CosmicBody():

    """
    Class desribes any CosmicBody. 
    Placed in (0,1000), has 10 kg mass and zero velosity, if doesn't mentioned another
    Mass [kg]
    P[m, m]
    V[m/s]
    """

    def __init__(self, mass: float = 10, vec_p=(0, 1000), vec_v=(0, 0)):
        self.exist = 1
        self.mass = mass
        self.vec_p = np.array(vec_p)
        self.vec_v = np.array(vec_v)

    def destroy(self):
        self.exist = 0
        self.mass = 0
        self.vec_v = (0, 0)

    def distance(self, body):
        return (body.vec_p - self.vec_p)

    def move(self, dt):
        """
        
        
        Parameters
        ----------
        dt : промежуток времени

        Returns
        -------
        Меняет координату как если бы тело двигалось равномерно прямолинейно в течении dt

        """

        self.vec_p = self.vec_p + self.vec_v * dt * self.exist

    def grav(self, dt, *bodies):
        """
        

        Parameters
        ----------
        dt : промежуток времени.
        *bodies : кортеж взаимодействующих тел.
        M : масса конкретного тела.

        Returns
        -------
        Изменяет скорость self по второму закону ньютона для self и *bodies

        """

        for i in range(len(bodies)):
            M = bodies[i].getMass()
            r = self.distance(bodies[i])
            self.vec_v = self.vec_v + G * M * r * dt / \
                (Norm(r) ** 3) * self.exist * bodies[i].exist

    def getMass(self):
        return self.mass

    def getPosition(self):
        return self.vec_p

    def getVelosity(self):
        return self.vec_v


if __name__ == '__main__':

    def test_2DMotion_StarBody():
        """
        Discrete motion. During time dt moves with constant velosity. 
        Momentum changes according to the Gravitational Law.
        Crirical velosity is 3.65. 
        1 ellipse < v_critical
        2 apple < v_central
        3 hyperbola > v_critical
        
        Time [sec]
        Mass [kg]
        P[m]
        V[m/sec]
        """
        velosity = [(3.2, 0), (0.27, 0), (3.9, 0)]
        for i in (0, 1, 2):

            Sun = Star()
            Aster = CosmicBody(0.1, (0, 1000), velosity[i])
            x = [Aster.getPosition()[0]]
            y = [Aster.getPosition()[1]]
            dt = 0.2
            n = 28000
            t = 0
            while t < dt*n:
                Aster.grav(dt, Sun)
                Aster.move(dt)
                x.append(Aster.getPosition()[0])
                y.append(Aster.getPosition()[1])

                if (Norm(Aster.getPosition()) <= 45):
                    Aster.destroy()

                t = t + dt
            plt.plot(x, y, '--', linewidth=1)

    def test_2DMotion_2Bodies():
        velosity = [(2, 0), (-3, 0)]
        Ast1 = CosmicBody(12, (0, 0), velosity[0])
        Ast2 = CosmicBody(16, (0, 5), velosity[1])
        x1 = [Ast1.getPosition()[0]]
        y1 = [Ast1.getPosition()[1]]
        x2 = [Ast2.getPosition()[0]]
        y2 = [Ast2.getPosition()[1]]

        dt = 0.01
        n = 1000
        t = 0
        while t < dt*n:
            Ast1.grav(dt, Ast2)
            Ast2.grav(dt, Ast1)
            Ast1.move(dt)
            Ast2.move(dt)
            x1.append(Ast1.getPosition()[0])
            y1.append(Ast1.getPosition()[1])
            x2.append(Ast2.getPosition()[0])
            y2.append(Ast2.getPosition()[1])

            if (Norm(Ast1.distance(Ast2)) <= 0.2):
                Ast1.destroy()
                Ast2.destroy()

            t = t + dt
        plt.show()
        plt.plot(x1, y1, '--', linewidth=1)
        plt.plot(x2, y2, '--', linewidth=1)

    def test_destruction():
        velosity = [(0.1, 0), (-0.1, 0)]
        Ast1 = CosmicBody(12, (0, 0), velosity[0])
        Ast2 = CosmicBody(5, (0, 8), velosity[1])
        x1 = [Ast1.getPosition()[0]]
        y1 = [Ast1.getPosition()[1]]
        x2 = [Ast2.getPosition()[0]]
        y2 = [Ast2.getPosition()[1]]

        dt = 0.001
        n = 25036
        t = 0
        while t < dt*n:
            Ast1.grav(dt, Ast2)
            Ast2.grav(dt, Ast1)
            Ast1.move(dt)
            Ast2.move(dt)
            x1.append(Ast1.getPosition()[0])
            y1.append(Ast1.getPosition()[1])
            x2.append(Ast2.getPosition()[0])
            y2.append(Ast2.getPosition()[1])

            if (Norm(Ast1.distance(Ast2)) <= 0.05):
                Ast1.destroy()
                Ast2.destroy()

            t = t + dt
        plt.show()
        plt.plot(x1, y1, 'o--', linewidth=1)
        plt.plot(x2, y2, 'o--', linewidth=1)
        assert Ast1.exist == 0
        print("Test Destruction is Ok")

    def test_Star():
        Sun = Star(35, (3, 9))
        assert Sun.getMass() == 35
        # print (Sun.getPosition())
        # assert Sun.getPosition() == (3,9)
        print("test_Star is Ok")

    def test_Body():
        Tesla_Roadster = CosmicBody(1814, (0, 6400), (7.91, 0))
        assert Tesla_Roadster.getMass() == 1814
        # print (Tesla_Roadster.getPosition())
        # print (Tesla_Roadster.getVelosity())
        print("Test Tesla Motorspots is Ok")

    # test_2DMotion_StarBody()
    # test_2DMotion_2Bodies()
    # test_destruction()
    test_Star()
    test_Body()
