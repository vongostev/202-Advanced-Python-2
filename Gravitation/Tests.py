import unittest
import numpy as np
from cosmicbody import CosmicBody2D, CosmicBody3D, CosmicBody
from Star import Star



class Tests(unittest.TestCase):
    def test_orbit_type(self):
        G = 6.674e-11
        m = 5.9722e24
        v_2 = np.sqrt(G * m / 6371e3)
        earth = Star(m)
        body1 = CosmicBody3D(1e3, np.array([7.91e3 / np.sqrt(2), 7.91e3 / np.sqrt(2), 0]),
                             np.array([0, 0, 6371e3]))
        body2 = CosmicBody3D(1e3, np.array([v_2, v_2, 0]),
                             np.array([0, 0, 6371e3]))
        body3 = CosmicBody3D(1e3, np.array([16.7e3 / np.sqrt(2), 16.7e3 / np.sqrt(2), 0]),
                             np.array([0, 0, 6371e3]))
        self.assertEqual(body1.orbit_type(earth), 'ellipse')
        self.assertEqual(body2.orbit_type(earth), 'parabola')
        self.assertEqual(body3.orbit_type(earth), 'hyperbole')

    def test_destroy(self):
        star = Star(1.9e30, 6e4)
        t_end = 10
        body = CosmicBody2D(5.97e24, np.array([-5000, 0]), np.array([146e2, 0]))
        body.locus(star, t_end)
        self.assertEqual(body.destroy(star), 1)

