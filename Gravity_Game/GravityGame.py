import numpy as np
import numba as nb
from collections import deque

G = 0.1  # Гравитационная постоянная


class Body:
    def __init__(self,
                 position: np.ndarray,
                 velocity: np.ndarray,
                 mass: float,
                 radius: float,
                 tick_length: float,
                 log_path: str):
        self.position = position
        self.velocity = velocity
        self.mass = mass
        self.radius = radius
        self.tick_length = tick_length
        self.log_path = log_path
        self.does_exist = 1
        # if tick > 0:
        #     log_file = open(log_path, 'a')
        #     for i in range(tick):
        #         log_file.write(f'0 0.0 0.0 {np.zeros_like(position)}')
        #     log_file.close()

    def move(self,
             acceleration: np.ndarray,):
        log_file = open(self.log_path, 'a')
        log_file.write(
            f'{self.does_exist} {self.mass} {self.radius} {self.position}')
        log_file.close()
        self.position += self.velocity * self.tick_length + \
            0.5 * acceleration * self.tick_length
        self.velocity += acceleration * self.tick_length

    def destroy(self):
        self.mass = 0
        self.radius = 0
        self.does_exist = 0
        self.position = np.zeros_like(self.position)
        self.velocity = np.zeros_like(self.position)

    def merge(self,
              other_body):
        self.position = (self.position * self.mass + other_body.position *
                         other_body.mass) / (self.mass + other_body.mass)
        self.velocity = (self.velocity * self.mass + other_body.velocity *
                         other_body.mass) / (self.mass + other_body.mass)
        self.mass += other_body.mass
        self.radius = (self.radius**3 + other_body.radius**3)**(1/3)
        other_body.destroy()


@nb.njit('float64[:](float64, float64[:], float64[:])', cache=True, nogil=False, fastmath=True, parallel=True)
def calculate_acceleration(attractor_mass: float,
                           attractor_position: np.ndarray,
                           body_position: np.ndarray):
    distance = attractor_position - body_position
    return G * attractor_mass * distance / np.linalg.norm(distance)**3


def gravitate(bodies: list):
    accelerations = deque()
    for body in bodies:
        accelerations.append(0.)
        if body.does_exist:
            for attractor in [a for a in bodies if a != body]:
                accelerations[-1] += calculate_acceleration(attractor.mass,
                                                            attractor.position,
                                                            body.position)
    for b in range(len(bodies)):
        bodies[b].move(accelerations.popleft())
        for c in range(b):
            if bodies[c].radius + bodies[b].radius <= np.linalg.norm(bodies[b].position - bodies[c].position) and bodies[b].does_exist*bodies[c].does_exist:
                bodies[b].merge(bodies[c])