import numpy as np


class Body:
    def __init__(self,
                 position: np.ndarray,
                 velocity: np.ndarray,
                 mass: float,
                 radius: float,
                 log_path: str,
                 tick: int = 0):
        self.position = position
        self.velocity = velocity
        self.mass = mass
        self.radius = radius
        self.log_path = log_path
        self.does_exist = 1
        if tick > 0:
            log_file = open(log_path, 'a')
            for i in range(tick):
                log_file.write('0 0 0 [0 0 0]')
            log_file.close()

    def move(self,
             acceleration: np.ndarray,
             tick_length: float):
        log_file = open(self.log_path, 'a')
        log_file.write(
            f'{self.does_exist} {self.mass} {self.radius} {self.position}')
        log_file.close()
        self.position += self.velocity*tick_length + 0.5*acceleration*tick_length
        self.velocity += acceleration*tick_length

    def destroy(self):
        self.mass = 0
        self.radius = 0
        self.does_exist = 0
        self.position = np.array([0, 0, 0])
        self.velocity = np.array([0, 0, 0])

    def merge(self,
              other_body):
        self.position = (self.position * self.mass + other_body.position *
                         other_body.mass)/(self.mass + other_body.mass)
        self.velocity = (self.velocity * self.mass + other_body.velocity *
                         other_body.mass)/(self.mass + other_body.mass)
        self.mass += other_body.mass
        self.radius = (self.radius**3 + other_body.radius**3)**(1/3)
        other_body.destroy()
