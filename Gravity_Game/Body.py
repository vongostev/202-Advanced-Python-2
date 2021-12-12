import numpy as np

class Body:
    def __init__(
            self,
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
        if tick > 0:
            log_file = open(log_path, 'a')
            for i in range(tick):
                log_file.write(f'{i} False 0 0 0 0 0')
            log_file.close()