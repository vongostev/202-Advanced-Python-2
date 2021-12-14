import numpy as np
import numba as nb
from collections import deque

G = 0.1  # Гравитационная постоянная
tick_length = 1e-2  # Длительность тика - "кванта времени"


class Body:
    def __init__(self,
                 position: np.ndarray,
                 velocity: np.ndarray,
                 mass: float,
                 radius: float,
                 log_path: str):
        """
        Конструктор тела, создающий тело с заданными координатой, скоростью, 
        массой и радиусом. Телу передаётся путь к файлу, куда сохраняется 
        траектория. Файл будет очищен при создании тела.

        Parameters
        ----------
        position : np.ndarray
            Координаты тела в пространстве.
        velocity : np.ndarray
            Скорость тела.
        mass : float
            Масса тела.
        radius : float
            Радиус тела.
        log_path : str
            Путь к файлу, в который будет сохранена траектория тела.

        Returns
        -------
        None.

        """
        self.position = position
        self.velocity = velocity
        self.mass = mass
        self.radius = radius
        self.log_path = log_path
        f = open(log_path, 'w')
        f.close()
        self.does_exist = 1
        # if tick > 0:
        #     log_file = open(log_path, 'a')
        #     for i in range(tick):
        #         log_file.write(f'0 0.0 0.0 {np.zeros_like(position)}')
        #     log_file.close()

    def move(self,
             acceleration: np.ndarray,):
        """
        Записывает в файл траектории строку с характеристиками тела в формате:
        'Флаг_существования Радиус [Координаты через пробел]'
        Затем изменяет координаты и скорость тела, как если бы тело двигалось
        равноускорено с заданным ускорением в течение одного тика.

        Parameters
        ----------
        acceleration : np.ndarray
            Ускорение, действующее на тело.

        Returns
        -------
        None.

        """
        log_file = open(self.log_path, 'a')
        log_file.write(
            f'{self.does_exist} {self.radius} {self.position}')
        log_file.close()
        self.position += self.velocity * tick_length \
            + 0.5 * acceleration * tick_length**2
        self.velocity += acceleration * tick_length

    def destroy(self):
        """
        Уничтожает тело: меняет значение флага does_exist на 0 и обнуляет все
        характеристики тела. Траектория при этом не удаляется и записывается далее.

        Returns
        -------
        None.

        """
        self.mass = 0
        self.radius = 0
        self.does_exist = 0
        self.position = np.zeros_like(self.position)
        self.velocity = np.zeros_like(self.position)

    def merge(self,
              other_body):
        """
        Выполняет слияние указанных тел. За координату и скорость принимаются 
        координата и скорость центра масс двух тел. Радиусы тел складываются,
        как если бы тела были однородными шарами, и их объём сохранился при
        слиянии. Новое тело занимает место объекта, из которого вызван метод,
        другое тело уничтожается.

        Parameters
        ----------
        other_body : Body
            Тело, с которым происходит слияние.

        Returns
        -------
        None.

        """
        self.position = (self.position * self.mass + other_body.position
                         * other_body.mass) / (self.mass + other_body.mass)
        self.velocity = (self.velocity * self.mass + other_body.velocity
                         * other_body.mass) / (self.mass + other_body.mass)
        self.mass += other_body.mass
        self.radius = (self.radius**3 + other_body.radius**3)**(1/3)
        other_body.destroy()


@nb.njit('float64[:](float64, float64[:], float64[:])', cache=True, nogil=False,
         fastmath=True, parallel=True)
def calculate_acceleration(attractor_mass: float,
                           attractor_position: np.ndarray,
                           body_position: np.ndarray) -> np.ndarray:
    """
    Вычисляет и возвращает ускорение, вызванное притягивающим телом (аттрактором)
    у притягиваемого тела.

    Parameters
    ----------
    attractor_mass : float
        Масса аттрактора.
    attractor_position : np.ndarray
        Координаты аттрактора.
    body_position : np.ndarray
        Положение притягиваемого тела.

    Returns
    -------
    acceleration : np.ndarray
        Ускорение притягиваемого тела, вызванное аттрактором.

    """
    distance = attractor_position - body_position
    return attractor_mass * G * distance / np.linalg.norm(distance)**3


def gravitate(bodies: list):
    """
    Выполняет логику одного тика движения тел. Вычисляются ускорения для всех
    тел из списка (кроме несуществующих), затем тела перемещаются под действием
    вычисленных ускорений. После каждого перемещения производится проверка
    столкновения тела с уже перемещёнными ранее (итерация по индексу используется,
    чтобы не затрагивать ещё не перемещённые тела). Если столкновение произошло,
    происходит слияние столкнувшихся тел.

    Parameters
    ----------
    bodies : list
        Список тел, участвующих в гравитационном взаимодействии.

    Returns
    -------
    None.

    """
    # Вычисление ускорений. Порядковые номера ускорения и соответствующего тела совпадают
    accelerations = deque()
    for body in bodies:
        accelerations.append(0.)
        if body.does_exist:
            for attractor in [a for a in bodies if a != body]:
                accelerations[-1] += calculate_acceleration(attractor.mass,
                                                            attractor.position,
                                                            body.position)

    for b in range(len(bodies)):
        # Перемещение тел
        bodies[b].move(accelerations.popleft())
        for c in range(b):
            # Проверка на столкновение (не учитываются несуществующие тела)
            if (bodies[b].does_exist*bodies[c].does_exist
                and bodies[c].radius + bodies[b].radius
                    <= np.linalg.norm(bodies[b].position - bodies[c].position)):
                bodies[b].merge(bodies[c])
