import numpy as np
import numba as nb
from collections import deque
from time import gmtime, strftime, time
from os import mkdir, path
import unittest
import matplotlib.pyplot as plt
from tqdm import tqdm


def read_trajectory(file_path: str) -> list:
    """
    Генератор информации о траектории из получаемого файла. Конвертирует
    траекторию по шаблону в список характеристик тела.

    Parameters
    ----------
    file_path : str
        Путь к файлу с траекторией.

    Yields
    ------
    list
        Свойства тела в конкретный тик в формате [does_exist, radius, position].

    """
    with open(file_path) as trajectory:
        for moment in trajectory:
            moment = deque(moment.split(' ['))
            moment[1] = np.array([float(a)
                                 for a in moment[1].rstrip(']\n').split()])
            moment.appendleft(int(moment[0].split()[0]))
            moment[1] = float(moment[1].split()[1])
            yield list(moment)


# def update_plot(num, trajectories, bodies_images):
#     for body_image, trajectory in zip(bodies_images, trajectories):
#         if trajectory[num][0]:
#             body_image.set_data(trajectory[num][2][:2])
#             body_image.set_3d_properties(trajectory[num][2][2])
#     return bodies_images


@nb.njit('float64[:](float64, float64[:], float64[:], float64)', cache=True,
         nogil=False, fastmath=True, parallel=True)
def calculate_acceleration(attractor_mass: float,
                           attractor_position: np.ndarray,
                           body_position: np.ndarray,
                           G: float) -> np.ndarray:
    """
    Вычисляет и возвращает ускорение, вызванное притягивающим телом 
    (аттрактором) у притягиваемого тела.

    Parameters
    ----------
    attractor_mass : float
        Масса аттрактора.
    attractor_position : np.ndarray
        Координаты аттрактора.
    body_position : np.ndarray
        Положение притягиваемого тела.
    G : float
        Гравитационная постоянная

    Returns
    -------
    acceleration : np.ndarray
        Ускорение притягиваемого тела, вызванное аттрактором.

    """
    distance = attractor_position - body_position
    return attractor_mass * G * distance / np.linalg.norm(distance)**3


@nb.njit('float64[:](float64[:], float64[:], float64)', cache=True,
         nogil=False, fastmath=True, parallel=True)
def calculate_position_change(velocity: np.ndarray,
                              acceleration: np.ndarray,
                              tick_length: float) -> np.ndarray:
    """
    Вычисляет изменение координат тела по просшествии одного тика

    Parameters
    ----------
    velocity : np.ndarray
        Скорость движения.
    acceleration : np.ndarray
        Ускорение.
    tick_length : float
        Длительность тика.

    Returns
    -------
    np.ndarray
        Изменение положения тела.

    """
    return velocity * tick_length + 0.5 * acceleration * tick_length**2


@nb.njit('float64[:](float64[:], float64)', cache=True,
         nogil=False, fastmath=True, parallel=True)
def calculate_velocity_change(acceleration: np.ndarray,
                              tick_length: float) -> np.ndarray:
    """
    Вычисляет изменение скорости по просшествии одного тика

    Parameters
    ----------
    acceleration : np.ndarray
        Ускорение.
    tick_length : float
        Длительность тика.

    Returns
    -------
    np.ndarray
        Изменение скорости.

    """
    return acceleration*tick_length


class Gravitation:
    def __init__(self,
                 G: float = 1,
                 tick_length: float = 1e-2,
                 log_dir: str = 'GravityGame_logs'):
        """
        Создание движка гравитации с конкретными параметрами.

        Parameters
        ----------
        G : float, optional
            Гравитационная постоянная. The default is 1.
        tick_length : float, optional
            Длина тика - "кванта" времени в секундах. The default is 1e-3.
        log_dir : str, optional
            Название директории с траекториями. Директория будет расположена
            по тому же адресу, что и скрипт. The default is 'GravityGame_logs'.

        Returns
        -------
        None.

        """
        self.G = G
        self.tick_length = tick_length
        self.log_dir = log_dir
        if not path.exists(log_dir):
            mkdir(log_dir)

    def gravitate(self):
        """
        Выполняет логику одного тика движения тел. Вычисляются ускорения для 
        всех тел (кроме несуществующих), затем тела перемещаются под действием
        вычисленных ускорений. После каждого перемещения производится проверка
        столкновения тела с уже перемещёнными ранее (итерация по индексу 
        используется, чтобы не затрагивать ещё не перемещённые тела). Если 
        столкновение произошло, происходит слияние столкнувшихся тел.

        Parameters
        ----------
        None.

        Returns
        -------
        None.

        """
        # Вычисление ускорений. Порядковые номера ускорения и соответствующего тела совпадают
        accelerations = deque()
        for body in self.bodies:
            accelerations.append(0.)
            if body.does_exist:
                for attractor in self.bodies:
                    if attractor != body:
                        accelerations[-1] += calculate_acceleration(attractor.mass,
                                                                    attractor.position,
                                                                    body.position,
                                                                    self.G)

        for b in range(self.bodies_number):
            # Перемещение тел
            self.bodies[b].move(accelerations.popleft(), self.tick_length)
            for c in range(b):
                # Проверка на столкновение (не учитываются несуществующие тела)
                if (self.bodies[b].does_exist * self.bodies[c].does_exist
                    and self.bodies[c].radius + self.bodies[b].radius
                        >= np.linalg.norm(self.bodies[b].position - self.bodies[c].position)):
                    self.bodies[b].merge(self.bodies[c])

    def create_bodies(self,
                      properties: list):
        """
        Создаёт внутри объекта гравитации список тел с указанными свойствами.
        Формат списка свойств:
            properties[i][0] : np.ndarray - начальное положение тела (массив из float64)
            properties[i][1] : np.ndarray - начальная скорость тела (массив из float64)
            properties[i][2] : float - масса тела
            properties[i][3] : float - радиус тела

        Parameters
        ----------
        properties : list
            Список свойств тел. Формат указан в описании функции!

        Returns
        -------
        None.

        """
        self.bodies = []
        self.paths = []
        for body_props in properties:
            self.paths.append(
                f'{self.log_dir}/m{body_props[2]}_r{body_props[3]}_{strftime("%d_%m_%y_%H_%M_%S", gmtime())}_{time()%1e-3*1e9:.0f}.txt')
            self.bodies.append(Body(body_props[0],
                                    body_props[1],
                                    body_props[2],
                                    body_props[3],
                                    self.paths[-1]))
        self.bodies_number = len(self.bodies)

    def simulate_trajectories(self,
                              simulation_time: float):
        """
        Запускает симулицию в течение указанного времени (если укладывается 
        нецелое число тиков, округление вниз)

        Parameters
        ----------
        simulation_time : float
            Время симуляции.

        Returns
        -------
        None.

        """
        self.simulation_time = simulation_time
        self.N = int(simulation_time / self.tick_length)
        for t in range(self.N):
            self.gravitate()

    def animate_trajectories(self):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        fig.show()

        bodies_trajectories_generators = [
            read_trajectory(log_path) for log_path in self.paths]

        bodies_trajectories = [[] for i in self.bodies]

        for generator, trajectory in zip(bodies_trajectories_generators, bodies_trajectories):
            for moment in generator:
                trajectory.append(moment)

        for t in np.arange(self.N):
            ax.set_xlim3d(-10, 10)
            ax.set_ylim3d(-10, 10)
            ax.set_zlim3d(-8, 8)
            for trajectory in bodies_trajectories:
                if trajectory[t][0]:
                    ax.plot(trajectory[t][2][0],
                            trajectory[t][2][1],
                            trajectory[t][2][2])
            fig.canvas.draw()


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
        self.position = position.astype(float)
        self.velocity = velocity.astype(float)
        self.mass = float(mass)
        self.radius = float(radius)
        self.log_path = log_path
        f = open(log_path, 'w')
        f.close()
        self.does_exist = 1

    def move(self,
             acceleration: np.ndarray,
             tick_length: float):
        """
        Записывает в файл траектории строку с характеристиками тела в формате:
        'Флаг_существования Радиус [Координаты через пробел]'
        Затем изменяет координаты и скорость тела, как если бы тело двигалось
        равноускорено с заданным ускорением в течение одного тика.

        Parameters
        ----------
        acceleration : np.ndarray
            Ускорение, действующее на тело.
        tick_length : float
            Длина тика в секундах.

        Returns
        -------
        None.

        """
        with open(self.log_path, 'a') as log_file:
            log_file.write(
                f'{self.does_exist} {self.radius} {self.position}\n')

        self.position += calculate_position_change(
            self.velocity, acceleration, tick_length)
        self.velocity += calculate_velocity_change(acceleration, tick_length)

    def destroy(self):
        """
        Уничтожает тело: меняет значение флага does_exist на 0 и обнуляет массу,
        радиус и скорость тела. Траектория при этом не удаляется.

        Returns
        -------
        None.

        """
        self.mass = 0
        self.radius = 0
        self.does_exist = 0
        self.velocity = np.zeros_like(self.velocity)

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


### Тесты ###

class TestBodyMethods(unittest.TestCase):
    def test_init(self):
        body = Body(np.arange(3), np.arange(3), 3, 2, 'test_body.txt')
        self.assertEqual(body.position.tolist(), [0, 1, 2])
        self.assertEqual(body.velocity.tolist(), [0, 1, 2])
        self.assertEqual(body.mass, 3)
        self.assertEqual(body.radius, 2)
        self.assertEqual(body.log_path, 'test_body.txt')
        self.assertEqual(body.does_exist, 1)

    def test_move(self):
        body = Body(np.arange(3),
                    np.arange(3)+1,
                    3, 2, 'test_body.txt')
        body.move(np.zeros(3), 1.)
        self.assertEqual(body.position.tolist(), [1, 3, 5])
        self.assertEqual(body.velocity.tolist(), [1, 2, 3])
        body.move(np.array([-1., 0., 1.]), 1.)
        self.assertEqual(body.position.tolist(), [1.5, 5, 8.5])
        self.assertEqual(body.velocity.tolist(), [0, 2, 4])
        with open('test_body.txt', 'r') as log_file:
            self.assertEqual(log_file.readline(), '1 2.0 [0. 1. 2.]\n')
            self.assertEqual(log_file.readline(), '1 2.0 [1. 3. 5.]\n')

    def test_destroy(self):
        body = Body(np.arange(3),
                    np.arange(3)+1,
                    3, 2, 'test_body.txt')
        body.destroy()
        self.assertEqual(body.position.tolist(), [0, 1, 2])
        self.assertEqual(body.velocity.tolist(), [0, 0, 0])
        self.assertEqual(body.mass, 0)
        self.assertEqual(body.radius, 0)

    def test_merge(self):
        body_1 = Body(np.arange(3),
                      np.arange(3)+1,
                      3, 2, 'test_body_1.txt')
        body_2 = Body(np.array([1, 2, 3]),
                      np.array([-1, -2, -3]),
                      3, 2, 'test_body_2.txt')
        body_1.merge(body_2)
        # body_2 is destroyed
        self.assertEqual(body_2.position.tolist(), [1, 2, 3])
        self.assertEqual(body_2.velocity.tolist(), [0, 0, 0])
        self.assertEqual(body_2.mass, 0)
        self.assertEqual(body_2.radius, 0)
        # body_1 is merged
        self.assertEqual(body_1.position.tolist(), [0.5, 1.5, 2.5])
        self.assertEqual(body_1.velocity.tolist(), [0, 0, 0])
        self.assertEqual(body_1.mass, 6)
        self.assertAlmostEqual(body_1.radius, 2.52, delta=0.001)


class TestFunctions(unittest.TestCase):
    def test_read_trajectory(self):
        body = Body(np.arange(3),
                    np.arange(3)+1,
                    3, 2, 'test_body.txt')
        body.move(np.zeros(3), 1.)
        body.move(np.array([-1., 0., 1.]), 1.)
        trajectory = read_trajectory('test_body.txt')
        positions = [[0, 1, 2], [1, 3, 5]]
        for i, moment in enumerate(trajectory):
            self.assertEqual(moment[0:2], [1, 2])
            self.assertEqual(moment[2].tolist(), positions[i])

    def test_calculate_acceleration(self):
        self.assertTrue(np.allclose(
            calculate_acceleration(7., np.array([1, 1, 1]).astype(float),
                                   np.array([3, 4, 7]).astype(float),
                                   0.49),  np.array([-0.02, -0.03, -0.06])))


if __name__ == '__main__':
    #unittest.main()
    g = Gravitation(tick_length=0.1)
    g.create_bodies([[np.arange(3), np.zeros(3), 3, 1],
                     [np.array([3, 5, 7]), np.array([-1, -1, 1]), 2, 2]])
    g.simulate_trajectories(2.)
    g.animate_trajectories()
