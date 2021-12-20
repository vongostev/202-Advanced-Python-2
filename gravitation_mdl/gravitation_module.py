import matplotlib.pyplot as plt 
import numpy as np 
from numpy.linalg import norm 
import random as rd
import unittest
import numba as nb

# Магическая строка для корректного отображения анимации
#matplotlib notebook

G = 4 * np.pi**2 # в системе единиц "годы, массы Солнца, а.е."
dt = 0.005
crash_distance = 0.001
#total_
total_time = 10 #в годах
dim = 3 #размерность задачи
destroy_flag = 0 
vmax = 3 #максимальная скорость моделируемых тел
radius = 2 #максимальное расстояние от начала координат для моделируемых тел
size = 5 #radius * 3
max_dist = 10

fig = plt.figure(figsize=(10, 10))
# Добавляем на нее объект Axes3D для отображения графиков в 3D
ax = fig.add_subplot(111, projection='3d')
time = 0
        
#единственная функция, где на мой взгляд стоит использовать numba, в остальных, 
#чтобы к атрибуту класса вектор, кажется, нет смысла писать отдельную функцию и ускорять ее
@nb.njit('float64[:](float64,float64[:],float64[:])', cache=False,
         nogil=False, fastmath=True, parallel=False) 
def calc_acceleration(mass, r_0, r):
    """
    Вычисляет и возвращает ускорение притягиваемого тела, 
    вызванное действием силы притягивающего тела
    Parameters
    ----------
    mass : float
        Масса притягивающего тела
    r_0 : np.ndarray
        Радиус-вектор притягивающего тела
    r : np.ndarray
        Радиус-вектор притягиваемого тела
    Returns
    -------
    calc_acceleration : np.ndarray
        Ускорение притягиваемого тела, вызванное действием силы притягивающего тела
    """
    return G * mass * (r_0 - r) / (norm(r - r_0)**3)

class Star():
    def __init__(self, mass = None, vec_r = None, vec_v = None): #, dv=0):
        """
        Создает объект класса Звезда
        Parameters
        ----------
        mass : float
            Масса звезды
        vec_r : np.ndarray
            Радиус-вектор звезды
        vec_v : np.ndarray
            Вектор скорости звезды
        Returns
        -------
        None
        """
        self.mass = mass if mass is not None else rd.uniform(1,10)
        self.vec_r = np.array(vec_r) if vec_r is not None else np.array([rd.uniform(-radius,radius) for i in range(dim)] )#, size=dim) )
        self.vec_v = np.array(vec_v) if vec_v is not None else np.array([rd.uniform(-vmax,vmax) for i in range(dim)] )
         #self.dv = dv #np.array(dv)
    
    def show(self):
        """
        Выводит в строку информацию об объекте класса
        Parameters
        ----------
        None
        Returns
        -------
        None
        """
        print(f'm = {self.mass}, vec_r = {self.vec_r}, vec_v = {self.vec_v}')

    def is_body(self, body):
        """
        Проверяет, совпадают ли объекты
        Parameters
        ----------
        body : Star 
            Звезда, с которым сравнивается объект класса
        Returns
        -------
        is_body : bool
            1, если совпадают; 0, если нет
        """
        res = 1
        for i in range(dim):
            res = res & (self.vec_r[i] == body.vec_r[i]) & (self.vec_v[i] == body.vec_v[i]) 
            #print(self.vec_r[i], body.vec_r[i], self.vec_v[i], body.vec_r[i])
        res = (self.mass == body.mass) & res #& (self.dv ==  body.dv)
        #print(res)
        return res

    def destroy(self, body):
        """
        Осуществляет столковение тела со звездой в предположении абсолютно неупругого удара, изменяя массу и скорость звезды
        Parameters
        ----------
        body : Star
            тело, налетающее на звезду
        Returns
        -------
        None
        """
        if (self.mass / body.mass < 1e4): #если масса body меньше, то влияние на звезду оказывается меньше погрешности метода
            self.vec_r = (self.vec_r + body.vec_r)/2
            self.vec_v = (self.mass * self.vec_v + body.mass * body.vec_r) / (self.mass + body.mass)       
            self.mass = self.mass + body.mass

    def move_system(self, bodies, time):
        """
        Моделирует движение на заданном промежутке времени и возвращает список столкнувшихся тел 
        Parameters
        ----------
        bodies : list of CosmicBodies
            список тел, движение которых моделируется в поле звезды и друг друга 
        Returns
        -------
        move_system : list of CosmicBodies
            список столкнувшихся тел на заданном промежутке времени 
        """
        if ((norm(self.vec_v) > 0)or(norm(self.vec_r)>0)): #проверка того, перешли ли мы уже в СО звезды
            for body in bodies:
                body.vec_v = body.vec_v - self.vec_v # переход в СО звезды
                body.vec_r = body.vec_r - self.vec_r
            self.vec_v = self.vec_v * 0
            self.vec_r = self.vec_r * 0
            
        self_destroy_flag = 0
        destroy_flag = 0
        
        #найдем ускорение самой звезды 
        self_dv = self.vec_v * 0  #для соблюдения размерности
        
        for obj in bodies:
            if (norm(obj.vec_r-self.vec_r) > crash_distance):
                if (self.mass / obj.mass < 1e4):
                    self_dv = self_dv + calc_acceleration(obj.mass, obj.vec_r, self.vec_r) #G * obj.mass * obj.vec_r  / (norm(obj.vec_r)**3)
            else:
                self_destroy_flag = 1
                destroyed_obj = obj
                break

        if (self_destroy_flag == 1):
            #new_list_of_bodies = [elem for elem in bodies if (not elem.is_body(destroyed_obj))] 
            self.destroy(destroyed_obj)
            #print(f'1 body crashed on a star, modeling {len(new_list_of_bodies)} bodies')
            #self.gravitate(new_list_of_bodies, time)
            return [destroyed_obj]
        
            #Найдем вклад каждого из тел в ускорение
        for body in bodies: 
            #delta_v = body.vec_v * 0 #для того, чтобы сделать delta_v такой же размерности
            delta_v = - calc_acceleration(self.mass, body.vec_r, self.vec_r)# G * self.mass * (body.vec_r) / (norm(body.vec_r))**3
            for obj in bodies:
                if (norm(obj.vec_r-body.vec_r) > crash_distance):
                    delta_v = delta_v + calc_acceleration(obj.mass, obj.vec_r, body.vec_r)#G * obj.mass * (obj.vec_r - body.vec_r) / (norm(obj.vec_r - body.vec_r))**3      
                elif (not obj.is_body(body)):
                    '''new_list_of_bodies = [elem for elem in bodies if ((not elem.is_body(obj)) or (not elem.is_body(body)))] 
                        #destroy() реализованный через моделирование тех же тел, только без столкнувшихся
                        """new_body = CosmicBody(obj.mass+body.mass, (obj.mass*obj.vec_v + body.mass*body.vec_v), (obj.vec_r+body.vec_r)/2, 0)
                        #Добавим новое тело, получившееся в результате абсолютно неупрогого столкновения
                        new_list_of_bodies.append(new_body)"""
                    self.gravitate(new_list_of_bodies) '''
                    destroy_flag = 1
                    destroyed_obj = obj
                    break #выход из цикла obj
                
            if (destroy_flag == 1):
                destroy_flag = 0
                '''new_list_of_bodies = [elem for elem in bodies if ((not elem.is_body(destroyed_obj)) & (not elem.is_body(body)))] 
                #destroy() реализованный через моделирование тех же тел, только без столкнувшихся
                """new_body = CosmicBody(obj.mass+body.mass, (obj.mass*obj.vec_v + body.mass*body.vec_v), (obj.vec_r+body.vec_r)/2, 0)
                    #Добавим новое тело, получившееся в результате абсолютно неупрогого столкновения
                    new_list_of_bodies.append(new_body)"""'''
                #print(f'2 bodies crashed, modeling {len(new_list_of_bodies)} bodies')
                #self.gravitate(new_list_of_bodies, time) #запуск моделирования с новым списком тел. Получается рекурсия, но не предполагается моделирование для > 50 тел, поэтому исключается возможность ошибки
                return [body, destroyed_obj] #выход из цикла для моделирования с прошлым списком тел 
            
            body.move(dt/2)            #так как при равноускоренном движении dr = dt/2*(v1+v2)                
            body.accelerate(delta_v - self_dv) # так как a_rel = a_abs - a_trans (a_trans = - a_star)
            body.move(dt/2)
            body.trajectory.append(body.vec_r)

        #result = [body.vec_r for body in bodies]
        return []

    def gravitate(self, bodies, time): # юпитер ругается на вызов до инициализации, если использовать глобальное время, поэтому придется его передавать
        """
        Делает анимацию движения системы с помощью matplotlib,
        моделируя движение системы и выводя 3D-график траекторий 
        частиц и их текущее положение
        Parameters
        ----------
        bodies : list of CosmicBodies
            список тел, движение которых моделируется в поле звезды и друг друга 
        Returns
        -------
        None
        """
        if (len(bodies) == 0):
            time = total_time #заканчивает процесс моделирования, если все тела столкнулись
        # Задаем границы
        ax.axes.set_xlim3d(-size, size)
        ax.axes.set_ylim3d(-size, size)
        ax.axes.set_zlim3d(-size, size)    

        # Отображаем пустую фигуру
        fig.show()
        fig.canvas.draw()

        while ((time < total_time) and (max_r(bodies) < max_dist) ): 
            destroyed_bodies = self.move_system(bodies, time) 
            if (len(destroyed_bodies) != 0):
                new_list_of_bodies = [elem for elem in bodies if (not (elem in destroyed_bodies))] #.is_body(destroyed_obj)) & (not elem.is_body(body)))]
                self.gravitate(new_list_of_bodies, time)
                return
            
            #анимация для bodies
            #angle = 60 + 60 * time / total_time

            # Удаляем графики, отображенные на предыдущем шаге
            ax.clear()

            # Задаем границы, потому что они скидываются на дефолтные
            ax.axes.set_xlim3d(-size, size)
            ax.axes.set_ylim3d(-size, size)
            ax.axes.set_zlim3d(-size, size)
            
            # Отображаем движение точек
            ax.scatter([0],[0],[0], label = f'time = {round(time, 3)}')

            for body in bodies: 
                traj = np.array(body.trajectory).T
                ax.plot(traj[0], traj[1], traj[2])
                ax.scatter(body.vec_r[0], body.vec_r[1], body.vec_r[2])
            # Изменяем угол отображения графика
            #ax.view_init(30 - angle * 0.2, angle)
            ax.legend()
            # Перерисовываем фигуру
            fig.canvas.draw()

            time += dt


class CosmicBody(Star): # idk, на кой тут наследование, но кажется логичным, что СosmicBody можно gravitate друг относительно друга, но если у нас конкретная звезда, она можен не обладать доп методами

    def __init__(self, mass = None, vec_r = None, vec_v = None):
        """
        Создает объект класса Космическое тело, дочернего к классу Звезда
        с добавлением тректории - списка радиус-векторов тела на каждом шаге моделирования
        Parameters
        ----------
        mass : float
            Масса тела
        vec_r : np.ndarray
            Радиус-вектор тела
        vec_v : np.ndarray
            Вектор скорости тела
        Returns
        -------
        None
        """       
        self.trajectory = []        
        super().__init__(mass if mass is not None else rd.uniform(0,1),
                         np.array(vec_r) if (vec_r) is not None else np.array([rd.uniform(-radius,radius) for i in range(dim)]), 
                         np.array(vec_v) if (vec_v) is not None else np.array([rd.uniform(-vmax,vmax) for i in range(dim)]))
    

    def accelerate(self, Dv): 
        """
        Осуществляет изменение скорости тела под действием ускорения Dv за время dt
        Parameters
        ----------
        Dv : float
            ускорение тела
        Returns
        -------
        None
        """ 
        self.vec_v = self.vec_v + Dv * dt

    def move(self, time):
        """
        Изменяет положение тела при линейном движении при заданном времени 
        Parameters
        ----------
        time : float 
            время перемещения
        Returns
        -------
        None
        """ 
        self.vec_r = self.vec_r + time * self.vec_v

def max_r(bodies):
    """
    Находит модуль расстояния до наиболее удаленного объекта системы
    Parameters
    ----------
    bodies : list of CosmicBodies
        список моделируемых тел
    Returns
    -------
    max_r : float
        модуль расстояния до наиболее удаленного тела
    """
    r = 0
    for body in bodies:
        r = norm(body.vec_r) * (norm(body.vec_r) > r)
    return r

def max_norm(array):
    r = 0
    for i in array:
        r = norm(i) * (norm(i) > r)
    return r
    
                         

