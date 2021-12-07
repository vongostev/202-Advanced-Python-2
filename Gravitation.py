import matplotlib.pyplot as plt 
import numpy as np 
from numpy.linalg import norm 
import random as rd
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

G=1 # в системе единиц "годы, массы Солнца, а.е."
dt = 0.01
#total_
total_time = 15
dim = 3 #размерность задачи
destroy_flag = 0 
vmax = 3 #максимальная скорость моделируемых тел
radius = 2 #максимальное расстояние от начала координат для моделируемых тел

class Star():
    def __init__(self, mass = rd.uniform(1, 10), vec_r = np.array([rd.uniform(-radius,radius) for i in range (dim)]), vec_v = np.array([rd.uniform(-vmax,vmax) for i in range (dim)])): #, dv=0):
        self.mass = mass
        self.vec_v = np.array(vec_v)
        self.vec_r = np.array(vec_r)
        #self.dv = dv #np.array(dv)
    
    def show(self):
        print(f'm = {self.mass}, vec_r = {self.vec_r}, vec_v = {self.vec_v}')

    def is_body(self, body):
        res = 1
        for i in range(dim):
            res = res & (self.vec_r[i] == body.vec_r[i]) & (self.vec_v[i] == body.vec_v[i]) 
            #print(self.vec_r[i], body.vec_r[i], self.vec_v[i], body.vec_r[i])
        res = (self.mass == body.mass) & res #& (self.dv ==  body.dv)
        #print(res)
        return res

    def destroy(self, body):
        self.vec_r = (self.vec_r + body.vec_r)/2
        self.vec_v = (self.mass * self.vec_v + body.mass * body.vec_r) / (self.mass + body.mass)       
        self.mass = self.mass + body.mass

    def gravitate(self, bodies):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        if ((norm(self.vec_v) > 0)or(norm(self.vec_r)>0)): #проверка того, перешли ли мы уже в СО звезды
            for body in bodies:
                body.vec_v = body.vec_v - self.vec_v # переход в СО звезды
                body.vec_r = body.vec_r - self.vec_r
            self.vec_v = self.vec_v * 0
            self.vec_r = self.vec_r * 0
        time = 0
        while((max_r(bodies) < 50) and (time < total_time)): #моделируем, пока что-то не улетит за пределы СС
        #while(max(norm(bodies).vec_r) < 50):

            self_destroy_flag = 0
            destroy_flag = 0

            #найдем ускорение самой звезды 
            self_dv = self.vec_v * 0  #для соблюдения размерности
            for obj in bodies:
                if (norm(obj.vec_r-self.vec_r) > 0.001):
                    self_dv = self_dv + G * obj.mass * obj.vec_r  / (norm(obj.vec_r)**3)
                else:
                    self_destroy_flag = 1
                    destroyed_obj = obj
                    break

            if (self_destroy_flag == 1):
                new_list_of_bodies = [elem for elem in bodies if (not elem.is_body(destroyed_obj))] 
                self.destroy(destroyed_obj)
                self.gravitate(new_list_of_bodies)
                return

            #Найдем вклад каждого из тел в ускорение
            for body in bodies: 
                delta_v = body.vec_v * 0 #для того, чтобы сделать delta_v такой же размерности
                body.move(dt/2)            #так как при равноускоренном движении dr = dt/2*(v1+v2)

                for obj in bodies:
                    if (norm(obj.vec_r-body.vec_r) > 0.001):
                        delta_v = delta_v + G * obj.mass * (obj.vec_r - body.vec_r) / (norm(obj.vec_r - body.vec_r))**3      
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
                
                body.accelerate(delta_v - self_dv) # так как a_rel = a_abs - a_trans

                if (destroy_flag == 1):
                    destroy_flag = 0
                    new_list_of_bodies = [elem for elem in bodies if ((not elem.is_body(destroyed_obj)) & (not elem.is_body(destroyed_obj)))] 
                    #destroy() реализованный через моделирование тех же тел, только без столкнувшихся
                    """new_body = CosmicBody(obj.mass+body.mass, (obj.mass*obj.vec_v + body.mass*body.vec_v), (obj.vec_r+body.vec_r)/2, 0)
                    #Добавим новое тело, получившееся в результате абсолютно неупрогого столкновения
                    new_list_of_bodies.append(new_body)"""
                    self.gravitate(new_list_of_bodies) #запуск моделирования с новым списком тел. Получается рекурсия, но не предполагается моделирование для > 50 тел, поэтому исключается возможность ошибки
                    return #выход из цикла для моделирования с прошлым списком тел 

                body.move(dt/2)            #так как при равноускоренном движении dr = dt/2*(v1+v2), это даст точность на порядок выше, чем dr = v dt 
                
            time = time + dt
            #анимация для body

            
"""
ax.set_xlim3d([-1.0, 1.0])
ax.set_xlabel('X')
ax.set_ylim3d([-1.0, 1.0])
ax.set_ylabel('Y')
ax.set_zlim3d([0.0, 10.0])
ax.set_zlabel('Z')
ani = animation.FuncAnimation(fig, update, N, fargs=(data, line), interval=10000/N, blit=False)
#ani.save('matplot003.gif', writer='imagemagick')
plt.show()"""
           

class CosmicBody(Star): 

    def __init__(self, mass = rd.uniform(0.1, 1), vec_r = np.array([rd.uniform(-radius,radius) for i in range (dim)]), vec_v = np.array([rd.uniform(-vmax, vmax) for i in range (dim)])): #, dv=0):
        super().__init__(mass, vec_v, vec_r)#, dv)
    
    '''def change_dv(self, changed_dv):
        self.dv = changed_dv'''

    def accelerate(self, Dv): 
        self.vec_v = self.vec_v + Dv
        #self.dv = Dv #на самом деле оно нам и не нужно, если есть delta_v в цикле

    def move(self, time):
        self.vec_r = self.vec_r + time * self.vec_v

def max_r(bodies):
    r = 0
    for body in bodies:
        if (norm(body.vec_r) > r):
            r = norm(body.vec_r)
    return r
'''
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.set_xlim3d([-1.0, 1.0])
ax.set_xlabel('X')

ax.set_ylim3d([-1.0, 1.0])
ax.set_ylabel('Y')

ax.set_zlim3d([0.0, 10.0])
ax.set_zlabel('Z')

ani = animation.FuncAnimation(fig, update, N, fargs=(data, line), interval=10000/N, blit=False)
#ani.save('matplot003.gif', writer='imagemagick')
plt.show()
'''

def test_functions():
    a = CosmicBody()
    b = Star()
    a.show()
    b.show()
    b.gravitate([a])

test_functions()