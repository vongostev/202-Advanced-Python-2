import pygame
import numpy as np
import matplotlib.pyplot as plt
import math
import sys
import numba as nb
from numba import jit

pygame.font.init()
print('😀')

G = 4*math.pi**2/332950

class Planet:
    
    
    def __init__(self, mass=0, r=None, v=None, R=0):
        if r is None: r = np.zeros(3)
        if v is None: v = np.zeros(3)
        self.m, self.r, self.v, self.R = float(mass), np.array(r), np.array(v), R
        self.last_coor = []
        self.last_coor.append(np.zeros(3))
        self.exist = 1

    def __str__(self):
        
        return f"Mass: {round(self.m, 2)}; Coordinates: {np.round(self.r, 2)}; Velocity: {np.round(self.v, 2)}"

    def Ek(self):
        
        return self.m * np.dot(self.v, self.v) / 2

    def if_destroyed(self, R):
        
        if np.linalg.norm(self.r) <= R:
            self.m, self.v, self.exist = 0, np.zeros(len(self.r)), 0


def Ep_2b(mom1, mom2) -> np.float64:
    return -G*mom1.m*mom2.m/np.linalg.norm(mom2.r-mom1.r)

def totalE(mom) -> np.float64:
    print('######')
    total = 0
    for i in range(len(mom)):
        total += mom[i].Ek()
        i_ = i+1
        while(i_ < len(mom)):
            total += Ep_2b(mom[i], mom[i_])
            print(i, i_)
            i_ += 1
            
    print('#######')
    return total 

def totalEForOne(mom, bodies) -> np.float64:
    total = mom.Ek()
    for i in bodies:
        if (i != mom):
            total += Ep_2b(mom, i)
    
    return total

def typeOrbit(mom):
    E = totalEForOne(mom, planets)
    if (E>0):
        return 'Гипербола'
    elif (E<0):
        return 'Эллипс'
    else:
        return 'Парабола'

##@jit(nopython = False)
def step(sun, bodies, method = 0):
    G = 4*math.pi**2/332950 
    dt = 5e-3
    if method == 0:
        #acs = [G*sun.m*body.m*(body.r - sun.r)/np.dot((body.r - sun.r), (body.r - sun.r))**1.5 for body in bodies]
        k1 = [dt*body.v for body in bodies]
        kv1 = [-G*sun.m*(bodies[j].r - sun.r)/np.dot((bodies[j].r - sun.r), (bodies[j].r - sun.r))**1.5*dt for j in range(len(bodies))]
        k2 = [dt*(bodies[j].v + kv1[j]/2) for j in range(len(bodies))]
        kv2 = [-G*sun.m*((bodies[j].r + k1[j]/2) - sun.r)/np.dot(((bodies[j].r + k1[j]/2) - sun.r), ((bodies[j].r + k1[j]/2) - sun.r))**1.5*dt for j in range(len(bodies))]
        k3 = [dt*(bodies[j].v + kv2[j]/2) for j in range(len(bodies))] 
        kv3 = [-G*sun.m*((bodies[j].r + k2[j]/2) - sun.r)/np.dot(((bodies[j].r + k2[j]/2) - sun.r), ((bodies[j].r + k2[j]/2) - sun.r))**1.5*dt for j in range(len(bodies))]
        k4 = [dt*(bodies[j].v + kv3[j]) for j in range(len(bodies))] 
        kv4 = [-G*sun.m*((bodies[j].r + k3[j]) - sun.r)/np.dot(((bodies[j].r + k3[j]) - sun.r), ((bodies[j].r + k3[j]) - sun.r))**1.5*dt for j in range(len(bodies))]
        for j in range(len(bodies)):
            #print(type(bodies[j].r), type(k1[j]))
            bodies[j].r += (k1[j] + 2*k2[j] + 2*k3[j] + k4[j])/6
            bodies[j].v += (kv1[j] + 2*kv2[j] +2*kv3[j] + kv4[j])/6
        #print('k1 m0', k1, '\n\n\n')
    
    if method == 1:
        bodies.append(sun)
        lb = len(bodies)
        #acs = [G*sun.m*body.m*(body.r - sun.r)/np.dot((body.r - sun.r), (body.r - sun.r))**1.5 for body in bodies]
        
        kv1 = np.zeros_like(bodies)
        kv2 = np.zeros_like(bodies)
        kv3 = np.zeros_like(bodies)
        kv4 = np.zeros_like(bodies)
        b = bodies
        
        k1 = [dt*body.v for body in b]
        for j in nb.prange(lb):
            for i in nb.prange(lb):
                if i != j:
                    kv1[j] += -G*b[i].m*(b[j].r - b[i].r)/np.dot((b[j].r - b[i].r), (b[j].r - b[i].r))**1.5*dt
        #print('k1 m1', k1[0], ' b ', b[0].r)
        
        k2 = [dt*(b[j].v + kv1[j]/2) for j in nb.prange(lb)]
        
        for j in nb.prange(lb):
            for i in nb.prange(lb):
                if i != j:
                    #print(b[i].m, len(b), len(k1), len(b[i].r))
                    dr = ((b[j].r + k1[j]/2) - b[i].r)
                    #print(dr, kv2)
                    kv2[j] += -dt*G*b[i].m*dr/np.linalg.norm(dr)**3
                    
        
        k3 = [dt*(b[j].v + kv2[j]/2) for j in nb.prange(lb)]
        
        for j in nb.prange(lb):
            for i in nb.prange(lb):
                if i != j:
                    kv3[j] += -G*b[i].m*((b[j].r + k2[j]/2) - b[i].r)/np.dot(((b[j].r + k2[j]/2) - b[i].r), ((b[j].r + k2[j]/2) - b[i].r))**1.5*dt 
        
        k4 = [dt*(b[j].v + kv3[j]) for j in nb.prange(lb)] 
        
        for j in nb.prange(lb):
            for i in nb.prange(lb):
                if i != j:
                    kv4[j] += -G*b[i].m*((b[j].r + k3[j]) - b[i].r)/np.dot(((b[j].r + k3[j]) - b[i].r), ((b[j].r + k3[j]) - b[i].r))**1.5*dt 
        
        for j in range(lb):
            #print(type(bodies[j].r), type(k1[j]))
            b[j].r += (k1[j] + 2*k2[j] + 2*k3[j] + k4[j])/6
            b[j].v += (kv1[j] + 2*kv2[j] + 2*kv3[j] + kv4[j])/6
        sun = b.pop()
        #print(sun.m, b)
    return 0


WHITE = (255, 255, 255)
RED = (225, 0, 50)
GREEN = (0, 225, 0)
BLUE = (0, 0, 225)

pygame.init()
sc = pygame.display.set_mode((1000, 1000))
sc.fill(WHITE)
#pygame.draw.rect(sc, GREEN, (101, 101, 19, 19))
pygame.display.update()
Sun = Planet(332950, [0., 0., 0.], [0., 0., 0.])
planets = []
planets.append(Planet(1, [1., 0., 0.], [0., 6.28/1.4, 0]))
#planets.append(Planet(100000, [-1., 0., 0.], [0., -3.14, -3.14]))
#print(A)
phi = 0
theta = 0
exit_flag = True
w_go = 1
ppae = 100

def set_axes(sc):
    
    Mz = np.array([[math.cos(phi), -math.sin(phi), 0],
                       [math.sin(phi), math.cos(phi), 0],
                           [0, 0, 1]])
    Mx = np.array([[1, 0, 0],
                       [0, math.cos(theta), -math.sin(theta)],
                       [0, math.sin(theta), math.cos(theta)]])
    xmin = np.array([-4, 0, 0])
    xmax = np.array([4, 0, 0])
    ymin = np.array([0, -4, 0])
    ymax = np.array([0, 4, 0])
    zmin = np.array([0, 0, -4])
    zmax = np.array([0, 0, 4])
    xmin_ = np.dot(np.dot(Mx, Mz), xmin.reshape(3, 1))
    xmax_ = np.dot(np.dot(Mx, Mz), xmax.reshape(3, 1))
    ymax_ = np.dot(np.dot(Mx, Mz), ymax.reshape(3, 1))
    ymin_ = np.dot(np.dot(Mx, Mz), ymin.reshape(3, 1))
    zmin_ = np.dot(np.dot(Mx, Mz), zmin.reshape(3, 1))
    zmax_ = np.dot(np.dot(Mx, Mz), zmax.reshape(3, 1))
    #pygame.draw.circle(sc, (250, 250, 250), (int(r_[0]*100)+500, int(r_[1]*100)+500), 3, 3)
    #pygame.draw.circle(sc, (250, 250, 250), (int(r_[0]*100)+500, int(r_[1]*100)+500), 3, 3)
    #pygame.draw.circle(sc, (250, 250, 250), (int(r_[0]*100)+500, int(r_[1]*100)+500), 3, 3)
    
    ##x
    pygame.draw.aaline(sc, (0, 0, 0), [int(xmin_[0]*100)+500, int(xmin_[1]*100)+500], [int(xmax_[0]*100)+500, int(xmax_[1]*100)+500])  
    f1 = pygame.font.Font(None, 36)
    text1 = f1.render('X', True, (180, 0, 0))
    sc.blit(text1, [int(xmax_[0]*100)+510, int(xmax_[1]*100)+510])
    i = 1
    ex = xmax_/np.linalg.norm(xmax_)
    while i < 5:
        pygame.draw.circle(sc, (0, 0, 0), (int(ex[0]*i*100)+500, int(ex[1]*i*100)+500), 2, 2)
        pygame.draw.circle(sc, (0, 0, 0), (int(-ex[0]*i*100)+500, int(-ex[1]*i*100)+500), 2, 2)
        i += 1
        
    ##y
    pygame.draw.aaline(sc, (0, 0, 0), [int(ymin_[0]*100)+500, int(ymin_[1]*100)+500], [int(ymax_[0]*100)+500, int(ymax_[1]*100)+500])
    text2 = f1.render('Y', True, (180, 0, 0))
    sc.blit(text2, [int(ymax_[0]*100)+510, int(ymax_[1]*100)+510])    
    i = 1
    ex = ymax_/np.linalg.norm(ymax_)
    while i < 5:
        pygame.draw.circle(sc, (0, 0, 0), (int(ex[0]*i*100)+500, int(ex[1]*i*100)+500), 2, 2)
        pygame.draw.circle(sc, (0, 0, 0), (int(-ex[0]*i*100)+500, int(-ex[1]*i*100)+500), 2, 2) 
        i += 1
    
    ##z
    pygame.draw.aaline(sc, (0, 0, 0), [int(zmin_[0]*100)+500, int(zmin_[1]*100)+500], [int(zmax_[0]*100)+500, int(zmax_[1]*100)+500])    
    text3 = f1.render('Z', True, (180, 0, 0))
    sc.blit(text3, [int(zmax_[0]*100)+510, int(zmax_[1]*100)+510])    
    i = 1
    ex = zmax_/np.linalg.norm(zmax_)
    while i < 5:
        pygame.draw.circle(sc, (0, 0, 0), (int(ex[0]*i*100)+500, int(ex[1]*i*100)+500), 2, 2)
        pygame.draw.circle(sc, (0, 0, 0), (int(-ex[0]*i*100)+500, int(-ex[1]*i*100)+500), 2, 2)  
        i += 1
    

E0 = totalE(planets) + totalEForOne(Sun, planets)
counter = 1
total_e_list = [E0]

for i in planets:
    print(i)
    print("Orbit type: " + typeOrbit(i))


while exit_flag:
    
    draw_time_scale = int(20)
     
    keys = pygame.key.get_pressed()
    if keys[pygame.K_RIGHT] and keys[pygame.K_r]:
        phi += math.pi/180
    if keys[pygame.K_LEFT] and keys[pygame.K_r]:
        phi -= math.pi/180
    if keys[pygame.K_UP] and keys[pygame.K_r]:
        theta += math.pi/180
    if keys[pygame.K_DOWN] and keys[pygame.K_r]:
        theta -= math.pi/180
    
    if keys[pygame.K_UP] and keys[pygame.K_s]:
        ppae = ppae*1.1
    if keys[pygame.K_DOWN] and keys[pygame.K_s]:
        ppae = ppae/1.1
        
    if keys[pygame.K_s] and keys[pygame.K_t] and keys[pygame.K_o] and keys[pygame.K_p]:
        exit_flag = False   
    if keys[pygame.K_w]:
        w_go = (w_go + 1) % 2    
    if w_go:
        step(Sun, planets, 1)
        
        counter += 1
        if (counter%10 == 0):
            E = totalE(planets) + totalEForOne(Sun, planets)
            counter = 1
            total_e_list.append(E)
            #assert np.isclose(E, E0, rtotal = 0.05)              
        
    
    sc.fill(WHITE)
    f1 = pygame.font.Font(None, 32)
    text1 = f1.render('phi = ' + str(round(phi*180/math.pi, )) + ' theta = ' + str(round(theta*180/math.pi, 1)), True, (180, 0, 0))
    sc.blit(text1, (20, 20))
    text2 = f1.render('ax scale = ' + str("{:.3e}".format(100/ppae)) + 'ae', True, (180, 0, 0))
    #print(ppae)
    sc.blit(text2, (20, 40))
        
        
    Mz = np.array([[math.cos(phi), -math.sin(phi), 0],
                           [math.sin(phi), math.cos(phi), 0],
                           [0, 0, 1]])
    Mx = np.array([[1, 0, 0],
                       [0, math.cos(theta), -math.sin(theta)],
                       [0, math.sin(theta), math.cos(theta)]])
        
    set_axes(sc)
    
    pygame.draw.circle(sc, (240, 240, 22), (int(np.dot(np.dot(Mx, Mz), Sun.r.reshape(3, 1))[0]*ppae)+500, int(np.dot(np.dot(Mx, Mz), Sun.r.reshape(3, 1))[1]*ppae)+500), 10, 10)
    #print('iter')
    for p in planets:
            
           
        r_ = np.dot(np.dot(Mx, Mz), p.r.reshape(3, 1))
        pygame.draw.circle(sc, (0, 0, 220), (int(r_[0]*ppae)+500, int(r_[1]*ppae)+500), 3, 3)
        #print(p.last_coor)
        
        p.last_coor.append(np.array([p.r[i] for i in range(3)]))
        trace =  [[int(np.dot(np.dot(Mx, Mz), v.reshape(3, 1))[0]*ppae + 500), int(np.dot(np.dot(Mx, Mz), v.reshape(3, 1))[1]*ppae + 500)] for v in p.last_coor]  
        for i in range(2, len(trace) - 1):
            pygame.draw.aaline(sc, (0, 0, 220), trace[i-1], trace[i])
        #print([[int(np.dot(np.dot(Mx, Mz), v.reshape(3, 1))[0]*ppae + 500), int(np.dot(np.dot(Mx, Mz), v.reshape(3, 1))[1]*ppae + 500) ] for v in p.last_coor])
       # p.last_coor.append(p.r)
    pygame.display.update()
    for i in pygame.event.get():
        if i.type == pygame.QUIT:
            sys.exit()
        #elif i.type == pygame.KEYDOWN:



  
    pygame.time.delay(draw_time_scale)
    
plt.plot(total_e_list)
plt.show()