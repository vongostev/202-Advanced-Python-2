# Задания по Python. Ноябрь - декабрь
## Игра "Астероиды и кометы"
Игра направлена на закрепление навыков работы с ООП, numpy, оптимизацией, отображением графической информации, автоматическим тестированием.
### Условия игры
Моделируется звездная система в собственной системе координат звезды. Изначально в системе существует только звезда, у которой задан один параметр, масса <img src="https://render.githubusercontent.com/render/math?math=M">, и которая находится в начале координат. Из точки <img src="https://render.githubusercontent.com/render/math?math=P"> запускается космический объект массы <img src="https://render.githubusercontent.com/render/math?math=m">, движущийся равномерно со скоростью <img src="https://render.githubusercontent.com/render/math?math=\vec{v}">. Космический объект попадает в поле тяготения звезды, и на него действует сила тяготения

<img style="float: center;" src="https://render.githubusercontent.com/render/math?math=F=G\dfrac{mM}{r_P^2},\quad\quad(1)">

ускоряющая космический объект. Здесь <img src="https://render.githubusercontent.com/render/math?math=r_P"> -- модуль радиус-вектора точки <img src="https://render.githubusercontent.com/render/math?math=P">. Изменение импульса объекта <img src="https://render.githubusercontent.com/render/math?math=p"> описывается вторым законом Ньютона:

<img style="float: center;" src="https://render.githubusercontent.com/render/math?math=\dfrac{dp}{dt}=F.\quad\quad(2)">

Изменение скорости приводит к искривлению траектории космического объекта. Существует четыре варианта:
1. Объект падает на звезду
2. Объект выходит на стабильную эллиптическую орбиту и становится астероидом
3. Объект выходит на стабильную орбиту в виде конического сечения и становится кометой
4. Объект уходит из звездной системы
Необходимо смоделировать формирование такой звездной системы при случайной инициализации параметров космических объектов. При наличии других тел в системе необходимо учитывать и их влияние на космический объект. Ниже предлагается одна из схем дискретизации задачи.
### Дискретизация задачи
Введем дискретную шкалу времени <img src="https://render.githubusercontent.com/render/math?math=t_i">. В момент времени <img src="https://render.githubusercontent.com/render/math?math=t_0"> объект со скоростью <img src="https://render.githubusercontent.com/render/math?math=\vec{v}"> и массой <img src="https://render.githubusercontent.com/render/math?math=m"> генерируется в точке пространства <img src="https://render.githubusercontent.com/render/math?math=P">. В каждый момент времени на него действует сила тяготения (1), в связи с чем его импульс изменяется согласно (2). Учитывая, что <img src="https://render.githubusercontent.com/render/math?math=p=m\vec{v}">, в дискретной шкале времени для произвольно момента времени <img src="https://render.githubusercontent.com/render/math?math=t_i"> можно приблизительно написать:

<img style="float: center;" src="https://render.githubusercontent.com/render/math?math=m\dfrac{\vec{v}_{i+1} - \vec{v}_i}{\Delta t}=G\dfrac{mM}{r_i^2},\quad\quad(3)">

Где <img src="https://render.githubusercontent.com/render/math?math=\Delta t=t_{i+1} - t_i">. Отсюда изменение скорости

<img style="float: center;" src="https://render.githubusercontent.com/render/math?math=\Delta\vec{v}_i=\vec{v}_{i+1}  - \vec{v}_i= G\dfrac{M\Delta t}{r_i^2}.\quad\quad(4)">

Если объект находится в поле тяготения нескольких тел, то формулу (4) можно обобщить:

<img style="float: center;" src="https://render.githubusercontent.com/render/math?math=\Delta\vec{v}_i= G\sum_k\dfrac{m_k\Delta t}{r_{ik}^2},\quad\quad(5)">

где индекс <img src="https://render.githubusercontent.com/render/math?math=k"> соответствует отдельному телу в системе.
Если на каждом шагу рассчитывать изменение скорости, то можно рассчитать новый радиус-вектор:

<img style="float: center;" src="https://render.githubusercontent.com/render/math?math=\vec{r}_{i+1}=\vec{r}_i+\Delta\vec{v}_i\Delta t.">

Таким образом, в дискретном виде возможно смоделировать систему многих тел.
### Задания
1. Описать базовые классы звезды и космического объекта `Star(mass: float)`, `CosmicBody(mass: float, vec_v: numpy.ndarray, vec_P: numpy.ndarray)` в двухмерном пространстве
2. Реализовать метод гравитационного взаимодействия между телом и другими телами, например `.gravitate(bodys: list)`
3. Реализовать тесты корректного создания и взаимодействия объектов
4. Реализовать цикл симуляции для одного объекта: в начальный момент времени генерируется объект со случайными характеристиками, -- необходимо найти его траекторию в течение определенного числа шагов.
5. Реализовать цикл симуляции для нескольких объектов, которые появляются в разные моменты времени.
6. Реализовать отображение траекторий тел по окончании цикла симуляции.
7. Реализовать детектирование четырех типов траекторий тел
8. Реализовать тесты детектирования траекторий
9. Реализовать оптимизационный алгоритм поиска оптимальной скорости для объекта с заданной массой для скорейшего выхода на стабильную траекторию.
10. Обобщить код на трехмерный случай, включая тесты
11. Реализовать отображение трехмерных траекторий объектов
12. Реализовать анимацию движений объектов в системе
