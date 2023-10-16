# Экстраполяция Ричардсона
## Входные данные
* Отрезок [a,b];
* n - количество уравнений в системе;
* k - количество элементов в разбиении [a,b]
* y<sub>10</sub> y<sub>20</sub> y<sub>30</sub> ... y<sub>n0</sub> - условие Коши для данной системы уравнений;
* e - заданная точность решения

## Структура проекта:
* funk_NN.c - содержит функции вида y'<sub>i</sub> = f<sub>i</sub>(x, y<sub>1</sub>, y<sub>2</sub>,..., y<sub>n</sub>), i = <span style="text-decoration: overline;">1,n</span>
* solve_ODE.c - содержит алгоритм решения задачи.
* main_NN.c - содержит описание входных данных и вызов функции решения solve_ODE.c
## Описание алгоритма
1. Отрезок [a,b] делится на отрезки [x<sub>i</sub>, x<sub>i+1</sub>], i = <span style="text-decoration: overline;">0,k</span>, |x<sub>i+1</sub> - x<sub>i</sub>| = h, где h - диаметр разбиения k.
2. Находим методом Рунге-Кутты s решений в отрезке [x<sub>i</sub>, x<sub>i+1</sub>] для разбиений h, 2h,..., 2sh.
3. Строим систему уравнений вида y<sub>m</sub>(h<sub>m</sub>, x<sub>i+1</sub>) = A<sub>m</sub>e<sup>h<sub>m</sub></sup>+A<sub>2</sub>e<sup>2h<sub>m</sub></sup>+A<sub>3</sub>e<sup>4h<sub>m</sub></sup>+...+A<sub>s</sub>e<sup>2sh<sub>m</sub></sup>, m = <span style="text-decoration: overline;">1,s</span>
4. Подставляем найденные коэффициенты в уравнение вида y(h) = A<sub>1</sub>e<sup>h</sup>+A<sub>2</sub>e<sup>2h</sup>+A<sub>3</sub>e<sup>4h</sup>+...+A<sub>s</sub>e<sup>2sh</sup>
5. Повторяем для h' = h/2, пока не будет выполнено неравенство:
$$
\sqrt{\Sigma(y_i - y_{2i})^2}<\epsilon
$$
6. Выводим результат для каждого x из разбиения h:
$$
\begin{matrix}
x_0 & y_{10} & y_{20} &...& y_{n0} \\
x_1 & y_{11} & y_{21} &...& y_{n1} \\
x_2 & y_{12} & y_{22} &...& y_{n2} \\
&&.......& \\
x_k & y_{1k} & y_{2k} &...& y_{nk}

\end{matrix}
$$


