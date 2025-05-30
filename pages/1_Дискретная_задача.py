import streamlit as st
 
menu = st.sidebar.radio('***',
    (
    "Вариационная задача", 
    "Конечно-элементное решение",    
    "Метод Галеркина", 
    "Система линейных уравнений",    
    )
)
  
if menu == "Вариационная задача":
    r"""
##### Вариационная задача

Пространства

* $V ~-~$ пространство достаточно гладких функций
* $u, v \in V ~-~$ естественные краевые условия

Вариационная задача: найти  $u \in V$ такую, что

$\begin{aligned}
 a(u,v) = l(v) ,
 \quad \forall v \in V
\end{aligned}$

Билинейная форма 

$\begin{aligned}
 & a(\sigma u_1 + u_2,v) = \sigma a(u_1,v) +  a(u_2,v) \\
 & a(u, \sigma v_1 + v_2) = \sigma a(u,v_1) +  a(u,v_2) 
\end{aligned}$

Линейная форма 

$\begin{aligned}
 l(\sigma v_1 + v_2) = \sigma l(v_1) +  l(v_2) 
\end{aligned}$

    """    
if menu == "Конечно-элементное решение":
    r"""
##### Конечно-элементное решение

Базис (пробные функции) в пространстве конечных элементов $V_h \ (V_h \subset V)$ 

$\varphi_i(x) , \quad i = 1,2, \ldots, n$

$\begin{aligned}
 \varphi_i(x_j) = \left \{ \begin{array}{cc}
  1 ,  &  i = j \\
  0 ,  &  i \neq j \\
\end{array}
\right .
\end{aligned}$

Приближенное решение

$\begin{aligned}
   y(x) = \sum_{i=1}^{n} y_i \varphi_i(x) ,
   \quad y_i = y(x_i),
   \quad i = 1,2, \ldots, n
\end{aligned}$

Вариационная задача: найти  $y \in V_h$ такую, что

$\begin{aligned}
 a(y,v) = l(v),
 \quad \forall v \in V_h
\end{aligned}$

    """
    
if menu == "Метод Галеркина":
    r"""
##### Метод Галеркина

В вариационном уравнении

$\begin{aligned}
 a(y,v) = l(v)
\end{aligned}$

положим

$\begin{aligned}
   y(x) = \sum_{j=1}^{n} y_j \varphi_j(x) ,
\end{aligned}$   
$\begin{aligned}   
   v = \varphi_i(x),
   \quad i = 1,2, \ldots, n
\end{aligned}$

Определение приближенного решения из

$\begin{aligned}
 \sum_{j=1}^{n} a(\varphi_j, \varphi_i) y_j = l(\varphi_i),
 \quad i = 1,2, \ldots , n
\end{aligned}$

Вычисление (квадратурные формулы)

$\begin{aligned}
  a(\varphi_j, \varphi_i) ,
  \quad l(\varphi_i) ,
  \quad i, j = 1,2, \ldots , n
\end{aligned}$


    """
    
if menu == "Система линейных уравнений":
    r"""
##### Система линейных уравнений

Для решения $y=\{y_i\}$

$\begin{aligned}
  A y = g
\end{aligned}$

Коэффициенты уравнения

$\begin{aligned}
  A = \{a_{ij}\}, 
  \quad a_{ij} =  a(\varphi_j, \varphi_i) ,
  \quad i, j = 1,2, \ldots , n
\end{aligned}$

$\begin{aligned}
  g = \{g_i\}, 
  \quad g_i =  l(\varphi_i) ,
  \quad i = 1,2, \ldots , n
\end{aligned}$

Для симметричной билинейной формы

$a(u,v) = a(v,u)$

система уравнений с симметричной матрицей

$\begin{aligned}
  A = \{a_{ij}\}, 
  \quad a_{ij} =  a_{ji} ,
  \quad i, j = 1,2, \ldots , n
\end{aligned}$

    """
    



























