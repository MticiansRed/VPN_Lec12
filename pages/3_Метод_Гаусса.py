import streamlit as st
 
menu = st.sidebar.radio('***',
    (
    "Исключение неизвестных", 
    "Решение системы уравнений",    
    "Компактная схема метода Гаусса", 
    "Расчетные формулы",    
    "Метод Гаусса с выбором главного элемента",     
    "Метод квадратного корня",   
    "Разреженные матрицы",       
    )
)
  
if menu == "Исключение неизвестных":
    r"""
##### Исключение неизвестных
Система уравнений

$\begin{aligned}
  a_{11} y_1 +  a_{12} y_2 + \cdots +  a_{1n} y_n  & = g_1  \\
  a_{21} y_1 +  a_{22} y_2 + \cdots +  a_{2n} y_n  & = g_2  \\
    \cdots   &  \\
  a_{n1} y_1 +  a_{n2} y_2 + \cdots +  a_{nn} y_n  & = g_n  \\
\end{aligned}$

Из первого уравнения

$\begin{aligned}
 y_1 +  u_{12} y_2 + \cdots +  u_{1n} y_n = v_1
\end{aligned}$

$\begin{aligned}
 u_{12} = \frac{a_{1j}}{a_{11}}, \quad j = 2,3, \ldots, n, 
 \quad v_1 = \frac{g_1}{a_{11}}   
\end{aligned}$

Подстановка $y_1$ в другие уравнения системы дает

$\begin{aligned}
  y_1 +  u_{12} y_2 + \cdots +  u_{1n} y_n  & = v_1  \\
  a_{22}^{(1)} y_2 + \cdots +  a_{2n}^{(1)} y_n  & = g_2^{(1)}  \\
    \cdots   &  \\
  a_{n2}^{(1)} y_2 + \cdots +  a_{nn}^{(1)} y_n  & = g_n^{(1)}  \\
\end{aligned}$

где

$\begin{aligned}
 a_{ij}^{(1)} = a_{ij} - u_{1j} a_{i1} ,
 \quad g_{i}^{(1)} = g_{i} - v_{1} a_{i1} ,
 \quad i,j = 2,3, \ldots, n  
\end{aligned}$

Далее последовательное исключение $y_2, y_3, \ldots, y_n$ 
    """

if menu == "Решение системы уравнений":
    r"""
##### Решение системы уравнений

Прямой ход метода Гаусса

$\begin{aligned}
  y_1 +  u_{12} y_2 + \cdots +  u_{1n} y_n  & = v_1  \\
  y_2 + \cdots +  u_{2n} y_n  & = v_2  \\
    \cdots   &  \\
  y_n  & = v_n \\
\end{aligned}$

Обратный ход

$\begin{aligned}
 y_n & = v_n \\
 y_i & = v_i - \sum_{j=i+1}^{n} u_{ij} y_j ,
 \quad i = n-1, n-2, \ldots, 1 
\end{aligned}$
    """

if menu == "Компактная схема метода Гаусса":
    r"""
##### Компактная схема метода Гаусса

$LU$-факторизация матрицы

$\begin{aligned}
 A = LU
\end{aligned}$

Нижняя ($L$) и верхняя ($U$) треугольные матрицы

$\begin{aligned}
 L = \begin{pmatrix}
  l_{11} & 0      & \cdots &  0       \\
  l_{21} & l_{22} & \cdots &  0       \\
  \cdots & \cdots & \cdots &  \cdots  \\
  l_{n1} & l_{n2} & \cdots &  l_{nn}  \\
\end{pmatrix} 
\end{aligned}$

$\begin{aligned}
  U = \begin{pmatrix}
  u_{11} & u_{12} & \cdots &  u_{1n}  \\
  0      & u_{22} & \cdots &  u_{2n}  \\
  \cdots & \cdots & \cdots &  \cdots  \\
  0      & 0      & \cdots &  u_{nn}  \\
\end{pmatrix}  
\end{aligned}$

при задании $l_{ii}$ или  $u_{ii}$ для $i = 1,2, \ldots, n$ 

Решение системы уравнений $A y = g$ 

$\begin{aligned}
 L v = g ,
 \quad U y = v 
\end{aligned}$
    """

if menu == "Расчетные формулы":
    r"""
##### Расчетные формулы

$\begin{aligned}
 A = LU,
 \quad a_{ij} = \sum_{s=1}^{n} l_{is} u_{sj} = \sum_{s=1}^{\min(i,j)} l_{is} u_{sj}
\end{aligned}$

Строки $1, 2, \ldots, k-1$ для $u$ и столбцы $1, 2, \ldots, k-1$ для $L$ вычислены

При $k=1$ 

$\begin{aligned}
 a_{11} = {\color{red} l_{11} u_{11}},
 \quad a_{1j} = l_{11} {\color{red} u_{1j}}, \quad j = 2,3, \ldots, n ,
 \quad a_{i1} = {\color{red} l_{i1}} u_{11}, \quad i = 2,3, \ldots, n  
\end{aligned}$

Для нового $k = 2, 3, \ldots, n$

$\begin{aligned}
 a_{kk} = {\color{red} l_{kk} u_{kk}} + \sum_{s=1}^{k-1} l_{ks} u_{sj}
\end{aligned}$

$\begin{aligned}
 a_{kj} = l_{kk} {\color{red} u_{kj}} + \sum_{s=1}^{k-1} l_{ks} u_{sj}, \quad j = k+1, k+2, \ldots, n 
\end{aligned}$

$\begin{aligned}
 a_{ik} = {\color{red} l_{ik}} u_{kk} + \sum_{s=1}^{k-1} l_{ks} u_{sj}, \quad i = k+1, k+2, \ldots, n  
\end{aligned}$
    """
    
if menu == "Метод Гаусса с выбором главного элемента":
    r"""
##### Метод Гаусса с выбором главного элемента

Метод Гаусса применим, когда все главные миноры $А$ отличны от нуля

$\begin{aligned}
  a_{11} \neq 0,
  \quad  \det \begin{pmatrix}
 a_{11} & a_{12} \\
 a_{21} & a_{22} \\
\end{pmatrix}  
\neq 0,
\ \cdots ,
\quad \det (A) \neq 0
\end{aligned}$

Выбор неизвестной для исключения

$\begin{aligned}
  a_{11} y_1 +  {\color{red} a_{12} y_2} + \cdots +  a_{1n} y_n  & = f_1  \\
  a_{21} y_1 +  a_{22} y_2 + \cdots +  a_{2n} y_n  & = f_2  \\
    \cdots   &  \\
  a_{n1} y_1 +  a_{n2} y_2 + \cdots +  a_{nn} y_n  & = f_n  \\
\end{aligned}$


Максимальный по модулю коэффициент в уравнении на каждом этапе прямого хода метода Гаусса
 
$\begin{aligned}
 |a_{12}| \geq |a_{1j}| , \quad j = 1,2, \ldots, n 
\end{aligned}$

Применимость: $\det (A) \neq 0$ 

    """

if menu == "Метод квадратного корня":
    r"""
##### Метод квадратного корня

Транспонированная матрица $A^*$: $\ A^* = \{a_{ji}\}$ 

Для симметричной 
вещественной матрицей $A = A^*$ используется 

$\begin{aligned}
  A = S^{*} D S
\end{aligned}$

$S ~-~$  верхняя треугольная матрица с $s_{ii} > 0, \ i = 1,2, \ldots, n$ 

$D ~-~$  диагональная матрица с элементами
$d_i, \ i =1,2,\ldots,n$, равными $\pm 1$

Расчетные формулы метода квадратного корня (метода Холецкого)

$\begin{aligned}
  d_{1} = {\rm sign} \ a_{11},
  \quad s_{11} = |a_{11}|^{1/2},
  \quad s_{1j} = a_{1j} / s_{11},
  \quad j = 2,3,\ldots,n 
\end{aligned}$

$\begin{aligned}
  d_{i} = {\rm sign} \left ( a_{ii} - \sum_{k=1}^{i-1}
  |s_{ki}|^{2} d_k \right ),
\quad 
  s_{ii} = \left | a_{ii} - \sum_{k=1}^{i-1} |s_{ki}|^{2} d_{k}
  \right |^{1/2} 
\end{aligned}$

$\begin{aligned}
  s_{ij} = \frac{1}{s_{ii} d_{i}}
  \left ( a_{ij} - \sum_{k=1}^{i-1} s_{ki} s_{kj} d_{k} \right ) 
\end{aligned}$

$\begin{aligned}
  i = 2,3,\ldots,n,
  \quad j = i+1,i+2,\ldots,n 
\end{aligned}$
    """
    
if menu == "Разреженные матрицы":
    r"""
##### Разреженные матрицы

* Вычислительные затраты метода Гаусса (плотная матрица) $\mathcal{O}(n^3)$ 

* Матрица с преимущественно нулевыми элементами ~-~ разреженная матрица

* Типичный случай: число ненулевых элементов $\mathcal{O}(n)$ 
    """

    

    
