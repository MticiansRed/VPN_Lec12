import streamlit as st
 
menu = st.sidebar.radio('***',
    (
    "Система линейных уравнений", 
    "Квадратные матрицы",    
    "Решение систем линейных уравнений", 
    "Простейшие задачи",    
    )
)
  
if menu == "Система линейных уравнений":
    r"""
##### Система линейных уравнений

Система из $n$ уравнений для $n$ неизвестных $y_1,y_2,\ldots,y_n$

$\begin{aligned}
  a_{11} y_1 +  a_{12} y_2 + \cdots +  a_{1n} y_n  & = g_1  \\
  a_{21} y_1 +  a_{22} y_2 + \cdots +  a_{2n} y_n  & = g_2  \\
    \cdots   &  \\
  a_{n1} y_1 +  a_{n2} y_2 + \cdots +  a_{nn} y_n  & = g_n  \\
\end{aligned}$

Вектор неизвестных и вектор правых частей

$\begin{aligned}
 y = \{y_i\} = \begin{pmatrix}
   y_1  \\
   y_2  \\
   \cdots  \\
   y_n \\
\end{pmatrix} ,
\quad  
 g = \{g_i\} = \begin{pmatrix}
   g_1  \\
   g_2  \\
   \cdots  \\
   g_n \\
\end{pmatrix} 
\end{aligned}$

Матрица коэффициентов

$\begin{aligned}
A = \{a_{ij}\} = \begin{pmatrix}
  a_{11} & a_{12} & \cdots &  a_{1n}  \\
  a_{21} & a_{22} & \cdots &  a_{2n}  \\
  \cdots & \cdots & \cdots &  \cdots  \\
  a_{n1} & a_{n2} & \cdots &  a_{nn}  \\
\end{pmatrix} 
\end{aligned}$
    """

if menu == "Квадратные матрицы":
    r"""
##### Квадратные матрицы

Умножение матрицы на вектор

$\begin{aligned}
 (A)_{ij} = a_{ij}, \ i, j = 1,2, \ldots, n, \quad (y)_{i} = y_i, \ i = 1,2, \ldots, n
\end{aligned}$

$\begin{aligned}
 (Ay)_i = \sum_{j=1}^{n} a_{ij} y_j, \quad i = 1,2, \ldots, n
\end{aligned}$

Умножение матрицы на матрицу

$\begin{aligned}
 (A)_{ij} = a _{ij}, \quad (B)_{ij} = b_{ij},\quad  \ i, j = 1,2, \ldots, n
\end{aligned}$

$\begin{aligned}
 (AB)_{ij} = \sum_{k=1}^{n} a_{ik} b_{kj} ,
 \quad i, j = 1,2, \ldots, n
\end{aligned}$
    """
    
if menu == "Решение систем линейных уравнений":
    r"""
##### Решение систем линейных уравнений

$A^{-1} ~-~$  обратная матрица: $A A^{-1} = I$ 

Обратная матрица существует при $ {\rm det}(A) \neq 0$

Уравнение

$\begin{aligned}
 A y = g
\end{aligned}$

Решение

$\begin{aligned}
 y = A^{-1} g
\end{aligned}$
    """
    
if menu == "Простейшие задачи":
    r"""
##### Простейшие задачи

Диагональная матрица

$\begin{aligned}
 A = {\rm diag}\{d_{1},d_{2},\ldots,d_{n}\} =  \begin{pmatrix}
  d_1 & 0 & \cdots &  0  \\
  0   & d_2 & \cdots &  0 \\
  \cdots & \cdots & \cdots &  \cdots  \\
  0 & 0 & \cdots &  d_n  \\
\end{pmatrix} 
\end{aligned}$

при $d_i \neq 0, \ i = 1,2, \ldots, n$ 

Нижняя $\big( a_{ij} = 0, i < j \big)$ треугольная матрица

$\begin{aligned}
 A = \{a_{ij}\} = \begin{pmatrix}
  a_{11} & 0      & \cdots &  0       \\
  a_{21} & a_{22} & \cdots &  0       \\
  \cdots & \cdots & \cdots &  \cdots  \\
  a_{n1} & a_{n2} & \cdots &  a_{nn}  \\
\end{pmatrix} 
\end{aligned}$

Верхняя $\big ( a_{ij} = 0, i > j \big )$ треугольная матрица

$\begin{aligned}
  A = \{a_{ij}\} = \begin{pmatrix}
  a_{11} & a_{12} & \cdots &  a_{1n}  \\
  0      & a_{22} & \cdots &  a_{2n}  \\
  \cdots & \cdots & \cdots &  \cdots  \\
  0      & 0      & \cdots &  a_{nn}  \\
\end{pmatrix}  
\end{aligned}$
 
при $a_{ii} \neq 0, \ i = 1,2, \ldots, n$ 

    """
    

    

    
