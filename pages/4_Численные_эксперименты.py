import streamlit as st
 
menu = st.sidebar.radio('***',
    (
    "Тестовая задача", 
    "Линейные решатели",     
    "Ключевые фрагменты кода (FEniCS)", 
    "Численные результаты", 
    )
)

if menu == "Тестовая задача":
    r"""
##### Тестовая задача

**Область**

$\Omega = \{ x \ | \ x = (x_1, x_2), \ 0 < x_1 < 1, \ 0 < x_2 < 1\} $

**Уравнение**

$\begin{aligned} - 
 \operatorname{div} \operatorname{grad} u + c u = x_1 x_2 ,
 \quad x = (x_1, x_2) \in \Omega
\end{aligned}$

с постоянным коэффициентом $c$

**Граничные условия**

$\begin{aligned}
  \frac{\partial u}{\partial n} = 0,
  \quad x \in \partial\Omega
\end{aligned}$

    """     
if menu == "Линейные решатели":
    r"""
##### Линейные решатели 

Прямой метод библиотеки NumPy

* метод Гаусса для плотных матриц
* решатель numpy.linalg.solve

Прямой метод библиотеки SciPy

* метод Гаусса для разреженных матриц
* решатель scipy.sparse.linalg.spsolve

Прямой метод FEniCS (библиотека PETSc: UMFPACK)

* метод Гаусса для разреженных матриц
* решатель solve
    
    """   
    
if menu == "Ключевые фрагменты кода (FEniCS)":
    r"""
##### Ключевые фрагменты кода (FEniCS)
     
**Импорт библиотек**
    """    
    code = """  
import numpy as np 
from fenics import *
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
import time  
    """ 
    st.code(code, language="python")   
    r"""         
**Вариационная формулировка задачи**
    """    
    code = """  
f = Expression("x[0]*x[1]", degree=p+2)
q = Expression("0", degree=p+2)    
cc = Constant(c)
a = dot(grad(u), grad(v))*dx + cc*u*v*dx
L = f*v*dx 
    """ 
    st.code(code, language="python")   
    r"""        
**Матрицы**
    """    
    code = """  
A, b = assemble_system(a, L)
mat = as_backend_type(A).mat()
As = csr_matrix(mat.getValuesCSR()[::-1], shape=mat.size)
An = As.toarray() 
    """ 
    st.code(code, language="python")  
        
    r"""        
**Решение задачи (NumPy)**
    """    
    code = """  
start_time = time.clock()
wn = Function(V)
yn = np.linalg.solve(An, b) 
wn.vector().set_local(yn)
tn = time.clock() - start_time   
    """ 
    st.code(code, language="python") 
             
    r"""           
**Решение задачи (SciPy)**
    """    
    code = """  
start_time = time.clock()
ws = Function(V)
ys = spsolve(As, b)
ws.vector().set_local(ys)
tn = time.clock() - start_time   
    """ 
    st.code(code, language="python")  
         
    r"""            
**Решение задачи (FEniCS)**
    """    
    code = """  
start_time = time.clock()
w = Function(V)
solve(a == L, w, solver_parameters={"linear_solver": "default", "preconditioner":"default"})
tf = time.clock() - start_time
    """ 
    st.code(code, language="python")           
                   
if menu == "Численные результаты":
    r"""
##### Численные результаты

    """

    import matplotlib.pyplot as plt
    import numpy as np 
    from fenics import *
    from scipy.sparse import csr_matrix
    from scipy.sparse.linalg import spsolve
    import pandas as pd
    import time
    time.clock = time.time
    
    c1, cc, c2 = st.columns([5,1,5])
    c1.write("$~~$")
    c1.write("Число ячеек $m$")        
    m = c2.slider("", 10, 50, 20)
    c1, cc, c2 = st.columns([5,1,5])
    c1.write("$~~$")
    c1.write("Порядок полинома $p$")        
    p = c2.slider("", 1, 3, 1)  
    c1, cc, c2 = st.columns([5,1,5])
    c1.write("$~~$")
    c1.write("Коэффициент $c$")        
    c = c2.slider("", -10, 10, 1, 1)     
      
    mesh = UnitSquareMesh(m, m)
    xm = mesh.coordinates()
    ym = np.zeros((m+1), "float") 
    
    V = FunctionSpace(mesh, "CG", p)
    n = V.dim()-1
    
    u = TrialFunction(V)
    v = TestFunction(V)
    
    f = Expression("x[0]*x[1]", degree=p+2)
    q = Expression("0", degree=p+2)    
    cc = Constant(c)
    a = dot(grad(u), grad(v))*dx + cc*u*v*dx
    L = f*v*dx 
    
    A, b = assemble_system(a, L)
    mat = as_backend_type(A).mat()
    As = csr_matrix(mat.getValuesCSR()[::-1], shape=mat.size)
    An = As.toarray()
    
    r"""
**Портрет матрицы**
    """    
          
    fig1 = plt.figure(1)
    plt.spy(As)        
 
    c1, c2, = st.columns([3,1]) 
    c1.pyplot(fig1) 
    
    start_time = time.clock()
    wn = Function(V)
    yn = np.linalg.solve(An, b) 
    wn.vector().set_local(yn)
    tn = time.clock() - start_time
    
    start_time = time.clock()
    ws = Function(V)
    ys = spsolve(As, b)
    ws.vector().set_local(ys)
    ts = time.clock() - start_time
    
    start_time = time.clock()
    w = Function(V)
    solve(a == L, w, solver_parameters={"linear_solver": "default", "preconditioner":"default"})
    tf = time.clock() - start_time

    r"""
**Решение**
    """    
    
    N = 200
    x = np.linspace(0,1,N)
    y = np.linspace(0,1,N)
    yy  = np.zeros((N,N)) 
    ye  = np.zeros((N,N)) 
    tk = np.linspace(0,1,m+1)
    
    for i in range(0, N): 
        for j in range(0, N): 
            pp = Point(x[i],y[j])
            yy[i,j] = w(pp)
              
    fig2 = plt.figure(2)
    ss = "$m = $" + str(m) + "$, \ p = $" + str(p) + "$, \ c = $" + str(c)
    plt.title(ss)
    plt.contourf(x,y,yy) 
    plt.gca().set_aspect("equal")
    plt.colorbar()
    
    c1, c2, = st.columns([3,1]) 
    c1.pyplot(fig2)  

    r"""
**Время решения**
    """   
    
    df = pd.DataFrame(
       [
       {"Метод": "NumPy", "Время (с)": tn},
       {"Метод": "SciPy", "Время (с)": ts},
       {"Метод": "FEniCS", "Время (с)": tf},
       ])

    st.table(df) 
    
 
    



























