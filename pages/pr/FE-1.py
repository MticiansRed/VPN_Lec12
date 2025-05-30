import matplotlib.pyplot as plt
import numpy as np 
from fenics import *
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
import time
import pandas as pd

m = 50
p = 1
c = 1

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
fig1 = plt.figure(1)
plt.spy(As) 

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


plt.show()
