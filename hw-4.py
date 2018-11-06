
# coding: utf-8

# In[ ]:


get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import numpy as np


# In[ ]:


def func(x):
    return np.exp(-2*x)*np.cos(10*x)


# In[ ]:


def trapezoid_core(f,x,h):
    return h*(f(x+h) + f(x))/2


# In[ ]:


def trapezoid_method(f,a,b,N):
    x = np.linspace(a,b,N)
    h = x[1] - x[0] 
    Fint = 0.0
    for i in range(0,len(x)-1,1):
        Fint += trapezoid_core(f,x[i],h)
    return Fint


# In[ ]:


print(trapezoid_method(func,0,np.pi,100))


# In[ ]:


def simpson_core(f,x,h):
    return h*( f(x) + 4*f(x+h) + f(x+2*h))/3.


# In[ ]:


def simpsons_method(f,a,b,N):
    x = np.linspace(a,b,N)
    h = x[1]-x[0]
    Fint = 0.0
    for i in range(0,len(x)-2,2):
        Fint += simpson_core(f,x[i],h)
    if((N%2)==0):
        Fint += simpson_core(f,x[-2],0.5*h)
    return Fint


# In[ ]:


print(simpsons_method(func,0,np.pi,100))


# In[ ]:


def romberg_core(f,a,b,i):
    h = b-a
    dh = h/2.**(i)
    K  = h/2.**(i+1)
    M = 0.0
    for j in range(2**i):
        M += f(a + 0.5*dh + j*dh)
    return K*M


# In[ ]:


def romberg_int(f,a,b,tol):
    i = 0
    imax = 1000
    delta = 100.0*np.fabs(tol)
    I = np.zeros(imax,dtype=float)
    I[0] = 0.5*(b-a)*(f(a) + f(b))
    i += 1
    while(delta>tol):
        I[i] = 0.5*I[i-1] + romberg_core(f,a,b,i)
        delta = np.fabs((I[i]-I[i-1])/I[i])
        print(i,I[i],I[i-1],delta)
        if(delta>tol):
            i += 1
            if(i>imax):
                print("Max iterations reached.")
                raise StopIteration('Stopping iterations after ',i)
    return I[i]


# In[ ]:


tolerance = 1.0e-6
print(romberg_int(func,0,np.pi,tolerance))

