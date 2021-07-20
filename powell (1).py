#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from goldensearch import *
import math
def powell(F,x,h=0.1,tol=1.0e-6):
    def f(s): return F(x + s*v) # F in direction of v
    n = len(x) # Number of design variables
    df = np.zeros(n) # Decreases of F stored here
    u = np. identity(n) # Vectors v stored here by rows
    for j in range(100): # Allow for 30 cycles:
        xOld = x.copy() # Save starting point
        fOld = F(xOld)
# First n line searches record decreases of F
        for i in range(n):
            v = u[i]
            a,b = bracket(f,0.0,h)
            s,fMin = search(f,a,b)
            df[i] = fOld - fMin
            fOld = fMin
            x = x + s*v
# Last line search in the cycle
        v = x - xOld
        a,b = bracket(f,0.0,h)
        s,fLast = search(f,a,b)
        x = x + s*v
# Check for convergence
        if math.sqrt(np.dot(x-xOld,x-xOld)/n) < tol: return x,j+1
# Identify biggest decrease & update search directions
        iMax = np.argmax(df)
        for i in range(iMax,n-1):
            u[i] = u[i+1]
        u[n-1] = v
    print("Powell did not converge")
    


# In[ ]:




