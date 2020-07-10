import numpy as np
import scipy
from scipy.sparse import diags
import math

#putting the bunch of linear differential equations in a matrix form. As we know, this matrix
#will be in the tridiagonal form.
#For the first case, we keep both the atom and the photon contribution to the coefficients. Then, we drop
#the photon contribution to look for some symmetry in the matrix.
#all constants set to 1
#getting the S value as the input
S=float(input())
m=S
diag1=np.array([])
while(m>-S):
    ele=math.sqrt(S*(S+1)-m*(m-1))*math.sqrt(S-m+1)
    diag1=np.append(diag1,ele)
    m=m-1

diag2=np.array([])
p=S-1
while(p>=-S):
    ele=math.sqrt(S*(S+1)-p*(p+1))*math.sqrt(S-p)
    diag2=np.append(diag2,ele)
    p=p-1

diagonals = [diag2,diag1] 
A =scipy.sparse.diags(diagonals, [-1,1], format='csc')
B=A.toarray()

print(B)

m=S
diag1=np.array([])
while(m>-S):
    ele=math.sqrt(S*(S+1)-m*(m-1))
    diag1=np.append(diag1,ele)
    m=m-1

diag2=np.array([])
p=S-1
while(p>=-S):
    ele=math.sqrt(S*(S+1)-p*(p+1))
    diag2=np.append(diag2,ele)
    p=p-1

diagonals = [diag2,diag1] 
A =scipy.sparse.diags(diagonals, [-1,1], format='csc')
C=A.toarray()

print(C)