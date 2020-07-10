import numpy as np 
import math
from numpy.linalg import lstsq
from scipy.linalg import orth
import cmath
import scipy
from sympy import Matrix
from scipy.linalg import qr
from scipy.linalg import null_space

def decimalToBinary(n):  
    return bin(n).replace("0b", "")  

# Function to convert Binary number 
# to Decimal number 
  
def binDec(array):
    sum=0
    for i in range(len(array)):
        sum=sum+2**(len(array)-i-1)*array[i]
    return sum
    

def GSO(m):
    mystr=decimalToBinary(m)
    arr=(np.array(list(mystr), dtype=int))
    
    #for i in range(int(2**N)):
#GSO(7)
#v=np.array([])
#v=np.append(v,)
#ini_array=np.array([[0,0,0],[1,1,1]])
#esult = np.vstack((ini_array, np.array([1,2,3]))) 
#print(result[1:,:])
#lis=np.array([1,2])
#print(lis)
m=np.array([1,2,3,4])
print(m[-4])