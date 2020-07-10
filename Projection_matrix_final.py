#header files
import numpy as np 
import math
from numpy.linalg import lstsq
from scipy.linalg import orth

#getting the number of atoms as the input
global N
N=int(input()) 

#routine to get the S values. S goes from |U-D|/2,......,(U+D)/2, where U and D are the number of upspins and downspins respectively
#when both U and D are equal (if odd, then differ by 1), we would get the full spectrum of values for S
global S
if(N%2==0):
    S_arr=np.arange((N/2)+1)
else:
    S_arr=np.arange(1/2,(N/2)+1,1)

#we start with the (S_max,m_max) state, since we know this state in the (m1,m2...) basis
global S_max
S_max=np.max(S_arr)

# Function to convert Decimal number  
# to Binary number  
def decimalToBinary(n):  
    return bin(n).replace("0b", "")  

# Function to convert Binary number 
# to Decimal number 
  
def binDec(array):
    sum=0
    for i in range(len(array)):
        sum=sum+2**(len(array)-i-1)*array[i]
    return sum
    

#function which performs the S_minus operation
def S_minus(arr):
    #print(arr)
    s=arr[0]
    m=arr[1]
    n=math.sqrt((s+m)*(s-m+1))
    #psi=np.zeros((c-2))
    psi=arr[2:]
    #print(psi)
    out=np.zeros(len(psi))
    for i in range(len(psi)):
        if(psi[i]==0):
            continue
        else:
            cop=psi[i]
            bino=decimalToBinary(i)
            bin_state=np.array(list(bino), dtype=int)
            #print(np.zeros(N+1-len(bin_state)))
            if(N-len(bin_state)==0):
                #print(N-len(bin_state))
                state=bin_state
            else: 
                #print(N-len(bin_state))
                #print(bin_state)
                state=np.concatenate([np.zeros(N-len(bin_state)),bin_state])
                #print(state)
            for l in range(N):
                test=state.copy()
                if(test[l]==1):
                    test[l]=0
                    num=int(binDec(test))
                    out[num]=out[num]+cop
                else:
                    continue
                
    return np.true_divide(out,n)


#generating the complete M matrix for other values of S also. Call this matrix P.
#initialising it with zeros of required size
#number of rows of the matrix
r=2**N
#number of columns
c=(2**N)+2
P=np.zeros(int(r),int(c))

i=0
count_1=0
count=0
#initialising the matrix with the first two columns as s and m.
while(S_max-i>=0):
    count_1=count_1+count
    count=0
    S=S_max-i
    for j in range(count_1,count_1+(int(2*S)+1)
        count=count+1
        P[j][0]=S
        P[j][1]=S-j+count_1
    i=i+1
    
#filling in the matrix with the states

P[0][int(N**2+1)]=1
for i in range(1,int(P.shape[0])):
    if(P[i][0]==P[i][1]):
    else:
        P[i,2:]=S_minus(P[i-1,:])
    

#getting the projection matrices
sin=input()
if(math.modf(sin)[0]!=0):
    bin=1
else:
    bin=2

for i in range(int(P.shape[0])):
    if(P[i][0]==sin):
        vmatrix=P[i:i+int(2*sin+bin)][:]
        break

dim=N**2
proj=np.zeros((dim,dim))
i=0
while(vmatrix[i][0]==sin):
    proj=proj+np.outer(vmatrix[i,2:],vmatrix[i,2:])
    i=i+1
#the final answer!
print(proj)