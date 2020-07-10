#header files
import numpy as np 
import math

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


#printing all the states for a given S value
S_in=float(input())
if(S_in>0):
    m_in=np.arange(-S_in,S_in+1)

#this routine will output the corresponding projection matrix given the user input on S.
#ones represent up spin and zeros down spin in the (m1,m2...) basis

#we start with the (S_max,m_max) state, since we know this state in the (m1,m2...) basis
global S_max
S_max=np.max(S_arr)


#number of rows of the matrix
r=2**N
#number of columns
c=(2**N)+2
#the matrix which contains all the states and their corresponding representations in the (m1,m2...) basis
M=np.zeros((r,c))

#filling the first two columns of the matrix with the S and M values. The other columns will contain the projection of this S,m vector along each of the (m1,m2...) bases.
for i in range(2*S_max+1):
    M[i][0]=S_max
    M[i][1]=S_max-i

#filling the state corresponding to the S_max, m_max state
M[0][c-1]=1

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
    s=arr[0]
    m=arr[1]
    n=math.sqrt((s+m)*(s-m+1))
    #psi=np.zeros((c-2))
    psi=arr[2:]
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
               state=bin_state
            else: 
               state=np.concatenate(np.zeros(N-len(bin_state)),bin_state)
            for l in range(N):
                test=state.copy()
                if(test[l]==1):
                    test[l]=0
                else:
                    continue
                num=binDec(test)
                out[num]=out[num]+cop
    return np.true_divide(out,n)
                
                
#filling the other rows by calling a routine which takes in the previous row as its input
for i in range(1,2*S_max+1):
    M[i,2:]=S_minus(M[i-1,:])

print(M)

