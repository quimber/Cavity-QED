#header files
import numpy as np 
import math
from numpy.linalg import lstsq
from scipy.linalg import orth
import cmath
import scipy
from sympy import Matrix
from scipy.linalg import qr
from scipy.linalg import null_space

#algo
#getting inputs (N and S)
#a 3d matrix with the first matrix containing the vectors for the maximum S and so on.
#find the first matrix using the S- method.
#find the second matrix using the nth roots method.
#start the loop for higher S. in each run of the loop, find the degeneracy for the maximum element (S,S) using an array initialised with all m values and incemented as and when the m value is added. From the previous matrices, find the elements with the m value as S
#send this vector for GSO. Then, apply S_minus to get the matrix elements. 

#pseudocode  
#input (N and S)
#if(N%2==0):
#    S_arr=np.arange((N/2)+1)
#else:
#    S_arr=np.arange(1/2,(N/2)+1,1)

#we start with the (S_max,m_max) state, since we know this state in the (m1,m2...) basis
#global S_max
#S_max=np.max(S_arr)
#all zeros, S_max in height, N**2 rows, N**2+2 columns
#S[0][0][0]=S[0][0][1]=S_max
# #S[0][0][N**2+1]=1 
#for i in range(1,2*S_max+1):
#S[0][i][0]=S_max
#S[0][i][1]=s_max-i
#S[0,i,2:]=S_minus(S[0,i-1,:])
#m=ones(length of (2*Smax+1)/2) if integer, length of S_max:
#m_max=zeros(length of (2*Smax+1)/2) if integer, length of S_max
#m_max[i]=N C (N/2-i) if N is even
#if N is odd, m_max[i]=N C int (N/2-(2*i+1)/2)

#omega=cos(2*pi/5)+i sin(2*pi/5)
#l=0
#p=m_max[N-1]-m[N-1]
#while(m[N-1]<=m_max[N-1]):
#S[1][l][0]=S[1][l][1]=S_max-1
#l=l+1
#omega=omega**l
#S[1,l,2:]=S[0,1,2:]*omega**0 till N-1
#m[N-1]=m[N-1]+1
#all mi's except the last =m[N-1]
#complete this matrix for now

#for later matrices:
#S=S_max-2
#i=2
#while(S>=0):
#degen=m_max[N-i]-m[N-i]

#fn that returns all the vectors as a matrix, with the same m
#for i in range(len(S_array)):
#vec_mat=zeros of specified length
#while(S[i][p][0]==S[i][0][0] and S[i][p][0]!=np.zeros of specified length):
#if(S[i][p][1]==m):
#add S[i,p,2:] as a new row to vec_mat

#p=p+1
#
#p=0
#while(m_max[N-i]-m[N-i]>=0):
#S[i][p][0]=S[i][p][1]=S
#S[i,p,2:]=GSO(that matrix i mentioned above)
#m=S.copy()
#j=1

#while(S-j>=-S):
#S[i][p+j*degen][0]=S
#S[i][p+j*degen][1]=S-j
#S[i,p+j*degen,2:]=S_minus(S[i,p+(j-1)*degen,:])
#j=j+1


#GSO(append this)
#p=p+1
#m[N-i]=m[N-i]+1

#
#S=S-1
#all m is same as m[N-i] 
#i=i+1

#for i in range(len(S_array)):
#if(S[i][0][0]==S):
#mat=zeros of specified size (N**2 * N**2)
#p=0
#while(S[i][p][0]=S and S[i,p.2:]!=np.zeros):
#mat=mat+np.outer(S[i,p,2:],S[i,p,2:])
#p=p+1


#final matrix: print(mat)

#FINAL CODE
#functions
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
    

def GSO(m):
    pos=np.array([])
    for i in range(int(2**N)):
        mystr=decimalToBinary(i)
        arr=np.array(list(mystr), dtype=int)
        if(N-len(arr)==0):
            #print(N-len(bin_state))
            state=arr
        else: 
        #print(N-len(bin_state))
        #print(bin_state)
            state=np.concatenate([np.zeros(N-len(arr)),arr])
        zeros=0
        ones=0
        for j in range(len(state)):
            if(state[j]==1):
                ones=ones+1
            if(state[j]==0):
                zeros=zeros+1
        if((ones-zeros)*0.5==m):
            pos=np.append(pos,[i])
    O=np.zeros(int(2**N))
    for j in range(len(S_arr)):
        for i in range(int(2**N)):
            if((S[j][i][0]==0) and S[j][i][1]==-1):
                break
            else:
                if(S[j][i][1]==m):
                    O=np.vstack((O, S[j,i,2:])) 
    pos=pos.astype(int)
    pos1=pos.tolist()
    O=O[1:,:]
    #print(pos1)
    O=O[:,pos1]
    #print(O)
    O=O[np.all(O!=0, axis=1)]
    #print(O)
    ns=null_space(O)
    vec=ns[:,0]
    vector=np.zeros(int(2**N),dtype=np.complex_)
    l=0
    for i in range(len(pos)):
        vector[int(pos[i])]=vec[l]
        l=l+1
    return(vector)


#function which performs the S_minus operation
def S_minus(arr):
    #print(arr)
    s=arr[0]
    m=arr[1]
    if((s+m)*(s-m+1)<0):
        print('error')
    n=math.sqrt((s+m)*(s-m+1))
    #psi=np.zeros((c-2))
    psi=arr[2:]
    #print(psi)
    out=np.zeros(len(psi),dtype=np.complex_)
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
                
    #print(np.true_divide(out,n))
    return np.true_divide(out,n)
def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)

global N
global S_val
N=int(input())
S_val=float(input())

global S_arr
if(N%2==0):
    S_arr=np.arange((N/2)+1)
else:
    S_arr=np.arange(1/2,(N/2)+1,1)

#print(S_arr)
global S_max
S_max=np.max(S_arr)

#3D matrix
global S
S=np.zeros((len(S_arr),int(2**N),int(2**N+2)),dtype=np.complex_)
for j in range(len(S_arr)):
    for i in range(int(2**N)):
        S[j][i][1]=-1
#initialising the max state
S[0][0][0]=S_max
S[0][0][1]=S_max
S[0][0][int(2**N+1)]=1

#first matrix
for i in range(1,int(2*S_max+1)):
    S[0][i][0]=S_max
    S[0][i][1]=S_max-i
    S[0,i,2:]=S_minus(S[0,i-1,:])


#m number counter
m=np.ones(len(S_arr))
m_max=np.zeros(len(S_arr))
if(N%2==0):
    for i in range(len(m_max)):
        m_max[i]=nCr(int(N),int(N/2-i))
else:
    for i in range(len(m_max)):
        m_max[i]=nCr(int(N),int(N/2-((2*i+1)/2)))

pie=math.pi
omega=math.cos(2*(pie/N))+1j*math.sin(2*(pie/N))
#print(omega)
l=0
#print(S[0,1,:])
#print(m[-2])
#print(m_max[-2])
deg=(m_max[-2]-m[-2]).copy()
while(m[-2]<m_max[-2]):
    m[-2]=m[-2]+1
    S[1][l][0]=S_max-1
    S[1][l][1]=S_max-1
    om=omega**(l+1)
    S[1,l,2:]=S[0,1,2:].copy()
    count=0
    for j in range(2,len(S[1,l,:])):
        if(S[1,l,j]!=0):
            S[1,l,j]=S[1,l,j]*(om**count)
            count=count+1
    S1=S_max-2
    j1=1
    while(S1>=-S_max+1):
        c=int(l+j1*deg)
        S[1][c][0]=S_max-1
        S[1][c][1]=S1
        #print(S[1][int(l+(j1-1)*deg)][0])
        #print(S[1][int(l+(j1-1)*deg)][1])
        S[1,c,2:]=S_minus(S[1,int(l+(j1-1)*deg),:])
        S1=S1-1
        j1=j1+1


    l=l+1

val=m[-2].copy()
for i in range(len(m)-1):
    m[i]=val

#second matrix over
#print(S[0:2,:,:])
#print(m)
#print(m_max)
#for later matrices
i=2
S2=S_max-i
while(S2>=0):
    #print(S2)
    l=0
    deg=(m_max[-i-1]-m[-i-1]).copy()
    while(m[-i-1]<m_max[-i-1]):
        m[-i-1]=m[-i-1]+1
        S[i][l][0]=S2
        S[i][l][1]=S2
        S[i,l,2:]=GSO(S2)
        S1=S_max-i-1
        j1=1
        #print(S[i][0][0])
        #print(S[i][0][1])
        while(S1>=-S_max+i):
            #print(deg)
            c=int(l+j1*deg)
            S[i][c][0]=S_max-i
            S[i][c][1]=S1
            #print(int(l+(j1-1)*deg))
            #print(S[i][int(l+(j1-1)*deg)][0])
            #print(S[i][int(l+(j1-1)*deg)][1])
            S[i,c,2:]=S_minus(S[i,int(l+(j1-1)*deg),:])
            S1=S1-1
            j1=j1+1
        l=l+1
    val=m[-i-1].copy()
    for j in range(len(m)-i):
        m[j]=val
    i=i+1
    S2=S_max-i
        

#print(m)
#print(m_max)
#the final matrix
print(S[2,:,:])

#projection_matrix
dimension=int(2**N)
mat=np.zeros((dimension,dimension),dtype=np.complex_)
for i in range(len(S_arr)):

    if(S[i][0][0]==S_val):
        p=0
        while((S[i][p][1]!=-1)):
            mat=mat+np.outer(S[i,p,2:],S[i,p,2:])
            p=p+1
        break
            
#print(mat)

