
import numpy as np
import random

#By Muhamad Fahim Aqil bin Muhamad Sahlan
#import numpy for matrix multplication
#function HammingG
#input: a number r
#output: G, the generator matrix of the (2^r-1,2^r-r-1) Hamming code
def hammingGeneratorMatrix(r):
    n = 2**r-1
    
    #construct permutation pi
    pi = []
    for i in range(r):
        pi.append(2**(r-i-1))
    for j in range(1,r):
        for k in range(2**j+1,2**(j+1)):
            pi.append(k)

    #construct rho = pi^(-1)
    rho = []
    for i in range(n):
        rho.append(pi.index(i+1))

    #construct H'
    H = []
    for i in range(r,n):
        H.append(decimalToVector(pi[i],r))

    #construct G'
    GG = [list(i) for i in zip(*H)]
    for i in range(n-r):
        GG.append(decimalToVector(2**(n-r-i-1),n-r))

    #apply rho to get Gtranpose
    G = []
    for i in range(n):
        G.append(GG[rho[i]])

    #transpose    
    G = [list(i) for i in zip(*G)]

    return G

#generate the message of the correct length
def randomMessage(r):
    k = 2**r -r -1
    #make array containing the message
    v = []
    for i in range(k):
        v.append(random.randint(0,1))
    return v


    


#encode
def encoder(m):
#find r
    found = False
    r=0
    while(found == False):
        if( len(m) == 2**r -r - 1):
            
            found = True
        else:
            r = r + 1
#create M of size v  
    v=2**r-1
    matrix=[[0 for x in range(v)] for y in range(v)]
    for i in range(0,v):
        matrix[i][i] = 1
    for row in matrix:
        row.append(1)
#G'=GM
    result= np.dot(hammingGeneratorMatrix(r),matrix)
    result=np.mod(result,2)
    Gprime = np.mod(result,2)
#Encode m with G'   

    encode=np.matmul(m,Gprime)
    

    return np.mod(encode,2)
 
 
def BSC(c,p):
    y=c
#Create a new list L for BSC
    L=[]
#Creating BSC
    for item in y:
        number=random.random()
        if number<=p:
            if item== 0:
                L.append(1)
            else:
                L.append(0) 
        else:
             L.append(item)
    return L
            
      
            
 
        
                
def syndrome(v):
#Find r
    found=False
    r=0
    while(found==False):
        r+=1
        if 2**r==len(v):
            found=True
#Create H
    n=len(v)-1
    N = n+1
    a = np.arange(N, dtype=int)[np.newaxis,:]

    l = int(np.log2(N))
    b = np.arange(l, dtype=int)[::-1,np.newaxis]

    matrix=np.array(a & 2**b > 0, dtype=int)
    m=np.delete(matrix, np.s_[0], axis=1)
    b=np.insert(m, m.shape[1], 0, axis=1)
    L=[]
    for i in range(n+1):
        L.append(1)
        
    M = np.vstack([b,L])
#Create Htranspose    
    c=np.transpose(M)
    bsc=v
#Htranspose * bsc
    result=np.matmul(bsc,c)
    result2=np.mod(result,2)
    print("* Destination * \nDecoding by Syndrome \nSyndrome:",result2)
#Finding value of syndrome
    value=int(''.join(map(str,result2)))
    value=int(str(value), 2)
    print ("syndrome:",value)
    failure=False
#if false which is syndrome==even number, failure will occured,estimated codework is None
#if syndrome==0,it will return BSC
#if syndrome==odd,it will retrieve the estimated codework
    if value==0:
        print ("i:-")
        return bsc
        
    else:
        if value%2==0:
            print("i:-")
            failure=True
            return None
           
        else:
            syndrome=(value-1)
            syndrome=syndrome/2
            print("i:",syndrome)
            s=int(syndrome-1)
            L=[]
            count=0
            for item in bsc:
                if count==s:
                    if bsc[s]==0:
                        L.append(1)
                        count+=1
                    else:
                        L.append(0)
                        count+=1
                else:
                    L.append(item)
                    count+=1
            return L
 

def retrieveMessage(c):
#if not None, it will decode the codework
    if c != None:
#find r
        Y=c
        found=False
        r=0
        while(found==False):
            r+=1
            if 2**r==len(c):
                found=True
        L=[]
        for x in Y:
           L.append(x)
         
        for i in range(r+1):
            d=r-i
            L.pop((2**d)-1)
        return L
    else:
        return None



 

    
#function decimalToVector
#input: numbers n and r (0 <= n<2**r)
#output: a string v of r bits representing n
def decimalToVector(n,r): 
    v = []
    for s in range(r):
        v.insert(0,n%2)
        n //= 2
    return v

def simulation(r, N, p):
    errors=0
    failures=0
    success=0
    
    for i in range(N):
        print("********* Experiment ",i+1," of ", N," *********")
        print(" ")
        m = randomMessage(r)
 
 
        print("*Source* \n Message \n m=",m)
        print(" ")

        c=encoder(m)
        print("Codework \n c=",c)
        print(" ")

        v=BSC(c,p)
        print("* Channel * \n Received Vector \n v=",v)
        print(" ")

        cprime=syndrome(v)
        if cprime==None:
            print("*Failure#*")
            print(" ")
            failures+=1

        if cprime:
            print("Estimated Codework:",cprime)
            print(" ")

        mprime=retrieveMessage(cprime)
#if m' and m is not the same, it will return error else it's a success. For None,it will be failure
        if mprime==None:
            print(mprime)
            print(" ")

        if mprime!=m and mprime !=None:
            print("MessageEstimate:",mprime)
            print(" ")
            print("Error")
            errors+=1
            print(" ")

        if mprime!= None and mprime==m:
            print("MessageEstimate:",mprime)
            print(" ")
            print("Success!")
            success+=1
    

    print("***** End of Experiments ****** \n ")
#Calculate errors,failures and success
    print("Errors:",errors)
    print("Failures:",failures)
    print("Success:",success)
#Count DEP
    DEP=errors/N
    print("DEP :",DEP)

    
          
#syndrome(BSC(encoder(randomMessage(3)),0.2))


