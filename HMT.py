#!/usr/bin/env python
# coding: utf-8

# written by Davide Lofano, November 1, 2021

# Davide Lofano, and Frank H. Lutz.
# Hadamard Matrix Torsion
# Preprint, 14 pages, 2021; arXiv:2109.13052. 


# # Program to create triangulations of the Hadamard matrix torsion simplicial complexes HMT(n)


# Import library to use Hadamard matrices

from scipy.linalg import hadamard


# ### Compute a valid sequence


# Choose m as a power of 2 and a corresponding starting valid sequence for the induction procedure

m=1;
StartSequence=[[1]]


# Choose n as a power of 2 to get the corresponding valid sequence and simplicial complex HMT(n).
# (If the chosen number is not a power of 2 the program will fail.)

n=4


# Recursive function to compute a valid sequence for n from a StartSequence for m

def ValSequence(order):
    if (order==m):
        return StartSequence
    else:
        n=int(order/2)
        Sequence=ValSequence(n)
        NewSequence=[]
        SequencePlusN = [[x+n for x in perm]for perm in Sequence]
                
        for i in range(order):
            if i<n:
                NewPermutation=Sequence[i]+SequencePlusN[i]
                NewSequence.append(NewPermutation)  
            else:
                NewPermutation= [None]*order
                for j in range(n):
                    if (j%2==0):
                        NewPermutation[j]=Sequence[i-n][j]
                        NewPermutation[j+n]=SequencePlusN[i-n][j]
                    else:
                        NewPermutation[j]=SequencePlusN[i-n][j]
                        NewPermutation[j+n]=Sequence[i-n][j]                        
                NewSequence.append(NewPermutation)
        
        
        return NewSequence
    

# Run the recursive function and get a valid sequence for n

FinalSequence=ValSequence(n)

print(FinalSequence)


# ### Use the resulting sequence to triangulate HMT(n)


H=hadamard(n)
SC=[]


# Create n-1 digons

for i in range(2,n+1):
    SC.extend(([0,4*i-3,4*i],[0,4*i-2,4*i-1],[4*i-3,4*i-2,4*i-1],[4*i-3,4*i-1,4*i]))


# Create the polygonal discs correspondinng to the valid sequence

for i in range(n):
    for j in range(n):
        if (H[i][FinalSequence[i][j]-1]>0):
            w=4*FinalSequence[i][j]
        else:
            w=4*FinalSequence[i][j]-2
        if (H[i][(FinalSequence[i][(j+1)%n])-1]>0):
            z=4*FinalSequence[i][(j+1)%n]-1
        else:
            z=4*FinalSequence[i][(j+1)%n]-3  
            
        SC.extend(([w-1,w,4*n+i+1],[w,z,4*n+i+1],[0,z,w]))       
    

# Change labeling of the vertices

SC= [[x-2 if x>2 else x for x in tri] for tri in SC]

    
# Print the complex

print(SC)
