'''
Date        : 31/10/2018
Description : FFT Library for micropython devices (nodeMCU,STM32F405xx,...)
Autor       : Yassine OUAISSA
'''

import math as m
import socket as sock

FEN_TYPES = ["HANN","HAMM","RECTONGLE","TRIANGLE","FLAT_TOP","WELCH"] 

def swap(a,b):
    a,b = b,a
    return a,b  

def complexToMagnitude(x,y):
    if len(x) != len(y):
        raise(" ERROR lenth ")
    for i in range(len(x)):
        x[i] = m.sqrt(m.pow(x[i],2) + m.pow(y[i],2))

def Windowing(data,windwoType):
    samplesMinusOne = (float(len(data))-1.0)
    
    for i in range(len(data)>>1):
        indexMinusOne = float(i)
        ratio = (indexMinusOne/samplesMinusOne)
        weighingFactor = 1.0
        
        if windwoType == "HAMM":
            weighingFactor = 0.54 -(0.46 * m.cos(2*m.pi * ratio))
        if windwoType == "HANN":
            weighingFactor = 0.54 -(1.0 * m.cos(2*m.pi * ratio)) 
        if windwoType == "RECTONGLE":
            weighingFactor = 1.0
        if windwoType == "TRIANGLE":
            weighingFactor = 1.0 - ((2.0 * m.fabs(indexMinusOne - (samplesMinusOne/2.0)))/samplesMinusOne)       
        if windwoType == "BLACKMAN":
            weighingFactor = 0.4232 - (0.49755 * (m.cos(2*m.pi *ratio))) + (0.07922 * (m.cos(4*m.pi * ratio)))
        if windwoType == "FLAT_TOP":
            weighingFactor = 0.2810639 - (0.5208972 * m.cos(2*m.pi *ratio)) + (0.1980399 * m.cos(4*m.pi *ratio)) 
        if windwoType == "WELCH":
            weighingFactor = 1.0 - m.pow(( indexMinusOne - samplesMinusOne/2.0)/(samplesMinusOne /2.0),2)
        
        data[i] *= weighingFactor
        data[len(data) -(i+1)] *= weighingFactor
        
    
def compute(Vreal,Vimag):
    j=0
    power = m.log(len(Vreal),2)
   #print(power)
    for i in range(len(Vreal)-1):
        if i<j:
            #print(Vreal[i],Vreal[j])
            Vreal[i],Vreal[j] = swap(Vreal[i],Vreal[j])
           # print(Vreal[i],Vreal[j])
            #print("----------------------------------------------------------------------")
        k = len(Vreal) >> 1
        
        while k <=j :
            j-=k
            k = k >> 1
            #print(j)
        j += k
    c1 = -1.0
    c2 = 0.0
    l2 = 1
    for i in range(int(power) -1):
        l1 = l2
        l2 = l2 << 1
        u1 = 1.0
        u2 = 0.0
        for j in range(l1):
            for i in range(j,len(Vreal),l2):
                
                i1 = i + l1
                #print(Vreal[i1],"Debut ")
                t1 = u1 * Vreal[i1] - u2 * Vimag[i1]
                t2 = u1 * Vimag[i1] + u2 *Vreal[i1]
                
                Vreal[i1] = Vreal[i] - t1
                Vimag[i1] = Vimag[i] - t2
                
                Vreal[i] += t1
                Vimag[i] += t2
                #print(Vreal[i1],"fin ")
            z = ((u1*c1) -(u2 *c2))
            u2 = ((u1*c2)+(u2 *c1))
            u1 = z
            
        c2 = m.sqrt((1.0-c1)/2.0)
        c2 = -c2
        c1 =m.sqrt((1.0 + c1)/2.0)
    #print(Vreal)
    return Vreal, Vimag

def FFT_Pack(rDATA):
    
    iDATA = []
    for i in range(len(rDATA)):
        iDATA.append(0.0)
    compute(rDATA,iDATA)
    complexToMagnitude(rDATA,iDATA)
    for fen in FEN_TYPES:
        Windowing(rDATA,fen)


    
