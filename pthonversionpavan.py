#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 13:47:18 2023

@author: titouan
"""
from math import*
import numpy as np
import matplotlib.pyplot as plt
M_soleil=1.89*10**30
M_terre=5.972*10**24
M_mars= 6.39*10**23

h_orbite_mars=400*10**3
h_orbite_geo=42157*10**3

d_terre_soleil=149597870.7*10**3
d_mars_soleil=227939200*10**3

periode_Terre=365
periode_Mars=687
vt=30000

G=6.67*10**(-11)

def posTerre(t):                                               #position de la Terre au temps t 
    x_T=d_terre_soleil*cos((2*(np.pi))/(365*3600*24)*t)
    y_T=d_terre_soleil*sin((2*(np.pi))/(365*3600*24)*t)
    
    return(x_T,y_T)

def posMars(t,Theta):                                           #position de Mars au temps t avec un angle Theta varible par rapport à la Terre
    x_M=d_mars_soleil*cos(((2*(np.pi))/(687*3600*24)*t)+Theta)
    y_M=d_mars_soleil*sin(((2*(np.pi))/(687*3600*24)*t)+Theta)
    
    return(x_M,y_M)

def accel(x,y,t,theta):                                         #accélération de la fusée à une position (x,y) dans l'espace au temps t avec un déphasage de Théta entre la Terre et Mars

     ax= -G*M_soleil*x/(x**2+y**2)**(3/2)-G*M_terre*(x-posTerre(t)[0])/((x-posTerre(t)[0])**2+(y-posTerre(t)[1])**2)**(3/2)-G*M_mars*(x-posMars(t,theta)[0])/((x-posMars(t,theta)[0])**2+(y-posMars(t,theta)[1])**2)**(3/2)
     ay= -G*M_soleil*y/(x**2+y**2)**(3/2)-G*M_terre*(y-posTerre(t)[1])/((x-posTerre(t)[0])**2+(y-posTerre(t)[1])**2)**(3/2)-G*M_mars*(y-posMars(t,theta)[1])/((x-posMars(t,theta)[0])**2+(y-posMars(t,theta)[1])**2)**(3/2)
    
     return(ax,ay)

def F(U,t,theta):
    
    return np.array((U[2],U[3],accel(U[0],U[1],t,theta)[0],accel(U[0],U[1],t,theta)[1])  )
  

def primitive(U,t,dt,theta):                                    #fonction permettant d'utiliser la méthode d'Euler pour obtenir les primitives
    return U+F(U,t,theta)*dt
       
    



class orbit():
    
    def __init__(self,x0,y0,Vx0,Vy0):
        
        self.U0=np.array((x0,y0,Vx0,Vy0))
        self.x0=x0
        self.y0=y0
        
    def trajectoire(self,dt,theta):
        import matplotlib.pyplot as plt
        au=1.5e11
        
        u=self.U0
        x=self.x0
        y=self.y0
        N=0
        X=[self.x0]
        Y=[self.y0]
        T=[0]
       
        while(N<1000):
            
            X=X+[primitive(u,N*dt,dt,theta)[0]]
            Y=Y+[primitive(u,N*dt,dt,theta)[1]]
            u=primitive(u,N*dt,dt,theta)
            x=primitive(u,N*dt,dt,theta)[0]
            y=primitive(u,N*dt,dt,theta)[1]
           
            T=T+[N*dt]
            if  ((x**2+y**2)**(1/2)>d_mars_soleil+100000*10**3): #si la distance fusée-soleil est plus grande que la distance Mars-Soleil + 100 000km
                XX=u[0]/au
                YY=u[1]/au
                #plt.scatter(X,Y)
                #plt.scatter(T,X)
                #plt.scatter(T,Y)
                return(["air ball"])
            if ((posMars(N*dt,theta)[0]-x)**2+(posMars(N*dt,theta)[1]-y)**2)**(1/2)<500*10**3 and ((posMars(N*dt,theta)[0]-x)**2+(posMars(N*dt,theta)[1]-y)**2)**(1/2)>300000: #si la fusée est entre 300 et 500 km de Mars 
                
                XX=u[0]/au
                YY=u[1]/au
                #plt.scatter(XX,YY)
                return(["3 pts",X,Y,N*dt])
            N=N+1
        XX=u[0]/au
        YY=u[1]/au
        #plt.scatter(XX,YY) 
        #plt.show()
        return(["more testing needed"])

def vitesseinitiale(n,m):             #fonction nous permettant de trouver pour quelles vitesses et quels angles de déphasages la fusée s'approche de Mars 
    Dv=np.linspace(11000,7*11000,n)   #on teste les vitesses entre 11 km/s et 77 km/s avec un pas de n      
    Theta=np.linspace(0,2*np.pi,m)    #on teste les angles entre 0 et 2π avec un pas de m
    y=[]
    x=[]
    d=d_terre_soleil
    aimsu=[]
    for i in (Dv):                                     #on teste chaque angle pour chaque vitesse 
        for j in (Theta):
            a=orbit(d_terre_soleil,42000*10**3,0,vt+i) #on appelle la classe "orbit" en fasant partir la fusée de l'orbite géostationnaire avec une vitesse initilale purement tangentielle 
            b=a.trajectoire(86400,j)                   #on calcule la trajectoire avec dt=1jour et pour chaque angle theta
            if len(b)==4:                              # si cette condition est vérifiée, notre fusée est suffisament près de Mars
                x=x+a.trajectoire(86400,j)[1]
                y=y+a.trajectoire(86400,j)[2]
                aimsu=aimsu+[[i,j]]                    #on note alors dans aimsu le couple [vitesse, angle] 
    print(aimsu)
    return(x,y,aimsu)
"""

def calcul_dv_Mars(v):
    dv_mars=abs((G*M_mars/3389*10**3)**(1/2)-v) #calcul du Dv sur Mars à une altitude de 400km
    return(dv_mars)

def calcul_masse_carburant(v1,v2):
    m1 =m*(1-exp(-v1/(g*500)))
    m2 =m*(1-exp(-v2/(g*500))) #avec ISP=500 pour les moteurs de fusée ; source: Wikipédia
    
    return(m1+m2)
    
    
def trajet_econome(aimsu):
    a=len(aimsu)
    b=calcul_masse_carburant()
    for i in range (a):
        aimsu[i][0]
    
"""
    
for i in range(0,1000,1):

    plt.plot(posTerre(i), color = 'r' )
    plt.plot(posMars(i , 0,0),color ='b')

plt.show()    
    
    
