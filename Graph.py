import numpy as np
from numpy.linalg import inv
from numpy.linalg import det
from math import pi, floor
import scipy.optimize as so
import matplotlib.pyplot as plt
import os

## Classe modélisant un graphe :

class Graph:
    
    def __init__(self):
        self.vertices = []
        self.edges = {}
        self.E = 0
        self.L = 0
    
    # Rajoute un nouveau sommet au graphe :
    def addNewVertice(self):
        label = len(self.vertices)
        self.vertices.append(label)
        self.edges[label] = []
        return label
    
    # Rajoute une arrête au graphe :
    def addEdge(self,v1,v2,d=1):
        if v1 in self.vertices and v2 in self.vertices:
            if v1==v2:
                self.edges[v1].append((v2,d))
                self.E += 1
                self.L += d
                return
            self.edges[v1].append((v2,d))
            self.edges[v2].append((v1,d))
            self.E += 1
            self.L += d
        elif v1 not in self.vertices:
            print(str(v1)+" not found in the graph!")
        else:
            print(str(v2)+" not found in the graph!")
            
    # Normalise à 1 la longueur totale du graphe :
    def normalize(self):
        n = len(self.vertices)
        L = sum([sum([w[1] for w in self.edges[v]]) for v in self.vertices])/2.0
        for v in self.vertices:
            for itt in range(len(self.edges[v])):
                w = self.edges[v].pop(0)
                self.edges[v].append((w[0],w[1]/L))
        self.L = 1
    
    # La méthode suivante résout le graphe au sens où elle rajoute des sommets simples sur les arrêtes multiples ou les boucles aux sommets :
    def resolve(self):
        resolved = Graph()
        n = len(self.vertices)
        
        # Ajout des sommets du graphe original
        for v in range(n):
            resolved.addNewVertice()
            
        for x in range(n):
            # Traitement du cas des arrêtes multiples
            for y in range(x+1,n):
                Exy = [e for e in self.edges[x] if e[0] == y]
                Nxy = len(Exy)
                if Nxy>0:
                    for itt in range(1,Nxy):
                        e = resolved.addNewVertice()
                        d = Exy[itt][1]
                        resolved.addEdge(x,e,d/2)
                        resolved.addEdge(e,y,d/2)
                    d = Exy[0][1]
                    resolved.addEdge(x,y,d)
            # Traitement du cas des sommets à boucles
            Exx = [e for e in self.edges[x] if e[0] == x]
            Nxx = len(Exx)
            if Nxx>0:
                for itt in range(Nxx):
                    e = resolved.addNewVertice()
                    f = resolved.addNewVertice()
                    d = Exx[itt][1]
                    resolved.addEdge(x,e,d/3)
                    resolved.addEdge(e,f,d/3)
                    resolved.addEdge(f,x,d/3)
                    
        return resolved
        
    # Matrices de continuité aux sommets du graphe :
    def A(self,a):
        resolved = self.resolve()
        A = {}
        for v in resolved.vertices:
            d = len(resolved.edges[v])
            tmp=np.zeros((d,d))
            tmp[-1][0]=np.sin(a)
            tmp[-1][-1]=1
            A[v]=np.eye(d)-np.diagflat((d-1)*[1],1)-tmp
        return A
    def B(self,a):
        resolved = self.resolve()
        B = {}
        for v in resolved.vertices:
            d = len(resolved.edges[v])
            tmp=np.zeros((d,d))
            tmp[-1][0]=np.sin(a)
            tmp[-1][-1]=1
            tmp=np.cos(a)*np.array(d*[1])
            B[v]=np.zeros((d,d))
            B[v][-1]=tmp
        return B
    def sigma(self,x,a):
        A = self.A(a)
        B = self.B(a)
        sigma = {}
        i=complex(0,1)
        for v in A.keys():
            if det(A[v]+i*x*B[v]) != 0:
                sigma[v] = -inv(A[v]+i*x*B[v]).dot(A[v]-i*x*B[v])
            else:
                sigma[v] = 0*A[v]
        return sigma
        
    # Retourne l'ensemble des arrêtes orientées du graphe
    def bonds(self):
        resolved = self.resolve()
        bonds = []
        for v in resolved.vertices:
            for w in resolved.edges[v]:
                bonds.append((v,w[0],w[1]))
        return bonds
        
    # Matrice de diffusion du graphe :
    def S(self,x=pi,a=0):
        resolved = self.resolve()
        bonds = self.bonds()
        B = len(bonds)
        S = np.full((B,B),complex(0,0))
        sigma = self.sigma(x,a)
        for b1 in range(B):
            for b2 in range(B):
                if bonds[b1][1] == bonds[b2][0]:
                    v = bonds[b1][1]
                    c = resolved.edges[v].index((bonds[b1][0],bonds[b1][2]))
                    l = resolved.edges[v].index((bonds[b2][1],bonds[b2][2]))
                    S[b1][b2] = sigma[v][c][l]
        return S
        
    def M(self,x=pi,a=0):
        i=complex(0,1)
        bonds = self.bonds()
        B = len(bonds)
        L = [b[2] for b in bonds]
        ExpL=np.diag(np.exp(-i*x*np.array(L)))
        I = np.eye(B)
        return I-self.S(x,a).dot(ExpL)
        
    # Membre gauche de l'équation séculaire du graphe :
    def secEq(self,x=pi,a=0):
        return det(self.M(x,a))
        
    # Première valeur propre non nulle du graphe :
    def k1(self,a=0):
        start = 0
        end = 10*pi*self.E/self.L
        def absSecEq(x):
            return abs(self.secEq(x,a))
        ret = [float("inf")]
        f = 20
        for k in range(1,int(f)):
            retk=so.fminbound(absSecEq,start+k*(end-start)/f,start+(k+1)*(end-start)/f,full_output=True)
            if retk[1]<1e-4:
                ret = retk
                break
        return ret[0]
        
    # Optimisation de la géométrie du graphe : minimisation de la première valeur propre
    def minimizeOverGeo(self):
        print("...")
        def f(b):
            G = Graph()
            for v in self.vertices:
                G.addNewVertice()
            i = 0
            for v in self.vertices:
                for w in self.edges[v]:
                    if v <= w[0]:
                        if b[i] < 0:
                            return float("inf")
                        G.addEdge(v,w[0],b[i])
                        i += 1
            self.vertices = G.vertices
            self.edges = G.edges
            self.normalize()
            return self.k1()
        b0 = []
        for v in self.vertices:
            for w in self.edges[v]:
                if v <= w[0]:
                    b0.append(w[1])
        b0 = np.array(b0)
        ret = so.minimize(f, b0, method='nelder-mead', options={'xatol': 1e-8, 'disp': True}) 
        return self.k1()
        
    # Optimisation de la géométrie du graphe : maximisation de la première valeur propre
    def maximizeOverGeo(self):
        print("...")
        def f(b):
            G = Graph()
            for v in self.vertices:
                G.addNewVertice()
            i = 0
            for v in self.vertices:
                for w in self.edges[v]:
                    if v <= w[0]:
                        if b[i] < 0:
                            return float("inf")
                        G.addEdge(v,w[0],b[i])
                        i += 1
            self.vertices = G.vertices
            self.edges = G.edges
            self.normalize()
            return -self.k1()
        b0 = []
        for v in self.vertices:
            for w in self.edges[v]:
                if v <= w[0]:
                    b0.append(w[1])
        b0 = np.array(b0)
        ret = so.minimize(f, b0, method='nelder-mead', options={'xatol': 1e-12, 'disp': True}) 
        return self.k1()
        
    # Module du membre gauche de l'équation séculaire en fonction du paramètre x :
    def plotNormOfLeft(self,a=0):
        def absSecEq(x):
            return abs(self.secEq(pi*x,a))
        start = 0
        end = 10*self.E*pi/self.L
        X = np.linspace(start,end,num=floor((end-start)*1000))
        Y = [absSecEq(x) for x in X]
        plt.clf()
        plt.plot(X,Y)
        plt.show()