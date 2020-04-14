from math import pi , sin
from copy import deepcopy
import numpy as np
from numpy.linalg import norm
import os
path = 'C:\\Users\\zouha\\Desktop\\Rapport final PSC\\PSC_CODE\\'
os.chdir(path)
from Graph import *

## Classe modélisant une fonction sur un graphe :

class GraphFunction:
    
    def __init__(self,G,epsilon=10**(-4),n=1,a=0):
        self.G = deepcopy(G).resolve()
        self.epsilon = epsilon
        self.n = n
        self.a = a

        self.bonds = self.G.bonds()
        self.B = len(self.bonds)
        self.f = np.full((self.B,1),complex(0,0))

        self.x = self.n*pi
        self.S = self.G.S(self.x,self.a)
        self.I = np.eye(self.B)
        self.index = {}
        for b in range(self.B):
            self.index[(self.bonds[b][0],self.bonds[b][1])] = b
    
    # Transfert des longueurs entre deux arrêtes du graphe :
    def transfertLengths(self,b1,b2):
        v1b1 = self.bonds[b1][0]
        v2b1 = self.bonds[b1][1]
        v1b2 = self.bonds[b2][0]
        v2b2 = self.bonds[b2][1]
        
        found = False
        x = None
        for itt in range(len(self.G.edges[v1b1])):
            w = self.G.edges[v1b1].pop(0)
            if not found and w[0] == v2b1:
                x = w
                found = True
                continue
            self.G.edges[v1b1].append(w)
        self.G.edges[v1b1].append((x[0],x[1]-self.epsilon))
        
        found = False
        x = None
        for itt in range(len(self.G.edges[v2b1])):
            w = self.G.edges[v2b1].pop(0)
            if not found and w[0] == v1b1:
                x = w
                found = True
                continue
            self.G.edges[v2b1].append(w)
        self.G.edges[v2b1].append((x[0],x[1]-self.epsilon))
        
        found = False
        x = None
        for itt in range(len(self.G.edges[v1b2])):
            w = self.G.edges[v1b2].pop(0)
            if not found and w[0] == v2b2:
                x = w
                found = True
                continue
            self.G.edges[v1b2].append(w)
        self.G.edges[v1b2].append((x[0],x[1]+self.epsilon))
        
        found = False
        x = None
        for itt in range(len(self.G.edges[v2b2])):
            w = self.G.edges[v2b2].pop(0)
            if not found and w[0] == v1b2:
                x = w
                found = True
                continue
            self.G.edges[v2b2].append(w)
        self.G.edges[v2b2].append((x[0],x[1]+self.epsilon))
        
        self.bonds = self.G.bonds()
    
        return norm((self.G.M()).dot(self.f))
    
    # Incrémentation des coefficients complexes d'une arrête :
    def update_f(self,b,u,v):
        b_bar = self.index[(self.bonds[b][1],self.bonds[b][0])]
        self.f[b] += u*self.epsilon
        self.f[b_bar] += v*self.epsilon
        
        return norm((self.G.M()).dot(self.f))
        
    # Incrémentation du paramètre x :
    def update_x(self,e):
        self.x += e*self.epsilon
        
        return norm((self.G.M()).dot(self.f))
    
    # Calcul du rapport de Rayleigh :
    def R(self):
        N = 0
        D = 0
        for b in self.bonds:
            L_b = b[2]
            a_b = self.f[self.index[(b[0],b[1])]][0]
            a_b_bar = self.f[self.index[(b[1],b[0])]][0]
            N += ((abs(a_b)**2)*L_b*self.x)-((a_b*np.conjugate(a_b_bar))*sin(L_b*self.x))
            D += ((abs(a_b)**2)*L_b*self.x)+((a_b*np.conjugate(a_b_bar))*sin(L_b*self.x))
        if D != 0:
            return abs(N/D)*(self.x)**2
        else:
            return float("-inf")
    
    """
    # Jacobien de la fonction qui associe à une géométrie une valeur propre :
    def jac(self):
        N = 0
        for b in self.bonds:
            L_b = b[2]
            a_b = self.f[self.index[(b[0],b[1])]][0]
            a_b_bar = self.f[self.index[(b[1],b[0])]][0]
            N += ((abs(a_b)**2)*L_b*self.x)-((a_b*np.conjugate(a_b_bar))*sin(L_b*self.x))
        return -abs(N)*(self.x)**2
    """