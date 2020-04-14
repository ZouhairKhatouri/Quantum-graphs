import os
path = 'C:\\Users\\zouha\\Desktop\\Rapport final PSC\\PSC_CODE\\'
os.chdir(path)
from random import uniform
from math import log , exp, floor
from Graph import *
from GraphFunction import *

## Algorithme de recuit simule pour optimisation de la geometrie d'un graphe : 

def h(u):
    return u/(1+u)
    
# q doit etre le plus grand parmis ces trois parametres afin de permettre a l'energie de converger vers la valeur propre de la configuration geometrique "courante"
# Doivent avoir pour somme 1! 
p = 0.3
q = 0.5
r = 0.2
# Parametre de la temperature du recuit simule
K = 10
# Seuil de nullité 
epsilon = 10**(-4)

# Graphe mandarine a deux boucles
G = Graph()
G.addNewVertice()
G.addNewVertice()
G.addEdge(0,1)                    
G.addEdge(0,1) 
G.addEdge(0,1) 
G.addEdge(0,1) 
G.addEdge(0,1) 
G.addEdge(0,0)
G.addEdge(1,1)
G.normalize()
F = GraphFunction(G)
B = F.B

""" 
Nomenclature :
    - Tn : temperature a l'instant n
    - E : energie de l'etat courant
    - E0 : plus grande valeur propre parmis celles deja trouvees
    - w : alea utilise dans les simulations en recuit simule
"""

E0 = 0

for n in range(2,10002):
    print("epoch = "+str(n-1))
    
    Tn = K/log(n)
    
    U = uniform(0,1)
    
    # Premier cas : on incremente le parametre F.x
    if U < p:
        # Transformation de l'etat :
        e = 2*uniform(0,1)-1
        Q = F.update_x(e)
        # Alea :
        w = uniform(0,1)
        # Energie :
        E = F.R()
        # Test de conformite du nouveau etat :
        if Q >= epsilon or E-E0 < epsilon or w > h(exp(-E/Tn)):
            # Retour vers l'etat precedent (de plus faible energie) :
            F.update_x(-e)
            
    # Deuxieme cas : on translate le vecteur F.f suivant un vecteur proportionnel à l'un des vecteurs de la base canonique (choisis uniformement parmis tous les vecteurs de la base canonique)
    elif p <= U and U < p+q:
        # Transformation de l'etat :
        V = uniform(0,1)
        b = floor(B*V)
        u_Re = 2*uniform(0,1)-1
        u_Im = 2*uniform(0,1)-1
        v_Re = 2*uniform(0,1)-1
        v_Im = 2*uniform(0,1)-1
        u = complex(u_Re,u_Im)
        v = complex(v_Re,v_Im)
        Q = F.update_f(b,u,v)
        # Alea :
        w = uniform(0,1)
        # Energie :
        E = F.R()
        # Test de conformite du nouveau etat :
        if Q >= epsilon or E-E0 < epsilon or w > h(exp(-E/Tn)):
            # Retour vers l'etat precedent (de plus faible energie) :
            F.update_f(b,-u, -v)
            
    # Troisieme cas : on transfert un fraction de longueur d'une arete vers l'autre (choisies uniformement et independament parmis toutes les aretes du graphe)
    elif p+q <= U :
        # Transformation de l'etat :
        V1 = uniform(0,1)
        V2 = uniform(0,1)
        b1 = floor(B*V1)
        b2 = floor(B*V2)
        Q = F.transfertLengths(b1,b2)
        # Alea :
        w = uniform(0,1)
        # Energie :
        E = F.R()
        # Test de conformite du nouveau etat :
        if Q >= epsilon or E-E0 < epsilon or w > h(exp(-F.R()/Tn)):
            # Retour vers l'etat precedent (de plus faible energie) :
            F.transfertLengths(b2,b1)
            
    else:
        print("Une erreur a survenue!")
        break
        
# Resultat de la simulation :
print("Résultat de la simulation :")
print("Valeur propre trouvée :"+str(F.x))
print("Energie finale (rapport de Rayleigh) :"+str(F.R()))