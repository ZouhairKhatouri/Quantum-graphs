import os
path = 'C:\\Users\\zouha\\Desktop\\Rapport final PSC\\PSC_CODE\\'
os.chdir(path)
from Graph import *

## Optimisation directe de la geometrie d'un graphe :

# Graphe mandarine à 5 arêtes et à deux boucles :
G = Graph()
G.addNewVertice()
G.addNewVertice()
G.addEdge(0,1) 
G.addEdge(0,1) 
G.addEdge(0,1) 
G.addEdge(0,1) 
G.addEdge(0,1)
G.addEdge(0,0) 
G.normalize()

print("Valeur propre initiale trouvée :"+str(G.k1()/pi)+"pi")

# Résultat de la simulation :
print("Valeur propre optimale trouvée :"+str(G.maximizeOverGeo()/pi)+"pi")