from Bio import pairwise2
import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random

from Bio import SeqIO

# Remplacez 'votre_fichier.fasta' par le nom de votre fichier FASTA
fasta_file = "sequence1.fasta"

sequences = []
sequence_names = []  # Pour stocker les noms des séquences

for record in SeqIO.parse(fasta_file, "fasta"):
    sequences.append(record.seq)
    # Ajoutez le nom à partir de la description (en supprimant tout après le premier espace)
    sequence_name = record.description.split(" ")[0]
    sequence_names.append(sequence_name)

alignments = []

for i in range(len(sequences)):
    for j in range(i + 1, len(sequences)):
        align = pairwise2.align.globalxx(
            sequences[i], sequences[j], one_alignment_only=True)[0]
        alignments.append(align)

# Créer un graphique NetworkX
G = nx.Graph()

# Ajouter des nœuds pour chaque séquence avec leur nom
for i, sequence_name in enumerate(sequence_names):
    G.add_node(i)  # Utilisez un indice entier comme nœud
    G.nodes[i]['name'] = sequence_name  # Stockez le nom du nœud comme attribut

# Ajouter des arêtes avec les scores d'alignement comme poids
for i in range(len(sequences)):
    for j in range(i + 1, len(sequences)):
        score = alignments.pop(0).score
        G.add_edge(i, j, weight=score)

# Créer une figure 3D
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')  # Utiliser la projection 3D

# Distribuer les nœuds en utilisant un algorithme de disposition
pos = nx.spring_layout(G, dim=3)  # Spécifier dim=3 pour une disposition en 3D

# Obtenez les positions des nœuds en 3D
pos_3d = {i: pos[i] for i in G.nodes()}

# Obtenez les poids d'arête comme labels
edge_labels = nx.get_edge_attributes(G, 'weight')

# Dessinez les nœuds
nx.draw_networkx_nodes(G, pos_3d, node_size=1000, node_color='skyblue')

# Ajoutez les arêtes une par une
for (u, v, data) in G.edges(data=True):
    x = [pos_3d[u][0], pos_3d[v][0]]
    y = [pos_3d[u][1], pos_3d[v][1]]
    z = [pos_3d[u][2], pos_3d[v][2]]
    ax.plot(x, y, z, c='gray', alpha=0.5, linewidth=data['weight'] * 0.01)

# Affichez les labels des nœuds
for node, (x, y, z) in pos_3d.items():
    ax.text(x, y, z, G.nodes[node]['name'], fontsize=8)

# Ajoutez des unités aux axes x, y et z
ax.set_xlabel('X (unité)')
ax.set_ylabel('Y (unité)')
ax.set_zlabel('Z (unité)')

# Titre et affichage
plt.title('Graphe 3D des scores d\'alignement entre les séquences')
plt.show()
