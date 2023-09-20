from Bio import pairwise2
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as cm
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
    G.add_node(sequence_name)

# Ajouter des arêtes avec les scores d'alignement comme poids
for i in range(len(sequences)):
    for j in range(i + 1, len(sequences)):
        score = alignments.pop(0).score
        G.add_edge(sequence_names[i], sequence_names[j], weight=score)

# Dessiner le graphe
# Distribuer les nœuds en utilisant un algorithme de disposition
pos = nx.spring_layout(G)
# Obtenez les poids d'arête comme labels
edge_labels = nx.get_edge_attributes(G, 'weight')

# Créer une liste de couleurs en fonction des scores d'arête
edge_colors = [d['weight'] for u, v, d in G.edges(data=True)]

plt.figure(figsize=(12, 8))  # Réglez la taille de la figure selon vos besoins

# Dessinez les nœuds avec une couleur fixe
nx.draw_networkx_nodes(
    G, pos, node_size=3000, node_color='lightblue')

nx.draw_networkx_labels(G, pos, font_size=10, font_color='black')
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=8)
# Dessinez les arêtes avec des couleurs basées sur les scores d'arête
nx.draw_networkx_edges(G, pos, edge_color=edge_colors, edge_cmap=cm.RdYlBu)

# Créez une barre de couleur pour afficher les scores
cbar = plt.colorbar(plt.cm.ScalarMappable(
    cmap=cm.RdYlBu), label='Score d\'alignement')

plt.axis('off')
# Titre et affichage
plt.title('Graphe des scores d\'alignement entre les séquences avec couleur')
plt.show()
