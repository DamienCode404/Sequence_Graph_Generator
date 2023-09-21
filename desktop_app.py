import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
from Bio import SeqIO
from Bio import pairwise2
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import customtkinter
import seaborn as sns


# Global variables
fasta_file_path = ""
sequences = []
sequence_names = []
alignments = []
sequence_checkboxes = []


def open_fasta_file():
    global fasta_file_path
    fasta_file_path = filedialog.askopenfilename(
        filetypes=[("FASTA files", "*.fasta")])
    if fasta_file_path:
        fasta_file_label.configure(text=f"FASTA File: {fasta_file_path}")
        load_sequences()


def reset_checkboxes():
    for checkbox, var in sequence_checkboxes:
        checkbox.destroy()
    sequence_checkboxes.clear()


def load_sequences():
    global sequences, sequence_names  # Utilisez les variables globales

    sequences = []  # Réinitialisez la liste des séquences
    sequence_names = []  # Réinitialisez la liste des noms de séquence

    try:
        fasta_sequences = SeqIO.parse(fasta_file_path, "fasta")
        for sequence in fasta_sequences:
            sequences.append(sequence.seq)
            sequence_name = sequence.description.split(
                "|")[0]  # Obtenez le nom avant le caractère "|"
            sequence_names.append(sequence_name)

        for i in range(len(sequences)):
            for j in range(i + 1, len(sequences)):
                align = pairwise2.align.globalxx(
                    sequences[i], sequences[j], one_alignment_only=True)[0]
                alignments.append(align)

        update_sequence_checkboxes()
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {str(e)}")


def update_sequence_checkboxes():
    for checkbox in sequence_checkboxes:
        checkbox.destroy()
    sequence_checkboxes.clear()

    for i, sequence_name in enumerate(sequence_names):
        checkbox_var = tk.IntVar()
        checkbox = customtkinter.CTkCheckBox(
            master=frame, text=sequence_name, variable=checkbox_var)
        checkbox.pack(pady=12, padx=10, anchor=tk.W)
        sequence_checkboxes.append((checkbox, checkbox_var))

    # Enable or disable the "Generate Graph" button based on checkbox selection
    generate_button.configure(
        state=tk.NORMAL if sequence_checkboxes else tk.DISABLED)

    # Enable or disable the "Generate Heatmap" button based on whether sequences are loaded
    heatmap_button.configure(
        state=tk.NORMAL if sequence_checkboxes else tk.DISABLED)


def generate_graph():
    if not alignments:
        messagebox.showinfo(
            "Info", "Please load or reload sequences from a FASTA file.")
        return

    selected_sequence_indices = [i for i, (_, var) in enumerate(
        sequence_checkboxes) if var.get() == 1]

    if not selected_sequence_indices:
        messagebox.showinfo(
            "Info", "Please select at least one sequence to generate the graph.")
        return

    # Create a NetworkX graph and add nodes
    G = nx.Graph()
    for i in selected_sequence_indices:
        sequence_name = sequence_names[i]
        G.add_node(sequence_name)

    # Add edges with alignment scores as weights
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            score = alignments.pop(0).score
            if i in selected_sequence_indices and j in selected_sequence_indices:
                sequence_name_i = sequence_names[i]
                sequence_name_j = sequence_names[j]
                G.add_edge(sequence_name_i, sequence_name_j, weight=score)

    # Calculate the maximum alignment score for scaling the colorbar
    max_alignment_score = max([d['weight'] for u, v, d in G.edges(data=True)])

    # Plot the graph
    pos = nx.spring_layout(G)
    edge_labels = {(u, v): f"{d['weight']:.2f}" for u,
                   v, d in G.edges(data=True)}
    edge_colors = [d['weight'] for u, v, d in G.edges(data=True)]

    plt.figure(figsize=(12, 8))
    nx.draw_networkx_nodes(
        G, pos, node_size=3000, node_color='white')
    nx.draw_networkx_labels(G, pos, font_size=10, font_color='black')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=12)
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, edge_cmap=cm.RdYlBu)
    # Specify vmin and vmax for the colorbar
    cbar = plt.colorbar(plt.cm.ScalarMappable(
        cmap=cm.RdYlBu, norm=plt.Normalize(vmin=0, vmax=max_alignment_score)))

    plt.axis('off')
    plt.title('Alignment Score Graph')
    plt.show()


def generate_heatmap():
    if not alignments:
        messagebox.showinfo(
            "Info", "Please load or reload sequences from a FASTA file.")
        return

    # Obtenez les indices des séquences sélectionnées à partir des cases à cocher
    selected_sequence_indices = [i for i, (_, var) in enumerate(
        sequence_checkboxes) if var.get() == 1]

    if not selected_sequence_indices:
        messagebox.showinfo(
            "Info", "Please select at least one sequence to generate the heatmap.")
        return

    # Créez une liste pour stocker les séquences sélectionnées
    selected_sequences = [sequences[i] for i in selected_sequence_indices]
    # Créez une liste pour stocker les noms des séquences sélectionnées
    selected_sequence_names = [sequence_names[i]
                               for i in selected_sequence_indices]

    # Créez une matrice carrée pour stocker les scores d'alignement entre les séquences sélectionnées
    num_selected_sequences = len(selected_sequences)
    scores_matrix = [
        [0] * num_selected_sequences for _ in range(num_selected_sequences)]

    for i in range(num_selected_sequences):
        for j in range(i + 1, num_selected_sequences):
            # Calculez les scores d'alignement entre les séquences sélectionnées
            align = pairwise2.align.globalxx(
                selected_sequences[i], selected_sequences[j], one_alignment_only=True)[0]
            score = align.score
            scores_matrix[i][j] = score
            # Assurez-vous de remplir les deux côtés de la matrice
            scores_matrix[j][i] = score

    # Utilisez Seaborn pour créer un heatmap
    sns.set()
    plt.figure(figsize=(10, 8))
    sns.heatmap(scores_matrix, annot=True, fmt=".2f",
                xticklabels=selected_sequence_names, yticklabels=selected_sequence_names)

    plt.title("Alignment Score Heatmap")
    plt.show()


customtkinter.set_appearance_mode("System")
customtkinter.set_default_color_theme("blue")

# Create the main window
root = customtkinter.CTk()
root.geometry("800x600")

frame = customtkinter.CTkFrame(master=root)
frame.pack(pady=20, padx=60, fill="both", expand=True)

label = customtkinter.CTkLabel(
    master=frame, text="Sequence Graph Generator", font=("Roboto", 24))
label.pack(pady=12, padx=10)

# Create a frame for buttons
buttons_frame = customtkinter.CTkFrame(master=root)
buttons_frame.pack(pady=12, padx=10)

# Create a label to display the selected FASTA file path
fasta_file_label = customtkinter.CTkLabel(
    master=frame, text="FASTA File: None", font=("Roboto", 18))
fasta_file_label.pack(pady=12, padx=10)

# Create buttons for opening the FASTA file and generating the graph
open_button = customtkinter.CTkButton(
    master=buttons_frame, text="Open FASTA File", command=open_fasta_file)
generate_button = customtkinter.CTkButton(
    master=buttons_frame, text="Generate Graph", command=generate_graph, state=tk.DISABLED)
heatmap_button = customtkinter.CTkButton(
    master=buttons_frame, text="Generate Heatmap", command=generate_heatmap, state=tk.DISABLED)


# Utilisez l'option side pour les mettre côte à côte
open_button.pack(pady=12, padx=10, fill=tk.BOTH, side=tk.LEFT)
generate_button.pack(pady=12, padx=10, fill=tk.BOTH, side=tk.LEFT)
heatmap_button.pack(pady=12, padx=10, fill=tk.BOTH, side=tk.LEFT)

# Start the main loop
root.mainloop()
