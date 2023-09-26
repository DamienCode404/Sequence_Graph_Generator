import tkinter as tk
from tkinter import filedialog, messagebox
from Bio import SeqIO, pairwise2
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import customtkinter
import seaborn as sns

# Constants
FASTA_FILE_TYPES = [("FASTA files", "*.fasta")]
NO_FASTA_FILE_SELECTED = "FASTA File: None"


class SequenceGraphGenerator:
    def __init__(self):
        self.fasta_file_path = ""
        self.sequences = []
        self.sequence_names = []
        self.alignments = []
        self.sequence_checkboxes = []
        # Initialisez le seuil ici (à ajuster selon vos besoins)
        self.seuil = 0.5

        customtkinter.set_appearance_mode("System")
        customtkinter.set_default_color_theme("blue")

        self.root = customtkinter.CTk()
        self.root.geometry("800x600")

        frame = customtkinter.CTkFrame(master=self.root)
        frame.pack(pady=20, padx=60, fill="both", expand=True)

        self.frame = frame

        label = customtkinter.CTkLabel(
            master=frame, text="Sequence Graph Generator", font=("Roboto", 24))
        label.pack(pady=12, padx=10)

        buttons_frame = customtkinter.CTkFrame(master=self.root)
        buttons_frame.pack(pady=12, padx=10)

        self.fasta_file_label = customtkinter.CTkLabel(
            master=frame, text=NO_FASTA_FILE_SELECTED, font=("Roboto", 18))
        self.fasta_file_label.pack(pady=12, padx=10)

        open_button = customtkinter.CTkButton(
            master=buttons_frame, text="Open FASTA File", command=self.open_fasta_file)
        self.heatmap_button = customtkinter.CTkButton(
            master=buttons_frame, text="Generate Heatmap", command=self.generate_heatmap, state=tk.DISABLED)
        self.generate_button = customtkinter.CTkButton(
            master=buttons_frame, text="Generate Graph", command=self.generate_graph, state=tk.DISABLED)

        open_button.pack(pady=12, padx=10, fill=tk.BOTH, side=tk.LEFT)
        self.heatmap_button.pack(pady=12, padx=10, fill=tk.BOTH, side=tk.LEFT)
        self.generate_button.pack(pady=12, padx=10, fill=tk.BOTH, side=tk.LEFT)

        self.root.mainloop()

    def open_fasta_file(self):
        self.fasta_file_path = filedialog.askopenfilename(
            filetypes=FASTA_FILE_TYPES)
        if self.fasta_file_path:
            self.fasta_file_label.configure(
                text=f"FASTA File: {self.fasta_file_path}")
            self.load_sequences()

    def reset_checkboxes(self):
        for checkbox, var in self.sequence_checkboxes:
            checkbox.destroy()
        self.sequence_checkboxes.clear()

    def load_sequences(self):
        self.sequences = []
        self.sequence_names = []
        self.alignments = []
        try:
            fasta_sequences = SeqIO.parse(self.fasta_file_path, "fasta")
            for sequence in fasta_sequences:
                self.sequences.append(sequence.seq)
                sequence_name = sequence.description.split("|")[0]
                self.sequence_names.append(sequence_name)

            self.calculate_alignments()
            self.update_sequence_checkboxes()
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {str(e)}")

    def calculate_alignments(self):
        for i in range(len(self.sequences)):
            for j in range(i + 1, len(self.sequences)):
                align = pairwise2.align.globalxx(
                    self.sequences[i], self.sequences[j], one_alignment_only=True)[0]
                self.alignments.append(align)
                # print(align)  # Debugging: Check alignment scores

    def update_sequence_checkboxes(self):
        self.reset_checkboxes()
        for i, sequence_name in enumerate(self.sequence_names):
            checkbox_var = tk.IntVar()
            checkbox = customtkinter.CTkCheckBox(
                master=self.frame, text=sequence_name, variable=checkbox_var)
            checkbox.pack(pady=12, padx=10, anchor=tk.W)
            self.sequence_checkboxes.append((checkbox, checkbox_var))

        self.generate_button.configure(
            state=tk.NORMAL if self.sequence_checkboxes else tk.DISABLED)
        self.heatmap_button.configure(
            state=tk.NORMAL if self.sequence_checkboxes else tk.DISABLED)

    def generate_graph(self):
        if not self.alignments:
            messagebox.showinfo(
                "Info", "Please load or reload sequences from a FASTA file.")
            return

        selected_sequence_indices = [i for i, (_, var) in enumerate(
            self.sequence_checkboxes) if var.get() == 1]

        if not selected_sequence_indices:
            messagebox.showinfo(
                "Info", "Please select at least one sequence to generate the graph.")
            return

        # Create a NetworkX graph and add nodes
        G = nx.Graph()
        for i in selected_sequence_indices:
            sequence_name = self.sequence_names[i]
            G.add_node(sequence_name)

        # Add edges with alignment scores as weights
        for i in range(len(self.sequences)):
            for j in range(i + 1, len(self.sequences)):
                score = self.alignments.pop(0).score
                if i in selected_sequence_indices and j in selected_sequence_indices:
                    sequence_name_i = self.sequence_names[i]
                    sequence_name_j = self.sequence_names[j]
                    G.add_edge(sequence_name_i, sequence_name_j, weight=score)

        # Calculate the maximum alignment score for scaling the colorbar
        max_alignment_score = max([d['weight']
                                  for u, v, d in G.edges(data=True)])

        # Plot the graph
        pos = nx.spring_layout(G)
        edge_labels = {(u, v): f"{d['weight']:.2f}" for u,
                       v, d in G.edges(data=True)}
        edge_colors = [d['weight'] for u, v, d in G.edges(data=True)]

        plt.figure(figsize=(12, 8))
        nx.draw_networkx_nodes(
            G, pos, node_size=3000, node_color='white')
        nx.draw_networkx_labels(G, pos, font_size=10, font_color='black')
        nx.draw_networkx_edge_labels(
            G, pos, edge_labels=edge_labels, font_size=12)
        nx.draw_networkx_edges(
            G, pos, edge_color=edge_colors, edge_cmap=cm.RdYlBu)
        # Specify vmin and vmax for the colorbar
        cbar = plt.colorbar(plt.cm.ScalarMappable(
            cmap=cm.RdYlBu, norm=plt.Normalize(vmin=0, vmax=max_alignment_score)))

        plt.axis('off')
        plt.title('Alignment Score Graph')
        plt.show()

    def generate_heatmap(self):
        if not self.alignments:
            messagebox.showinfo(
                "Info", "Please load or reload sequences from a FASTA file.")
            return

        selected_sequence_indices = [i for i, (_, var) in enumerate(
            self.sequence_checkboxes) if var.get() == 1]

        if not selected_sequence_indices:
            messagebox.showinfo(
                "Info", "Please select at least one sequence to generate the heatmap.")
            return

        selected_sequences = [self.sequences[i]
                              for i in selected_sequence_indices]
        selected_sequence_names = [self.sequence_names[i]
                                   for i in selected_sequence_indices]

        num_selected_sequences = len(selected_sequences)
        scores_matrix = [
            [0] * num_selected_sequences for _ in range(num_selected_sequences)]

        for i in range(num_selected_sequences):
            for j in range(i + 1, num_selected_sequences):
                align = pairwise2.align.globalxx(
                    selected_sequences[i], selected_sequences[j], one_alignment_only=True)[0]
                score = align.score
                scores_matrix[i][j] = score
                scores_matrix[j][i] = score

        sns.set()
        plt.figure(figsize=(10, 8))
        sns.heatmap(scores_matrix, annot=True, fmt=".2f",
                    xticklabels=selected_sequence_names, yticklabels=selected_sequence_names)

        plt.title("Alignment Score Heatmap")
        plt.show()


if __name__ == "__main__":
    app = SequenceGraphGenerator()
