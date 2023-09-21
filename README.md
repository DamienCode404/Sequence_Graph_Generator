# Sequence Graph Generator from FASTA Files

This project is a sequence graph generator for biological sequences from FASTA files. It allows users to visualize similarities between different sequences and generate graphs and heatmaps for comparative analysis.

## Features

- **FASTA File Loading**: The program allows users to load FASTA files containing biological sequences.

- **Similarity Analysis**: It utilizes the Biopython library to compute similarities between sequences using global alignment.

- **Graph Visualization**: Similarities between sequences are represented as graphs, where sequences are nodes and alignment scores are edge weights.

- **Alignment Score Heatmap**: Alignment scores between selected sequences can be visualized as a heatmap.

## How to Use

1. Select a FASTA file using the "Open FASTA File" button.

2. Sequences will be loaded and displayed. You can select the sequences you want to include in the graph or heatmap.

3. Use the "Generate Graph" button to create a graph of similarities between the selected sequences.

4. Use the "Generate Heatmap" button to create a heatmap of alignment scores between the selected sequences.

## Dependencies

This project uses the following libraries:

- Biopython for FASTA file processing and alignment calculations.
- NetworkX for graph creation and visualization.
- Matplotlib for generating graphs.
- CustomTkinter for a customized user interface.
- Seaborn for heatmap creation.
