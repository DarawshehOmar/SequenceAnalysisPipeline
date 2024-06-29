#!/usr/bin/env python

import tkinter as tk
from tkinter import messagebox, scrolledtext, filedialog
from Bio import Entrez, AlignIO, Phylo
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import os
import threading

def fetch_sequences():
    """Fetch sequences from NCBI databases using Entrez."""
    Entrez.email = email_entry.get()
    query = query_entry.get()
    if not Entrez.email or not query:
        messagebox.showerror("Error", "Email and query must be provided!")
        return

    try:
        handle = Entrez.esearch(db='nucleotide', term=query)
        record = Entrez.read(handle)
        ids = record['IdList']
        seq_handle = Entrez.efetch(db='nucleotide', id=','.join(ids), rettype='fasta')
        sequences = seq_handle.read()
        seq_handle.close()
        output_text.delete('1.0', tk.END)
        output_text.insert(tk.END, sequences)
        messagebox.showinfo("Success", "Sequences fetched successfully.")
    except Exception as e:
        messagebox.showerror("Error", str(e))

def save_sequences():
    """Save fetched sequences to a file."""
    sequences = output_text.get("1.0", tk.END)
    filename = filedialog.asksaveasfilename(
        defaultextension=".fasta",
        filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")],
    )
    if not filename:
        messagebox.showerror("Error", "Name of the file must be provided!")
        return
    with open(filename, 'w') as file:
        file.write(sequences)
    messagebox.showinfo("Success", "Sequences saved to file.")

def perform_alignment():
    """Perform sequence alignment using Clustal Omega and construct a phylogenetic tree."""
    progress_label.config(text="Alignment in progress. This might take a few minutes.")
    root.update_idletasks()
    
    sequences = output_text.get("1.0", tk.END)
    if not sequences.strip():
        messagebox.showerror("Error", "No sequences to align.")
        progress_label.config(text="")
        return

    temp_filename = 'temp_sequences.fasta'
    with open(temp_filename, 'w') as file:
        file.write(sequences)

    clustalo_exe = "clustalo"  # Ensure clustalo is in your PATH, or provide full path if necessary
    clustalo_cline = ClustalOmegaCommandline(clustalo_exe, infile=temp_filename, outfile="temp_sequences.aln", verbose=True, auto=True, force=True)

    def run_alignment():
        try:
            stdout, stderr = clustalo_cline()
            if stderr:
                raise Exception(f"Clustal Omega execution failed: {stderr}")

            if os.path.exists('temp_sequences.aln') and os.path.getsize('temp_sequences.aln') > 0:
                alignment = AlignIO.read('temp_sequences.aln', 'fasta')
                create_phylogenetic_tree(alignment)
            else:
                messagebox.showerror("Error", "Alignment failed or output file is empty.")
        except Exception as e:
            messagebox.showerror("Error", str(e))
        finally:
            progress_label.config(text="")

    threading.Thread(target=run_alignment).start()

def create_phylogenetic_tree(alignment):
    """Construct a phylogenetic tree from the alignment."""
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(alignment)
    Phylo.draw(tree)

root = tk.Tk()
root.title("Sequence Retrieval and Analysis Pipeline")

tk.Label(root, text="Email:").grid(row=0, column=0, sticky="e")
email_entry = tk.Entry(root, width=30)
email_entry.grid(row=0, column=1, padx=5, pady=5)

tk.Label(root, text="Query:").grid(row=1, column=0, sticky="e")
query_entry = tk.Entry(root, width=30)
query_entry.grid(row=1, column=1, padx=5, pady=5)

tk.Button(root, text="Fetch Sequences", command=fetch_sequences).grid(row=2, column=0, columnspan=2, pady=5)
tk.Button(root, text="Save Sequences", command=save_sequences).grid(row=3, column=0, columnspan=2, pady=5)
tk.Button(root, text="Perform Alignment", command=perform_alignment).grid(row=4, column=0, columnspan=2, pady=5)

output_text = scrolledtext.ScrolledText(root, width=60, height=10)
output_text.grid(row=5, column=0, columnspan=2, pady=5)

progress_label = tk.Label(root, text="")
progress_label.grid(row=6, column=0, columnspan=2, pady=5)

root.mainloop()
