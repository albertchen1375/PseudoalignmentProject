import numpy as np
import pandas as pd
import gzip
import networkx as nx
from collections import defaultdict
from itertools import product


#functions to gather data from our input files

def read_fasta_gz(file_path):
    with gzip.open(file_path, 'rt') as file:
        seqs = {}
        seq_id = None
        seq = []
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    seqs[seq_id] = ''.join(seq)
                seq_id = line[1:]
                seq = []
            else:
                seq.append(line)
        if seq_id:
            seqs[seq_id] = ''.join(seq)
    return seqs

def parse_fasta(file_path):
    transcript_sequences = {}
    with open(file_path, 'r') as file:
        transcript_id = None
        sequence = []
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if transcript_id:
                    transcript_sequences[transcript_id] = ''.join(sequence)
                transcript_id = line.split()[0][1:]  
                sequence = []
            else:
                sequence.append(line)
        if transcript_id:
            transcript_sequences[transcript_id] = ''.join(sequence)
    return transcript_sequences

#reverse complement function
def get_reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(seq))



#Build a De Bruijn Graph Using Transcriptome

def build_colored_de_bruijn_graph(isoforms, k):
    G = nx.DiGraph()
    for transcript_id, sequence in isoforms.items():
        for strand_sequence in [sequence]:
            for i in range(len(strand_sequence) - k + 1):
                kmer = strand_sequence[i:i+k]
                prefix = kmer[:-1]
                suffix = kmer[1:]
                if G.has_edge(prefix, suffix):
                    G[prefix][suffix]['isoforms'].add(transcript_id)
                else:
                    G.add_edge(prefix, suffix, isoforms={transcript_id})
    return G


#get the equivalence class of a read

def get_equivalence_class(kmers, de_bruijn_graph):
    possible_isoforms = set()

    for idx, kmer in enumerate(kmers):
        prefix = kmer[:-1]
        suffix = kmer[1:]

        if 'N' in kmer:  # Check if 'N' is present in the k-mer
            n_count = kmer.count('N')
            if n_count > 0:
                updated_possible_isoforms = set()
                for replacement_bases in product('ACGT', repeat=n_count):
                    new_bases = ''.join(replacement_bases)
                    new_kmer = ''.join(new_bases if base == 'N' else base for base in kmer)
                    prefix = new_kmer[:-1]
                    suffix = new_kmer[1:]
                    if de_bruijn_graph.has_edge(prefix, suffix):
                        if idx == 0:
                            updated_possible_isoforms.update(de_bruijn_graph[prefix][suffix]['isoforms'])
                        else:
                            possible_isoforms.intersection_update(de_bruijn_graph[prefix][suffix]['isoforms'])
                if idx == 0:
                    possible_isoforms.update(updated_possible_isoforms)

        elif de_bruijn_graph.has_edge(prefix, suffix):
            if idx == 0:
                possible_isoforms.update(de_bruijn_graph[prefix][suffix]['isoforms'])
            else:
                possible_isoforms.intersection_update(de_bruijn_graph[prefix][suffix]['isoforms'])

        else:
            possible_isoforms = set()
            break

    return possible_isoforms

#Pseudoalignment

def pseudoalign_reads_with_colored_graph(rna_seq_reads, de_bruijn_graph, k):
    equivalence_classes = defaultdict(int)

    for key in rna_seq_reads:
        read = rna_seq_reads[key]
        reverse_read = get_reverse_complement(read)

        kmers = [read[i:i+k] for i in range(len(read) - k + 1)]
        reverse_kmers = [reverse_read[i:i+k] for i in range(len(reverse_read) - k + 1)]
        
        forward_possible_isoforms = set()
        reverse_possible_isoforms = set()


        forward_possible_isoforms = get_equivalence_class(kmers, de_bruijn_graph)
        reverse_possible_isoforms = get_equivalence_class(reverse_kmers, de_bruijn_graph)
        

        possible_isoforms = set()
        possible_isoforms= possible_isoforms.union(forward_possible_isoforms, reverse_possible_isoforms)
        equivalence_classes[frozenset(possible_isoforms)] += 1

    return equivalence_classes


#Creating the Table

def format_equivalence_class_counts(equivalence_classes):
    formatted_counts = []
    for eq_class, count in equivalence_classes.items():
        if len(eq_class) == 0:
            formatted_counts.append((count, 0, "NA"))
        else:
            formatted_counts.append((count, len(eq_class), ','.join(eq_class)))

    # Sort the list of tuples by the count in descending order
    formatted_counts.sort(key=lambda x: x[0], reverse=True)

    # Convert to a DataFrame
    df = pd.DataFrame(formatted_counts, columns=['Count', 'Number of items in equivalence class', 'Isoforms In equivalence Class'])
    return df


# Main

def main():
    isoforms = parse_fasta("chr11_transcriptome.fasta")
    rna_seq_reads = read_fasta_gz("reads.fasta.gz")


    k = 31
    
    colored_de_bruijn_graph = build_colored_de_bruijn_graph(isoforms, k)

    equivalence_classes = pseudoalign_reads_with_colored_graph(rna_seq_reads, colored_de_bruijn_graph, k)
    data_frame = format_equivalence_class_counts(equivalence_classes)

    data_frame.to_csv('projectdata.csv', index=False)
    print(data_frame)

if __name__ == "__main__":
    main()