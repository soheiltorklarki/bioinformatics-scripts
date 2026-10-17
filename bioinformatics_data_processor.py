# Title: Genomic Sequence & Gene Expression Data Processor
# Author: Soheil Tork Larki
# Description: Python script for parsing FASTA sequences and processing expression metrics

import pandas as pd
import numpy as np

def calculate_gc_content(sequence: str) -> float:
    """Calculates GC content percentage of a DNA sequence."""
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0.0

def process_gene_expression_data(df: pd.DataFrame) -> pd.DataFrame:
    """Filters high-expression genes and calculates log2 fold change."""
    df['Log2_FoldChange'] = np.log2(df['Treated_TPM'] / df['Control_TPM'] + 1e-5)
    return df[df['Log2_FoldChange'].abs() > 1.0]

if __name__ == "__main__":
    # Test DNA GC Content
    sample_seq = "ATGCGATCGATCGATCGATCGATCGATC"
    print(f"Sample Sequence: {sample_seq}")
    print(f"GC Content: {calculate_gc_content(sample_seq):.2f}%")

    # Sample Gene Expression Processing
    data = {
        'Gene_ID': ['GENE1', 'GENE2', 'GENE3', 'GENE4'],
        'Control_TPM': [10.5, 2.0, 45.1, 0.5],
        'Treated_TPM': [50.2, 2.1, 12.0, 8.4]
    }
    df = pd.DataFrame(data)
    filtered_df = process_gene_expression_data(df)
    print("\nDifferentially Expressed Genes:")
    print(filtered_df)