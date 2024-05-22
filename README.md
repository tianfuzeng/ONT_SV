# ONT_SV: Identifying Somatic Structural Variations with Long-Read Nanopore Sequencing

A robust method for identifying **somatic structural variations (SVs)** and **HBV integration events**, facilitating the inference of SV types induced by HBV integration.

## Workflow for Detecting HBV Integration Events

To detect HBV integration using long-read sequencing, we have developed a comprehensive workflow with the following steps:

1. **Extraction of Chimeric Reads**: Identify reads that contain sequences from both HBV and human genomes.
2. **Extraction of Insertion (INS) and Clipped Sequences**: Obtain INS sequences and clipped sequences from all long reads.
3. **Re-alignment using BLAST**: Re-align INS and clipped sequences to a custom reference genome using BLAST.
4. **Determination of HBV Integration Breakpoints**: Pinpoint the exact breakpoints where HBV DNA integrates into the human genome.

## Identifying Somatic SVs

To identify somatic SVs, we first merge the alignment results of each tumor sample with their paired blood samples. Using the widely recognized variant caller Sniffles, which clusters long reads to identify those supporting the same SV, we proceed as follows:

- **SV Calling**: A somatic SV candidate is identified if the supporting reads are absent in blood samples.
- **Quality Control and Manual Inspection**: After thorough quality control and manual inspection, we identify five types of SVs:
  - **Insertions (INS)**
  - **Deletions (DEL)**
  - **Tandem Duplications (DUP)**
  - **Inversions (INV)**
  - **Translocations (TRA)**

This approach ensures accurate identification and classification of SVs, providing valuable insights into the genomic alterations induced by HBV integration.
