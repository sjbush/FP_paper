# FP_paper

This repository contains the 6 benchmarking datasets used in the paper "Generalisable characteristics of false positive bacterial variant calls".

Each dataset comprises a triplet of files: a fasta, a VCF, and a BED. The associated fastqs are sourced from the [FDA-ARGOS reference collection](https://www.nature.com/articles/s41467-019-11306-6) and hosted on the SRA, as detailed below.

**Background**

To benchmark variant calling pipelines, we first required a truth set of variant calls against which their output could be compared. To do so, we can obtain sequenced reads from one sample with a closed genome (**_fastq_**), map them to another closed genome (**_fasta_**), call variants from these alignments, and then compare that set of calls to the set of calls made using pairwise whole genome alignment, which for our purpose we consider the truth set (**_VCF_**). The pipeline VCF is then compared to the truth set VCF using the haplotype comparison tool [hap.py](https://github.com/Illumina/hap.py), with analysis restricted to a list of higher-confidence regions using parameter -f (**_BED_**).

This repository contains fasta, VCF and BED files created for use with six fastqs sourced from the [FDA-ARGOS reference collection](https://www.nature.com/articles/s41467-019-11306-6), a public database of microbial genomes for diagnostic use: _Enterococcus faecalis_ (accession FDAARGOS_338), _Escherichia coli_ (FDAARGOS_536), _Francisella tularensis_ (FDAARGOS_598), _Salmonella enterica_ (FDAARGOS_687), _Bacillus anthracis_ (FDAARGOS_700) and _Mycobacterium tuberculosis_ (FDAARGOS_751). For further details, see the [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/231221).

Note that one sample, FDAARGOS_700, does not have an associated BED file. This is because all positions in its associated fasta could be called; none required exclusion.

**Location of sequencing data**

For the aforementioned samples, we require the original FDA-ARGOS fastqs: 150bp paired-end Illumina HiSeq4000 sequencing reads. These have corresponding SRA run accession IDs [SRR5448651](https://www.ebi.ac.uk/ena/browser/view/SRR5448651), [SRR8180486](https://www.ebi.ac.uk/ena/browser/view/SRR8180486), [SRR8283296](https://www.ebi.ac.uk/ena/browser/view/SRR8283296), [SRR9163323](https://www.ebi.ac.uk/ena/browser/view/SRR9163323), [SRR9171533](https://www.ebi.ac.uk/ena/browser/view/SRR9171533), and [SRR9176751](https://www.ebi.ac.uk/ena/browser/view/SRR9176751), respectively. For the paper, each fastq was randomly down-sampled to 1,000,000 reads using [seqtk ‘sample’ v1.3](https://github.com/lh3/seqtk) with seed 42.
