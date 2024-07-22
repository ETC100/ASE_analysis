# ASE_analysis
The allelic specific expression for Muscovy ducks  

## What is ASE (allele specific expression) ?
Allelic specific expression (ASE) refers to the phenomenon where alleles of a gene are expressed at different levels in a diploid organism. In a diploid organism, each individual inherits two alleles (one from each parent) for each gene. ASE occurs when one allele is expressed at a higher level than the other allele, leading to allele-specific differences in gene expression.  

ASE can be influenced by various factors, including genetic variations such as single nucleotide polymorphisms (SNPs), epigenetic modifications, transcription factor binding, and chromatin structure. It is an important aspect of gene regulation and can contribute to phenotypic diversity and disease susceptibility.  

ASE can be studied using techniques such as RNA sequencing (RNA-seq) that can differentiate between alleles based on their sequence variations. By analyzing the expression levels of alleles within an individual, researchers can uncover allele-specific patterns of gene expression and gain insights into the regulatory mechanisms underlying gene expression variability.  

## How to identify ASE?
I recommend this reference: Perspectives on Allele-Specific Expression, it provides detailed information for causes and identification of ASE.
The first step is to find suitable method to obtain paternal-maternal origin read counts table. Here is my sharing of experiences.
1. Use hisat2/star directly
When species with distant genetic distances are used for hybridization, this method can be utilized. For example, when studying mules (hybrids of horses and donkeys), one can directly merge the genomes of both species and align them using an RNA-seq aligner. Subsequently, one can remove suboptimal alignments and proceed with counting.
2. Split bam according to SNP
This method is currently the most widely used method in situations involving the same species but different strains. For instance, when studying mice, one can first align the RNA-seq data using an RNA-seq aligner to obtain a BAM file, and then use strain-specific SNP information to phase the BAM file.


The second step is to find a viable statistical method for ASE identification.
The main methods includes:
1. Binomial distribution test
We assumed that there

2. Bayes test
  
3. Proportional distribution test
This method has been used in Qllelic for ASE identification based on techinical replicates. You can refer to the following link. Noted that, Qllelic can also be used for biological replicates, however, the author dont recommend to do so.
https://github.com/gimelbrantlab/Qllelic/issues
