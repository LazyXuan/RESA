Support for the RNA Element Selection Assay (Yartseva et al, 2017)

## Quick start:
1. Prepare the following:
  * A BED12 file with the annotations for your (UTR) regions of interest
  * A corresponding FASTA file (in the correct strand orientation) of these regions (can be produced using `bedtools fastaFromBed`)
  * RESA RNA-Seq fastq files or equivalent, for the experimental and the control conditions (which we will name "expt" and "ctrl")
  
2. Map your sequencing reads to the **_whole genome_** using your favorite aligner (e.g., Tophat) and output as bamfiles

3. Calculate coverage using `resa_coverage.py`, which will produce a pickled data file (in this example, it will be called `my_resa.coverage.pkl`):

    `python3 resa_coverage.py coverage --samples_as_bam -u <utr_file> -f <fasta_file> -l expt,ctrl -o my_resa <expt_bamfile> <ctrl_bamfile>`

4. Identify responsive regions using `resa_discovery.py`, which will produce a data file of regions and corresponding bed and fasta files

    `python3 resa_discovery.py regions -e expt -c ctrl -d my_resa.coverage.pkl -o my_resa`
