# Whisper2

[![GitHub downloads](https://img.shields.io/github/downloads/refresh-bio/whisper/total.svg?style=flag&label=GitHub%20downloads)](https://github.com/refresh-bio/Whisper/releases)
[![GitHub Actions CI](../../actions/workflows/main.yml/badge.svg)](../../actions/workflows/main.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


## Quick start

```bash
# download and unpack E.coli str. K-12 substr. MG1655 reference genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/904/425/475/GCF_904425475.1_MG1655/GCF_904425475.1_MG1655_genomic.fna.gz
gzip -k -d GCF_904425475.1_MG1655_genomic.fna.gz 

# download NovaSeq 6000 reads from E.coli MG1655 IR-10-94 population sequencing
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/030/SRR10051130/SRR10051130_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/030/SRR10051130/SRR10051130_2.fastq.gz

# clone and build Whisper
git clone https://github.com/refresh-bio/Whisper
make -C Whisper/src

# build an index named ecoli for the reference genome
mkdir index
mkdir temp
./Whisper/src/whisper-index ecoli GCF_904425475.1_MG1655_genomic.fna ./index ./temp

# map reads to the reference genome and store the results in mappings.sam
./Whisper/src/whisper -rp -out mappings -temp ./temp/ ./index/ecoli  SRR10051130_1.fastq.gz SRR10051130_2.fastq.gz
```

Note, that Whisper was optimized for processing data from sequencing large samples with high coverage (e.g. human). 

## Installation and configuration

Whisper comes with a set of precompiled binaries for Windows and Linux. They can be found under *Releases* tab.

The software can be also built from the sources distributed as:
* Visual Studio 2017 solution for Windows,
* MAKE project (G++ 7.2 required) for Linux.

NOTE

Linux systems limit number of files that can be opened by the process. Make sure this limit is sufficient for Whisper, which requires:

`num_files = num_bins (384 by default) + num_threads + 2 * total_size_of_FASTQ_in_GB` 

E.g., if sample reads in uncompressed FASTQ format have 100GB and the processing is done by 16 CPU threads, Whisper opens 600 files. To change the limit use `ulimit` Linux command:

`ulimit -n 600`

## Usage

The preliminary step of the analysis, performed only once for a given reference genome, is construction of an index. The index may be then used for mapping reads from different samples to the reference. 

### Indexing reference genome

Indexing can be executed in two wariants, depending on the representation of the reference (single versus multiple FASTA files):

`whisper-index <index_name> <ref_seq_file_name> <dest_dir> <temp_dir>`

`whisper-index <index_name> @<ref_seq_files_name> <dest_dir> <temp_dir>`

Parameters:

* `index_name` - name of index (output),
* `ref_seq_file_name` - reference sequence in the FASTA format (input)
* `ref_seq_files_name` - file with a list of multiple reference FASTA files (input)
* `dest_dir` - existing destination directory where index files should be placed (output)
* `temp_dir`- existing temporary directory (output)

Examples: 
* `whisper-index hg38-chr20 chr20.fa index-dir temp-dir` 

Generates index named `hg38-chr20` for `chr20.fa` reference sequence and places it in `index-dir` directory.

* `whisper-index hg38 @hg38.list index-dir temp-dir` 

Generates index named `hg38` for all FASTA files listed in `hg38.list` file and places it in `index-dir` directory.

### Mapping reads

   `whisper [options] <index_name> @<files>`
   
   `whisper [options] <index_name> file_se`
   
   `whisper [options] <index_name> file_pe_1 file_pe_2`
   
Parameters:
 
 * `index_name`   - name of the index (as created by whisper-index)
 * `files`        - name of the file containing list of FASTQ files with seq. reads
 * `file_se`      - FASTQ file (single-end)
 * `file_pe_[1|2]` - FASTQ files (paired-end)
 
Options:

 * `-b <value>` - no. of temporary files (minimum: 100, default: 384)
 * `-clipping-distance <value>` - no. of sigmas for max. additional distance in clipping (default: 14)
 * `-d[fr/ff/rf]` - mapping orientation (default: -dfr (forward - reverse)
 * `-dist_paired <value>` - max. distance for paired read (default: 1000)
 * ` -e <value>` - max. no of errors (default: auto)
 * `-e-paired <value>` - max. fraction of errors in paired read (default: 0.09)
 * `-enable-boundary-clipping <value>` - enable clipping at boundaries when a lot of mismatches appears (default: 1)
 * `-enable-mapping_indels <value>` - enable looking for long indels during mapping stages (default: 1)
 * `-enable-short-indel-refinement <value>` - enable short indel refinement after mapping (default: 1)
 * `-enable-short-reads <value>` - enable reads shorter than 90% of the longest reads (default: 0)
 * `-filter <value>` - store only mappings for given chromosome (default: )
 * `-gap-del-open <value>` - score for gap (del) open (default: -5)
 * `-gap-del-extend <value>` - score for gap (del) extend (default: -0.4)
 * `-gap-ins-open <value>` - score for gap (ins) open (default: -5)
 * `-gap-ins-extend <value>` - score for gap (ins) extend (default: -0.4)
 * `-gzipped-SAM <value>` - gzip compression level of SAM/BAM, 0 - no compression (default: 0)
 * `-high-confidence-sigmas <value>` - (default: 4)
 * `-hit-merging-threshold <value>` - minimal distance between different mappings (default: 12)
 * `-hit-merging-wrt-first <value>` - calculate distance in marged group w.r.t. first (default: 1)
 * `-m[f/s/a]` - mode: first stratum/second stratum/all strata (default: first stratum)
 * `-mask-lqb <value>` - mask bases of quality lower than value (default: 0)
 * `-max-indel-len <value>` - max. indel length (default: 50)
 * `-min-clipped-factor <value>` - mask bases of quality lower than value (default: 1)
 * `-out <name>` - name of the output file (default: whisper)
 * `-penalty-saturation <value>` - no. of sigmas for max. penalty in matching pairs (default: 7)
 * `-rg <read_group>` - complete read group header line, '\t' character will be converted to a TAB in the output SAM while the read group ID will be attached to every read (example line: '@RG     ID:foo  SM:bar')
 * `-r[s|p]` - single or paired-end reads (default: single)
 * `-score-discretization-threshold` (default: 0.5)
 * `-score-clipping <value>` score for clipping (default: -6)
 * `-score-match <value>` - score for matching symbol (default: 1)
 * `-score-mismatch <value>` - score for mismatching symbol (default: -5)
 * `-sens <value>` - turn on/off sensitive mode (default: 1)
 * `-sens-factor <value>` - sensitivity factor (default: 2.5)
 * `-stdout` - use stdout to store the output (default: 0)
 * `-store-BAM` - turn on saving in BAM (default: 0)
 * `-t <value>` - no. of threads (0-adjust to hardware) (default: 0)
 * `-temp <name>` - prefix for temporary files (default: ./whisper_temp_)
 * `-x` - load complete suffix arrays in main memory (default: 0)
  
  
Examples:

* `whisper -out result.sam -rp -t 12 hg38 reads_1.fq reads_2.fq`

Maps paired-end reads from `reads_1.fq` and `reads_2.fq` FASTQ files using `hg38` index. Computations are distributed over 12 threads,  results are stored in `result.sam` file.

* `whisper hg38 @reads_se.list`

Maps single-end reads from FASTQ files listed in `reads_se.list` file using `hg38` index. The example contents of `reads_se.list` file:
```
readsA
readsB
readsC
...
```

* `whisper -rp hg38 @reads_pe.list`

Maps paired-end reads from FASTQ files listed in `reads_pe.list` file using `hg38` index. The example contents of `reads_pe.list` file:
```
readsA_1;readsA_2
readsB_1;readsB_2
readsC_1;readsC_2
...
```


## Citing

[Deorowicz, S., Debudaj-Grabysz, A., Gudyś, A., Grabowski, S. (2018) Whisper: Read sorting allows robust mapping of sequencing data, Bioinformatics, 35(12):2043–2050](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/bty927/5165374?redirectedFrom=fulltext)


[Deorowicz, S., Gudyś, A. (2021) Whisper 2: Indel-sensitive short read mapping, SoftwareX, 14:100692](https://doi.org/10.1016/j.softx.2021.100692) 

 
