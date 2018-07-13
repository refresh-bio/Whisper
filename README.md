# Whisper

## Installation and configuration

Whisper comes with a set of precompiled binaries for Windows and Linux. They can be found under *Releases* tab.

The software can be also built from the sources distributed as:
* Visual Studio 2015 solution for Windows,
* MAKE project (G++ 6.2 required) for Linux.

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
 * `-d[fr/ff/rf]` - mapping orientation (default: -dfr (forward - reverse)
 * `-dist_paired <value>` - max. distance for paired read (default: 1000)
 * `-e <value>` - max. fraction of errors in % (default: 0.04, max: 0.05)
 * `-e-paired <value>` - max. fraction of errors in paired read (default: 0.06)
 * `-enable-boundary-clipping <value>` - enable clipping at boundaries when a lot of mismatches appears (default: 0)
 * `-filter <value>` - store only mappings for given chromosome (default: )
 * `-gap-open <value>` - score for gap open (default: -6)
 * `-gap-extend <value>` - score for gap extend (default: -1)
 * `-gzipped-SAM-level <value>` - gzip compression level of SAM, 0 - no compression (default: 0)
 * `-hit-merging-threshold <value>` - minimal distance between different mappings (default: 12)
 * `-high-confidence-sigmas <value>` - (default: 4)
 * `-hit-merging-wrt-first <value>` - calculate distance in marged group w.r.t. first (default: 1)
 * `-m[f/s/a]` - mode: first stratum/second stratum/all strata (default: first stratum)
 * `-mask-lqb <value>` - mask bases of quality lower than value (default: 0)
 * `-out <name>` - name of the output file (default: whisper)
 * `-penalty-saturation <value>` - no. of sigmas for max. penalty in matching pairs (default: 7)
 * `-rg <value>` - complete read group header line, (example: '@RG\\tID:foo\\tSM:bar'), '\t' character will be converted to a TAB in the output SAM/BAM, while the read group ID will be attached to every read
 * `-r[s|p]` - single of paired-end reads (default: single)
 * `-score-discretization-threshold <value>` - (default: 0.5)
 * `-score-match <value>` - score for matching symbol (default: 1)
 * `-score-clipping <value>` score for clipping (default: -10)
 * `-score-mismatch <value>` - score for mismatching symbol (default: -5)
 * `-sens <value>` - turn on/off sensitive mode (default: 1)
 * `-sens-factor <value>` - sensitivity factor (default: 3)
 * `-stdout` - use stdout to store the output
 * `-store-BAM` - turn on saving in BAM
 * `-t <value>` - no. of threads (0-adjust to hardware) (default: 0)
 * `-temp <name>` - prefix for temporary files (default: ./whisper_temp_)
 * `-x <value>` - load complete suffix arrays in main memory (default: 0)`
  
  
Examples:

* `whisper -out result.sam -t 12 hg38 reads_1.fq reads_2.fq`

Maps paired-end reads from `reads_1.fq` and `reads_2.fq` FASTQ files using `hg38` index. Computations are distributed over 12 threads,  results are stored in `result.sam` file.

* `whisper hg38 @reads.list`

Maps reads from FASTQ files listed in `reads.list` file using `hg38` index. File `reads.list` may contain either single-end reads:

```
readsA
readsB
readsC
...
```

or paired-end reads:
```
readsA_1;readsA_2
readsB_2;readsB_2
readsC_3;readsC_2
...
```


## Citing

[Deorowicz, S., Debudaj-Grabysz, A., Gudyś, A., Grabowski, S. (2017) Whisper: Read sorting allows robust mapping of sequencing data, bioRxiv](https://www.biorxiv.org/content/early/2017/12/28/240358) DOI: 10.1101/240358  



 
