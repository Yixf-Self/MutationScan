## MutationScan
MutationScan is a project implemented by pure Rust language for searching Muatation point directly in fastq files (support single-end or pair-end)

## Install
> 1. Clone this project. 
`git clone https://github.com/OpenGene/MutationScan.git`  
2. Install rust compiler.  
On Unix/Linux, run 
`curl https://sh.rustup.rs -sSf | sh`  
On Windows, download and run the [rustup-init.exe for the i686-pc-windows-gnu target](https://static.rust-lang.org/rustup/dist/i686-pc-windows-gnu/rustup-init.exe)  
3. Update Rust nightly version.  
`rustup update nightly`  
`rustup default nightly`  
4. Build the project.  
`cd MutationScan`  
`cargo build --release`  
**PS: Ingore the built Warning**

## Config
1. Prepare the fastq file(s).
2. Set snpfile.
```
T790M(0) 25
CTCACCTCCACCGTGCAGCTCATCACGCAGCTCATGCCCTTCGGCTGCCTC
T790M(1) 25
CTCACCTCCACCGTGCAGCTCATCATGCAGCTCATGCCCTTCGGCTGCCTC
```
The first line is a snp header label, and `25` is the position of mutation (0-based from left). 
The second line is `ATCG` sequence

ps: If you have more than one mutations in the sequence, your snp header label can be set like `T790M(0) 25 27 40`.

## Examples
***1. Run with pair-end fastq files(R1.fastq and R2.fastq).***

**`MutationScan r1_file_path r2_file_path snpfile_path [threshold_value] [fusion_num]`**

`./MutationScan/target/release/MutationScan ./R1.fastq ./R2.fastq ./snpfile` 

threshold_value = 2, fusion_num =4 by default.

`./MutationScan/target/release/MutationScan ./R1.fastq ./R2.fastq ./snpfile 3`

`./MutationScan/target/release/MutationScan ./R1.fastq ./R2.fastq ./snpfile 3 6`

ps: On Windows, you must add .exe with `MutationScan`

***2. Run with one fastq file(R.fastq).***

**`MutationScan r_file_path snpfile_path [threshold_value]`**

`./MutationScan/target/release/MutationScan R.fastq ./snpfile` 

threshold_value = 2 by default.

`./MutationScan/target/release/MutationScan R.fastq ./snpfile 3`
