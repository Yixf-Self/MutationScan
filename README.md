## MutationScan
MutationScan is a project implemented by pure Rust language for searching Muatation point direct in Pair of Fastq(r1.fastq and r2.fastq) files or one Fastq file.

## Install
> 1. Clone this project. `git clone https://github.com/OpenGene/MutationScan.git`  
2. Install rust compiler.  
On Unix/Linux, run `curl https://sh.rustup.rs -sSf | sh`  
On Windows, download and run the [rustup-init.exe for the i686-pc-windows-gnu target](https://static.rust-lang.org/rustup/dist/i686-pc-windows-gnu/rustup-init.exe)  
3. Update Rust nightly version.  
`rustup update nightly`  
`rustup default nightly`  
4. Build the project.  
`cd MutationScan`  
`cargo build --release`  

**PS: Ingore the built Warning**

## Config

1. Prepare the fastq file(pair or one is ok).
2. Set snpfile.
```
T790M(0) 25
CTCACCTCCACCGTGCAGCTCATCACGCAGCTCATGCCCTTCGGCTGCCTC
T790M(1) 25
CTCACCTCCACCGTGCAGCTCATCATGCAGCTCATGCCCTTCGGCTGCCTC
```
The first line is a snp header label, and `25` is a Mutation loci in the sequence. 
The count of location is start with zero.
CTCACCTCCACCGTGCAGCTCATCACGCAGCTCATGCCCTTCGGCTGCCTC  
0123456789  
The second line is `ATCG` sequence

ps: If you have more the one Mutation loci in the sequence, your snp header label can set like this `T790M(0) 25 27 40`.

## Example
1. Run with pait fastq files(R1.fastq and R2.fastq).

**`MutationScan r1_file_path r2_file_path snpfile_path [threshold_value] [fusion_num]`**

`./MutationScan/target/release/MutationScan ./R1.fastq ./R2.fastq snpfile`
Set threshold_value = 2, fusion_num =4 by default.

`./MutationScan/target/release/MutationScan ./R1.fastq ./R2.fastq snpfile 3`

`./MutationScan/target/release/MutationScan ./R1.fastq ./R2.fastq snpfile 3 6`

ps: On window, you must add .exe with `MutationScan`

2. Run with one fastq file(R.fastq).

**`MutationScan r_file_path snpfile_path [threshold_value]`**

`./MutationScan/target/release/MutationScan R.fastq snpfile`
Set threshold_value = 2 by default.

`./MutationScan/target/release/MutationScan R.fastq snpfile 3`

## License

The MIT License (MIT)

Copyright (c) 2016 OpenGene

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
