# Reference Bias
[![Build Status](https://travis-ci.org/supernifty/reference-bias.svg?branch=master)](https://travis-ci.org/supernifty/mgsa)
[![Coverage Status](https://coveralls.io/repos/supernifty/reference-bias/badge.svg?branch=master&service=github)](https://coveralls.io/github/supernifty/reference-bias?branch=master)

Methods for quantifying reference bias in resequencing

## Installation

* git clone https://github.com/supernifty/reference-bias
* cd reference-bias
* sudo python setup.py install

## Dependencies

* bedtools
* bwa or bowtie2
* progressiveMauve
* samtools

## Usage

### Example
python ./bin/calculate_bias.py --donor ./data/e-coli-mg1655.fasta --reference ./data/phylogeny/NC_002695.1.fasta --tmpdir ./tmp ./data/SRR892241-ecoli-mg1655_1.fastq.gz --align bwa > results/calculate_bias.mg1655.NC_002695.1.bwa.txt

### Parameters

* donor: fasta file specifying the donor. Genome closely resembling the specified short read archive.
* reference: fasta file specifying the reference. The bias associated with this genome choice will be calculated.
* tmpdir: place to store generated files
* align: choose either bwa or bowtie2
* start: use this to resume a job from a given stage
* job: if resuming a job, specify the job id
* donorbam: specify a precalculated donor bam file
* donorsam: specify a precalculated donor sam file

### Generated output

* The headline number is at the bottom of what is written to stdout, the estimated bias.

