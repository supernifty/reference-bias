# Reference Bias
[![Build Status](https://travis-ci.org/supernifty/reference-bias.svg?branch=master)](https://travis-ci.org/supernifty/mgsa)
[![Coverage Status](https://coveralls.io/repos/supernifty/reference-bias/badge.svg?branch=master&service=github)](https://coveralls.io/github/supernifty/reference-bias?branch=master)

Methods for quantifying reference bias in resequencing

## Installation

* git clone https://github.com/supernifty/reference-bias
* cd reference-bias
* sudo python setup.py install

## Dependencies

* andi - (phylogeny only) - http://github.com/evolbioinf/andi/
* bedtools - http://bedtools.readthedocs.org/en/latest/
* bwa or bowtie2
* datamash - https://www.gnu.org/software/datamash/
* progressiveMauve - http://darlinglab.org/mauve/mauve.html
* samtools (v1.2+ recommended) - https://github.com/samtools

## Usage

### Bias for a single genome

#### Example
```
python ./bin/calculate_bias.py --donor ./data/e-coli-mg1655.fasta --reference ./data/phylogeny/NC_002695.1.fasta --tmpdir ./tmp ./data/SRR892241-ecoli-mg1655_1.fastq.gz --align bwa > results/calculate_bias.mg1655.NC_002695.1.bwa.txt
```

#### Parameters

* donor: fasta file specifying the donor. Genome closely resembling the specified short read archive.
* reference: fasta file specifying the reference. The bias associated with this genome choice will be calculated.
* tmpdir: place to store generated files
* align: choose either bwa or bowtie2
* start: use this to resume a job from a given stage
* job: if resuming a job, specify the job id
* donorbam: specify a precalculated donor bam file
* donorsam: specify a precalculated donor sam file

#### Generated output

* The headline number is at the bottom of what is written to stdout, the estimated bias.

### Multiple genomes with uncertainty

#### Example

`python ./bin/draw_bias.py ./data/ecoli.map ./results/calculate_bias.* > ../results/graph.pdf`

#### Parameters

* map file: comma separated list of fasta file, name, pathogenicity
* list of generated result files

#### Generated output

* A graph showing each measured genome and its estimated bias

### Phylogeny 

#### Example steps

```
andi *.fasta > ./results/ecoli.distances
python ./bin/fix_andi_names.py ./data/ecoli65.list < ./results/ecoli.distances > ./results/ecoli.distances.fixed
python ./bin/draw_phylogeny.py ./data/ecoli.map ./results/ecoli.distances.fixed ./results/calculate_bias.* > maketree.R
```

* run the resulting R code to generate the phylogeny
* gradient.pdf is also generated which represents the colour scale

### Generate bias matrix

#### Example

```
python ./bin/calculate_bias_matrix.py ./results/ecoli_matrix.out < ./data/ecoli_matrix.txt
python ./bin/draw_matrix.py ./data/ecoli_matrix.map < ./results/ecoli_matrix.out > matrix.pdf
```

#### Parameters
* calculate_bias_matrix takes a set of readsets and genomes and calculates bias for each pairwise combination, then writes this to the first argument
* draw_matrix takes this generated output as stdin and uses the names from the matrix file to generate a graph.

