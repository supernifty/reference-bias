# Reference Bias
[![Build Status](https://travis-ci.org/supernifty/reference-bias.svg?branch=master)](https://travis-ci.org/supernifty/mgsa)
[![Coverage Status](https://coveralls.io/repos/supernifty/reference-bias/badge.svg?branch=master&service=github)](https://coveralls.io/github/supernifty/reference-bias?branch=master)

Methods for quantifying reference bias in resequencing

## Installation

* git clone https://github.com/supernifty/reference-bias
* cd reference-bias
* sudo python setup.py install

## Dependencies

* bedtools - http://bedtools.readthedocs.org/en/latest/
* bwa or bowtie2
* datamash - https://www.gnu.org/software/datamash/
* progressiveMauve - http://darlinglab.org/mauve/mauve.html
* samtools (v1.2+ recommended) - https://github.com/samtools

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

Detailed description:
* Donor reads total: total reads reported by samtools flagstat
* Donor reads mapped: mapped reads reported by samtools flagstat
* Donor reads % mapped: mapped / total
* Reference reads total: total reads reported by samtools flagstat
* Reference reads mapped: mapped reads reported by samtools flagstat
* Reference reads % mapped: mapped / total
* Remapped reads total: total reads reported in remapped alignment
* Remapped reads mapped: mapped reads reported in remapped alignment
* Remapped reads % mapped: mapped / total
* Mismatched reads total: 2499
* Mismatched reads mapped: 2431
* Mismatched reads % mapped: 97.3
* Notcovered reads total: 27857 (2.9%)
* Notcovered reads mapped: 27857 (2.9%)
* Notcovered reads % mapped: 100.0
* Almost correct reads total: 27
* Almost correct reads mapped: 27
* Almost correct reads % mapped: 100.0

* Mapped to correct location: 938065 (99.7%)
* Mapped correctly or within 50bp: 938092 (99.7%)
* Mapped incorrectly <50bp: 27 (0.0%)
* Mapped incorrectly >50bp: 2431 (0.3%)

* Donor not covered: 3105686045 (1246.01%)
* Donor not covered with mauve target: 3212456 (0.10%)
* Donor covered: -2856435424 (-1146.01%)
* Donor gaps: 74450
* Donor max gap: 249250621
* Donor mean coverage: 2.891958750152
* Donor max coverage: 53
* Reference not covered: 3056743057 (2284.6%)
* Reference covered: -2922945635 (-2184.6%)
* Reference gaps: 74833
* Reference max gap: 248956422
* Remapped not covered: 20706620 (8.3%)
* Remapped covered: 228544001 (91.7%)
* Remapped gaps: 71284
* Remapped max gap: 16865934
* Remapped mean coverage: 2.8927606170652
* Remapped max coverage: 53

* Bases affected by mismatch: 92699
* Max mismatch coverage: 13

* Off target bases: 987831
* Max coverage of off target: 19

* Mapped bases: 33723857 (13.5%)
* Not mapped bases: 215526764 (86.5%)
* Mapped blocks: 5
* Covered reads: 940523 (97.1)
* Covered partial reads: 941196 (97.1)
* Not mapped reads: 28530 (2.9)

* Donor not covered by direct alignment: 3105686045 (1246.01%)
* Best case loss from reference coverage: 3172196256 / 249250621: 1272.7%
* Best case loss from remapping: 215526764 / 249250621: 86.5%
* Loss after remap coverage: 20706620 / 249250621: 8.3%
* Loss due to remap: -3151489636 / 249250621: -1264.4%
* Potential mismatch impact: 92699 / 249250621: 0.0%
* Off target: 987831 / 249250621: 0.4%
* Donor not covered with mauve target: 3212456 (0.10%)
* ESTIMATED BIAS: 6.6 -> 7.0 -> 8.3
