#!/usr/bin/env python

import argparse
import re
import sys

TOLERANCE=50

def cigar_to_length(cigar):
    result = 0
    dropped = 0
    matches = re.findall(r'(\d+)(\w)', cigar)
    for match in matches:
        if match[1] in ('M', 'D', 'N', 'EQ', 'X', 'P'):
            result += int(match[0])
        else: # I, S
            dropped += int(match[0])
    return (result, dropped)

def analyze_reference(reference_file, remapped_file, mode='remapped'):
    sys.stderr.write("reading remapped...\n")
    remapped = {}
    remapped_min = 1e9
    remapped_max = 0
    remapped_reads = unmapped_reads = 0
    for idx, line in enumerate(open(remapped_file, 'r')):
        if idx % 100000 == 0:
            sys.stderr.write("read {} lines...\n".format(idx))
        if line.startswith('@'):
            continue
        fields = line.strip('\n').split('\t')

        if len(fields) > 9: # aligned read
            flag = int( fields[1] )
            if flag & 0x04 != 0 or fields[2] == '*': # unmapped
                unmapped_reads += 1
                continue
            else:
                start_pos = int(fields[3])
                read_len, dropped = cigar_to_length(fields[5]) #len(fields[9]) # not always correct, not considering cigar
                remapped[fields[0]] = (start_pos, read_len)
                remapped_min = min(remapped_min, start_pos)
                remapped_max = max(remapped_max, start_pos + read_len)
                remapped_reads += 1
    sys.stderr.write("reading remapped: done: range {} to {}. {} mapped reads. {} unmapped reads.\n".format(remapped_min, remapped_max, remapped_reads, unmapped_reads))

    written = 0
    for idx, line in enumerate(open(reference_file, 'r')):
        if idx % 100000 == 0:
            sys.stderr.write("read {} lines, wrote {}...\n".format(idx, written))
        if line.startswith('@'):
            sys.stdout.write(line)
        else:
            fields = line.strip('\n').split('\t')
            if len(fields) > 9: # aligned read
                if fields[0] in remapped: # read name isn't in remapped
                    if mode == 'remapped':
                        sys.stdout.write(line)
                        written += 1
                else:
                    if mode == 'untransformed':
                        sys.stdout.write(line)
                        written += 1
  
def analyze_donor(donor_file, remapped_file, mode='correct'):
    sys.stderr.write("reading donor...\n")
    donor = {}
    donor_min = 1e9
    donor_max = 0
    donor_actual = None
    donor_header = []
    remapped_reads = unmapped_reads = 0
    for idx, line in enumerate(open(donor_file, 'r')):
        if idx % 100000 == 0:
            sys.stderr.write("read {} lines...\n".format(idx))
        if line.startswith('@'):
            donor_header.append(line)
            if line.startswith('@SQ') and 'LN:' in line:
                # @SQ     SN:22   LN:50818468
                donor_actual = int(line.strip('\n').split('LN:')[1])
                sys.stderr.write("donor length is {}...\n".format(donor_actual))
            continue
        fields = line.strip('\n').split('\t')

        if len(fields) > 9: # aligned read
            flag = int( fields[1] )
            if flag & 0x04 != 0 or fields[2] == '*': # unmapped
                unmapped_reads += 1
                continue
            else:
                start_pos = int(fields[3])
                read_len, dropped = cigar_to_length(fields[5]) #len(fields[9]) # not always correct, not considering cigar
                donor[fields[0]] = (start_pos, read_len)
                donor_min = min(donor_min, start_pos)
                donor_max = max(donor_max, start_pos + read_len)
                remapped_reads += 1
    sys.stderr.write("reading donor: done: range {} to {}. {} mapped reads. {} unmapped reads.\n".format(donor_min, donor_max, remapped_reads, unmapped_reads))

    sys.stderr.write("reading remapped and comparing to donor...\n")
    correct_reads = incorrect_reads = unmapped_reads = not_transformed_reads = 0

    remapped_min = 1e9
    remapped_max = 0

    #correct = [0] * donor_max
    #total = [0] * donor_max

    #sys.stdout.write(''.join(donor_header))

    for idx, line in enumerate(open(remapped_file, 'r')):
        if idx % 100000 == 0:
            sys.stderr.write("read {} lines. {} correct...\n".format(idx, correct_reads))
        if line.startswith('@'):
            if line.startswith('@SQ') and 'LN:' in line:
                components = line.strip('\n').split('LN:')
                sys.stdout.write(''.join([components[0], 'LN:', str(donor_actual), '\n'])) # fix length
            else:
                sys.stdout.write(line) # include donor headers
            continue 
        fields = line.strip('\n').split('\t')

        if len(fields) > 9: # aligned read
            flag = int( fields[1] )
            if flag & 0x04 != 0 or fields[2] == '*': # unmapped
                unmapped_reads += 1
                if mode == 'unmapped':
                    sys.stdout.write(line)
                continue
            else:
                start_pos = int(fields[3]) # donor aligned position
                remapped_min = min(remapped_min, start_pos)
                remapped_max = max(remapped_max, start_pos + read_len)
                read_len, dropped = cigar_to_length(fields[5]) #len(fields[9]) # not always correct, not considering cigar
                if fields[0] in donor: # read name is in remapped
                    if abs(donor[fields[0]][0] - start_pos) < TOLERANCE: # correct
                        correct_reads += 1
                        #for pos in range(start_pos, min(remapped_max, start_pos + read_len)): # read length on donor mapping
                        #for pos in range(start_pos, min(donor_max, start_pos + donor[fields[0]][1])): # read length on remapped mapping
                        #    correct[pos] += 1
                        if mode == 'correct' or mode == 'any':
                            sys.stdout.write(line)
                    else: # incorrect
                        incorrect_reads += 1
                        if mode == 'incorrect' or mode == 'any':
                            sys.stdout.write(line)

                    #for pos in range(start_pos, min(donor_max, start_pos + read_len)):
                    #    total[pos] += 1
                    if donor_max < start_pos + read_len:
                        sys.stderr.write("WARNING: read extends past max size {}: start {} read len {}\n".format(remapped_max, start_pos, read_len))
                    del donor[fields[0]]
                else:
                    not_transformed_reads += 1
                    if mode == 'untransformed':
                        sys.stdout.write(line)
 
    sys.stderr.write("reading donor: done: correctly remapped reads {} incorrectly remapped reads {} not transformed {} unmapped reads {} total reads {}\n".format(correct_reads, incorrect_reads, not_transformed_reads, unmapped_reads, correct_reads + incorrect_reads + unmapped_reads + not_transformed_reads))

if __name__ == '__main__':
    sys.stderr.write("version {}\n".format(0.1))
    parser = argparse.ArgumentParser(description='Compare BAMs')
    parser.add_argument('--donor', help='sam files to extract reads')
    parser.add_argument('--remapped', help='sam files to extract reads')
    parser.add_argument('--reference', help='sam files to extract reads')
    parser.add_argument('--mode', help='what to write: correct, incorrect, unmapped, untransformed')
    args = parser.parse_args()
    if args.reference:
        analyze_reference(args.reference, args.remapped, args.mode)
    else:
        analyze_donor(args.donor, args.remapped, args.mode)

