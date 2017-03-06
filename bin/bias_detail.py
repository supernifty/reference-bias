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

def analyze(reference_file, donor_file, remapped_file):
    sys.stderr.write("reading remapped...\n")
    remapped = {}
    remapped_min = 1e9
    remapped_max = 0
    remapped_dropped = 0
    for idx, line in enumerate(open(remapped_file, 'r')):
        if idx % 100000 == 0:
            sys.stderr.write("read {} lines...\n".format(idx))
        if line.startswith('@'):
            continue
        fields = line.strip('\n').split('\t')

        if len(fields) > 9: # aligned read
            flag = int( fields[1] )
            if flag & 0x04 != 0 or fields[2] == '*': # unmapped
                continue
            else:
                start_pos = int(fields[3]) # position of remapped alignment
                read_len, dropped = cigar_to_length(fields[5]) #len(fields[9]) # not always correct, not considering cigar
                remapped_dropped += dropped
                remapped[fields[0]] = start_pos # record with name of read
                remapped_min = min(remapped_min, start_pos)
                remapped_max = max(remapped_max, start_pos + read_len)
    sys.stderr.write("reading remapped: done: range {} to {}. dropped due to cigar: {}\n".format(remapped_min, remapped_max, remapped_dropped))

    sys.stderr.write("reading reference...\n")
    reference = {}
    reference_min = 1e9
    reference_max = 0
    reference_actual = None
    for idx, line in enumerate(open(reference_file, 'r')):
        if idx % 100000 == 0:
            sys.stderr.write("read {} lines...\n".format(idx))
        if line.startswith('@'):
            if line.startswith('@SQ') and 'LN:' in line:
                # @SQ     SN:22   LN:50818468
                reference_actual = int(line.strip('\n').split('LN:')[1])
                sys.stderr.write("reading length is {}...\n".format(reference_actual))
            continue
        fields = line.strip('\n').split('\t')

        if len(fields) > 9: # aligned read
            flag = int( fields[1] )
            if flag & 0x04 != 0 or fields[2] == '*': # unmapped
                continue
            else:
                start_pos = int(fields[3]) # position of aligned read on reference
                read_len, dropped = cigar_to_length(fields[5]) #len(fields[9]) # not always correct, not considering cigar
                reference[fields[0]] = start_pos # record position with name of read
                reference_min = min(reference_min, start_pos)
                reference_max = max(reference_max, start_pos + read_len)
    sys.stderr.write("reading reference: done: range {} to {}\n".format(reference_min, reference_max))
 
    correct = [0] * reference_max
    total = [0] * reference_max

    sys.stderr.write("reading donor and comparing to remapped...\n")
    correct_reads = incorrect_reads = 0
    for idx, line in enumerate(open(donor_file, 'r')):
        if idx % 100000 == 0:
            sys.stderr.write("read {} lines...\n".format(idx))
        if line.startswith('@'):
            continue
        fields = line.strip('\n').split('\t')

        if len(fields) > 9: # aligned read
            flag = int( fields[1] )
            if flag & 0x04 != 0 or fields[2] == '*': # unmapped
                continue
            else:
                start_pos = int(fields[3]) # direct alignment pos
                read_len, dropped = cigar_to_length(fields[5]) #len(fields[9]) # not always correct, not considering cigar
                if fields[0] in remapped:
                    reference_start = reference[fields[0]]
                    if abs(remapped[fields[0]] - start_pos) < TOLERANCE: # correct
                        correct_reads += 1
                        for pos in range(reference_start, min(reference_max, reference_start + read_len)):
                            correct[pos] += 1
                    else: # incorrect
                        incorrect_reads += 1
                    for pos in range(reference_start, min(reference_max, reference_start + read_len)):
                        total[pos] += 1
                    if reference_max < reference_start + read_len:
                        sys.stderr.write("WARNING: read extends past reference_max {}: start {} length {}\n".format(reference_max, reference_start, read_len))
                    del remapped[fields[0]]
                    del reference[fields[0]]
 
    sys.stderr.write("reading donor: done: correctly remapped reads {} incorrectly remapped reads {}\n".format(correct_reads, incorrect_reads))

    sys.stderr.write("writing bias from {} to {}...\n".format(reference_min, reference_max))

    # write bias at every position
    sys.stdout.write('{},{},{}\n'.format('pos', 'correct', 'total'))
    not_covered = reference_min + max(0, reference_actual - reference_max)
    correct_ref = 0
    incorrect_ref = 0
    for pos in range(reference_min, reference_max):
        if pos % 1000000 == 0:
            sys.stderr.write("read {} lines...\n".format(pos))
        if total[pos] != 0:
            sys.stdout.write('{},{},{}\n'.format(pos, correct[pos], total[pos]))
            if correct[pos] > total[pos] - correct[pos]:
                correct_ref += 1
            else:
                incorrect_ref += 1
        else:
            not_covered += 1
    sys.stderr.write("writing bias from {} to {}: done\n".format(reference_min, reference_max))
    sys.stderr.write("reference stats: not_covered: {} correct majority: {} incorrect majority: {}\n".format(not_covered, correct_ref, incorrect_ref))
    sys.stderr.write("reference stats: not_covered: {} correct majority: {} incorrect majority: {}\n".format(100. * not_covered / reference_actual, 100. * correct_ref / reference_actual, 100. * incorrect_ref / reference_actual))

def analyze_donor(donor_file, remapped_file):
    sys.stderr.write("reading remapped...\n")
    remapped = {}
    remapped_min = 1e9
    remapped_max = 0
    remapped_reads = unmapped_reads = 0
    # record all the positions that are remapped
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

    sys.stderr.write("reading donor and comparing to remapped...\n")
    correct_reads = incorrect_reads = unmapped_reads = not_transformed_reads = 0
    donor_actual = None
    donor_min = 1e9
    donor_max = 0

    correct = [0] * remapped_max # correctly remapped
    total = [0] * remapped_max # correctly or incorrectly remapped

    for idx, line in enumerate(open(donor_file, 'r')):
        if idx % 100000 == 0:
            sys.stderr.write("read {} lines...\n".format(idx))
        if line.startswith('@'):
            if line.startswith('@SQ') and 'LN:' in line:
                # @SQ     SN:22   LN:50818468
                donor_actual = int(line.strip('\n').split('LN:')[1])
                donor_covered = [0] * donor_actual
                sys.stderr.write("donor length is {}...\n".format(donor_actual))
            continue
        fields = line.strip('\n').split('\t')

        if len(fields) > 9: # aligned read
            flag = int( fields[1] )
            if flag & 0x04 != 0 or fields[2] == '*': # unmapped
                unmapped_reads += 1
                continue
            else:
                start_pos = int(fields[3]) # donor aligned position
                donor_min = min(donor_min, start_pos)
                read_len, dropped = cigar_to_length(fields[5]) #len(fields[9]) # not always correct, not considering cigar
                donor_max = max(donor_max, start_pos + read_len)
                if fields[0] in remapped: # read name is in remapped
                    if abs(remapped[fields[0]][0] - start_pos) < TOLERANCE: # correct
                        correct_reads += 1
                        #for pos in range(start_pos, min(remapped_max, start_pos + read_len)): # read length on donor mapping
                        for pos in range(start_pos, min(remapped_max, start_pos + remapped[fields[0]][1])): # read length on remapped mapping
                            correct[pos] += 1
                    else: # incorrect
                        incorrect_reads += 1

                    for pos in range(start_pos, min(remapped_max, start_pos + read_len)):
                        total[pos] += 1
                    if remapped_max < start_pos + read_len:
                        sys.stderr.write("WARNING: read extends past max size {}: start {} read len {}\n".format(remapped_max, start_pos, read_len))
                    del remapped[fields[0]]
                else:
                    not_transformed_reads += 1

                # mark bases covered by donor
                for pos in range(start_pos, min(donor_max, start_pos + read_len)):
                    donor_covered[pos] += 1
 
    sys.stderr.write("reading donor: done: correctly remapped reads {} incorrectly remapped reads {} not transformed {} unmapped reads {} total reads {}\n".format(correct_reads, incorrect_reads, not_transformed_reads, unmapped_reads, correct_reads + incorrect_reads + unmapped_reads + not_transformed_reads))

    # write bias at every position which was covered by donor alignment
    sys.stderr.write("writing bias from {} to {}...\n".format(donor_min, donor_max))

    not_covered = 0
    correct_donor = 0
    incorrect_donor = 0
    sys.stdout.write('{},{},{}\n'.format('pos', 'correct', 'total'))
    for pos in range(donor_min, donor_max):
        if pos % 1000000 == 0:
            sys.stderr.write("read {} lines from remapped...\n".format(pos))
        #if pos < remapped_max and total[pos] != 0:
        if pos < donor_max and donor_covered[pos] != 0 and pos < remapped_max:
            sys.stdout.write('{},{},{}\n'.format(pos, correct[pos], total[pos]))
            if correct[pos] > total[pos] - correct[pos]:
                correct_donor += 1
            else:
                incorrect_donor += 1
        else:
            not_covered += 1
    sys.stderr.write("writing bias from {} to {}: done\n".format(donor_min, donor_max))

    # count locations with at least one correct
    sys.stderr.write("counting locations with >0 correct...\n")
    any_correct_donor = 0
    any_donor = 0
    for pos in range(donor_min, donor_max):
        if pos < remapped_max:
            if correct[pos] > 0:
                any_correct_donor += 1
            if total[pos] > 0:
                any_donor += 1

    sys.stderr.write("donor stats: not_covered: {} correct majority: {} incorrect majority: {} any correct: {} any: {}\n".format(not_covered, correct_donor, incorrect_donor, any_correct_donor, any_donor))
    sys.stderr.write("donor stats: not_covered: {} correct majority: {} incorrect majority: {} any correct: {} any: {}\n".format(100. * not_covered / donor_actual, 100. * correct_donor / donor_actual, 100. * incorrect_donor / donor_actual, 100. * any_correct_donor / donor_actual, 100. * any_donor / donor_actual ))

if __name__ == '__main__':
    sys.stderr.write("version {}\n".format(0.1))
    parser = argparse.ArgumentParser(description='Compare BAMs')
    parser.add_argument('--reference', help='sam files to extract reads')
    parser.add_argument('--donor', help='sam files to extract reads')
    parser.add_argument('--remapped', help='sam files to extract reads')
    args = parser.parse_args()
    if args.reference:
        analyze(args.reference, args.donor, args.remapped)
    else:
        analyze_donor(args.donor, args.remapped)

