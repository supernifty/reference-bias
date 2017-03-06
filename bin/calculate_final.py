#!/usr/bin/env python
import sys

correct_cov = int(sys.stdin.readline().strip('\n').split(':')[-1])
any_cov = int(sys.stdin.readline().strip('\n').split(':')[-1])
untransformed = int(sys.stdin.readline().strip('\n').split(':')[-1])
incorrect = int(sys.stdin.readline().strip('\n').split(':')[-1])
donor_length = int(sys.stdin.readline().strip('\n').split(':')[-1])
donor_cov = int(sys.stdin.readline().strip('\n').split(':')[-1])

sys.stdout.write('{:.1f} {:.1f} {:.1f}\n'.format( 100. * (donor_cov - (any_cov + untransformed)) / donor_cov, 100. * (donor_cov - any_cov) / donor_cov, 100. * (donor_cov - correct_cov) / donor_cov ))
