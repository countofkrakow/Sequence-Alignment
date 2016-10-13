import scoring
import argparse
import sys
import random

def smith_waterman(seq1, seq2, gap_penalty=-4, sf=scoring.ScoreFunction()):
    n = len(seq1)
    m = len(seq2)
    score_matrix = [[0] * (m + 1) for i in range(n + 1)]
    trace_matrix = [[None] * (m + 1) for i in range(n + 1)]
    global_max = 0
    global_max_i = 0
    global_max_j = 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score_matrix[i][j] = max(
                    score_matrix[i-1][j-1] + sf.getScore(seq1[i-1], seq2[j-1]),
                    score_matrix[i-1][j] + gap_penalty,
                    score_matrix[i][j-1] + gap_penalty,
                    0)

            if score_matrix[i-1][j-1] + sf.getScore(seq1[i-1], seq2[j-1]) == score_matrix[i][j]:
                trace_matrix[i][j] = 'ul'
            elif score_matrix[i-1][j] + gap_penalty == score_matrix[i][j]:
                trace_matrix[i][j] = 'u'
            elif score_matrix[i][j-1] + gap_penalty == score_matrix[i][j]:
                trace_matrix[i][j] = 'l'
            # start of sequence. Must be zero
            else:
                trace_matrix[i][j] = None

            if score_matrix[i][j] > global_max:
                global_max = score_matrix[i][j]
                global_max_i = i
                global_max_j = j

    return trace_matrix, global_max_i, global_max_j, global_max, score_matrix

def calculate_pval(seq1, seq2, score, n, algo=smith_waterman):
    k = 0
    for i in range(n):
        perm_seq2 = ''.join(random.sample(seq2, len(seq2)))
        _, _, _, max, _ = algo(seq1, perm_seq2)
        if max >= score:
            k += 1
    return (k+1)/(n+1)

def trace_alignment(trace_matrix, trace_max_i, trace_max_j, seq1, seq2, sf=scoring.ScoreFunction()):
    seq1_aligned = ''
    seq2_aligned = ''
    s1i = trace_max_i
    s2i = trace_max_j
    relation = ''
    while(trace_matrix[s1i][s2i] is not None):

        if trace_matrix[s1i][s2i] == 'u':
            seq1_aligned = seq1[s1i - 1] + seq1_aligned
            seq2_aligned = '-' + seq2_aligned
            relation = ' ' + relation
            s1i -= 1
        elif trace_matrix[s1i][s2i] == 'l':
            seq1_aligned = '-' + seq1_aligned
            seq2_aligned = seq2[s2i - 1] + seq2_aligned
            relation = ' ' + relation
            s2i -= 1
        else:
            seq1_aligned = seq1[s1i - 1] + seq1_aligned
            seq2_aligned = seq2[s2i - 1] + seq2_aligned
            if (seq1[s1i - 1] == seq2[s2i - 1]):
                relation = seq1[s1i - 1] + relation
            elif sf.getScore(seq1[s1i - 1], seq2[s2i - 1]) > 0:
                relation = '+' + relation
            else:
                relation = ' ' + relation

            s1i -= 1
            s2i -= 1
    return seq1_aligned, seq2_aligned, relation

def needleman_wunsch():
    pass

def parse_args():
    # argument parsing
    parser = argparse.ArgumentParser(description="performs sequence alignment on two genetic sequences.")
    parser.add_argument('-f', '--file', action='store_true')
    parser.add_argument('-am', '--alignment_matrix', action='store_true', default=False)
    parser.add_argument('--pval_trials', default=1000, type=int)
    algo = parser.add_mutually_exclusive_group()
    algo.add_argument('-sw', '--smith_waterman', action='store_true', default=True)
    algo.add_argument('-nw', '--needleman_wunsch', action='store_true')
    parser.add_argument('sequence1')
    parser.add_argument('sequence2')
    return parser.parse_args(sys.argv[1:])

if __name__ == "__main__":
    args = parse_args()

    # file names as arguments
    if args.file:
        fname1 = args.sequence1
        fname2 = args.sequence2
        with open(fname1, 'r') as file:
            f = file.read().split('\n')
            seqname1 = f[0].split('|')[1]
            seq1=''.join(f[1:])

        with open(fname2, 'r') as file:
            f = file.read().split('\n')
            seqname2 = f[0].split('|')[1]
            seq2=''.join(f[1:])

    # raw sequences as arguments
    else:
        seqname1 = 'Sequence 1'
        seqname2 = 'Sequence 2'
        seq1 = args.sequence1
        seq2 = args.sequence2

    if args.smith_waterman:
        trace_matrix, trace_max_i, trace_max_j, global_max, alignment_matrix = smith_waterman(seq1, seq2)
        s1, s2, relation = trace_alignment(trace_matrix, trace_max_i, trace_max_j, seq1, seq2)
        pval = calculate_pval(seq1, seq2, global_max, args.pval_trials)
    elif args.needleman_wunsch:
        pass
    else:
        print("No algorithm to align with")

    line1 = trace_max_i - len(s1) + s1.count('-') + 1
    line2 = trace_max_j - len(s2) + s2.count('-') + 1
    while s1:
        print('{0:15}{1:5d}  {2:60}'.format(seqname1 + ':', line1, s1[:60]))
        print('{0:20}  {1:60}'.format(' ',relation[:60]))
        print('{0:15}{1:5d}  {2:60}'.format(seqname2 + ':', line2, s2[:60]))
        print()

        line1 = line1 + 60 - s1[:60].count('-')
        line2 = line2 + 60 - s2[:60].count('-')

        s1 = s1[60:]
        s2 = s2[60:]
        relation = relation[60:]


    print(s1 + '\n' + relation + '\n' + s2 + '\npval: ' + str(pval) + '\noptimal score: ' + str(global_max))
    if args.alignment_matrix:
        for row in alignment_matrix:
            print('\t'.join([str(e) for e in row]))
