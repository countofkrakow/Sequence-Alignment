import align
fnames = ['P15172', 'P17542', 'P10085', 'P16075', 'P13904', 'Q90477', 'Q8IU24', 'P22816', 'Q10574', 'O95363']

for i, f1, in enumerate(fnames):
    fname1 = f1 + '.fasta'
    for f2 in fnames[i+1:]:
        fname2 = f2 + '.fasta'

        with open(fname1, 'r') as file:
            f = file.read().split('\n')
            seq1=''.join(f[1:])

        with open(fname2, 'r') as file:
            f = file.read().split('\n')
            seq2=''.join(f[1:])

        trace_matrix, trace_max_i, trace_max_j, global_max, alignment_matrix = align.smith_waterman(seq1, seq2)
        s1, s2, relation = align.trace_alignment(trace_matrix, trace_max_i, trace_max_j, seq1, seq2)
        line1 = trace_max_i - len(s1) + s1.count('-') + 1
        line2 = trace_max_j - len(s2) + s2.count('-') + 1
        print('%s vs %s:' % (f1, f2))
        print('-'*40)
        while s1:
            print('{0:15}{1:5d}  {2:60}'.format(f1 + ':', line1, s1[:60]))
            print('{0:20}  {1:60}'.format(' ',relation[:60]))
            print('{0:15}{1:5d}  {2:60}'.format(f2 + ':', line2, s2[:60]))
            print()

            line1 = line1 + 60 - s1[:60].count('-')
            line2 = line2 + 60 - s2[:60].count('-')

            s1 = s1[60:]
            s2 = s2[60:]
            relation = relation[60:]
        print('total score: ' + str(global_max))
        print()