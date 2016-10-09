import scoring
import argparse
import sys


def smith_waterman(seq1, seq2, gap_penalty=-4):
    sf = scoring.ScoreFunction()
    n = len(seq1)
    m = len(seq2)
    matrix = [[0] * (m + 1) for i in range(n + 1)]
    global_max = 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            matrix[i][j] = max(
                    matrix[i-1][j-1] + sf.getScore(seq1[i-1], seq2[j-1]),
                    matrix[i-1][j] + gap_penalty,
                    matrix[i][j-1] + gap_penalty,
                    0)
            global_max = matrix[i][j] if matrix[i][j] > global_max else global_max

    print(matrix)
    print(global_max)

def needleman_wunsch():
    pass

def parse_args():
    # argument parsing
    parser = argparse.ArgumentParser(description="performs sequence alignment on two genetic sequences.")
    parser.add_argument('-f', '--file', action='store_true')
    algo = parser.add_mutually_exclusive_group()
    algo.add_argument('-sw', '--smith_waterman', action='store_true')
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

    # raw sequences as arguments
    else:
        pass

    smith_waterman(args.sequence1, args.sequence2)
    if args.smith_waterman:
        print("b")
    else:
        print('t')