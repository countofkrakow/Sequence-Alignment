base_to_index = {
    "A": 0,
    "R": 1,
    "N": 2,
    "D": 3,
    "C": 4,
    "Q": 5,
    "E": 6,
    "G": 7,
    "H": 8,
    "I": 9,
    "L": 10,
    "K": 11,
    "M": 12,
    "F": 13,
    "P": 14,
    "S": 15,
    "T": 16,
    "W": 17,
    "Y": 18,
    "V": 19
}

NUM_BASE_PAIRS = 20

class ScoreFunction:
    def __init__(self, scoring_matrix_file="BLOSUM62.txt"):
        self.score_matrix = []
        for line in open(scoring_matrix_file):
            row = [int(num) for num in line.split('\t')]
            self.score_matrix.append(row)

    def getScore(self, base1, base2):
        b1i = base_to_index[base1.upper()]
        b2i = base_to_index[base2.upper()]
        return self.score_matrix[b1i][b2i]


