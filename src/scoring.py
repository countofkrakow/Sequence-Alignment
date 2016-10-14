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

    # INGREDIENTS
    # 4 tablespoons (1/2 stick) butter
    # 4 tablespoons minced chervil, parsley or dill
    # 1 salmon fillet, 1 1/2 to 2 pounds
    # Salt and freshly ground black pepper to taste
    # Lemon wedges

    # PREPARATION
    # 1
    # Preheat the oven to 475 degrees. Place the butter and half the herb in a roasting pan just large
    # enough to fit the salmon and place it in the oven. Heat about 5 minutes, until the butter melts and the herb begins to sizzle.

    # 2
    # Add the salmon to the pan, skin side up. Roast 4 minutes. Remove from the oven, then peel the skin off. (If the skin does not
    # lift right off, cook 2 minutes longer.) Sprinkle with salt and pepper and turn the fillet over. Sprinkle with salt and pepper again.

    # 3
    # Roast 3 to 5 minutes more, depending on the thickness of the fillet and the degree of doneness you prefer. Cut into serving portions,
    # spoon a little of the butter over each and garnish with the remaining herb. Serve with lemon wedges.
