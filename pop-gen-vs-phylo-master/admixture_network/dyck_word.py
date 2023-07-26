# http://www.cs.otago.ac.nz/staffpriv/mike/Papers/RandomGeneration/RandomBinaryTrees.pdf
# https://www.math.utah.edu/mathcircle/tessler-catalan-notes.pdf (slide 14)

import random

def _phi(w):
    """
    Implements the function phi as described in the algorithm "Random Bracket Sequence."
    """
    if len(w) == 0: # If n = 0, return the empty string
        return ""
    else:
        # Initialize variables that persist across loop iterations
        u = ""
        v = ""
        zigzag_height = 0
        u_defect = 0
        for i, char in enumerate(w):
            # Keep track of zigzag line height in order to determine if substring up to i is balanced
            if char == "1":
                zigzag_height += 1
            elif char == "2":
                zigzag_height -= 1

            # Partial summation procedure
            if i % 2 == 0 and zigzag_height < 0: # Checks for *even* i to compensate for zero-based indexing
                u_defect += 1
        
            if zigzag_height == 0: # First balanced substring will be irreductible
                u = w[:i+1] # start of string through character at index i
                v = w[i+1:] # character at index i+1 through end of string
                break

        if u_defect == 0: # Check if u is well-formed
            return u + _phi(v) # Recursive call
        else:
            return "1" + _phi(v) + "2" + "".join(["1" if c == "2" else "2" for c in u[1:-1]]) # Recursive call, final part strips off first and last characters of u and then inverts the remainder

def generate_dyck_word(n):
    """
    Randomly generates a Dyck word that is a member of W_n, following a uniform distribution.
    """
    L = set(random.sample(range(1, 2*n + 1), n)) # Step 1 of algorithm "Random Bracket Sequence"
    x = "".join(["1" if i in L else "2" for i in range(1, 2*n + 1)]) # Step 2 of algorithm "Random Bracket Sequence"
    return _phi(x) # Step 3 of algorithm "Random Bracket Sequence"

if __name__ == "__main__":
    print(generate_dyck_word(5))