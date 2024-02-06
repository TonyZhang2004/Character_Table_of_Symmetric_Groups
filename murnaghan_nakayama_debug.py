import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sympy.utilities.iterables import partitions


def list_to_string(arr: list):
    return "".join([str(elem) for elem in arr])


def string_to_list(s: str):
    return [int(char) for char in s]

memo = {}

def murnaghan_nakayama(n: int, lambda_p, sigma_p):
    print("very start")
    lambda_ = string_to_list(lambda_p)
    sigma = string_to_list(sigma_p)
    print(f"lambda: {lambda_}")
    print(f"sigma: {sigma}")

    # Base Cases
    if (n == 1 and len(lambda_) == 1 and len(sigma) == 1) or (n == 0 and len(lambda_) == 0 and len(sigma) == 0):
        print(f"returning: {1}")
        return 1

    if n < 1:
        print(f"returning: {0}")
        return 0

    # DP
    # if (n, lambda_p, sigma_p) in memo:
    # print(f"returning: {memo[(n, lambda_p, sigma_p)]}")
    # return memo[(n, lambda_p, sigma_p)]
    # Get our hooks
    hooks = [
        (x, y)
        for x in range(len(lambda_))
        for y in range(len(lambda_))
        if (lambda_[x] == 0 and lambda_[y] == 1)
    ]

    # Invalid diagram -- return 0
    if len(hooks) == 0:
        memo[(n, lambda_p, sigma_p)] = 0
        # print(f'returning: {0}')
        # return memo[(n, lambda_p, sigma_p)]

    cycle_lengths = (
        {}
    )  # dictionary of unique cycle lengths where key is cycle length and value is starting 0 and ending 1
    curr = 0
    j = 0
    for i in range(len(sigma)):
        if sigma[i] == 0:
            curr += 1
        if sigma[i] == 1:
            if curr not in cycle_lengths:
                cycle_lengths[curr] = (j, i)
            j = i + 1

    print(f"cycle lengths: {cycle_lengths}")
    sum = 0
    # Loop through our hooks
    print(f"hooks: {hooks}")
    for i, j in hooks:
        t = j - i
        print(f"t: {t}")
        if t in cycle_lengths:  # there exists rho
            # Now need to find rho and tau
            # Just now need to remove the actual row
            tau = sigma[: cycle_lengths[t][0] + 1] + sigma[cycle_lengths[t][1] + 1 :]
            # Get rid of trailing 0's
            print(f"tau pre-removal: {tau}")
            while len(tau) > 0:
                if tau[-1] != 1:
                    del tau[-1]
                else:
                    break

            print(f"tau : {tau}")
            height = lambda_[i:j].count(1)
            erased_bs = lambda_.copy()
            erased_bs[i] = lambda_[j]
            erased_bs[j] = lambda_[i]

            # Remove 1's from start of list
            print(f"erased bs : {erased_bs}")
            while len(erased_bs) > 0:
                if erased_bs[0] != 0:
                    del erased_bs[0]
                else:
                    break

            # Remove 0's from end of list
            print(f"erased bs: {erased_bs}")
            while len(erased_bs) > 0:
                if erased_bs[-1] != 1:
                    del erased_bs[-1]
                else:
                    break

            print(f"erased bs : {erased_bs}")
            print(f"height: {height}")
            sum += ((-1) ** height) * murnaghan_nakayama(
                n - t, list_to_string(erased_bs), list_to_string(tau)
            )

    # Done!
    memo[(n, lambda_p, sigma_p)] = sum
    print(f"returning: {sum}")
    return memo[(n, lambda_p, sigma_p)]


def get_character_table(N: int):
    parts = list(partitions(N))
    bit_strings = []
    char_table = []
    for lambda_ in parts:
        # Iterate through each partition

        # Make our bit string
        arr = []
        prev = 0
        lamda_prime = reversed(list(lambda_.keys()))
        for key in list(lamda_prime):
            for _ in range(key - prev):
                arr.append(0)
            for _ in range(lambda_[key]):
                arr.append(1)
            prev = key

        bit_strings.append(list_to_string(arr))
    i = 0
    for lambda_ in bit_strings:
        char_table.append([])
        for sigma in bit_strings:
            to_add = murnaghan_nakayama(N, lambda_, sigma)
            print(f"adding new char: {to_add}")
            char_table[i].append(to_add)
        i += 1

    return np.fliplr(np.matrix(char_table))

if __name__ == "__main__":
# Get list of partitions
    N = 3
    char_table = get_character_table(N)
    print(char_table)
