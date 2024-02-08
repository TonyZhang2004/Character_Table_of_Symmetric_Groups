# import pandas as pd
# import matplotlib.pyplot as plt
import numpy as np
from sympy.utilities.iterables import partitions


def make_bit_strings(parts: list) -> list[str]:
    bit_strings = []
    for lambda_ in parts:
        # Iterate through each partition
        # Make our bit string
        curr_bit_str = ""
        prev = 0
        lamda_prime = reversed(list(lambda_.keys()))
        for key in list(lamda_prime):
            for _ in range(key - prev):
                curr_bit_str += '0'
            for _ in range(lambda_[key]):
                curr_bit_str += '1'
            prev = key

        bit_strings.append(curr_bit_str)
    return bit_strings


def validate_bit_string(bit_string: str) -> str:
    # Remove 1's from start of bit_string
    while len(bit_string) > 0:
        if bit_string[0] == '1':
            bit_string = bit_string[1:]
        else:
            break
    # Remove 0's from end of bit_string
    while len(bit_string) > 0:
        if bit_string[-1] == '0':
            bit_string = bit_string[:-1]
        else:
            break
    return bit_string


def get_cycle_lengths(perm: list) -> dict:
    # key is cycle length and value is starting 0 and ending 1
    cycle_lengths = {}
    curr = 0
    j = 0
    for i, val in enumerate(perm):
        if val == '0':
            curr += 1
        if val == '1':
            if curr not in cycle_lengths:
                cycle_lengths[curr] = (j, i)
            j = i + 1
    return cycle_lengths


def get_erased_bs(lambda_: str, i: int, j: int) -> str:
    erased_bs = ''
    for idx, c in enumerate(lambda_):
        if idx == i:
            erased_bs += lambda_[j]
        elif idx == j:
            erased_bs += lambda_[i]
        else:
            erased_bs += c
    erased_bs = validate_bit_string(erased_bs)
    return erased_bs


def murnaghan_nakayama(n: int, lambda_: str, sigma: str):
    # TODO: I think n is not a required argument, we only care about the partitions
    # lambda_ = string_to_list(lambda_p)
    # sigma = string_to_list(sigma_p)

    # Base Cases
    if len(lambda_) == 0 and len(sigma) == 0: # X_0(0) = 1
        return 1
    # if lambda_ == [0,1] and len(sigma) == [0,1]:
    #     return 1
    # if n < 1:
    #     return 0

    # DP
    if (n, lambda_, sigma) in memo:
        return memo[(n, lambda_, sigma)]

    # Get our hooks
    hooks = [
        (x, y)
        for x in range(len(lambda_))
        for y in range(len(lambda_))
        if (lambda_[x] == '0' and lambda_[y] == '1')
    ]

    # Invalid diagram -- return 0
    if len(hooks) == 0:
        memo[(n, lambda_, sigma)] = 0
        return memo[(n, lambda_, sigma)]

    cycle_lengths = get_cycle_lengths(sigma)

    mn_sum = 0
    # Loop through our hooks
    for hook in hooks:
        i, j = hook
        t = j - i
        if t in cycle_lengths:
            # There exists rho (t-cycle)
            # Now need to find rho and tau
            # Just now need to remove the actual row
            tau = sigma[: cycle_lengths[t][0] + 1] + sigma[cycle_lengths[t][1] + 1 :]
            tau = validate_bit_string(tau)

            height = lambda_[i:j].count('1')

            erased_bs = get_erased_bs(lambda_, i, j)

            temp = ((-1) ** height) * murnaghan_nakayama(
                n - t, erased_bs, tau
            )
            mn_sum += temp

    # Done!
    memo[(n, lambda_, sigma)] = mn_sum
    return memo[(n, lambda_, sigma)]


def get_character_table(N: int):
    parts = list(partitions(N))  # Get list of partitions
    bit_strings = make_bit_strings(parts)
    char_table = []

    for lambda_ in bit_strings:
        curr_row = []
        for sigma in bit_strings:
            val = murnaghan_nakayama(N, lambda_, sigma)
            curr_row.append(val)
        char_table.append(curr_row)

    return np.fliplr(np.matrix(char_table))


memo = {}

if __name__ == "__main__":
    N = 3
    char_table = get_character_table(N)
    print(char_table)
