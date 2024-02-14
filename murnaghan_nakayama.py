# import matplotlib.pyplot as plt
from typing import Dict, List, Tuple
import pickle
import pandas as pd
import numpy as np
from sympy.utilities.iterables import partitions


def make_bit_strings(parts: List[Dict[int, int]]) -> List[str]:
    bit_strings = []
    for part in parts:
        # Iterate through each partition
        curr_bit_str = ""
        prev = 0
        lambda_prime = reversed(list(part.keys()))
        for key in list(lambda_prime):
            for _ in range(key - prev):
                curr_bit_str += "0"
            for _ in range(part[key]):
                curr_bit_str += "1"
            prev = key
        bit_strings.append(curr_bit_str)
    return bit_strings


def validate_bit_string(bit_string: str) -> str:
    # Remove 1's from start of bit_string
    while len(bit_string) > 0:
        if bit_string[0] == "1":
            bit_string = bit_string[1:]
        else:
            break
    # Remove 0's from end of bit_string
    while len(bit_string) > 0:
        if bit_string[-1] == "0":
            bit_string = bit_string[:-1]
        else:
            break
    return bit_string


def get_cycle_lengths(perm: str) -> Dict[int, Tuple[int, int]]:
    # key is cycle length and value is starting 0 and ending 1
    cycle_lengths = {}
    curr = 0
    j = 0
    for i, val in enumerate(perm):
        if val == "0":
            curr += 1
        if val == "1":
            if curr not in cycle_lengths:
                cycle_lengths[curr] = (j, i)
            j = i + 1
    return cycle_lengths


def get_tau(sigma: str, j: int) -> str:
    # j := idx of 1
    tau = sigma[:j] + sigma[(j + 1) :]  # just remove the 1
    tau = validate_bit_string(tau)
    return tau


def get_erased_bs(lambda_: str, i: int, j: int) -> str:
    erased_bs = ""
    for idx, c in enumerate(lambda_):
        if idx == i:
            erased_bs += lambda_[j]
        elif idx == j:
            erased_bs += lambda_[i]
        else:
            erased_bs += c
    erased_bs = validate_bit_string(erased_bs)
    return erased_bs


def murnaghan_nakayama(n: int, lambda_: str, sigma: str) -> int:
    # TODO: I think n is not a required argument, we only care about the partitions
    global memo

    # Base Cases
    if len(lambda_) == 0 and len(sigma) == 0:  # X_0(0) = 1
        return 1

    # DP
    if (n, lambda_, sigma) in memo:
        return memo[(n, lambda_, sigma)]

    # Get our hooks
    hooks = [
        (x, y)
        for x in range(len(lambda_))
        for y in range(len(lambda_))
        if (lambda_[x] == "0" and lambda_[y] == "1")
    ]

    # Invalid diagram -- return 0
    if len(hooks) == 0:
        memo[(n, lambda_, sigma)] = 0
        return memo[(n, lambda_, sigma)]

    cycle_lengths = get_cycle_lengths(sigma)

    max_cycle_length = max(list(cycle_lengths.keys()))
    tau = get_tau(sigma, cycle_lengths[max_cycle_length][1])

    mn_sum = 0
    # Loop through our hooks
    for hook in hooks:
        i, j = hook
        t = j - i  # formula for hook length
        if t == max_cycle_length:
            height = lambda_[i:j].count("1")
            erased_bs = get_erased_bs(lambda_, i, j)
            temp = ((-1) ** height) * murnaghan_nakayama(n - t, erased_bs, tau)
            mn_sum += temp

    memo[(n, lambda_, sigma)] = mn_sum
    return memo[(n, lambda_, sigma)]


def get_character_table(
    n: int, memo_file_name: str = "memo.txt", csv_file_name: str = ""
):
    parts = list(partitions(n))  # Get list of partitions
    bit_strings = make_bit_strings(parts)

    read_memo_from_file(memo_file_name)
    char_table = []

    for lambda_ in bit_strings:
        curr_row = []
        for sigma in bit_strings:
            val = murnaghan_nakayama(n, lambda_, sigma)
            curr_row.append(val)
        char_table.append(curr_row)

    write_memo_to_file(memo_file_name)
    char_table = np.fliplr(np.array(char_table, dtype=np.int64))

    if csv_file_name:
        np.savetxt(csv_file_name, char_table, delimiter=",", fmt="%d")

    return char_table


def read_memo_from_file(file_name: str = "memo.txt"):
    global memo
    try:
        with open(file_name, "rb") as memo_file:
            memo = pickle.load(memo_file)
    except IOError:
        return


def write_memo_to_file(file_name: str = "memo.txt"):
    with open(file_name, "wb") as memo_file:
        pickle.dump(memo, memo_file)


memo = {}

if __name__ == "__main__":
    N = 21
    char_table = get_character_table(N, csv_file_name="S21.csv")
    print(char_table)
