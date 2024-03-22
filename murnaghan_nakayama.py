from typing import Dict, List, Tuple
import pickle
import numpy as np
from sympy.utilities.iterables import partitions


MEMO = {}


def make_bit_strings(parts: List[Dict[int, int]]) -> List[int]:
    bit_strings = []
    for part in parts:
        # Iterate through each partition
        curr_bit_str = 0
        prev = 0
        lambda_prime = reversed(list(part.keys()))
        for key in list(lambda_prime):
            for _ in range(key - prev):
                curr_bit_str = curr_bit_str << 1
            for _ in range(part[key]):
                curr_bit_str = (curr_bit_str << 1) | 1
            prev = key
        bit_strings.append(curr_bit_str)
    return bit_strings


def reverse_bits(num: int, n: int) -> int:
    result = 0
    looped = 0
    while num:
        result = (result << 1) + (num & 1)
        num >>= 1
        looped += 1
    result <<= (n - looped + 1)
    return result


def validate_bit_string(bit_string: int, n: int) -> int:
    # Remove 0's from end of bit_string
    while bit_string > 0:
        if bit_string % 2 == 0:
            bit_string = bit_string >> 1
        else:
            break
    # Remove 1's from start of bit_string
    bit_string = reverse_bits(bit_string, n)
    while bit_string > 0:
        if bit_string % 2 == 1:
            bit_string = bit_string >> 1
        else:
            break
    bit_string = reverse_bits(bit_string, n)
    return bit_string


def get_hooks(bit_string: int, n: int) -> List[Tuple[int, int]]:
    hooks = []
    bit_string_rv = reverse_bits(bit_string, n)
    idx_0 = 0
    while bit_string_rv:
        if bit_string_rv % 2 == 0:
            idx_1 = idx_0    # idx_1 is relative to idx_0
            bit_string_rv2 = bit_string_rv
            while bit_string_rv2:
                if bit_string_rv2 % 2 == 1:
                    hooks.append((idx_0, idx_1))
                idx_1 += 1
                bit_string_rv2 = bit_string_rv2 >> 1
        idx_0 += 1
        bit_string_rv = bit_string_rv >> 1
    return hooks


def get_cycle_lengths(perm: int, n: int) -> Dict[int, Tuple[int, int]]:
    # key is cycle length and value is starting 0 and ending 1
    cycle_lengths = {}
    perm_rv = reverse_bits(perm, n)
    curr = idx_0 = idx_1 = 0
    while perm_rv > 0:
        if perm_rv % 2 == 0:
            curr += 1
        else:
            if curr not in cycle_lengths:
                cycle_lengths[curr] = (idx_0, idx_1)
            idx_0 = idx_1 + 1
        idx_1 += 1
        perm_rv = perm_rv >> 1
    return cycle_lengths


def get_tau(sigma: int, n: int) -> Tuple[int, int]:
    # # j := idx of 1
    # sigma_rv = reverse_bits(sigma, n)
    # # Set all first j bits to 0, right shift it by 1, then add all the j-1 bits back
    # tau = ((sigma_rv & (0 << j)) >> 1) | (sigma_rv & ((1 << j) - 1))
    # tau = validate_bit_string(reverse_bits(tau, n), n)
    tau = sigma >> 1
    # tau_temp = tau
    # while tau_temp:
    #     if tau_temp
    #     tau_temp >> 1
    tau = validate_bit_string(tau, n - cycle_length)
    return (tau, cycle_length)


def get_height(bit_string: int, n: int, idx_start: int = 0, idx_end: int = -1) -> int:
    bit_string_rv = reverse_bits(bit_string, n)
    bit_string_rv = bit_string_rv >> idx_start    # Remove the first i bits
    count = 0
    while bit_string_rv and idx_end:
        if bit_string_rv % 2 == 1:
            count += 1
        bit_string_rv = bit_string_rv >> 1
        idx_end -= 1
    return count


def get_erased_bs(lambda_: int, n: int, idx_0: int, idx_1: int) -> int:
    lambda_rv = reverse_bits(lambda_, n)
    erased_bs = lambda_rv | (1 << idx_0)
    erased_bs = lambda_rv & ~(1 << idx_1)
    erased_bs = validate_bit_string(reverse_bits(erased_bs, n), n)
    return erased_bs


def murnaghan_nakayama(n: int, lambda_: int, sigma: int) -> int:
    # TODO: I think n is not a required argument, we only care about the partitions
    global MEMO

    # Base Cases
    if lambda_ == 0 and sigma == 0:  # X_0(0) = 1
        return 1
    if lambda_ == 1 and sigma == 1:  # X_1(1) = 1
        return 1

    # DP
    if (n, lambda_, sigma) in MEMO:
        return MEMO[(n, lambda_, sigma)]

    # Get our hooks
    hooks = get_hooks(lambda_, n)

    # Invalid diagram -- return 0
    if len(hooks) == 0:
        MEMO[(n, lambda_, sigma)] = 0
        return MEMO[(n, lambda_, sigma)]

    # cycle_lengths = get_cycle_lengths(sigma, n)

    # max_cycle_length = max(list(cycle_lengths.keys()))
    tau, cycle_length = get_tau(sigma, n)

    mn_sum = 0
    # Loop through our hooks
    for hook in hooks:
        i, j = hook
        t = j - i  # formula for hook length
        if t == max_cycle_length:
            height = get_height(lambda_, n, i, j)
            erased_bs = get_erased_bs(lambda_, n, i, j)
            temp = ((-1) ** height) * murnaghan_nakayama(n - t, erased_bs, tau)
            mn_sum += temp

    MEMO[(n, lambda_, sigma)] = mn_sum
    return MEMO[(n, lambda_, sigma)]


def get_character_table(
    n: int, csv_file_name: str = "", memo_file_name: str = "memo.txt"
):
    parts = list(partitions(n))  # Get list of partitions
    bit_strings = make_bit_strings(parts)

    if memo_file_name:
        read_memo_from_file(memo_file_name)
    
    char_table = []
    for lambda_ in bit_strings:
        curr_row = []
        for sigma in bit_strings:
            val = murnaghan_nakayama(n, lambda_, sigma)
            curr_row.append(val)
        char_table.append(curr_row)

    if memo_file_name:
        write_memo_to_file(memo_file_name)
    char_table = np.fliplr(np.array(char_table, dtype=np.int64))

    if csv_file_name:
        np.savetxt(csv_file_name, char_table, delimiter=",", fmt="%d")

    return char_table


def get_character_value_of_column(
    n: int,
    csv_file_name: str = "",
    col_bit_string="",
    memo_file_name: str = "memo.txt"
):
    parts = list(partitions(n))  # Get list of partitions
    bit_strings = make_bit_strings(parts)

    if memo_file_name:
        read_memo_from_file(memo_file_name)
    char_table = []

    for lambda_ in bit_strings:
        curr_row = []
        sigma = col_bit_string
        val = murnaghan_nakayama(n, lambda_, sigma)
        curr_row.append(val)
        char_table.append(curr_row)

    if memo_file_name:
        write_memo_to_file(memo_file_name)
    char_table = np.fliplr(np.array(char_table, dtype=np.int64))

    if csv_file_name:
        np.savetxt(csv_file_name, char_table, delimiter=",", fmt="%d")

    return char_table


def read_memo_from_file(file_name: str = "memo.txt"):
    try:
        with open(file_name, "rb") as memo_file:
            global MEMO
            MEMO = pickle.load(memo_file)
    except IOError:
        return


def write_memo_to_file(file_name: str = "memo.txt"):
    with open(file_name, "wb") as memo_file:
        global MEMO
        pickle.dump(MEMO, memo_file)


if __name__ == "__main__":
    my_char_table = get_character_table(3, "S3_new.txt", "memo_new.txt")
