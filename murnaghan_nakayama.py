from typing import Dict, List
import pickle
import numpy as np
from sympy.utilities.iterables import partitions


MEMO = {}


def read_memo_from_file(file_name: str = "memo.txt"):
    """
    Read existing memo file.
    """
    try:
        with open(file_name, "rb") as memo_file:
            global MEMO
            MEMO = pickle.load(memo_file)
    except IOError:
        return


def write_memo_to_file(file_name: str = "memo.txt"):
    """
    Save memo to a file.
    """
    with open(file_name, "wb") as memo_file:
        global MEMO
        pickle.dump(MEMO, memo_file)


def make_bit_strings(parts: List[Dict[int, int]]) -> List[Tuple[int, int]]:
    """
    Returns a list of partitions in abaci format bit strings.
    Parameters:
    parts: List of partitions where each partition is given by
    a dictionary whose key is the number appearing in the
    partition and the value is the number of times it appears.
    """
    bit_strings = []
    for part in parts:
        # Iterate through each partition
        curr_bit_str = 0
        no_of_bits = 0
        prev = 0
        lambda_prime = reversed(list(part.keys()))
        for key in list(lambda_prime):
            for _ in range(key - prev):
                curr_bit_str = curr_bit_str << 1  # Appends a 0
                no_of_bits += 1
            for _ in range(part[key]):
                curr_bit_str = (curr_bit_str << 1) | 1  # Appends a 1
                no_of_bits += 1
            prev = key
        bit_strings.append((curr_bit_str, no_of_bits))
    return bit_strings


def reverse_bits(bit_string: Tuple[int, int]) -> int:
    """
    Returns the number corresponding to the reverse bits 
    in the bit string
    """
    perm, no_of_bits = bit_string
    result = 0
    looped = 0
    while perm > 0:
        result = (result << 1) + (perm & 1)  # Append digit to end
        perm >>= 1  # Remove that digit
        looped += 1
    result <<= no_of_bits - looped
    return result


def validate_bit_string(bit_string: Tuple[int, int]) -> Tuple[int, int]:
    """
    Returns validated abaci bit string which must start with
    a 0 and end with a 1.
    """
    perm, no_of_bits = bit_string
    # Remove 0's from end of the invalid permutation
    while perm > 0:
        if perm % 2 == 0:
            perm >>= 1
            no_of_bits -= 1
        else:
            break
    # Remove 1's from end of the invalid permutation
    perm = reverse_bits((perm, no_of_bits))
    while perm > 0:
        if perm % 2 == 1:
            perm >>= 1
            no_of_bits -= 1
        else:
            break
    perm = reverse_bits((perm, no_of_bits))
    return (perm, no_of_bits)


def get_hooks(bit_string: Tuple[int, int]) -> List[Tuple[int, int]]:
    """
    Returns all the hooks present in the given bit string's partition.
    """
    hooks = []
    bit_string_rv = reverse_bits(bit_string)
    idx_0 = 0
    while bit_string_rv:
        if bit_string_rv % 2 == 0:
            idx_1 = idx_0  # idx_1 is relative to idx_0
            bit_string_rv2 = bit_string_rv
            # Iterate through the bit_string again to find the 1's
            while bit_string_rv2:
                if bit_string_rv2 % 2 == 1:
                    hooks.append((idx_0, idx_1))
                idx_1 += 1
                bit_string_rv2 >>= 1
        idx_0 += 1
        bit_string_rv >>= 1
    return hooks


def get_tau(sigma: Tuple[int, int]) -> Tuple[Tuple[int, int], int]:
    """
    Returns the abaci bit string of a partition after removing
    the top-most row, i.e., the largest cycle type, and the
    largest cycle type removed.
    """
    perm, no_of_bits = sigma
    tau = perm >> 1
    no_of_bits_r = no_of_bits - 1
    tau = validate_bit_string((tau, no_of_bits_r))
    # Calculate the cycle length removed:
    cycle_len_rem = 0
    looped = 0
    while perm:
        if perm % 2 == 0:
            cycle_len_rem += 1
        perm >>= 1
        looped += 1
    cycle_len_rem += no_of_bits - looped
    return (tau, cycle_len_rem)


def get_height(
    bit_string: Tuple[int, int], idx_start: int = 0, idx_end: int = -1
) -> int:
    """
    Returns the height of the hook.
    """
    perm_rv = reverse_bits(bit_string)
    perm_rv >>= idx_start  # Remove the first bits till idx_start
    idx_end -= idx_start
    height = 0
    # No. of 1's till idx_end is the height
    while perm_rv and idx_end != 0:  # We do not count the 1 at index end
        if perm_rv % 2 == 1:
            height += 1
        perm_rv >>= 1
        idx_end -= 1
    return height


def get_erased_bs(lambda_: Tuple[int, int], idx_0: int, idx_1: int) -> Tuple[int, int]:
    """
    Returns the abaci bit string of a partition after removing
    the border strip corresponding to a hook starting at index
    i and ending at index j.
    """
    perm, no_of_bits = lambda_
    # Switch 0 and 1 and fix 0-indexing
    erased_bs = (perm | (1 << (no_of_bits - idx_0 - 1))) & ~(1 << (no_of_bits - idx_1 - 1))
    erased_bs = validate_bit_string((erased_bs, no_of_bits))
    return erased_bs


def murnaghan_nakayama(n: int, lambda_: Tuple[int, int], sigma: Tuple[int, int]) -> int:
    """
    Returns the character value corresponding to the partitions
    lambda_ (row) and sigma (column) using the recursive
    Murnaghan-Nakayama Rule.
    """
    global MEMO

    # Base Cases
    if lambda_[0] == 0 and sigma[0] == 0:  # X_0(0) = 1
        return 1
    if lambda_[0] == 1 and sigma[0] == 1:  # X_1(1) = 1
        return 1

    # DP
    if (lambda_, sigma) in MEMO:
        return MEMO[(lambda_, sigma)]

    # Get our hooks
    hooks = get_hooks(lambda_)

    # Invalid diagram -- return 0
    if len(hooks) == 0:
        MEMO[(n, lambda_, sigma)] = 0
        return MEMO[(n, lambda_, sigma)]

    tau, max_cycle_length = get_tau(sigma)

    mn_sum = 0
    # Loop through our hooks
    for hook in hooks:
        i, j = hook
        t = j - i  # formula for hook length
        if t == max_cycle_length:
            height = get_height(lambda_, i, j)
            erased_bs = get_erased_bs(lambda_, i, j)
            temp = ((-1) ** height) * murnaghan_nakayama(erased_bs, tau)
            mn_sum += temp

    MEMO[(lambda_, sigma)] = mn_sum
    return MEMO[(lambda_, sigma)]


def get_character_table(n: int, csv_file_name: str = "", memo_file_name: str = ""):
    """
    Returns the whole character table for n, i.e., the character
    table of Sn
    Parameters:
    n: The number for which character table is copmuted
    output_file_name: file to store character table
    memo_file_name: read/write from existing memo file
    """
    parts = list(partitions(n))  # Get list of partitions
    bit_strings = make_bit_strings(parts)

    if memo_file_name:
        read_memo_from_file(memo_file_name)

    char_table = []
    for lambda_ in bit_strings:
        curr_row = []
        for sigma in bit_strings:
            val = murnaghan_nakayama(lambda_, sigma)
            curr_row.append(val)
        char_table.append(curr_row)
    char_table = np.fliplr(np.array(char_table, dtype=np.int64))

    if memo_file_name:
        write_memo_to_file(memo_file_name)

    if output_file_name:
        np.savetxt(output_file_name, char_table, delimiter=",", fmt="%d")

    return char_table


def get_character_value_of_column(
    n: int, csv_file_name: str = "", col_bit_string="", memo_file_name: str = "memo.txt"
):
    """
    TODO: 
    Returns the character values for just for a column in Sn
    Parameters:
    n: The number for which character table is copmuted
    col_bit_string: the abaci bit string value of the column
    for which the characater values will be computed
    output_file_name: file to store character table
    memo_file_name: read/write from existing memo file
    """
    parts = list(partitions(n))  # Get list of partitions
    bit_strings = make_bit_strings(parts)

    if memo_file_name:
        read_memo_from_file(memo_file_name)

    char_table = []
    for lambda_ in bit_strings:
        curr_row = []
        sigma = col_bit_string
        val = murnaghan_nakayama(lambda_, sigma)
        curr_row.append(val)
        char_table.append(curr_row)
    char_table = np.fliplr(np.array(char_table, dtype=np.int64))

    if memo_file_name:
        write_memo_to_file(memo_file_name)

    if csv_file_name:
        np.savetxt(csv_file_name, char_table, delimiter=",", fmt="%d")

    return char_table


if __name__ == "__main__":
    print(get_character_table(8, "S8_new.csv"))
    # print(reverse_bits(11, 11))
