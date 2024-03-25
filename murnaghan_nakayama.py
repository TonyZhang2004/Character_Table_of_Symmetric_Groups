from typing import Dict, List
import pickle
import numpy as np
from sympy.utilities.iterables import partitions


MEMO = {}


def make_bit_strings(parts: List[Dict[int, int]]) -> List[str]:
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
    """
    Returns validated abaci bit string which must start with
    a 0 and end with a 1.
    """
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


def get_tau(sigma: str) -> str:
    """
    Returns the abaci bit string of a partition after removing
    the top-most row, i.e., the largest cycle type.
    """
    # j := idx of 1
    tau = sigma[:-1]  # just remove the 1
    tau = validate_bit_string(tau)
    return tau


def get_erased_bs(lambda_: str, i: int, j: int) -> str:
    """
    Returns the abaci bit string of a partition after removing
    the border strip corresponding to a hook starting at index
    i and ending at index j.
    """
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


def murnaghan_nakayama(lambda_: str, sigma: str) -> int:
    """
    Returns the character value corresponding to the partitions
    lambda_ (row) and sigma (column) using the recursive
    Murnaghan-Nakayama Rule.
    """
    global MEMO

    # Base Cases
    if len(lambda_) == 0 and len(sigma) == 0:  # X^0(0) = 1
        return 1
    if len(lambda_) == 1 and len(sigma) == 1:  # X^1(1) = 1
        return 1

    # DP
    if (lambda_, sigma) in MEMO:
        return MEMO[(lambda_, sigma)]

    hooks = [
        (x, y)
        for x in range(len(lambda_))
        for y in range(len(lambda_))
        if (lambda_[x] == "0" and lambda_[y] == "1")
    ]

    # Invalid diagram -- return 0
    if len(hooks) == 0:
        MEMO[(lambda_, sigma)] = 0
        return MEMO[(lambda_, sigma)]

    tau = get_tau(sigma)
    cycle_length_removed = sigma.count("0")

    mn_sum = 0
    # Loop through our hooks
    for hook in hooks:
        i, j = hook
        t = j - i  # formula for hook length
        if t == cycle_length_removed:
            height = lambda_[i:j].count("1")
            erased_bs = get_erased_bs(lambda_, i, j)
            temp = ((-1) ** height) * murnaghan_nakayama(erased_bs, tau)
            mn_sum += temp

    MEMO[(lambda_, sigma)] = mn_sum
    return MEMO[(lambda_, sigma)]


def get_character_table(
    n: int, output_file_name: str = "", memo_file_name: str = "memo.txt"
):
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

    if memo_file_name:
        write_memo_to_file(memo_file_name)
    char_table = np.fliplr(np.array(char_table, dtype=np.int64))

    if output_file_name:
        np.savetxt(output_file_name, char_table, delimiter=",", fmt="%d")

    return char_table


def get_character_value_of_column(
    n: int, col_bit_string="", csv_file_name: str = "", memo_file_name: str = "memo.txt"
):
    """
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

    if memo_file_name:
        write_memo_to_file(memo_file_name)
    char_table = np.fliplr(np.array(char_table, dtype=np.int64))

    if csv_file_name:
        np.savetxt(csv_file_name, char_table, delimiter=",", fmt="%d")

    return char_table


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


if __name__ == "__main__":
    N = 78
    i = 12
    get_character_value_of_column(
        N, "01" * i, f"S{N}_staircase.csv", "memo_staircase.txt"
    )
