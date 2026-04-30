from typing import Dict, Iterable, List, Tuple
import csv
import json
import logging
import math
import os
import pickle
import time
import numpy as np


MEMO = {}
RIM_HOOK_CACHE = {}
LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.NullHandler())


def configure_logging(log_file: str = "", level: int = logging.INFO) -> None:
    """
    Configure progress logging for long computations. Calling this more than
    once replaces handlers installed by this module.
    """
    LOGGER.setLevel(level)
    LOGGER.handlers = [
        handler
        for handler in LOGGER.handlers
        if not getattr(handler, "_mn_handler", False)
    ]

    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    stream_handler._mn_handler = True
    LOGGER.addHandler(stream_handler)

    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        file_handler._mn_handler = True
        LOGGER.addHandler(file_handler)


def partitions(n: int) -> Iterable[Dict[int, int]]:
    """
    Generate integer partitions of n as dictionaries compatible with
    sympy.utilities.iterables.partitions. Keys are inserted in decreasing
    order because make_bit_strings reverses the keys while constructing the
    abacus word.
    """

    def generate(remaining: int, max_part: int, prefix: List[int]):
        if remaining == 0:
            counts = {}
            for part in sorted(prefix, reverse=True):
                counts[part] = counts.get(part, 0) + 1
            yield counts
            return
        for part in range(min(remaining, max_part), 0, -1):
            yield from generate(remaining - part, part, prefix + [part])

    yield from generate(n, n, [])


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


def make_bit_string(part: Dict[int, int]) -> str:
    """
    Return the abacus bit-string encoding of one partition dictionary.
    """
    curr_bit_str = ""
    prev = 0
    lambda_prime = reversed(list(part.keys()))
    for key in list(lambda_prime):
        for _ in range(key - prev):
            curr_bit_str += "0"
        for _ in range(part[key]):
            curr_bit_str += "1"
        prev = key
    return curr_bit_str


def get_partition_bit_strings(n: int) -> List[str]:
    """
    Return the abaci bit-string encodings of all partitions of n, in the row
    order used by this module.
    """
    return make_bit_strings(list(partitions(n)))


def iter_partition_bit_strings(n: int) -> Iterable[str]:
    """
    Stream abacus bit-string encodings of all partitions of n in row order.
    """
    for partition in partitions(n):
        yield make_bit_string(partition)


def partition_number(n: int) -> int:
    """
    Return the number of integer partitions of n.
    """
    if n < 0:
        raise ValueError("n must be nonnegative")
    counts = [0] * (n + 1)
    counts[0] = 1
    for part in range(1, n + 1):
        for total in range(part, n + 1):
            counts[total] += counts[total - part]
    return counts[n]


def staircase_rank_from_n(n: int) -> int:
    """
    Return k when n is the triangular number k(k + 1) / 2.
    """
    if n <= 0:
        raise ValueError("n must be a positive triangular number")
    discriminant_root = math.isqrt(8 * n + 1)
    if discriminant_root * discriminant_root != 8 * n + 1:
        raise ValueError(f"n={n} is not triangular")
    if (discriminant_root - 1) % 2:
        raise ValueError(f"n={n} is not triangular")
    return (discriminant_root - 1) // 2


def staircase_bit_string(n: int) -> str:
    """
    Return the abacus bit string for the staircase partition of n.
    """
    return "01" * staircase_rank_from_n(n)


def staircase_cycle_lengths(n: int) -> Tuple[int, ...]:
    """
    Return the staircase cycle type (k, k - 1, ..., 1) for n = k(k + 1) / 2.
    """
    k = staircase_rank_from_n(n)
    return tuple(range(k, 0, -1))


def clear_caches() -> None:
    """
    Clear all in-memory dynamic-programming caches.
    """
    global MEMO, RIM_HOOK_CACHE
    MEMO = {}
    RIM_HOOK_CACHE = {}


def enforce_memo_limit(max_memo_entries: int) -> None:
    """
    Keep the recursive memo from growing without bound during large streamed
    computations. Clearing the memo preserves correctness but may trade memory
    for recomputation.
    """
    global MEMO
    if max_memo_entries and len(MEMO) > max_memo_entries:
        LOGGER.info("Clearing recursive memo after %d entries", len(MEMO))
        MEMO = {}


def default_progress_file(output_file_name: str) -> str:
    """
    Return the default progress sidecar path for an output file.
    """
    return f"{output_file_name}.progress.json"


def read_progress(progress_file: str) -> Dict:
    """
    Read a progress sidecar file. Missing or malformed files are treated as
    empty progress rather than trusted.
    """
    if not progress_file:
        return {}
    try:
        with open(progress_file, "r") as file:
            progress = json.load(file)
            return progress if isinstance(progress, dict) else {}
    except (IOError, json.JSONDecodeError):
        return {}


def write_progress(progress_file: str, progress: Dict) -> None:
    """
    Atomically write a progress sidecar file.
    """
    if not progress_file:
        return
    progress["updated_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
    temp_file = f"{progress_file}.tmp"
    with open(temp_file, "w") as file:
        json.dump(progress, file, indent=2, sort_keys=True)
        file.write("\n")
    os.replace(temp_file, progress_file)


def progress_matches(progress: Dict, n: int, table_size: int, output_kind: str) -> bool:
    """
    Check whether a progress file describes the table we are about to write.
    """
    return (
        progress.get("n") == n
        and progress.get("table_size") == table_size
        and progress.get("output_kind") == output_kind
    )


def make_progress(
    n: int, table_size: int, output_kind: str, output_file_name: str, completed_rows: int
) -> Dict:
    """
    Build progress metadata for a table computation.
    """
    return {
        "n": n,
        "table_size": table_size,
        "output_kind": output_kind,
        "output_file_name": output_file_name,
        "completed_rows": completed_rows,
    }


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


def cycle_lengths_from_bit_string(sigma: str) -> Tuple[int, ...]:
    """
    Convert an abaci bit string for a conjugacy class into the sequence of
    cycle lengths removed by the Murnaghan-Nakayama recursion.
    """
    cycle_lengths = []
    sigma = validate_bit_string(sigma)
    while sigma:
        cycle_lengths.append(sigma.count("0"))
        sigma = get_tau(sigma)
    return tuple(cycle_lengths)


def get_rim_hook_transitions(lambda_: str, cycle_length: int) -> List[Tuple[str, int]]:
    """
    Return all ways to remove a rim hook of the requested length from lambda_.
    Each transition is (new_lambda, sign). The sign convention matches the
    existing abacus implementation: sign = (-1) ** number_of_ones_between.
    """
    cache_key = (lambda_, cycle_length)
    if cache_key in RIM_HOOK_CACHE:
        return RIM_HOOK_CACHE[cache_key]

    transitions = []
    for i, char in enumerate(lambda_):
        if char != "0":
            continue
        j = i + cycle_length
        if j < len(lambda_) and lambda_[j] == "1":
            height = lambda_[i:j].count("1")
            sign = -1 if height % 2 else 1
            transitions.append((get_erased_bs(lambda_, i, j), sign))

    RIM_HOOK_CACHE[cache_key] = transitions
    return transitions


def murnaghan_nakayama_cycles(lambda_: str, cycle_lengths: Tuple[int, ...]) -> int:
    """
    Return the character value for row lambda_ and a conjugacy class described
    by cycle lengths, using cached rim-hook transitions.
    """
    global MEMO

    if len(lambda_) == 0 and len(cycle_lengths) == 0:
        return 1
    if len(lambda_) == 0 or len(cycle_lengths) == 0:
        return 0

    key = (lambda_, cycle_lengths)
    if key in MEMO:
        return MEMO[key]

    cycle_length = cycle_lengths[0]
    tau = cycle_lengths[1:]
    mn_sum = 0
    for erased_bs, sign in get_rim_hook_transitions(lambda_, cycle_length):
        mn_sum += sign * murnaghan_nakayama_cycles(erased_bs, tau)

    MEMO[key] = mn_sum
    return mn_sum


def murnaghan_nakayama(lambda_: str, sigma: str) -> int:
    """
    Returns the character value corresponding to the partitions
    lambda_ (row) and sigma (column) using the recursive
    Murnaghan-Nakayama Rule.
    """
    return murnaghan_nakayama_cycles(lambda_, cycle_lengths_from_bit_string(sigma))


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
    bit_strings = get_partition_bit_strings(n)

    if memo_file_name:
        read_memo_from_file(memo_file_name)
    char_table = []
    sigma_cycle_lengths = [cycle_lengths_from_bit_string(sigma) for sigma in bit_strings]

    for lambda_ in bit_strings:
        curr_row = []
        for cycle_lengths in sigma_cycle_lengths:
            val = murnaghan_nakayama_cycles(lambda_, cycle_lengths)
            curr_row.append(val)
        char_table.append(curr_row)

    if memo_file_name:
        write_memo_to_file(memo_file_name)
    char_table = np.fliplr(np.array(char_table, dtype=np.int64))

    if output_file_name:
        np.savetxt(output_file_name, char_table, delimiter=",", fmt="%d")

    return char_table


def write_character_table_csv(
    n: int,
    output_file_name: str,
    memo_file_name: str = "memo.txt",
    checkpoint_interval: int = 0,
    max_memo_entries: int = 0,
    progress_file: str = "",
    resume: bool = True,
    log_interval: int = 100,
) -> None:
    """
    Write the character table of Sn to CSV one row at a time.
    This avoids constructing the full p(n) by p(n) NumPy array in memory.
    """
    progress_file = progress_file or default_progress_file(output_file_name)
    bit_strings = get_partition_bit_strings(n)
    sigma_cycle_lengths = [
        cycle_lengths_from_bit_string(sigma) for sigma in reversed(bit_strings)
    ]
    table_size = len(bit_strings)
    progress = read_progress(progress_file) if resume else {}
    start_row = 0
    if (
        resume
        and os.path.exists(output_file_name)
        and progress_matches(progress, n, table_size, "csv")
    ):
        start_row = int(progress.get("completed_rows", 0))
        if start_row >= table_size:
            LOGGER.info("CSV table already complete: %s", output_file_name)
            return
        LOGGER.info("Resuming CSV table %s from row %d", output_file_name, start_row)
    else:
        progress = make_progress(n, table_size, "csv", output_file_name, 0)
        write_progress(progress_file, progress)
        LOGGER.info("Starting CSV table %s with %d rows", output_file_name, table_size)

    if memo_file_name:
        read_memo_from_file(memo_file_name)

    mode = "a" if start_row else "w"
    start_time = time.perf_counter()
    with open(output_file_name, mode, newline="") as csv_file:
        writer = csv.writer(csv_file)
        for row_index, lambda_ in enumerate(bit_strings[start_row:], start=start_row):
            row = [
                murnaghan_nakayama_cycles(lambda_, cycle_lengths)
                for cycle_lengths in sigma_cycle_lengths
            ]
            writer.writerow(row)
            completed_rows = row_index + 1
            progress["completed_rows"] = completed_rows
            write_progress(progress_file, progress)
            enforce_memo_limit(max_memo_entries)
            if (
                memo_file_name
                and checkpoint_interval
                and completed_rows % checkpoint_interval == 0
            ):
                write_memo_to_file(memo_file_name)
            if log_interval and completed_rows % log_interval == 0:
                elapsed = time.perf_counter() - start_time
                LOGGER.info(
                    "CSV progress: %d/%d rows complete (%.2f%%), %.1fs elapsed",
                    completed_rows,
                    table_size,
                    100 * completed_rows / table_size,
                    elapsed,
                )

    if memo_file_name:
        write_memo_to_file(memo_file_name)
    LOGGER.info("Finished CSV table %s", output_file_name)


def write_character_table_npy(
    n: int,
    output_file_name: str,
    memo_file_name: str = "memo.txt",
    checkpoint_interval: int = 0,
    max_memo_entries: int = 0,
    progress_file: str = "",
    resume: bool = True,
    log_interval: int = 100,
) -> None:
    """
    Write the character table of Sn to an int64 .npy file backed by a memory
    map. This is the preferred full-table output format for S40-scale runs.
    """
    progress_file = progress_file or default_progress_file(output_file_name)
    bit_strings = get_partition_bit_strings(n)
    sigma_cycle_lengths = [
        cycle_lengths_from_bit_string(sigma) for sigma in reversed(bit_strings)
    ]
    table_size = len(bit_strings)
    progress = read_progress(progress_file) if resume else {}
    start_row = 0
    mode = "w+"

    if (
        resume
        and os.path.exists(output_file_name)
        and progress_matches(progress, n, table_size, "npy")
    ):
        start_row = int(progress.get("completed_rows", 0))
        if start_row >= table_size:
            LOGGER.info("NPY table already complete: %s", output_file_name)
            return
        mode = "r+"
        LOGGER.info("Resuming NPY table %s from row %d", output_file_name, start_row)
    else:
        progress = make_progress(n, table_size, "npy", output_file_name, 0)
        write_progress(progress_file, progress)
        LOGGER.info("Starting NPY table %s with %d rows", output_file_name, table_size)

    if memo_file_name:
        read_memo_from_file(memo_file_name)

    char_table = np.lib.format.open_memmap(
        output_file_name,
        mode=mode,
        dtype=np.int64,
        shape=(table_size, table_size),
    )
    start_time = time.perf_counter()
    for row_index, lambda_ in enumerate(bit_strings[start_row:], start=start_row):
        char_table[row_index, :] = [
            murnaghan_nakayama_cycles(lambda_, cycle_lengths)
            for cycle_lengths in sigma_cycle_lengths
        ]
        completed_rows = row_index + 1
        progress["completed_rows"] = completed_rows
        write_progress(progress_file, progress)
        enforce_memo_limit(max_memo_entries)
        if (
            memo_file_name
            and checkpoint_interval
            and completed_rows % checkpoint_interval == 0
        ):
            char_table.flush()
            write_memo_to_file(memo_file_name)
        if log_interval and completed_rows % log_interval == 0:
            elapsed = time.perf_counter() - start_time
            LOGGER.info(
                "NPY progress: %d/%d rows complete (%.2f%%), %.1fs elapsed",
                completed_rows,
                table_size,
                100 * completed_rows / table_size,
                elapsed,
            )

    char_table.flush()
    if memo_file_name:
        write_memo_to_file(memo_file_name)
    LOGGER.info("Finished NPY table %s", output_file_name)


def write_partition_labels_csv(n: int, output_file_name: str) -> None:
    """
    Write the abaci bit-string row labels and displayed column labels for Sn.
    Columns in character tables are reversed to match the existing API.
    """
    bit_strings = get_partition_bit_strings(n)
    with open(output_file_name, "w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["axis", "index", "bit_string"])
        for index, bit_string in enumerate(bit_strings):
            writer.writerow(["row", index, bit_string])
        for index, bit_string in enumerate(reversed(bit_strings)):
            writer.writerow(["column", index, bit_string])


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
    bit_strings = get_partition_bit_strings(n)

    if memo_file_name:
        read_memo_from_file(memo_file_name)
    char_table = []
    cycle_lengths = cycle_lengths_from_bit_string(col_bit_string)

    for lambda_ in bit_strings:
        curr_row = []
        val = murnaghan_nakayama_cycles(lambda_, cycle_lengths)
        curr_row.append(val)
        char_table.append(curr_row)

    if memo_file_name:
        write_memo_to_file(memo_file_name)
    char_table = np.fliplr(np.array(char_table, dtype=np.int64))

    if csv_file_name:
        np.savetxt(csv_file_name, char_table, delimiter=",", fmt="%d")

    return char_table


def get_staircase_character_column(n: int, memo_file_name: str = "") -> List[int]:
    """
    Return all character values on the staircase conjugacy class for S_n.
    Here n must be triangular, n = k(k + 1) / 2, and the conjugacy class is
    (k, k - 1, ..., 1).
    """
    cycle_lengths = staircase_cycle_lengths(n)

    if memo_file_name:
        read_memo_from_file(memo_file_name)
    values = [
        murnaghan_nakayama_cycles(lambda_, cycle_lengths)
        for lambda_ in iter_partition_bit_strings(n)
    ]
    if memo_file_name:
        write_memo_to_file(memo_file_name)
    return values


def write_staircase_character_column_csv(
    n: int,
    output_file_name: str,
    memo_file_name: str = "",
    checkpoint_interval: int = 0,
    max_memo_entries: int = 0,
    progress_file: str = "",
    resume: bool = True,
    log_interval: int = 1000,
) -> None:
    """
    Stream the staircase conjugacy-class column of the character table of S_n
    to a one-column CSV. This avoids building a large NumPy array and preserves
    exact Python integers.
    """
    cycle_lengths = staircase_cycle_lengths(n)
    progress_file = progress_file or default_progress_file(output_file_name)
    table_size = partition_number(n)
    progress = read_progress(progress_file) if resume else {}
    start_row = 0

    if (
        resume
        and os.path.exists(output_file_name)
        and progress_matches(progress, n, table_size, "staircase_csv")
    ):
        start_row = int(progress.get("completed_rows", 0))
        if start_row >= table_size:
            LOGGER.info("Staircase column already complete: %s", output_file_name)
            return
        LOGGER.info(
            "Resuming staircase column %s from row %d", output_file_name, start_row
        )
    else:
        progress = make_progress(n, table_size, "staircase_csv", output_file_name, 0)
        progress["cycle_lengths"] = list(cycle_lengths)
        write_progress(progress_file, progress)
        LOGGER.info(
            "Starting staircase column %s for S_%d with %d rows",
            output_file_name,
            n,
            table_size,
        )

    if memo_file_name:
        read_memo_from_file(memo_file_name)

    mode = "a" if start_row else "w"
    start_time = time.perf_counter()
    with open(output_file_name, mode, newline="") as csv_file:
        writer = csv.writer(csv_file)
        for row_index, lambda_ in enumerate(iter_partition_bit_strings(n)):
            if row_index < start_row:
                continue
            value = murnaghan_nakayama_cycles(lambda_, cycle_lengths)
            writer.writerow([value])
            completed_rows = row_index + 1
            progress["completed_rows"] = completed_rows
            write_progress(progress_file, progress)
            enforce_memo_limit(max_memo_entries)
            if (
                memo_file_name
                and checkpoint_interval
                and completed_rows % checkpoint_interval == 0
            ):
                write_memo_to_file(memo_file_name)
            if log_interval and completed_rows % log_interval == 0:
                elapsed = time.perf_counter() - start_time
                LOGGER.info(
                    "Staircase progress: %d/%d rows complete (%.2f%%), %.1fs elapsed",
                    completed_rows,
                    table_size,
                    100 * completed_rows / table_size,
                    elapsed,
                )

    if memo_file_name:
        write_memo_to_file(memo_file_name)
    LOGGER.info("Finished staircase column %s", output_file_name)


def read_memo_from_file(file_name: str = "memo.txt"):
    """
    Read existing memo file.
    """
    try:
        with open(file_name, "rb") as memo_file:
            global MEMO
            memo = pickle.load(memo_file)
            MEMO = {
                key: value
                for key, value in memo.items()
                if (
                    isinstance(key, tuple)
                    and len(key) == 2
                    and isinstance(key[0], str)
                    and isinstance(key[1], tuple)
                )
            }
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
    # N = 78
    # write_staircase_character_column_csv(N, f"S{N}_staircase.csv")
    get_character_table(20, "S20.csv", "memo_6.txt")
