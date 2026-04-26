import argparse
import csv
import math
import sys
import time
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set


FULL_ORTHOGONALITY_LIMIT = 1000
MAX_REPORTED_ERRORS = 10


@dataclass
class CheckResult:
    name: str
    passed: bool
    details: List[str] = field(default_factory=list)


@dataclass
class VerificationResult:
    passed: bool
    rows_scanned: int
    elapsed_seconds: float
    checks: List[CheckResult]


class ProgressReporter:
    def __init__(self, total_rows: int, interval: int = 100):
        self.total_rows = total_rows
        self.interval = interval
        self.start = time.perf_counter()
        self.last_length = 0
        self.last_rows_scanned = -1

    def update(self, rows_scanned: int, force: bool = False) -> None:
        if self.interval <= 0:
            return
        if rows_scanned == self.last_rows_scanned:
            return
        if not force and rows_scanned % self.interval != 0:
            return

        elapsed = time.perf_counter() - self.start
        percent = 100 * rows_scanned / self.total_rows if self.total_rows else 100
        rate = rows_scanned / elapsed if elapsed else 0
        remaining_rows = max(self.total_rows - rows_scanned, 0)
        eta = remaining_rows / rate if rate else 0
        bar_width = 30
        filled = round(bar_width * min(rows_scanned, self.total_rows) / self.total_rows)
        bar = "#" * filled + "-" * (bar_width - filled)
        message = (
            f"Verifying [{bar}] {rows_scanned}/{self.total_rows} "
            f"({percent:5.1f}%) elapsed={elapsed:.1f}s eta={eta:.1f}s"
        )
        padding = " " * max(0, self.last_length - len(message))
        print(f"\r{message}{padding}", end="", file=sys.stderr, flush=True)
        self.last_length = len(message)
        self.last_rows_scanned = rows_scanned

    def finish(self, rows_scanned: int) -> None:
        if self.interval <= 0:
            return
        self.update(rows_scanned, force=True)
        print(file=sys.stderr, flush=True)


def partition_dict_to_list(partition: Dict[int, int]) -> List[int]:
    parts = []
    for part in sorted(partition, reverse=True):
        parts.extend([part] * partition[part])
    return parts


def partitions(n: int) -> Iterable[Dict[int, int]]:
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
    bit_strings = []
    for part in parts:
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


def get_partition_bit_strings(n: int) -> List[str]:
    return make_bit_strings(list(partitions(n)))


def hook_length_dimension(partition: Dict[int, int]) -> int:
    parts = partition_dict_to_list(partition)
    n = sum(parts)
    hook_product = 1
    for row, row_length in enumerate(parts):
        for col in range(row_length):
            below = sum(1 for later_row in parts[row + 1 :] if later_row > col)
            hook_product *= row_length - col + below
    return math.factorial(n) // hook_product


def conjugacy_class_size(partition: Dict[int, int]) -> int:
    n = sum(part * count for part, count in partition.items())
    denominator = 1
    for part, count in partition.items():
        denominator *= (part**count) * math.factorial(count)
    return math.factorial(n) // denominator


def sign_character_value(n: int, cycle_type: Dict[int, int]) -> int:
    number_of_cycles = sum(cycle_type.values())
    return -1 if (n - number_of_cycles) % 2 else 1


def selected_row_indices(table_size: int, sample_rows: int) -> Set[int]:
    if table_size <= 0:
        return set()

    indices = {0, table_size // 2, table_size - 1}
    if sample_rows <= 0:
        return indices
    if sample_rows == 1:
        indices.add(0)
        return indices

    for i in range(sample_rows):
        indices.add(round(i * (table_size - 1) / (sample_rows - 1)))
    return indices


def add_error(errors: List[str], message: str) -> None:
    if len(errors) < MAX_REPORTED_ERRORS:
        errors.append(message)


def parse_integer_row(raw_row: Sequence[str], row_index: int, errors: List[str]):
    parsed = []
    for col_index, value in enumerate(raw_row):
        try:
            parsed.append(int(value))
        except ValueError:
            add_error(
                errors,
                f"row {row_index}, column {col_index}: cannot parse integer {value!r}",
            )
            return None
    return parsed


def read_csv_rows(csv_file: str) -> Iterable[List[str]]:
    with open(csv_file, "r", newline="") as file:
        yield from csv.reader(file)


def verify_labels_file(n: int, labels_file: str, table_size: int) -> CheckResult:
    expected_bit_strings = get_partition_bit_strings(n)
    expected_rows = [["axis", "index", "bit_string"]]
    expected_rows.extend(
        ["row", str(index), bit_string]
        for index, bit_string in enumerate(expected_bit_strings)
    )
    expected_rows.extend(
        ["column", str(index), bit_string]
        for index, bit_string in enumerate(reversed(expected_bit_strings))
    )

    errors = []
    try:
        with open(labels_file, "r", newline="") as file:
            actual_rows = list(csv.reader(file))
    except OSError as exc:
        return CheckResult("labels file", False, [f"could not read labels file: {exc}"])

    if len(actual_rows) != 1 + 2 * table_size:
        add_error(
            errors,
            f"expected {1 + 2 * table_size} label rows, found {len(actual_rows)}",
        )

    for index, expected in enumerate(expected_rows):
        if index >= len(actual_rows):
            break
        if actual_rows[index] != expected:
            add_error(
                errors,
                f"label row {index}: expected {expected}, found {actual_rows[index]}",
            )

    return CheckResult("labels file", not errors, errors)


def verify_character_table(
    n: int,
    csv_file: str,
    labels_file: str = "",
    sample_rows: int = 20,
    full_orthogonality: bool = False,
    force: bool = False,
    progress_interval: int = 0,
) -> VerificationResult:
    start = time.perf_counter()
    row_partitions = list(partitions(n))
    column_partitions = list(reversed(row_partitions))
    table_size = len(row_partitions)
    n_factorial = math.factorial(n)
    dimensions = [hook_length_dimension(partition) for partition in row_partitions]
    class_sizes = [
        conjugacy_class_size(partition) for partition in column_partitions
    ]
    sign_values = [
        sign_character_value(n, partition) for partition in column_partitions
    ]
    sample_indices = selected_row_indices(table_size, sample_rows)

    checks = []
    rows_scanned = 0
    shape_errors: List[str] = []
    parse_errors: List[str] = []
    first_column_errors: List[str] = []
    trivial_row_errors: List[str] = []
    sign_row_errors: List[str] = []
    row_norm_errors: List[str] = []
    full_orthogonality_errors: List[str] = []
    dim_square_sum = 0
    full_rows: List[List[int]] = []
    progress: Optional[ProgressReporter] = None
    if progress_interval:
        progress = ProgressReporter(table_size, progress_interval)

    if full_orthogonality and table_size > FULL_ORTHOGONALITY_LIMIT and not force:
        checks.append(
            CheckResult(
                "full row orthogonality guard",
                False,
                [
                    "refusing full orthogonality for "
                    f"p(n)={table_size}; use --force to override"
                ],
            )
        )
        full_orthogonality = False

    try:
        row_iterator = read_csv_rows(csv_file)
        for row_index, raw_row in enumerate(row_iterator):
            rows_scanned += 1
            if row_index >= table_size:
                add_error(
                    shape_errors,
                    f"extra row {row_index}: expected exactly {table_size} rows",
                )
                continue
            if not raw_row:
                add_error(shape_errors, f"row {row_index}: blank row")
                continue
            if len(raw_row) != table_size:
                add_error(
                    shape_errors,
                    f"row {row_index}: expected {table_size} columns, "
                    f"found {len(raw_row)}",
                )
                continue

            row = parse_integer_row(raw_row, row_index, parse_errors)
            if row is None:
                continue

            first_value = row[0]
            expected_dimension = dimensions[row_index]
            dim_square_sum += first_value * first_value
            if first_value != expected_dimension:
                add_error(
                    first_column_errors,
                    f"row {row_index}, column 0: expected "
                    f"{expected_dimension}, found {first_value}",
                )

            if row_index == 0:
                for col_index, value in enumerate(row):
                    if value != 1:
                        add_error(
                            trivial_row_errors,
                            f"row 0, column {col_index}: expected 1, found {value}",
                        )

            if row_index == table_size - 1:
                for col_index, value in enumerate(row):
                    expected_sign = sign_values[col_index]
                    if value != expected_sign:
                        add_error(
                            sign_row_errors,
                            f"row {row_index}, column {col_index}: expected "
                            f"{expected_sign}, found {value}",
                        )

            if row_index in sample_indices:
                norm = sum(
                    class_size * value * value
                    for class_size, value in zip(class_sizes, row)
                )
                if norm != n_factorial:
                    add_error(
                        row_norm_errors,
                        f"row {row_index}: expected weighted norm "
                        f"{n_factorial}, found {norm}",
                    )

            if full_orthogonality:
                full_rows.append(row)

            if progress:
                progress.update(rows_scanned)
    except OSError as exc:
        checks.append(CheckResult("csv file readable", False, [str(exc)]))
    finally:
        if progress:
            progress.finish(rows_scanned)

    if rows_scanned < table_size:
        add_error(
            shape_errors,
            f"expected {table_size} rows, found {rows_scanned}",
        )

    checks.extend(
        [
            CheckResult("shape and row width", not shape_errors, shape_errors),
            CheckResult("integer parsing", not parse_errors, parse_errors),
            CheckResult(
                "first column hook-length dimensions",
                not first_column_errors,
                first_column_errors,
            ),
            CheckResult("trivial row", not trivial_row_errors, trivial_row_errors),
            CheckResult("sign row", not sign_row_errors, sign_row_errors),
            CheckResult(
                "identity-column dimension squares",
                dim_square_sum == n_factorial and rows_scanned == table_size,
                []
                if dim_square_sum == n_factorial and rows_scanned == table_size
                else [
                    f"expected sum(dim^2)={n_factorial}, "
                    f"found {dim_square_sum}"
                ],
            ),
            CheckResult("sampled row norms", not row_norm_errors, row_norm_errors),
        ]
    )

    if labels_file:
        checks.append(verify_labels_file(n, labels_file, table_size))

    if full_orthogonality and not full_orthogonality_errors:
        if len(full_rows) != table_size:
            add_error(
                full_orthogonality_errors,
                f"expected {table_size} parsed rows, found {len(full_rows)}",
            )
        else:
            for i, row_i in enumerate(full_rows):
                for j in range(i, table_size):
                    row_j = full_rows[j]
                    inner_product = sum(
                        class_size * value_i * value_j
                        for class_size, value_i, value_j in zip(
                            class_sizes, row_i, row_j
                        )
                    )
                    expected = n_factorial if i == j else 0
                    if inner_product != expected:
                        add_error(
                            full_orthogonality_errors,
                            f"rows {i},{j}: expected weighted inner product "
                            f"{expected}, found {inner_product}",
                        )
        checks.append(
            CheckResult(
                "full row orthogonality",
                not full_orthogonality_errors,
                full_orthogonality_errors,
            )
        )

    elapsed = time.perf_counter() - start
    return VerificationResult(
        passed=all(check.passed for check in checks),
        rows_scanned=rows_scanned,
        elapsed_seconds=elapsed,
        checks=checks,
    )


def print_result(result: VerificationResult) -> None:
    for check in result.checks:
        status = "PASS" if check.passed else "FAIL"
        print(f"{status} {check.name}")
        for detail in check.details:
            print(f"  {detail}")
    print(f"Rows scanned: {result.rows_scanned}")
    print(f"Elapsed seconds: {result.elapsed_seconds:.2f}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Verify a generated symmetric-group character table CSV."
    )
    parser.add_argument("--n", type=int, required=True, help="Verify a table for S_n.")
    parser.add_argument("--csv-file", required=True, help="Character table CSV file.")
    parser.add_argument(
        "--labels-file",
        default="",
        help="Optional labels CSV from write_partition_labels_csv.",
    )
    parser.add_argument(
        "--sample-rows",
        type=int,
        default=20,
        help="Number of evenly spaced rows to include in row-norm checks.",
    )
    parser.add_argument(
        "--full-orthogonality",
        action="store_true",
        help="Check all weighted row inner products. Intended for small n.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Allow expensive checks beyond default safety limits.",
    )
    parser.add_argument(
        "--progress-interval",
        type=int,
        default=100,
        help="Update the stderr progress bar every this many rows. Use 0 to disable.",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable progress output.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = verify_character_table(
        args.n,
        args.csv_file,
        labels_file=args.labels_file,
        sample_rows=args.sample_rows,
        full_orthogonality=args.full_orthogonality,
        force=args.force,
        progress_interval=0 if args.no_progress else args.progress_interval,
    )
    print_result(result)
    return 0 if result.passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
