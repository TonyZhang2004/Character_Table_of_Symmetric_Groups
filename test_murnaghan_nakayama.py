import csv
import json
import math
import os
import tempfile
import unittest

import numpy as np

import murnaghan_nakayama as mn


def partition_dict_to_list(partition):
    parts = []
    for part in sorted(partition, reverse=True):
        parts.extend([part] * partition[part])
    return parts


def hook_length_dimension(partition):
    parts = partition_dict_to_list(partition)
    n = sum(parts)
    hook_product = 1
    for row, row_length in enumerate(parts):
        for col in range(row_length):
            below = sum(1 for later_row in parts[row + 1 :] if later_row > col)
            hook_product *= row_length - col + below
    return math.factorial(n) // hook_product


def conjugacy_class_size(partition):
    n = sum(part * count for part, count in partition.items())
    denominator = 1
    for part, count in partition.items():
        denominator *= (part**count) * math.factorial(count)
    return math.factorial(n) // denominator


class MurnaghanNakayamaTests(unittest.TestCase):
    def setUp(self):
        mn.MEMO = {}
        mn.RIM_HOOK_CACHE = {}

    def test_identity_column_matches_hook_length_dimensions(self):
        for n in range(1, 9):
            parts = list(mn.partitions(n))
            table = mn.get_character_table(n, memo_file_name="")
            dimensions = np.array([hook_length_dimension(part) for part in parts])
            np.testing.assert_array_equal(table[:, 0], dimensions)

    def test_row_orthogonality(self):
        for n in range(1, 8):
            parts = list(mn.partitions(n))
            table = mn.get_character_table(n, memo_file_name="")
            class_sizes = np.array(
                [conjugacy_class_size(part) for part in reversed(parts)],
                dtype=np.int64,
            )
            gram = table @ np.diag(class_sizes) @ table.T
            np.testing.assert_array_equal(gram, math.factorial(n) * np.eye(len(parts)))

    def test_npy_writer_resumes_from_progress_file(self):
        n = 7
        completed_rows = 4
        expected = mn.get_character_table(n, memo_file_name="")
        with tempfile.TemporaryDirectory() as tmp_dir:
            output_file = os.path.join(tmp_dir, "S7.npy")
            progress_file = mn.default_progress_file(output_file)
            partial = np.lib.format.open_memmap(
                output_file,
                mode="w+",
                dtype=np.int64,
                shape=expected.shape,
            )
            partial[:completed_rows, :] = expected[:completed_rows, :]
            partial.flush()
            mn.write_progress(
                progress_file,
                mn.make_progress(n, expected.shape[0], "npy", output_file, completed_rows),
            )

            mn.clear_caches()
            mn.write_character_table_npy(n, output_file, memo_file_name="", log_interval=0)

            actual = np.load(output_file, mmap_mode="r")
            np.testing.assert_array_equal(actual, expected)
            with open(progress_file, "r") as file:
                progress = json.load(file)
            self.assertEqual(progress["completed_rows"], expected.shape[0])

    def test_csv_writer_resumes_from_progress_file(self):
        n = 6
        completed_rows = 3
        expected = mn.get_character_table(n, memo_file_name="")
        with tempfile.TemporaryDirectory() as tmp_dir:
            output_file = os.path.join(tmp_dir, "S6.csv")
            progress_file = mn.default_progress_file(output_file)
            with open(output_file, "w", newline="") as file:
                writer = csv.writer(file)
                writer.writerows(expected[:completed_rows, :])
            mn.write_progress(
                progress_file,
                mn.make_progress(n, expected.shape[0], "csv", output_file, completed_rows),
            )

            mn.clear_caches()
            mn.write_character_table_csv(n, output_file, memo_file_name="", log_interval=0)

            actual = np.loadtxt(output_file, delimiter=",", dtype=np.int64)
            np.testing.assert_array_equal(actual, expected)
            with open(progress_file, "r") as file:
                progress = json.load(file)
            self.assertEqual(progress["completed_rows"], expected.shape[0])

    def test_csv_writer_preserves_integers_larger_than_int64(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            output_file = os.path.join(tmp_dir, "large.csv")
            value = 2**63
            with open(output_file, "w", newline="") as file:
                writer = csv.writer(file)
                writer.writerow([value])
            with open(output_file, "r") as file:
                actual = int(next(csv.reader(file))[0])
            self.assertEqual(actual, value)


if __name__ == "__main__":
    unittest.main()
