import csv
import gzip
import os
import tempfile
import unittest

import murnaghan_nakayama as mn
import verify_character_table as verifier


def write_small_table(n: int, output_file: str) -> None:
    mn.clear_caches()
    table = mn.get_character_table(n, memo_file_name="")
    with open(output_file, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerows(table)


def read_rows(csv_file: str):
    with open(csv_file, "r", newline="") as file:
        return list(csv.reader(file))


class VerifyCharacterTableTests(unittest.TestCase):
    def test_small_tables_pass_streaming_checks(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            for n in range(1, 9):
                csv_file = os.path.join(tmp_dir, f"S{n}.csv")
                write_small_table(n, csv_file)
                result = verifier.verify_character_table(n, csv_file)
                self.assertTrue(
                    result.passed,
                    [(check.name, check.details) for check in result.checks],
                )

    def test_gzip_csv_passes_streaming_checks(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            csv_file = os.path.join(tmp_dir, "S6.csv")
            gzip_file = os.path.join(tmp_dir, "S6.csv.gz")
            write_small_table(6, csv_file)
            with open(csv_file, "rb") as source, gzip.open(gzip_file, "wb") as target:
                target.write(source.read())

            result = verifier.verify_character_table(6, gzip_file)
            self.assertTrue(
                result.passed,
                [(check.name, check.details) for check in result.checks],
            )

    def test_small_tables_pass_full_orthogonality(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            for n in range(1, 8):
                csv_file = os.path.join(tmp_dir, f"S{n}.csv")
                write_small_table(n, csv_file)
                result = verifier.verify_character_table(
                    n, csv_file, full_orthogonality=True
                )
                self.assertTrue(
                    result.passed,
                    [(check.name, check.details) for check in result.checks],
                )

    def test_corrupt_first_column_fails_hook_length_check(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            csv_file = os.path.join(tmp_dir, "S6.csv")
            write_small_table(6, csv_file)
            rows = read_rows(csv_file)
            rows[2][0] = str(int(rows[2][0]) + 1)
            with open(csv_file, "w", newline="") as file:
                csv.writer(file).writerows(rows)

            result = verifier.verify_character_table(6, csv_file)
            self.assertFalse(result.passed)
            self.assertFalse(
                next(
                    check
                    for check in result.checks
                    if check.name == "first column hook-length dimensions"
                ).passed
            )

    def test_missing_row_fails_shape_check(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            csv_file = os.path.join(tmp_dir, "S6.csv")
            write_small_table(6, csv_file)
            rows = read_rows(csv_file)
            with open(csv_file, "w", newline="") as file:
                csv.writer(file).writerows(rows[:-1])

            result = verifier.verify_character_table(6, csv_file)
            self.assertFalse(result.passed)
            self.assertFalse(
                next(
                    check for check in result.checks if check.name == "shape and row width"
                ).passed
            )

    def test_missing_column_fails_shape_check(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            csv_file = os.path.join(tmp_dir, "S6.csv")
            write_small_table(6, csv_file)
            rows = read_rows(csv_file)
            rows[1] = rows[1][:-1]
            with open(csv_file, "w", newline="") as file:
                csv.writer(file).writerows(rows)

            result = verifier.verify_character_table(6, csv_file)
            self.assertFalse(result.passed)
            self.assertFalse(
                next(
                    check for check in result.checks if check.name == "shape and row width"
                ).passed
            )

    def test_non_integer_token_fails_parse_check(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            csv_file = os.path.join(tmp_dir, "S6.csv")
            write_small_table(6, csv_file)
            rows = read_rows(csv_file)
            rows[1][1] = "not-an-int"
            with open(csv_file, "w", newline="") as file:
                csv.writer(file).writerows(rows)

            result = verifier.verify_character_table(6, csv_file)
            self.assertFalse(result.passed)
            self.assertFalse(
                next(
                    check for check in result.checks if check.name == "integer parsing"
                ).passed
            )

    def test_corrupt_trivial_row_fails_trivial_row_check(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            csv_file = os.path.join(tmp_dir, "S6.csv")
            write_small_table(6, csv_file)
            rows = read_rows(csv_file)
            rows[0][1] = "2"
            with open(csv_file, "w", newline="") as file:
                csv.writer(file).writerows(rows)

            result = verifier.verify_character_table(6, csv_file)
            self.assertFalse(result.passed)
            self.assertFalse(
                next(check for check in result.checks if check.name == "trivial row").passed
            )

    def test_labels_file_passes_and_corruption_fails(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            csv_file = os.path.join(tmp_dir, "S6.csv")
            labels_file = os.path.join(tmp_dir, "S6_labels.csv")
            write_small_table(6, csv_file)
            mn.write_partition_labels_csv(6, labels_file)

            result = verifier.verify_character_table(
                6, csv_file, labels_file=labels_file
            )
            self.assertTrue(result.passed)

            rows = read_rows(labels_file)
            rows[1][2] = "corrupt"
            with open(labels_file, "w", newline="") as file:
                csv.writer(file).writerows(rows)

            result = verifier.verify_character_table(
                6, csv_file, labels_file=labels_file
            )
            self.assertFalse(result.passed)
            self.assertFalse(
                next(check for check in result.checks if check.name == "labels file").passed
            )


if __name__ == "__main__":
    unittest.main()
