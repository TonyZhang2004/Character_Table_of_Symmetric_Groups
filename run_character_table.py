import argparse
import os
import time

import numpy as np

import murnaghan_nakayama as mn


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute a symmetric-group character table with resumable output."
    )
    parser.add_argument("--n", type=int, default=40, help="Compute the table for S_n.")
    parser.add_argument(
        "--output-dir",
        default="runs/S40",
        help="Directory for table, labels, logs, and progress metadata.",
    )
    parser.add_argument(
        "--format",
        choices=("npy", "csv"),
        default="npy",
        help="Output format. Use npy for large tables such as S40.",
    )
    parser.add_argument(
        "--max-memo-entries",
        type=int,
        default=5_000_000,
        help="Clear recursive memo after this many entries. Use 0 for no cap.",
    )
    parser.add_argument(
        "--log-interval",
        type=int,
        default=100,
        help="Log progress every this many completed rows. Use 0 to disable.",
    )
    parser.add_argument(
        "--checkpoint-interval",
        type=int,
        default=0,
        help="Write memo file every this many rows. Use 0 to disable.",
    )
    parser.add_argument(
        "--memo-file",
        default="",
        help="Optional pickle memo file. Empty default avoids large memo files.",
    )
    parser.add_argument(
        "--fresh",
        action="store_true",
        help="Start over instead of resuming from progress metadata.",
    )
    parser.add_argument(
        "--skip-labels",
        action="store_true",
        help="Do not write row/column label CSV.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    prefix = f"S{args.n}"
    log_file = os.path.join(args.output_dir, f"{prefix}.log")
    labels_file = os.path.join(args.output_dir, f"{prefix}_labels.csv")
    output_file = os.path.join(args.output_dir, f"{prefix}.{args.format}")

    mn.configure_logging(log_file)
    mn.clear_caches()

    start = time.perf_counter()
    bit_strings = mn.get_partition_bit_strings(args.n)
    table_size = len(bit_strings)
    raw_bytes = table_size * table_size * np.dtype(np.int64).itemsize
    mn.LOGGER.info(
        "Preparing S_%d table: p(n)=%d, entries=%d, raw int64 size=%.2f GiB",
        args.n,
        table_size,
        table_size * table_size,
        raw_bytes / 1024**3,
    )

    if not args.skip_labels:
        mn.write_partition_labels_csv(args.n, labels_file)
        mn.LOGGER.info("Wrote labels to %s", labels_file)

    writer = (
        mn.write_character_table_npy
        if args.format == "npy"
        else mn.write_character_table_csv
    )
    writer(
        args.n,
        output_file,
        memo_file_name=args.memo_file,
        checkpoint_interval=args.checkpoint_interval,
        max_memo_entries=args.max_memo_entries,
        resume=not args.fresh,
        log_interval=args.log_interval,
    )

    elapsed = time.perf_counter() - start
    mn.LOGGER.info("Finished S_%d in %.1f seconds", args.n, elapsed)


if __name__ == "__main__":
    main()
