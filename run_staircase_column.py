import argparse
import os
import time

import murnaghan_nakayama as mn


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute the staircase conjugacy-class character column."
    )
    parser.add_argument(
        "--n",
        type=int,
        required=True,
        help="Compute for S_n. Must be triangular: n = k(k + 1) / 2.",
    )
    parser.add_argument(
        "--output-dir",
        default="runs/staircase",
        help="Directory for the column CSV, logs, and progress metadata.",
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
        default=1000,
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
    return parser.parse_args()


def main():
    args = parse_args()
    k = mn.staircase_rank_from_n(args.n)
    os.makedirs(args.output_dir, exist_ok=True)

    prefix = f"S{args.n}_staircase"
    log_file = os.path.join(args.output_dir, f"{prefix}.log")
    output_file = os.path.join(args.output_dir, f"{prefix}.csv")

    mn.configure_logging(log_file)
    mn.clear_caches()

    start = time.perf_counter()
    table_size = mn.partition_number(args.n)
    mn.LOGGER.info(
        "Preparing staircase column for S_%d: k=%d, cycle type=%s, p(n)=%d",
        args.n,
        k,
        mn.staircase_cycle_lengths(args.n),
        table_size,
    )

    mn.write_staircase_character_column_csv(
        args.n,
        output_file,
        memo_file_name=args.memo_file,
        checkpoint_interval=args.checkpoint_interval,
        max_memo_entries=args.max_memo_entries,
        resume=not args.fresh,
        log_interval=args.log_interval,
    )

    elapsed = time.perf_counter() - start
    mn.LOGGER.info("Finished staircase column for S_%d in %.1f seconds", args.n, elapsed)


if __name__ == "__main__":
    main()
