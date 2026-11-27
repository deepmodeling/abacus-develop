#!/usr/bin/env python3
import argparse
import os
import re
import sys

import numpy as np


def sanitize_name(name):
    return re.sub(r"[^A-Za-z0-9_]", "_", name)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("descriptor_dir")
    args = parser.parse_args()

    descriptor_dir = args.descriptor_dir
    if not os.path.isdir(descriptor_dir):
        print(f"Descriptor dir not found: {descriptor_dir}", file=sys.stderr)
        return 1

    files = sorted(f for f in os.listdir(descriptor_dir) if f.endswith(".npy"))
    if not files:
        print(f"No .npy files in: {descriptor_dir}", file=sys.stderr)
        return 1

    for filename in files:
        path = os.path.join(descriptor_dir, filename)
        data = np.load(path, allow_pickle=False)
        if data.size == 0:
            mean_value = 0.0
        else:
            mean_value = float(np.mean(np.abs(data)))
        key = f"ml_desc_mean_{sanitize_name(filename)}"
        print(f"{key} {mean_value:.12f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
