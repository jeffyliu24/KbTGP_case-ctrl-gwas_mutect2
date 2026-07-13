#!/usr/bin/env python
import argparse
import gzip


def open_text(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--tags", default="case,1kgp,ctrl")
    args = parser.parse_args()

    tags = [x.strip() for x in args.tags.split(",") if x.strip()]

    with open_text(args.input, "r") as fin, open_text(args.output, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            variant_id = fields[2]
            ref = fields[3]
            alt = fields[4]

            if all(tag in variant_id for tag in tags) and ref and alt:
                fout.write(line)


if __name__ == "__main__":
    main()
