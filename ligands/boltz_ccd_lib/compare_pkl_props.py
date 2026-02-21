#!/usr/bin/env python3
import sys
import pickle


def main():
    ref = pickle.load(open(sys.argv[1], "rb"))
    new = pickle.load(open(sys.argv[2], "rb"))

    ref_props = set(ref.GetPropNames())
    new_props = set(new.GetPropNames())

    print("Missing in new:", sorted(ref_props - new_props))
    print("Extra in new:", sorted(new_props - ref_props))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: compare_pkl_props.py <ref.pkl> <new.pkl>")
        raise SystemExit(2)
    main()
