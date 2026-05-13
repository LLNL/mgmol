#!/usr/bin/env python3
"""
Read multiple XY datasets from one text file and plot them together.

Expected file format (repeat blocks):
    <dataset title line>
    x1  y1
    x2  y2
    ...
    (blank line separates datasets)

- Title: any non-empty line that is NOT two numeric columns.
- Data lines: at least two whitespace- or comma-separated numbers (x y).
- You can include comments starting with #; they are ignored.
"""

import re
import argparse
import matplotlib.pyplot as plt

_NUM = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
_DATA_RE = re.compile(rf"^\s*({_NUM})\s*[,;\s]\s*({_NUM})\s*(?:#.*)?$")

def read_blocks(path):
  datasets = []
  title = None
  xs, ys = [], []

  def flush():
    nonlocal title, xs, ys
    if xs and ys:
        datasets.append((title, xs, ys))
    title, xs, ys = None, [], []

  with open(path, "r", encoding="utf-8") as lines:
    for raw in lines:
      line = raw.strip()
      if not line:
        # blank line ends a dataset
        flush()
        continue

      data = _DATA_RE.match(line)
      if data:
        xs.append(float(data.group(1)))
        ys.append(float(data.group(2)))
      else:
        # treat as a title line; starting a new title flushes any prior dataset
        if xs:
          flush()
        title = line
        print("Dataset {}".format(title))

  flush()
  return datasets

def main():
  ap = argparse.ArgumentParser()
  ap.add_argument("file", help="Input text file containing multiple XY blocks")
  ap.add_argument("--out", default=None, help="Save figure to file (e.g. plot.png)")
  args = ap.parse_args()
  print(args.file)
  datasets = read_blocks(args.file)
  if not datasets:
    raise SystemExit("No datasets found.")

  for name, xs, ys in datasets:
    plt.plot(xs, ys, label=name)

  plt.xlabel("radius")
  plt.ylabel("potential (Ry)")
  plt.grid(True, alpha=0.3)
  plt.legend()

  if args.out:
    plt.savefig(args.out, dpi=200)
  else:
    plt.show()

if __name__ == "__main__":
  main()
