"""Calcultation of atomic distances."""

import math

# from pytrr import GroTrrReader

f = open("data.gro", "r")
# print(f.read(5))
# lines = f.readlines()
# print(lines[2])
count = 0
for line in f:
    if line.startswith("G"):
        count += 1
print(count)
