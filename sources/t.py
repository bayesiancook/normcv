import sys
import math
from scipy.stats import t

print("const double studentC[] = {", end='')
print("{0}, {1}".format(0,0), end='')
for n in range(2,1001):
    print(", {0}".format(t.isf(0.025, n)), end='');
print("};")
