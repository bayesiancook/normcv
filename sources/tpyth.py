import sys
import math
from scipy.stats import t

nmax = 1000

print("student_table = {")
for n in range(1,nmax):
    if n < nmax-1:
        print("\t{0}:{1},".format(n,t.isf(0.025, n)))
    else:
        print("\t{0}:{1}".format(n,t.isf(0.025, n)))
print("}")
