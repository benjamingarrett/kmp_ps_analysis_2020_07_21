# Usage: <folder>
# Assumes csize,misses
# Each csv contains two rows: COMMENTS & SORT_BY
# Investigating
# Conjecture: When the prefix-suffix function of our recursive version of KMP is memoized with LRU cache replacement on inputs of the form a^(n-1)b for all positive n, (cache size)(cache misses) is between (n^2)/2 and (n^2)/2 + n*sqrt(n) when the cache size < sqrt(n).


import os, sys, matplotlib.pyplot as plt
from math import sqrt
import numpy as np
from scipy import stats
from sklearn.metrics import r2_score


MAX_ITERATIONS = 5000000


def get_list(dr, fn):
  print('fn {}'.format(fn))
  return [ln.rstrip('\n').split(',') for ln in open(dr+'/'+fn)]


def get_sort_by(dr, fn):
  print('fn {}'.format(fn))
  lst = [ln.rstrip('\n').split(',') for ln in open(dr+'/'+fn)]
  j = 0
  while j < len(lst):
    if lst[j][0] == 'SORT_BY':
      return int(lst[j][1])
    j += 1
  print('SORT_BY field not found in {}/{}'.format(dr, fn))
  exit()


def bounds_above(pts, f, minx, maxx):
  j = 0
  while j < len(pts):
    if minx < pts[j][0] and pts[j][0] < maxx:
      if float(pts[j][1]) >= f(pts[j][0]):
        return False
    j += 1
  return True


def bounds_below(pts, f, minx, maxx):
  j = 0
  while j < len(pts):
    if minx < pts[j][0] and pts[j][0] < maxx:
      if float(pts[j][1]) <= f(pts[j][0]):
        return False
    j += 1
  return True


def points_from_csv_with_n(dr, fname, col=1, rows_to_skip=0):
  print('fn {}'.format(fname))
  n = get_sort_by(dr, fname)
  lst = get_list(dr, fname)[rows_to_skip:]
  return n, [(int(k[0]), int(k[col])) for k in lst]


print('args: {}'.format(sys.argv))
if len(sys.argv) < 2:
  print('Usage: <folder>')
  exit()
folder = sys.argv[1]
fnames = [fn for fn in os.listdir(folder) if 'csv' in fn]
point_sequences = {}

num_errors = 0
num_trials = 0
confirmed_list = []
for fn in fnames:
  n, lst = points_from_csv_with_n(folder, fn, rows_to_skip=2)
  minx = min([k[0] for k in lst])
  maxx = sqrt(n)
  if n in {9, 18, 37, 75, 150, 300}:
    pruned_lst = [(k[0], k[1]) for k in lst if k[0] < maxx and k[0] > minx]
    point_sequences[n] = pruned_lst
  lower_bound = lambda x: n**2 / (2*x)
  upper_bound = lambda x: n**2 / (2*x) + n*sqrt(n)/x
  num_trials += 1
  if (not bounds_below(lst, lower_bound, minx, maxx)) or (not bounds_above(lst, upper_bound, minx, maxx)):
    print('Found data not confirmed by conjecture in {}'.format(fn))
    num_errors += 1
  else:
    confirmed_list.append(n)

print('num_errors/trials: {}/{}'.format(num_errors, num_trials))
print('confirmed_list {}'.format(confirmed_list))

fig = plt.figure()
ax = fig.add_subplot(111)
for n, seq in point_sequences.items():
  minx = min([k[0] for k in seq])
  maxx = max([k[0] for k in seq])
  poly_pts_x = np.linspace(minx, maxx, 1000)
  ax.plot([k[0] for k in seq], [k[1] for k in seq], '*', label='data for n = {}'.format(n))
  ax.plot([k for k in poly_pts_x], [upper_bound(k) for k in poly_pts_x], label='upper bound for n = {}'.format(n))
  ax.plot([k for k in poly_pts_x], [lower_bound(k) for k in poly_pts_x], label='lower bound for n = {}'.format(n))
plt.xlabel('cache size')
plt.ylabel('cache misses')
t = 'KMP PS function, LRU, instance type: (a^(n-1))b'
t+= '\nconjecture: (cache misses)(cache size) is between n^2/2 and n^2/2 + n*sqrt(n) where cache size is between 1 and sqrt(n)'
ax.legend(loc='upper right')
plt.title(t)
plt.show()
