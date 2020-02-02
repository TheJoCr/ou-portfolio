# Subset sum problem

## Problem

My crytography professor liked to assign us class challenges where
we would all be given the same problem (with slightly different
numbers) and asked to find the 'closest' numerical solution. We
would then be graded on a curve based on how accurate our results
were. Here, we were assigned a variant of the approximate subset sum
problem: Given summands equal to the 10^100 times the cube root
of the `i`h prime number for i from 1 to 100, find the subset of
those whose sum approximates your student id number * 2 * 10^94 
(as computed in `target` in the code).

## Solution

My solution relied on 2 key strategies:

1) Itterative reductions in complexity: Approximate subset sum
problems are much easier to solve with smaller numbers. By dividing
through the entire problem by 10^r for some `r`, I can test many more
potential solutions and hopefully find a better one more quickly.
Starting with a large value for `r` and then slowly decreasing, I
can ensure that I find some solutions quickly, but that the longer I
run the program, the more likely I am to find a solution. 

2) Rewrite the problem as finding orthogonal basis vectors for a
special lattice and solve using BKZ (a variant of LLL). It is well
known that the exact subset sum problem can be solved by finding
nearly orthogonal basis vectors to a lattice constructed from the
summands. For summands `a1`, `a2`, ... and sum `M`, we construct the
matrix:

[ [     ] -a1 ]
[ [  I  ] -a2 ]
[ [     ] -a3 ]
    ...   ... ]
[ 0 0  ...  M ]

If the simplified form of the matrix has a row of all 1s and 0s
(i.e. the lattice has a basis vector with just 1s and 0s), then that
row vector represents a solution to the exact subset sum problem.
Because the BKZ algorithm gives us an extremly quick way to compute
approximatly orthogonal basis vectors, we can simply make many trial
runs: given some set of summands and a target value `M`, we try to
solve the exact subset sum problem for `M+i` for `0 <= i < S`, where
S is some large value (10^6 in this code). For very large `M`, this
produces very good approximates. 


## Results

Given summands as described above and a target sum of:
```
target = 2261705720000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
```

I computed for about 15 CPU hours across 20 cores, for a total of
about 45 minutes of run time. The best sum vector found was: 

```python
solution = [ 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1,
0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0,
0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1,
0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1,
0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1 ]
```

For a total sum of:
```python
total_sum = 2261705720000000097852663298125902444624958626831127967977815635183370083714921603072309350146251303412
```

This is only about 9.79*10^85 larger than the target value - the
first 17 (0f 102) digits are exactly the same. 
