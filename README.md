# EMST
[![Build Status](https://travis-ci.org/lithom/EMST.jl.svg?branch=master)](https://travis-ci.org/lithom/EMST.jl)

This package provides an implementation of the algorithm presented in [1] to compute euclidean minimum spanning trees.

Usage is very simple:
```
  edges = compute_emst(x)
```

returns the edges of the emst of the dataset x
- x is a (d x n) matrix
- edges is a (n-1 x 2) matrix


[1] March, William B., Parikshit Ram, and Alexander G. Gray. "Fast euclidean minimum spanning tree: algorithm, analysis, and applications." Proceedings of the 16th ACM SIGKDD international conference on Knowledge discovery and data mining. ACM, 2010.
http://www.mlpack.org/papers/emst.pdf
