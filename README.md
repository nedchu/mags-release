# Graph Summarization: Compactness Meets Efficiency

Paper Link: [TODO: link of paper]

This repository contains a reference implementation of our SIGMOD 2024 paper:

> Graph Summarization: Compactness Meets Efficiency. Deming Chu, Fan Zhang, Wenjie Zhang, Ying Zhang, Xuemin Lin. SIGMOD 2024.



Please cite our paper if you use the code:

[TODO: bixtex]

## Requirements

- G++
- CMake
- OpenMP



## Build


Build with the code below in the project folder.

``` shell
mkdir build
cd build
cmake ..
make -j
```

After that, you will get four executable programs.

- `mags`: serial implementation of *Mags*
- `mags_dm`: serial implementation of *Mags-DM*
- `pmags`: parallel implementation of *Mags*
- `pmags_dm`: parallel implementation of *Mags-DM*

## Run

The serial implementation (`mags` and `mags_dm`) requires one parameter.

1. the input graph file path.

The parallel implementation (`pmags` and `pmags_dm`) requires 1-2 parameters.

1. the input graph file path.
2. the number of thread (optional, 40 by default)

Input graph: undirected graph; no self-loop; one edge in each line.


## Code Structure


The source code of our algorithms is under `./src`:
- `util.h, util.cpp, global.h`: some helping functions and definitions
- `graph.h, graph.cpp`: graph data structures
- `gsum.h, gsum.cpp`: serial graph summarization algorithms
- `pgsum.h, pgsum.cpp`: parallel graph summarization algorithms
- `/parallel_hashmap`: a hash map implementation (see [greg7mdp/parallel-hashmap](https://github.com/greg7mdp/parallel-hashmap))

The main programs that runs the algorithms are under `./run`