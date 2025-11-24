# TeraTools

TeraTools is a software suite for constructing and analyzing compressed full-text indexes on terabase-scale pangenomic datasets.

TeraTools requires the run-length Burrows-Wheeler transform (RLBWT) of the text/dataset as input. Given the RLBWT, TeraTools enables efficient construction of RLBWT-based indexes and their supporting data structures in $O(r)$ space, where $r$ is the number of runs in the BWT. This allows scalable analysis of massive, repetitive biological datasets such as human pangenomes and large bacterial collections.

## Scientific Motivation

Compressed indexes, especially those based on the RLBWT, enable fast substring queries using space proportional to the nonredundant content of the data. However, constructing these indexes and their support structures (such as LCP arrays) efficiently has remained a major challenge. TeraTools addresses this by providing algorithms that build RLBWT-based indexes and all supporting structures in compressed space for large, repetitive datasets.

## Main Features

- **RLBWT-Based Index Construction**: Efficient algorithms for building compressed full-text index support structures (LF, $\psi$, $\phi$, $\phi^{-1}$, LCP, PLCP) in $O(r)$ space and $O(n)$ time for repetitive datasets.
- **LCP and PLCP Summaries**: First practical algorithms for computing LCP-related information in $O(r)$, greatly reducing memory requirements for large datasets.
- **Matching Statistics and MEMs**: Fast computation of matching statistics and maximal exact matches (MEMs) using compressed indexes.
- **Scalable to Terabase Datasets**: Demonstrated on datasets such as the Human Pangenome Reference Consortium and large bacterial collections, with significant reductions in memory usage compared to previous methods.
- **Parallelized Implementations**: Parallelized for efficient use of multi-core systems.

## Why TeraTools?

Traditional methods for computing LCP-related structures require prohibitive memory for large datasets. For example:

- **HPRC v2** (466 human haplotypes): Previous tools require ~2.1 TB RAM
- **TeraTools**: Reduces memory to **170 GiB** (12.6x reduction)

TeraTools scales to terabase-sized datasets like CommonBacteria and human472.

## Peak Memory Usage (GiB) for LCP Summary Tools

|                                                                | mtb152 | chr19.1000 | human100      | human472        | CommonBacteria  |
| --------------                                                 | ------ | ---------- | ------------- | --------------- | --------------- |
| [rlbwt2lcp](https://github.com/nicolaprezza/rlbwt2lcp)         | 0.92   | 76.00      | (est.) 460.18 | (est.) 1,947.93 | (est.) 9,901.64 |
| [pfp-thresholds](https://github.com/maxrossi91/pfp-thresholds) | 0.70   | 40.66      | -             | 2,134.90        | -               |
| [TeraLCP](src/TeraLCP)                                         | 0.16   | 2.83       | 136.99        | 172.16          | 613.65          |

## Project Structure

- `src/` — Core source code for TeraTools modules:
  - `TeraLCP/` — RLBWT support structure construction (LF, $\psi$, $\phi$, $\phi^{-1}$, LCP, PLCP)
  - `TeraMEM/` — Maximal Exact Match (MEM) enumeration
  - `TeraMS/` — Matching statistics computation
  - `TeraIndex/` — Efficient index construction and querying
  - `include/` — Header files for modular code organization
  - `thirdparty/` — External dependencies (e.g., kseq.h for FASTA parsing)
- `data/` — Datasets, synthetic data generation scripts, and testing files
- `less_efficient_but_accurate_algos/` — Reference and baseline implementations for accuracy

## Testing

To test any of the four current programs ([TeraLCP](src/TeraLCP), [TeraMEM](src/TeraMEM), [TeraMS](src/TeraMS), [TeraIndex](src/TeraIndex)), run `make test` in its directory.

For example for TeraLCP:

```sh
cd src/TeraLCP && make test
```

This will compute the irreducible LCP/PLCP values of the idx and mtb152 datasets.

### Test Data

The repository includes a set of test data files to help verify the functionality of TeraTools modules. `make test` will run programs on two datasets: `idx` and `mtb152`. `idx` is a very small toy example with 2 strings and 20 basepairs. `mtb152` is 152 M. tuberculosis genomes obtained from <https://doi.org/10.5281/zenodo.13147120> with 304 strings and ~1.5 billion basepairs. All tests should take between a few minutes and a few seconds on a modern machine and <4 400 MiB of RAM.

These files are located in `data/testing/` and include:

- Example `.fmd` files and their corresponding `.len` files for various crop and single datasets
- Index files (`idx.fmd`, `idx.input`, `idx.len`)
- Reference files for [MTB152](https://doi.org/10.5281/zenodo.13147120) and other datasets (`mtb152.fmd`, `mtb152.fmd.ssa`, `mtb152.fmr`)

These files allow you to quickly test the construction and analysis features of TeraTools without needing to generate large datasets. Use them to validate installation, run example workflows, and ensure correct output from each module.

## Installation

```sh
git clone https://github.com/ucfcbb/TeraTools.git
cd TeraTools
make
```

### Input Requirements

TeraTools requires the **RLBWT** of your input text. You can generate this using existing tools:

- [ropebwt3](https://github.com/lh3/ropebwt3)

Use the provided Makefiles to build modules:

## Reference

If you use TeraTools in your research, please cite: TBA.

## Dependencies

- C++ compiler supporting C++17
- OpenMP

## License

MIT License. See `LICENSE` for details.
