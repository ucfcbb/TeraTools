# TeraTools

TeraTools is a software suite for constructing and analyzing compressed full-text indexes on terabase-scale pangenomic datasets.

TeraTools enables efficient construction of run-length Burrows-Wheeler transform (RLBWT)-based indexes and their supporting data structures in $O(r)$ space and $O(n)$ time for repetitive datasets, where $r$ is the number of runs in the BWT and $n$ is the text length. This allows scalable analysis of massive, repetitive biological datasets such as human pangenomes and large bacterial collections.

## Scientific Motivation

Traditional full-text indexes become infeasible for terabase-scale pangenomic datasets due to memory requirements. Compressed indexes, especially those based on the RLBWT, enable fast substring queries using space proportional to the nonredundant content of the data. However, constructing these indexes and their support structures (such as LCP arrays) efficiently has remained a major challenge. TeraTools addresses this by providing algorithms that build RLBWT-based indexes and all supporting structures in optimal time and space for large, repetitive datasets.

## Main Features

- **RLBWT-Based Index Construction**: Efficient algorithms for building compressed full-text indexes and support structures (LF, $\\psi$, $\phi$, $\phi^{-1}$, LCP, PLCP) in $O(r)$ space and $O(n)$ time for repetitive datasets.
- **LCP and PLCP Summaries**: First practical algorithms for computing LCP-related information in optimal space and time, greatly reducing memory requirements for large datasets.
- **Matching Statistics and MEMs**: Fast computation of matching statistics and maximal exact matches (MEMs) using compressed indexes.
- **Scalable to Terabase Datasets**: Demonstrated on datasets such as the Human Pangenome Reference Consortium and large bacterial collections, with significant reductions in memory usage compared to previous methods.
- **Parallelized Implementations**: Modules are parallelized for efficient use of multi-core systems.

## Why TeraTools?

Traditional methods for computing LCP-related structures require prohibitive memory for large datasets. For example:

- **HPRC v2** (466 human haplotypes): Previous tools require ~2.1 TB RAM
- **TeraTools**: Reduces memory to **170 GiB** (12.6x reduction)

TeraTools scales to terabase-sized datasets like CommonBacteria and human472 that were previously computationally infeasible.

## Peak Memory Usage (GiB) for LCP Summary Tools

|                | mtb152 | chr19.1000 | human100      | human472        | CommonBacteria  |
| -------------- | ------ | ---------- | ------------- | --------------- | --------------- |
| rlbwt2lcp      | 0.92   | 76.00      | (est.) 460.18 | (est.) 1,947.93 | (est.) 9,901.64 |
| pfp-thresholds | 0.70   | 40.66      | -             | 2,134.90        | -               |
| TeraLCP        | 0.16   | 2.83       | 136.99        | 172.16          | 613.65          |

RLBWT generation with sequential grlbwt had a lower memory consumption that all LCP summary tools for mtb152 and chr19.1000. For human100 and CommonBacteria, rlbwt generation with ropebwt3 takes $\leq 83$~GiB. For human472, rlbwt generation with ropebwt3 likely takes $<170$ GiB according to the 99 GiB memory usage on human320. pfp-thresholds is currently running on human100.

## Project Structure

- `src/` — Core source code for TeraTools modules:
  - `TeraLCP/` — RLBWT support structure construction (LF, $\psi$, $\phi$, $\phi^{-1}$, LCP, PLCP)
  - `TeraMEM/` — Maximal exact match (MEM) finding
  - `TeraMS/` — Matching statistics computation
  - `TeraIndex/` — Efficient index construction and querying
  - `include/` — Header files for modular code organization
  - `thirdparty/` — External dependencies (e.g., kseq.h for FASTA parsing)
- `data/` — Datasets, synthetic data generation scripts, and testing files
- `less_efficient_but_accurate_algos/` — Reference and baseline implementations for accuracy

## Test Data

The repository includes a set of test data files to help verify the functionality of TeraTools modules. These files are located in `data/testing/` and include:

- Example `.fmd` files and their corresponding `.len` files for various crop and single datasets
- Index files (`idx.fmd`, `idx.input`, `idx.len`)
- Reference files for MTB152 and other datasets (`mtb152.fmd`, `mtb152.fmd.ssa`, `mtb152.fmr`)

These files allow you to quickly test the construction and analysis features of TeraTools without needing to generate large datasets. Use them to validate installation, run example workflows, and ensure correct output from each module.

## Installation

### Input Requirements

TeraTools requires the **RLBWT** of your input text. You can generate this using existing tools:

- [ropebwt3](https://github.com/lh3/ropebwt3)

For DNA sequences, include both forward and reverse complement strands.

Use the provided Makefiles to build modules:

```sh
git clone https://github.com/ucfcbb/TeraTools.git
cd TeraTools
make -C src
make -C less_efficient_but_accurate_algos
make -C data/generation
```

## Usage

Each module can be run independently. See module documentation and comments for usage instructions and examples. Typical workflows include:

- Constructing RLBWT and support structures for large genomic datasets
- Computing LCP summaries and thresholds in compressed space
- Fast matching statistics and MEM enumeration for pangenome analysis

## Reference

If you use TeraTools in your research, please cite: TBA.

## Dependencies

- C++ compiler (e.g., g++, clang++)

## License

MIT License. See `LICENSE` for details.
