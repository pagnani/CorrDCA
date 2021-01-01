# CorrDCA

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://pagnani.github.io/CorrDCA/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://pagnani.github.io/CorrDCA/dev)
[![Build Status](https://github.com/pagnani/CorrDCA/workflows/CI/badge.svg)](https://github.com/pagnani/CorrDCA/actions)
[![Coverage](https://codecov.io/gh/pagnani/CorrDCA/branch/main/graph/badge.svg)](https://codecov.io/gh/pagnani/CorrDCA)

Simple package that:

* Reads FASTA files and translate it in a numerical matrix

* Computes the reweighting score as described in ["Fast and accurate multivariate Gaussian modeling of protein families: Predicting residue contacts and protein-interaction partners"](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0092721) by Carlo Baldassi, Marco Zamparo, Christoph Feinauer, Andrea Procaccini, Riccardo Zecchina, Martin Weigt and Andrea Pagnani, (2014) PLoS ONE 9(3): e92721. doi:10.1371/journal.pone.0092721

* Compute empirical frequency counts (with and without pseudocount)

## Usage

The computation of the sequence weights is typically very expensive computationaly. A considerable speed-up can be achieved by exploiting parallel computation. To do so, just start julia with the ` -p nprocs` argument where `nprocs` is the number of workers available on your machine. Alternatively, from julia REPL, just do a:

```
julia> using Distributed; 
julia> addprocs(8) # put here the number of cores available 
julia> @everywhere using CorrDCA
```
## Credits

All methods available here are present also in the (so far) unregistered ["GaussDCA"](https://github.com/carlobaldassi/GaussDCA.jl) package.

The `remove_duplicate_sequences` has been basically copied by it.