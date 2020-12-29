module CorrDCA

using Distributed: @distributed
using SharedArrays: SharedArray, sdata
using FastaIO: FastaReader

export read_fasta_alignment, remove_duplicate_sequences,compute_weights, compute_frequencies,add_pseudocount

include("readfasta.jl")
include("reweighting.jl")
include("frequencies.jl")
end
