"""
    read_fasta_alignment(filename::AbstractString, max_gap_fraction::Real)

Return a `L × M` matrix of integers (`L` is the sequence length, and `M` is the number of sequences) of the 
multiple sequence alignment contained in the fasta file `filename` including all sequences with a fraction of gaps (`-`) ≤ 
`max_gap_fraction`.
"""
function read_fasta_alignment(filename::AbstractString, max_gap_fraction::Real)

    f = FastaReader(filename)

    max_gap_fraction = Float64(max_gap_fraction)

    # pass 1

    seqs = Int[]
    inds = Int[]
    fseqlen = 0

    for (name, seq) in f
        ngaps = 0
        if f.num_parsed == 1
            ls = length(seq)
            resize!(inds, ls)
            for i = 1:ls
                c = seq[i]
                if c != '.' && c == uppercase(c)
                    fseqlen += 1
                    inds[fseqlen] = i
                    c == '-' && (ngaps += 1)
                end
            end
        else
            ls = length(seq)
            ls == length(inds) || error("inputs are not aligned")
            tstfseqlen = 0
            for i = 1:ls
                c = seq[i]
                if c != '.' && c == uppercase(c)
                    tstfseqlen += 1
                    inds[tstfseqlen] == i || error("inconsistent inputs")
                    c == '-' && (ngaps += 1)
                end
            end
            tstfseqlen == fseqlen || error("inconsistent inputs")
        end
        ngaps / fseqlen <= max_gap_fraction && push!(seqs, f.num_parsed)
    end

    length(seqs) > 0 || error("Out of $(f.num_parsed) sequences, none passed the filter (max_gap_fraction=$max_gap_fraction)")

    # pass 2
    Z = Array{Int8}(undef, fseqlen, length(seqs))

    seqid = 1
    for (name, seq) in f
        seqs[end] < f.num_parsed && break
        seqs[seqid] == f.num_parsed || continue
        for i = 1:fseqlen
            c = seq[inds[i]]
            Z[i, seqid] = letter2num(c)
        end
        seqid += 1
    end
    @assert seqid == length(seqs) + 1

    close(f)

    return Z
end

let alphabet = [ 1,21, 2, 3, 4, 5, 6, 7, 8,21, 9,10,11,12,21,13,14,15,16,17,21,18,19,21,20]
               # A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y
    global letter2num
    function letter2num(c::Union{Char,UInt8})
        i = UInt8(c) - 0x40
        1 <= i <= 25 && return alphabet[i]
        return 21
    end
end

"""
    remove_duplicate_sequences(Z::Matrix{Ti}) where Ti<:Integer
Remove duplicate sequences (columns) in the alignment matrix `Z`

# Examples
```jldoctest
julia> Z = [1 2 3 1;
            1 3 2 1;]
2×4 Array{Int64,2}:
 1  2  3  1
 1  3  2  1

julia> remove_duplicate_sequences(Z)
removing duplicate sequences... done: 4 -> 3
([1 2 3; 1 3 2], [1, 2, 3])
```
"""
function remove_duplicate_sequences(Z::Matrix{Ti}) where Ti<:Integer

    N, M = size(Z)
    hZ = Array{UInt}(undef, M)
    @inbounds for i = 1:M
        hZ[i] = hash(Z[:,i])
    end
    print("removing duplicate sequences... ")
    ref_seq_ind = Array{Int}(undef, M)
    ref_seq = Dict{UInt,Int}()
    @inbounds for i = 1:M
        ref_seq_ind[i] = get!(ref_seq, hZ[i], i)
    end
    uniqueseqs = collect(values(ref_seq))
    # Check for collisions
    collided = falses(M)
    @inbounds for i in 1:M
        k = ref_seq_ind[i]
        k == i && continue
        for j = 1:N
            Z[j,i] != Z[j,k] && (collided[i] = true; break)
        end
    end
     if any(collided)
        nowcollided = BitArray(undef, M)
        while any(collided)
            # Collect index of first row for each collided hash
            empty!(ref_seq)
            @inbounds for i = 1:M
                collided[i] || continue
                ref_seq_ind[i] = get!(ref_seq, hZ[i], i)
            end
            for v in values(ref_seq)
                push!(uniqueseqs, v)
            end

            # Check for collisions
            fill!(nowcollided, false)
            @inbounds for i = 1:M
                k = ref_seq_ind[i]
                (!collided[i] || k == i) && continue
                for j = 1:N
                    Z[j,i] != Z[j,k] && (nowcollided[i] = true; break)
                end
            end
            collided, nowcollided = nowcollided, collided
        end
    end

    sort!(uniqueseqs)

    newM = length(uniqueseqs)
    newZ = Z[:,uniqueseqs]
    println("done: $M -> $newM")
    return newZ,uniqueseqs
end