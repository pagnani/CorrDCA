function coldist(Z,i,j,N)
    s = 0
    @inbounds @simd for k in 1:N
        s += Z[k,i] != Z[k,j]
    end
    s
end

function compute_theta(Z::Matrix{Ti}) where Ti <: Integer
    N,M = size(Z)
    meanfracid = @distributed (+) for i = 1:M-1
        s = 0.0
        for j = i+1:M            
            s += 1.0 - coldist(Z,i,j,N) / N
        end
        s
    end
    meanfracid /= 0.5 * M * (M-1)
    theta = min(0.5, 0.38 * 0.32 / meanfracid)
    return theta
end



"""
    compute_weights(Z::Matrix{Ti}, theta::Real) where Ti <: Integer
Compute the normalized counts of the number of sequences at hamming distance ≤ `theta` from any given sequence 
in `Z`.
"""
function compute_weights(Z::Matrix{Ti}, theta::Real) where Ti <: Integer
    N, M = size(Z)
    thresh = floor(N * theta)
    #W = ones(M) |> SharedArray
    W = zeros(M) |> SharedArray
    #sZ = SharedArray(Z)
    @sync @distributed for i in 1:M
        @simd for j in 1:M
            dij = coldist(Z,i,j,N)
            val = (dij < thresh)
            W[i] += val
        end
    end
    
    @sync @distributed for i in 1:M
        W[i] = 1.0 / W[i]
    end
    Meff = sum(W)
    return sdata(W)./Meff , Meff
end

"""
    compute_weights(Z::Matrix{Ti}, theta::Symbol) where Ti<:Integer
Compute the normalized counts of the number of sequences at hamming distance ≤ of a precomputed optimal threshold.
See `Fast and accurate multivariate Gaussian modeling of protein families: Predicting residue contacts and 
protein-interaction partners` [https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0092721]
and in particular the Supplementary Information at section `2 Reweighting Scheme` for details.
"""
function compute_weights(Z::Matrix{Ti}, theta::Symbol) where Ti<:Integer
    theta !== :auto && error("only theta=:auto is supported")
    thr = compute_theta(Z)
    return compute_weights(Z,thr)
end
