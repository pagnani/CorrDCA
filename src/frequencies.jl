function compute_frequencies(Z, W, q)
    N, M = size(Z)
    length(W) == M || error("incompatible length of W")
    Pij = zeros(q-1, q-1, N, N)
    for i in 1:N - 1
        for j in i + 1:N
            for s in 1:M
                coli = Z[i,s]
                colj = Z[j,s]
                (coli == q || colj == q) && continue
                Pij[Z[i,s],Z[j,s],i,j] += W[s]
                Pij[Z[j,s],Z[i,s],j,i] = Pij[Z[i,s],Z[j,s],i,j]
            end
        end
    end
    Pi = zeros(q-1,N)
    for i in 1:N
        for s in 1:M
            coli = Z[i,s]
            (coli == q) && continue
            Pi[coli,i] += W[s]
        end
    end
    return Pij,Pi
end

compute_frequencies(Z,W) = compute_frequencies(Z,W,maximum(Z))
compute_frequencies(Z) = compute_frequencies(Z,ones(size(Z,2))/size(Z,2),maximum(Z))


function add_pseudocount(Pi_true::Array{Float64,2}, Pij_true::Array{Float64,4},pc)
    0 <= pc <= 1 || error("pseudocount $pc should be ∈ [0,1]") 
    s,N = size(Pi_true)
    q = s+1
    Pij_matrix = reshape(Pij_true,N*s,N*s)
    Pi_vector = reshape(Pi_true, N*s)
    Pipc,Pijpc = add_pseudocount(Pi_vector,Pij_matrix,pc,N,q)
    return reshape(Pipc,s,N), reshape(Pijpc,s,s,N,N)
end

function add_pseudocount(Pi_true::Vector{Float64}, Pij_true::Matrix{Float64},pc,N,q)
    pcq = pc / q
    Pij = (1 - pc) * Pij_true .+ pcq / q
    Pi = (1 - pc) * Pi_true .+ pcq
    s = q - 1
    i0 = 0
    for i = 1:N
        xr::UnitRange{Int}= i0 .+ (1:s)
	    Pij[xr, xr] = (1 - pc) * Pij_true[xr, xr]
        for alpha = 1:s
            x = i0 + alpha
            Pij[x, x] += pcq
	    end
        i0 += s
    end
    return Pi, Pij
end