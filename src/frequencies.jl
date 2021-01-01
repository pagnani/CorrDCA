function compute_frequencies(Z, W, q)
    N, M = size(Z)
    length(W) == M || error("incompatible length of W")
    Pij = zeros(q-1, q-1, N, N)
    for i in 1:N
        for j in i:N
            for s in 1:M
                coli = Z[i,s]
                colj = Z[j,s]
                (coli == q || colj == q) && continue
                Pij[colj,coli,j,i] += W[s]
                Pij[coli,colj,i,j] = Pij[colj,coli,j,i]
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
    Pipc  = (1.0-pc)*Pi_true  .+ pc/q
    Pijpc = (1.0-pc)*Pij_true .+ pc/q^2
    for i in 1:N
        for a in 1:s
            Pijpc[a,a,i,i] = Pipc[a,i]
        end
    end
    Pijpc,Pipc
end