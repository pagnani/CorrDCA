function compute_frequencies(Z::Matrix, W::Vector, q)
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
        Pijpc[:,:,i,i] .= (1-pc) * Pij_true[:,:,i,i]
        for a in 1:s
            Pijpc[a,a,i,i] = Pipc[a,i]
        end
    end
    Pijpc,Pipc
end

function Cij(Pi,Pij)
    q,N = size(Pi)
    cij = Matrix{Float64}(undef, N*q,N*q)
    mi = Vector{Float64}(undef,N*q)
    for i in 1:N
        rowblk = ((i-1)*q+1):i*q
        vpi = @view Pi[:,i]
        mi[rowblk].= vpi
        for j in 1:N
            vpj = @view Pi[:,j]
            vpij = @view Pij[:,:,i,j]
            colblk = ((j-1)*q+1):j*q
            cij[rowblk,colblk] .= vpij - vpi*vpj' 
        end 
    end
    return cij,mi
end


"""
    covariance_matrix(Z, W; pc::Real=0)
Compute the covariance matrix from numerical alignment `Z` and weights `W`. The output is a `N(q-1) × N(q-1)` 
(skipping color `q`) symmetric matrix. 
## Keywords arguments:
*  `pc` in [0,1]: pseudocount [default = `0`]
end
"""
function covariance_matrix(Z, W; pc::Real=0)
    0<=pc<=1 || error("pc = $pc should be ∈ [0,1]") 
    Pij,Pi = compute_frequencies(Z,W)
    if pc > 0
        Pijpc,Pipc=add_pseudocount(Pi,Pij,pc)
        return Cij(Pipc,Pijpc)
            
    end
    return Cij(Pi,Pij)
end

covariance_matrix(Z; kwds...) = covariance_matrix(Z,1.0/size(Z,2); kwds...)