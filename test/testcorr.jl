module TestCorr
using Test,CorrDCA,DelimitedFiles

datadir = joinpath(dirname(@__FILE__),"data")
filePijtrue = joinpath(datadir,"Pij.dat")
filePitrue = joinpath(datadir,"Pi.dat")
fileal = joinpath(datadir,"pf14short.fasta.gz")
Z = read_fasta_alignment(fileal,0.8)[1:3,1:100]
Pijtrue = readdlm(filePijtrue)
Pitrue = readdlm(filePitrue)
W,Meff = compute_weights(Z,:auto)
Pij,Pi = Pij,Pi=compute_frequencies(Z,W)
q  = maximum(Z)
N  = size(Z,1)
s = q-1
@test sum(abs2,reshape(Pitrue,s,3)-Pi) < 1e-20
for i in 1:N
    for j in i:N
        idx =  ((i-1)*s+1):i*s
        jdx =  ((j-1)*s+1):j*s
        pmat = Pijtrue[idx, jdx] 
        vp = view(Pij,:,:,i,j)
        delta = sum(abs2, pmat-vp) 
        #println("i=$i j=$j delta=$delta")
        @test pmat ≈ vp
    end
end

Pijpc,Pipc=add_pseudocount(Pi,Pij,0.2)

@test permutedims(Pijpc,[2,1,4,3]) == Pijpc
for i in 1:N
    for a in 1:s
        @test Pijpc[a,a,i,i] ≈ Pipc[a,i]
    end
end

printstyled("All TestCorr passed!\n",color=:light_green,bold=true)
end # end module