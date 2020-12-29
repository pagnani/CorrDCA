module TestWeight
using Test,CorrDCA,DelimitedFiles

Z = [1 2 3 3 3;
    1 2 4 3 3;
    3 2 1 3 3;]

W,Meff = compute_weights(Z,0.9)

@test W == [0.25,0.25,0.25, 0.125,0.125]
@test Meff == 4.0
W,Meff = compute_weights(Z,:auto)
@test W == [0.25,0.25,0.25,0.125,0.125]
@test Meff == 4.0

datadir = joinpath(dirname(@__FILE__),"data")

fileal = joinpath(datadir,"pf14short.fasta.gz")
fileW = joinpath(datadir,"W.dat")
Z = read_fasta_alignment(fileal,0.8)
Wtrue = readdlm(fileW) |> vec
W,Meff=compute_weights(Z,:auto)
Mefftrue = sum(Wtrue)
@test Meff == Mefftrue
@test W .* Meff == Wtrue

printstyled("All TestWeight passed!\n",color=:light_green,bold=true)
end # end module