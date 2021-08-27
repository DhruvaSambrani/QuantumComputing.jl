module QuantumInformation
using LinearAlgebra
export entropy,trace_distance

using ..QuantumObjects

partial_shannon(x) = x ≈ 0 ? float(x) : x*log(x)
entropy(prob::Array{Float64})::Float64 = - sum((partial_shannon.(prob)))
entropy(ρ::DensityMatrix):: Float64 = entropy(eigvals(ρ.matrix))

trace_distance(prob1::Array{Float64},prob2::Array{Float64}) = sum(abs,(prob1-prob2))
trace_distance(ρ::DensityMatrix,σ::DensityMatrix) = sum(abs,(eigvals(ρ.matrix - σ.matrix)))

end

