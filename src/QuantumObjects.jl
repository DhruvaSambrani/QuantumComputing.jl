module QuantumObjects
import Base: *, ∘, kron, adjoint
import LinearAlgebra: Hermitian, tr, normalize
export QuantumObject, Ket, Bra, DensityMatrix, inner, dims

abstract type QuantumObject{N} end

(dims(::QuantumObject{N}) where N) = N

struct Ket{N} <: QuantumObject{N}
    coefficients::AbstractArray{ComplexF64,1}
    function Ket(a::AbstractArray{T, 1}) where T<:Number
        new{length(a)}(normalize(a))
    end
end

struct Bra{N} <: QuantumObject{N}
    k::Ket{N}
end

Bra(a::AbstractArray{Number, 1}) = Bra(Ket(conj.(a)))

Ket(b::Bra) = b.k

Base.adjoint(k::Ket) = Bra(k)
Base.adjoint(b::Bra) = Ket(b)

Base.:*(k1::Ket, k2::Ket) = inner(k1, k2)
Base.:*(b::Bra, k::Ket) = inner(Ket(b), k)

struct DensityMatrix{N} <: QuantumObject{N}
    matrix::Hermitian{ComplexF64, Matrix{ComplexF64}}
end

inner(ρ::DensityMatrix, σ::DensityMatrix) = 1/2 * tr(ρ.matrix * σ.matrix) 
inner(k1::Ket, k2::Ket) = k1.coefficients' * k2.coefficients


Base.:*(ρ::DensityMatrix, σ::DensityMatrix) = inner(ρ, σ) 

Base.kron(k1::Ket, k2::Ket) = Ket(k1.coefficients ⊗ k2.coefficients)
Base.kron(b1::Bra, b2::Bra) = Bra(b1.k ⊗ b2.k)
Base.kron(ρ1::DensityMatrix, ρ2::DensityMatrix) = DensityMatrix(ρ1.matrix ⊗ ρ2.matrix)

⊗ = Base.kron

end
