module QuantumObjects
import Base: *, ∘, kron, ⊗, adjoint

export QuantumObjects, Ket, Bra
abstract type QuantumObject end

struct Ket <: QuantumObject
    coefficients::AbstractArray{Complex,1}
end

Ket(b::Bra) = Ket(b.coefficients')
Base.adjoint(k::Ket) = Bra(k)

inner(k1::Ket, k2::Ket) = k1' * k2

struct Bra <: QuantumObject
    coefficients::AbstractMatrix{Complex,2}
end

Bra(k::Ket) = Bra(k.coefficients')
Base.adjoint(b::Bra) = Ket(b)

Base.:*(k1::Ket, k2::Ket) = inner(k1, k2)
Base.:*(b::Bra, k::Ket) = b.coefficients * k.coefficients

Base.kron(k1::Ket, k2::Ket) = Ket(k1.coefficients ⊗ k2.coefficients)
Base.kron(b1::Bra, b2::Bra) = Bra(b1.coefficients ⊗ b2.coefficients)
end
