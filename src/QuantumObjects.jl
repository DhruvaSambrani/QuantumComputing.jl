module QuantumObjects
import Base: *, ∘, kron, adjoint

export QuantumObject, Ket, Bra, inner

abstract type QuantumObject end

struct Ket <: QuantumObject
    coefficients::AbstractArray{ComplexF64,1}
end

inner(k1::Ket, k2::Ket) = k1' * k2

struct Bra <: QuantumObject
    k::Ket
end

Ket(b::Bra) = b.k

Base.adjoint(k::Ket) = Bra(k)
Base.adjoint(b::Bra) = Ket(b)

Base.:*(k1::Ket, k2::Ket) = inner(k1, k2)
Base.:*(b::Bra, k::Ket) = inner(Ket(b), k)

Base.kron(k1::Ket, k2::Ket) = Ket(k1.coefficients ⊗ k2.coefficients)
Base.kron(b1::Bra, b2::Bra) = Bra(b1.k ⊗ b2.k)

⊗ = Base.kron

end
