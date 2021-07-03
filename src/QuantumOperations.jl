module QuantumOperations
import Base: *, ∘, kron, adjoint
using ..QuantumObjects

export QuantumOperation, Operator, expectation

abstract type QuantumOperation end

struct Operator <: QuantumOperation
    trm::AbstractMatrix{ComplexF64}
end

Base.adjoint(o1::Operator) = Operator(o1.trm')

Base.:∘(o1::Operator, o2::Operator) = Operator(o.trm * o2.trm)

(o::Operator)(k::Ket)::Ket = Ket(o.trm * k.coefficients)
(o::Operator)(b::Bra)::Bra = Bra(o(b.k))
(o::Operator)(ρ::DensityMatrix)::DensityMatrix = o.trm * ρ.matrix * (o.trm)'

Base.:*(o::Operator, k::Ket) = o(k)
Base.:*(b::Bra, o::Operator) = o(b)
Base.:*(o1::Operator, o2::Operator) = o1(o2)

Base.kron(o1::Operator, o2::Operator) = Operator(o1.trm ⊗ o2.trm)

expectation(o::Operator, ψ::Ket) = ψ' * o * ψ

end
