module QuantumOperations
import Base: *, ∘, kron, ⊗

export QuantumOperation, Operator

struct Operator <: QuantumOperation
    trm::AbstractMatrix{Complex,2}
end

(o::Operator)(k::Ket)::Ket = Ket(o.trm * k.coefficients)
(o::Operator)(b::Bra)::Bra = Bra(b.coefficients * o.trm)
(o::Operator)(o2::Operator)::Operator = Operator(o.trm * o2.trm)

Base.:*(o::Operator, k::Ket) = o(k)
Base.:*(b::Bra, o::Operator) = o(b)
Base.:*(o1::Operator, o2::Operator) = o1(o2)

Base.:∘(o1::Operator, o2::Operator) = o1(o2)
Base.kron(o1::Operator, o2::Operator) = Operator(o1.trm ⊗ o2.trm)

expectation(o::Operator, ψ::Ket) = inner(ψ, o(ψ))
end
