module QuantumOperations
import Base: *, ∘, kron, adjoint
import LinearAlgebra: tr, Hermitian
using ..QuantumObjects

export QuantumOperation, PureOperator, expectation, Hermitian

abstract type QuantumOperation{N} <: QuantumObject{N} end

struct PureOperator{N} <: QuantumOperation{N}
    trm::AbstractMatrix{ComplexF64}
    function PureOperator(a::Hermitian{ComplexF64, Matrix{ComplexF64}})
        reduce(==, size(a)) || error("Must Be Square")
        tr(a) <= 1 || error("Must be trace non-increasing")
        new{size(a, 1)}(a)
    end
end

Base.adjoint(o1::PureOperator) = Operator(o1.trm')

Base.:∘(o1::PureOperator{N}, o2::PureOperator{N}) where N = Operator(o.trm * o2.trm)

(o::PureOperator{N})(k::Ket{N}) where N = Ket(o.trm * k.coefficients)
(o::PureOperator{N})(b::Bra{N}) where N = Bra(o(b.k))
(o::PureOperator{N})(ρ::DensityMatrix{N}) where N = DensityMatrix(o.trm * ρ.matrix * (o.trm)')
(o::PureOperator{N})(o2::PureOperator{N}) where N = o∘o2

Base.:*(o::PureOperator{N}, k::Ket{N}) where N = o(k)
Base.:*(b::Bra{N}, o::PureOperator{N}) where N = o(b)

Base.kron(o1::PureOperator, o2::PureOperator) = Operator(o1.trm ⊗ o2.trm)

expectation(o::PureOperator{N}, ψ::Ket{N}) where N = ψ' * o * ψ
expectation(o::PureOperator{N}, ρ::DensityMatrix{N}) where N = tr(o*ρ)

end
