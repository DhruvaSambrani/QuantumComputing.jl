module QuantumObjects
import Base: *, ∘, kron, adjoint
import LinearAlgebra: Hermitian, tr
export QuantumObject, Ket, Bra, DensityMatrix, inner,normalization

abstract type QuantumObject end

normalization(coefficients::Array{ComplexF64,1})  = sqrt(sum((abs.(coefficients)).^2))

struct Ket <: QuantumObject
    coefficients::AbstractArray{ComplexF64,1}
    function Ket(coefficients)
        new(coefficients ./ normalization(coefficients))
    end  
end

struct Bra <: QuantumObject
    k::Ket
end

normalization(k1:: Ket) = normalization(k1.coefficients)
normalization(k1:: Bra) = normalization(Ket(k1))

Ket(b::Bra) = b.k

Base.adjoint(k::Ket) = Bra(k)
Base.adjoint(b::Bra) = Ket(b)

Base.:*(k1::Ket, k2::Ket) = inner(k1, k2)
Base.:*(b::Bra, k::Ket) = inner(Ket(b), k)

struct DensityMatrix <: QuantumObject
    matrix::Hermitian{ComplexF64, Matrix{ComplexF64}}
    function DensityMatrix(rho::Matrix{ComplexF64}, to_normalize=true:: Bool)
        if to_normalize
            rho = rho ./ tr(rho)
        end
        if rho != rho'
            println("DensityMatrix has to be Hermitian, so Hermitian(Array) is done automatically!")
        end 
            new(Hermitian(rho))
    end
end

inner(ρ::DensityMatrix, σ::DensityMatrix) = 1/2 * tr(ρ.matrix * σ.matrix) 
inner(k1::Ket, k2::Ket) = k1.coefficients' * k2.coefficients
Base.:*(ρ::DensityMatrix, σ::DensityMatrix) = inner(ρ, σ) 

Base.kron(k1::Ket, k2::Ket) = Ket(k1.coefficients ⊗ k2.coefficients)
Base.kron(b1::Bra, b2::Bra) = Bra(b1.k ⊗ b2.k)
Base.kron(ρ1::DensityMatrix, ρ2::DensityMatrix) = DensityMatrix(ρ1.matrix ⊗ ρ2.matrix)

⊗ = Base.kron

end
