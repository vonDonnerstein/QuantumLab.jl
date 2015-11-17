module IntegralsModule
export Overlap, normalize!, doublefactorial, computeMatrixOverlap, computeMatrixKinetic
using ..BaseModule
using ..BasisModule

doublefactorial(n::Int) = prod(n:-2:1)

function GaussianIntegral1D(mqn::Int,exponent::Float64)
  # I_x = Integrate[x^m Exp[-ζ x^2], {x,-∞,∞}] (acc. to Fundament. of Mol. Integr. Eval. by Fermann, Valeev)
  m = mqn
  ζ = exponent

  if (mqn%2!=0) # integral over uneven function
    return 0 
  else
    return (doublefactorial(m-1)*sqrt(π)) / ((2ζ)^(m/2)*sqrt(ζ))
  end
end


function GaussianIntegral1D(mqn::Int,exponent::Float64)
  # I_x = Integrate[x^m Exp[-ζ x^2], {x,-∞,∞}] (acc. to Mathematica 9)
  m = mqn
  ζ = exponent

  if (mqn%2!=0) # integral over uneven function
    return 0 
  else
    t=(m+1)/2
    return ζ^(-t) * gamma(t)
  end
end

#function FundamentalIntegral(
#  pgb1::PrimitiveGaussianBasisFunction,
#  pgb2::PrimitiveGaussianBasisFunction,
#  func::Function,
#  pgb3::PrimitiveGaussianBasisFunction,
#  pgb4::PrimitiveGaussianBasisFunction)
#  # (pgb1(x1) pgb2(x1) | func(x1,x2) | pgb3(x2) pgb4(x2) )
#  # pgb1..4 are considered s-Functions
#  # acc. to Molecular Integrals of Gaussian Basis Functions by Gill
#
#  # I = (π G_AB G_CD) / ( 2(ζ η)^(3/2) R^3 ) Integrate[u Sin[u] Exp[-u^2/(4T)] FT[f][u/R],{u,0,∞}]
#end

function GaussProductFundamental(
  pgb1::PrimitiveGaussianBasisFunction,
  pgb2::PrimitiveGaussianBasisFunction)
  # Fundamental means primitives are taken as s-Functions
  # K Exp[-γ r^2_P] = pgb1 * pgb2 (acc. to Fundament. of Mol. Integr. Eval. by Fermann, Valeev)

  α1 = pgb1.exponent
  α2 = pgb2.exponent
  A  = pgb1.center
  B  = pgb2.center

  γ = α1 + α2
  P = (α1*A + α2*B) / γ
  AB = distance(A,B)
  K = exp(-α1*α2*AB^2/γ)

  return (K,P,γ)
end

type PolynomialFactors
  x::Array{(Real,Int),1} # Real*x^Int
  y::Array{(Real,Int),1} # Real*y^Int
  z::Array{(Real,Int),1} # Real*z^Int
end

function GaussProductPolynomialFactor(
  pgb1::PrimitiveGaussianBasisFunction,
  pgb2::PrimitiveGaussianBasisFunction)
  # _Sum_K[f_k r^k]_ K Exp[-γ r^2_P] = pgb1 * pgb2 (acc. to Fundament. of Mol. Integr. Eval. by Fermann, Valeev)

  A  = pgb1.center
  B  = pgb2.center
  (K,P,γ) = GaussProductFundamental(pgb1,pgb2)

  factors = PolynomialFactors([],[],[])
  for (xyz in [:x,:y,:z])
    l1 = pgb1.mqn.(xyz)
    l2 = pgb2.mqn.(xyz)
    for k in 0:l1+l2
      ij = [((k+q) ÷ 2,(k-q) ÷ 2) for q in max(-k,k-2l2):2:min(k,2l1-k)]
      f = sum([binomial(l1,i) * binomial(l2,j) * (P-A).x^(l1-i) * (P-B).x^(l2-j) for(i,j) in ij])
      append!(factors.(xyz),[(f,k)])
    end
  end

  return factors
end

function Overlap(
  pgb1::PrimitiveGaussianBasisFunction,
  pgb2::PrimitiveGaussianBasisFunction)
  # S_ij = Integrate[phi_i(r) phi_j(r),{r el R}] (acc. to Fundament. of Mol. Integr. Eval. by Fermann, Valeev)

  (K,P,γ) = GaussProductFundamental(pgb1,pgb2)
  factors = GaussProductPolynomialFactor(pgb1,pgb2)

  Ix = sum([f*GaussianIntegral1D(i,γ) for (f,i) in factors.x])
  Iy = sum([f*GaussianIntegral1D(i,γ) for (f,i) in factors.y])
  Iz = sum([f*GaussianIntegral1D(i,γ) for (f,i) in factors.z])
  return K*Ix*Iy*Iz
end

function Overlap(
  cgb1::ContractedGaussianBasisFunction,
  cgb2::ContractedGaussianBasisFunction)

  integral = 0
  for (coeff1,pgb1) in zip(cgb1.coefficients,cgb1.primitiveBFs)
    for (coeff2,pgb2) in zip(cgb2.coefficients,cgb2.primitiveBFs)
      integral += coeff1*coeff2*Overlap(pgb1,pgb2)
    end
  end
  return integral
end

function normalize!(cgb::ContractedGaussianBasisFunction)
  N = Overlap(cgb,cgb)
  scale!(cgb.coefficients,1/sqrt(N))
end

function normalize!(basis::GaussianBasis)
  for (cgb in basis.contractedBFs)
    normalize!(cgb)
  end
end

function incrmqn(mqn::MQuantumNumber,xyz::Symbol,inc::Int)
  mqn_t = [mqn.x,mqn.y,mqn.z]
  xyznum = {:x => 1, :y => 2, :z => 3}
  num = xyznum[xyz]
  mqn_t[num] += inc
  return MQuantumNumber(mqn_t[1],mqn_t[2],mqn_t[3])
end

function KineticIntegral(
  pgb1::PrimitiveGaussianBasisFunction,
  pgb2::PrimitiveGaussianBasisFunction)
  #Ix + Iy + Iz = (pgb1 | 1/2 ∇^2 | pgb2) 
  #Ix = 1/2 l1 l2 <-1|-1> + 2 α1 α2 <+1|+1> - α1 l2 <+1|-1> - α2 l1 <-1|+1>
  #(acc. to Fundamentals of Mol. Integr. Eval. by Fermann, Valeev (eq. 4.1 + 4.13))

  integral = 0

  for (xyz in (:x,:y,:z))
    pgb1decr = pgb1
    pgb1decr.mqn = incrmqn(pgb1.mqn,xyz,-1)
    pgb1incr = pgb1
    pgb1incr.mqn = incrmqn(pgb1.mqn,xyz,1)
    pgb2decr = pgb2
    pgb2decr.mqn = incrmqn(pgb2.mqn,xyz,-1)
    pgb2incr = pgb2
    pgb2incr.mqn = incrmqn(pgb2.mqn,xyz,1)

    l1 = pgb1.mqn.(xyz)
    l2 = pgb2.mqn.(xyz)
    α1 = pgb1.exponent
    α2 = pgb2.exponent

    integral +=
  	(1/2*l1*l2) * Overlap(pgb1decr,pgb2decr) +
	(2*α1*α2)   * Overlap(pgb1incr,pgb2incr) +
	(-α1*l2)    * Overlap(pgb1incr,pgb2decr) +
	(-α2*l1)    * Overlap(pgb1decr,pgb2incr)
  end

  return integral
end

function KineticIntegral(
  cgb1::ContractedGaussianBasisFunction,
  cgb2::ContractedGaussianBasisFunction)

  integral = 0
  for (coeff1,pgb1) in zip(cgb1.coefficients,cgb1.primitiveBFs)
    for (coeff2,pgb2) in zip(cgb2.coefficients,cgb2.primitiveBFs)
      integral += coeff1*coeff2*KineticIntegral(pgb1,pgb2)
    end
  end
  return integral
end

function computeMatrixKinetic(basis::GaussianBasis)
  return [KineticIntegral(cgb1,cgb2) for cgb1 in basis.contractedBFs, cgb2 in basis.contractedBFs]
end

function computeMatrixOverlap(basis::GaussianBasis)
return [Overlap(cgb1,cgb2) for cgb1 in basis.contractedBFs, cgb2 in basis.contractedBFs]
end

function OverlapFundamental(
  pgb1::PrimitiveGaussianBasisFunction,
  pgb2::PrimitiveGaussianBasisFunction)
  (K,P,γ) = GaussProductFundamental(pgb1,pgb2)
  # Overlap(s-type pgb1, s-type pgb2) (acc. to. Mathematica)
  return K*(π/γ)^(3/2)
end

end # module
