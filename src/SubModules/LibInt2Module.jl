# LibInt2 installation:
# sudo apt-get install libgmp3-dev libeigen3-dev
# ../libint/configure CXX=g++ CXXFLAGS="-I/usr/include/eigen3 -std=c++11" --enable-shared
# dos2unix -f libtool
# make -j8
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.

module LibInt2Module
export LibInt2Shell, destroy!, LibInt2Engine, LibInt2EngineCoulomb, LibInt2OneBodyEngine, LibInt2EngineOverlap, libInt2Initialize, libInt2Finalize
export computeMatrixOverlap, computeElectronRepulsionIntegral, computeBasisShellsLibInt2
import Base.convert
import Base.display, ..IntegralsModule.computeElectronRepulsionIntegral
using TensorOperations
using ..BaseModule
using ..GeometryModule
using ..BasisSetModule
using ..ShellModule


## Bitstypes

## LibInt2Shell
#  Type Declaration
bitstype 64 LibInt2Shell64
bitstype 32 LibInt2Shell32
if is(Int,Int32)
  typealias LibInt2Shell LibInt2Shell32
else
  typealias LibInt2Shell LibInt2Shell64
end

#  Constructors
"""
The libint2 library expects the coefficients to be input the same way as with basis set definition files. It then renormalizes the coefficients accordingly.
This behavior is the default for efficiency reasons. To directly enter the coefficients set the renorm flag to false.
When converting a Shell into a LibInt2Shell (e.g. with convert()), this is taken care of automatically.
"""
function LibInt2Shell(origin::Vector{Float64},lqn::Int,nprim::Int,exponents::Vector{Float64},coefficients::Vector{Float64}; renorm::Bool=true)
  if renorm==true
    return reinterpret(LibInt2Shell,ccall((:_Z11createShellPdiiS_S_,"libint2/libint2jl.so"),Ptr{Void},(Ptr{Float64},Int,Int,Ptr{Float64},Ptr{Float64}),origin,lqn,nprim,exponents,coefficients))
  else
    return convert(LibInt2Shell,Shell(LQuantumNumber(lqn),Position(origin...),exponents,coefficients))
  end
end

function LibInt2Shell(origin::Position,lqn::LQuantumNumber,exponents::Vector{Float64},coefficients::Vector{Float64}; renorm::Bool=true)
  orig = convert(Vector{Float64},origin)
  nprim = length(coefficients)
  LibInt2Shell(orig,lqn.exponent,nprim,exponents,coefficients;renorm=renorm)
end

function convert(::Type{LibInt2Shell},sh::Shell)
  scaledcoefficients = Float64[]
  for (coeff,expon) in zip(sh.coefficients,sh.exponents)
    pgb = PrimitiveGaussianBasisFunction(origin,expon,MQuantumNumber(sh.lqn.exponent,0,0))
    push!(scaledcoefficients,coeff * sqrt(IntegralsModule.Overlap(pgb,pgb)))
  end

  LibInt2Shell([sh.center.x,sh.center.y,sh.center.z],sh.lqn.exponent,length(sh.exponents),sh.exponents,scaledcoefficients)
end

#  Destructor
function destroy!(shell::LibInt2Shell)
  ccall((:_Z12destroyShellPN7libint25ShellE,"libint2/libint2jl.so"),Void,(LibInt2Shell,),shell)
end

#  Further Functions
function display(sh::LibInt2Shell)
  ccall((:_Z10printShellPN7libint25ShellE,"libint2/libint2jl.so"),Void,(LibInt2Shell,),sh)
end

function lqn(l2sh::LibInt2Shell)
  # cmp. convert(::Type{Shell},::LibInt2Shell)
  sh_ptrptrptr = reinterpret(Ptr{Ptr{Ptr{Cdouble}}},l2sh)
  contract_ptrptr = unsafe_load(sh_ptrptrptr,4)
  contract_ptr = reinterpret(Ptr{Cint},contract_ptrptr)
  return Int(unsafe_load(contract_ptr,1))
end

function nprims(l2sh::LibInt2Shell)
  # cmp. convert(::Type{Shell},::LibInt2Shell)
  sh_ptrptr    = reinterpret(Ptr{Ptr{Cdouble}},l2sh)
  alpha_ptr = unsafe_load(sh_ptrptr,1)
  alphaEnd_ptr = unsafe_load(sh_ptrptr,2)
  return Int(div(alphaEnd_ptr-alpha_ptr,8))
end


function convert(::Type{Shell},l2sh::LibInt2Shell)
  # LibInt2Shell objects reside in memory as
  # class Shell {
  #   std::vector<real_t> alpha;
  #   std::vector<Contraction> contr;
  #   std::array<real_t, 3> O;   		//!< this is inlined in memory
  #   std::vector<real_t> max_ln_coeff;
  # }
  # where
  # struct Contraction {
  #   int l;
  #   bool pure;
  #   std::vector<real_t> coeff;
  # }
  # A std::vector is of the same size as 3 ptrs = 24 bytes on 64bit with gnu compiler.
  # The first pointer gives the start of the allocated space, the second pointer gives the pointer
  # beyond the last vector element and the third pointer beyond the allocated space of the vector 
  # (cmp. std::vector.reserve()).

  # I assume one contraction per shell -> no shared-sp

  sh_ptrptrptr = reinterpret(Ptr{Ptr{Ptr{Cdouble}}},l2sh)
  sh_ptrptr    = reinterpret(Ptr{Ptr{Cdouble}},l2sh)
  sh_ptr       = reinterpret(Ptr{Cdouble},l2sh)

  alpha_ptr = unsafe_load(sh_ptrptr,1)
  alphaEnd_ptr = unsafe_load(sh_ptrptr,2)
  primitivecount = div(alphaEnd_ptr-alpha_ptr,8)
  expons = pointer_to_array(alpha_ptr,primitivecount)
  
  center_ptr = sh_ptr + 24*2
  center = Position(pointer_to_array(center_ptr,3)...)

  contract_ptrptr = unsafe_load(sh_ptrptrptr,4)
  coeff_ptr = unsafe_load(contract_ptrptr,2) # the int and bool together take the space of one pointer
  coeffs = pointer_to_array(coeff_ptr,primitivecount)

  contract_ptr = reinterpret(Ptr{Cint},contract_ptrptr)
  lqn = LQuantumNumber(Int(unsafe_load(contract_ptr,1)))

  return Shell(lqn,center,expons,coeffs)
end


## LibInt2Engine
#  Type Declaration
bitstype 64 LibInt2Engine64
bitstype 32 LibInt2Engine32
if is(Int,Int32)
  typealias LibInt2Engine LibInt2Engine32
else
  typealias LibInt2Engine LibInt2Engine64
end

#  Constructors
function LibInt2EngineCoulomb(maxNumberPrimitives::Int,maxLQN::LQuantumNumber)
  reinterpret(LibInt2Engine,ccall((:_Z19createEngineCoulombii,"libint2/libint2jl.so"),Ptr{Void},(Cint,Cint),maxNumberPrimitives,maxLQN.exponent))
end

#  Destructor
function destroy!(engine::LibInt2Engine)
  ccall((:_Z13destroyEnginePN7libint26EngineE,"libint2/libint2jl.so"),Void,(LibInt2Engine,),engine)
end

## LibInt2OneBodyEngine
#  Type Declaration
bitstype 64 LibInt2OneBodyEngine64
bitstype 32 LibInt2OneBodyEngine32
if is(Int,Int32)
  typealias LibInt2OneBodyEngine LibInt2OneBodyEngine32
else
  typealias LibInt2OneBodyEngine LibInt2OneBodyEngine64
end

#  Constructors
function LibInt2EngineOverlap(maxNumberPrimitives::Int,maxLQN::LQuantumNumber)
  reinterpret(LibInt2OneBodyEngine,ccall((:_Z19createEngineOverlapii,"libint2/libint2jl.so"),Ptr{Void},(Cint,Cint),maxNumberPrimitives,maxLQN.exponent))
end

#  Destructor
function destroy!(engine::LibInt2OneBodyEngine)
  ccall((:_Z20destroyOneBodyEnginePN7libint213OneBodyEngineE,"libint2/libint2jl.so"),Void,(LibInt2OneBodyEngine,),engine)
end


## library functionality
function libint2Initialize()
  ccall((:_ZN7libint210initializeEv,"libint2/libint2jl.so"),Void,())
end

function libint2Finalize()
  ccall((:_ZN7libint28finalizeEv,"libint2/libint2jl.so"),Void,())
end

function computeMatrixBlockOverlap(engine::LibInt2OneBodyEngine, μ::LibInt2Shell, ν::LibInt2Shell)
  μmqns = div((lqn(μ)+1)^2+(lqn(μ)+1),2)
  νmqns = div((lqn(ν)+1)^2+(lqn(ν)+1),2)
  buf = ccall((:_Z14computeOverlapPN7libint213OneBodyEngineEPNS_5ShellES3_,"libint2/libint2jl.so"),Ptr{Cdouble},(LibInt2OneBodyEngine,LibInt2Shell,LibInt2Shell), engine, μ,ν)
  return reshape(pointer_to_array(buf,μmqns*νmqns),(μmqns,νmqns))
end

function computeMatrixBlockOverlap(μlib::LibInt2Shell,νlib::LibInt2Shell)
  (μ, ν) = map(sh->convert(Shell,sh), (μlib, νlib))
  maxprims = max(length(μ.coefficients), length(ν.coefficients))
  maxlqn   = max(μ.lqn, ν.lqn)
  engine = LibInt2EngineOverlap(maxprims,maxlqn)

  result = copy(computeMatrixBlockOverlap(engine,μlib,νlib))
  
  destroy!(engine)
  return result
end

function computeElectronRepulsionIntegral(engine::LibInt2Engine, μ::LibInt2Shell, ν::LibInt2Shell, λ::LibInt2Shell, σ::LibInt2Shell)
  μmqns = div((lqn(μ)+1)^2+(lqn(μ)+1),2)
  νmqns = div((lqn(ν)+1)^2+(lqn(ν)+1),2)
  λmqns = div((lqn(λ)+1)^2+(lqn(λ)+1),2)
  σmqns = div((lqn(σ)+1)^2+(lqn(σ)+1),2)
  buf = ccall((:_Z10computeERIPN7libint26EngineEPNS_5ShellES3_S3_S3_,"libint2/libint2jl.so"),Ptr{Cdouble},(LibInt2Engine,LibInt2Shell,LibInt2Shell,LibInt2Shell,LibInt2Shell), engine, μ,ν,λ,σ)
  return reshape(pointer_to_array(buf,μmqns*νmqns*λmqns*σmqns),(μmqns,νmqns,λmqns,σmqns))
end

function computeElectronRepulsionIntegral(μlib::LibInt2Shell, νlib::LibInt2Shell, λlib::LibInt2Shell, σlib::LibInt2Shell)
  (μ, ν, λ, σ) = map(sh->convert(Shell,sh), (μlib, νlib, λlib, σlib))
  maxprims = max(length(μ.coefficients), length(ν.coefficients), length(λ.coefficients), length(σ.coefficients))
  maxlqn   = max(μ.lqn, ν.lqn, λ.lqn, σ.lqn)
  engine = LibInt2EngineCoulomb(maxprims,maxlqn)

  result = copy(computeElectronRepulsionIntegral(engine, μlib,νlib,λlib,σlib))

  destroy!(engine)
  return result
end

function computeBasisShellsLibInt2(basSet::BasisSet,geo::Geometry)
  shells = Vector{LibInt2Shell}()
  for atom in geo.atoms
    for contractedDefinition in basSet.definitions[atom.element]
      exponents = [prim.exponent for prim in contractedDefinition.primitives]
      coefficients = [prim.prefactor for prim in contractedDefinition.primitives]
      sh = LibInt2Shell(atom.position,contractedDefinition.lQuantumNumber,exponents,coefficients)
      push!(shells,sh)
    end
  end
  return shells
end

end # Module

LibInt2Module.libint2Initialize()
