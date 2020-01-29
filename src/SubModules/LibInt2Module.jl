module LibInt2Module
using Libdl

  """
  The maximal LQuantumNumber to which libint2 has been compiled. Compare :/deps/usr/src/Makefile: MAX_AM.
  Creating shells with higher LQuantumNumber causes segfaults within the libint2.
  """
  const MAX_AM = 7


export libInt2Initialize, libInt2Finalize
export LibInt2Shell
export LibInt2Engine, LibInt2EngineCoulomb, LibInt2EngineOverlap, LibInt2EngineKinetic, LibInt2EngineNuclearAttraction
export destroy!, lqn, nprims
export computeBasisShellsLibInt2
import ..IntegralsModule.computeMatrixBlockOverlap
import ..IntegralsModule.computeMatrixBlockKinetic
import ..IntegralsModule.computeMatrixBlockNuclearAttraction
import ..IntegralsModule.computeTensorBlockElectronRepulsionIntegrals
import ..ShellModule.computeFactorsNonAxial



import ...QuantumLab.libint2_available

if (libint2_available) # the normal case

  import ...QuantumLab.lib_path
  libint2jl = Ptr{Nothing}(0)
  function __init__()
    push!(Libdl.DL_LOAD_PATH,lib_path)
    try
      Libdl.dlopen("libint2-QuantumLab.so")
      global libint2jl = Libdl.dlopen("libint2jl.so")
      LibInt2Module.libInt2Initialize()
    catch e
      @warn "QuantumLab was compiled with libint2_QuantumLab.so. Maybe it has disappeared?"
      error(e)
    end
  end

  
  
  
  import Base.convert
  import Base.show
  import ..ShellModule.nbf
  using TensorOperations
  using ..BaseModule
  using ..GeometryModule
  using ..BasisSetModule
  using ..ShellModule
  using ..BasisFunctionsModule
  using ..IntegralsModule
  using ..AtomModule


  ## Bitstypes

  ## LibInt2Shell
  #  Type Declaration
  primitive type LibInt2Shell64 <: AbstractShell 64 end
  primitive type LibInt2Shell32 <: AbstractShell 32 end
  if Int === Int32
    const LibInt2Shell = LibInt2Shell32
  else
    const LibInt2Shell = LibInt2Shell64
  end

  #  Constructors
  """
      LibInt2Shell(origin::Vector{Float64}, lqn::Int, nprim::Int, exponents::Vector{Float64}, coefficients::Vector{Float64}; renorm::Bool=false)
  The libint2 library expects the coefficients to be input the same way as with basis set definition files. It then renormalizes the coefficients accordingly.
  This behavior is the default for efficiency reasons. To directly enter the coefficients set the renorm flag to false.
  When converting a Shell into a LibInt2Shell (e.g. with convert()), this is taken care of automatically.
  Note, that coefficients are generally only specified up to a global scaling factor - only relative factors are handled by renormalization.
  """
  function LibInt2Shell(origin::Vector{Float64},lqn::Int,nprim::Int,exponents::Vector{Float64},coefficients::Vector{Float64}; renorm::Bool=true)
    if lqn > MAX_AM
      error("libint2 has only been compiled to MAX_AM $MAX_AM.
      Can't create LibInt2Shell with lqn larger than that ($lqn).
      Consider recompiling to higher MAX_AM, or use `Shell`s instead of `LibInt2Shell`s.")
    end
    if renorm==true
      return reinterpret(LibInt2Shell,ccall(dlsym(libint2jl,:_Z11createShellPdiiS_S_),Ptr{Nothing},(Ptr{Float64},Int,Int,Ptr{Float64},Ptr{Float64}),origin,lqn,nprim,exponents,coefficients))
    else
      return convert(LibInt2Shell,Shell(Position(origin...),LQuantumNumber(lqn),exponents,coefficients;renorm=false))
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
      push!(scaledcoefficients,coeff * sqrt(computeIntegralOverlap(pgb,pgb)))
    end

    LibInt2Shell([sh.center.x,sh.center.y,sh.center.z],sh.lqn.exponent,length(sh.exponents),sh.exponents,scaledcoefficients)
  end

  #  Destructor
  function destroy!(shell::LibInt2Shell)
    ccall(dlsym(libint2jl,:_Z12destroyShellPN7libint25ShellE),Nothing,(LibInt2Shell,),shell)
  end

  #  Further Functions
  function show(io::IO,::MIME"text/plain",sh::LibInt2Shell)
    ccall(dlsym(libint2jl,:_Z10printShellPN7libint25ShellE),Nothing,(LibInt2Shell,),sh)
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

  nbf(shell::LibInt2Shell) = BaseModule.numberMQNsCartesian(LQuantumNumber(lqn(shell)))

  function computeDimensions(shells::Vector{LibInt2Shell})
    totaldim = 0
    maxprims = 0
    maxlqn = 0
    for sh in shells
      totaldim += nbf(sh)
      maxprims = maxprims < nprims(sh) ? nprims(sh) : maxprims
      maxlqn = maxlqn < lqn(sh) ? lqn(sh) : maxlqn
    end
    return (totaldim, maxprims, LQuantumNumber(maxlqn))
  end

  function convert(::Type{Shell},l2sh::LibInt2Shell)
    # LibInt2Shell objects reside in memory as
    # class Shell {
    #   std::vector<real_t> alpha;
    #   std::vector<Contraction> contr;
    #   std::array<real_t, 3> O;        //!< this is inlined in memory
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
    expons = unsafe_wrap(Array,alpha_ptr,primitivecount)

    center_ptr = sh_ptr + 24*2
    center = Position(unsafe_wrap(Array,center_ptr,3)...)

    contract_ptrptr = unsafe_load(sh_ptrptrptr,4)
    coeff_ptr = unsafe_load(contract_ptrptr,2) # the int and bool together take the space of one pointer
    coeffs = unsafe_wrap(Array,coeff_ptr,primitivecount)

    contract_ptr = reinterpret(Ptr{Cint},contract_ptrptr)
    lqn = LQuantumNumber(Int(unsafe_load(contract_ptr,1)))

    return Shell(center,lqn,expons,coeffs;renorm=false)
  end


  ## LibInt2Engine
  #  Type Declaration
  primitive type LibInt2Engine64 64 end
  primitive type LibInt2Engine32 32 end
  if Int === Int32
    const LibInt2Engine = LibInt2Engine32
  else
    const LibInt2Engine = LibInt2Engine64
  end

  #  Constructors
  function LibInt2EngineCoulomb(maxNumberPrimitives::Int,maxLQN::LQuantumNumber)
    reinterpret(LibInt2Engine,ccall(dlsym(libint2jl,:_Z19createEngineCoulombii),Ptr{Nothing},(Cint,Cint),maxNumberPrimitives,maxLQN.exponent))
  end
  function LibInt2EngineKinetic(maxNumberPrimitives::Int,maxLQN::LQuantumNumber)
    reinterpret(LibInt2Engine,ccall(dlsym(libint2jl,:_Z19createEngineKineticii),Ptr{Nothing},(Cint,Cint),maxNumberPrimitives,maxLQN.exponent))
  end
  function LibInt2EngineNuclearAttraction(maxNumberPrimitives::Int,maxLQN::LQuantumNumber,atom::Atom)
    l2atoms = [(atom.element.atomicNumber,atom.position.x,atom.position.y,atom.position.z)]
    reinterpret(LibInt2Engine,ccall(dlsym(libint2jl,:_Z29createEngineNuclearAttractioniiPN7libint24AtomEi),Ptr{Nothing},(Cint,Cint,Ptr{Nothing},Cint),maxNumberPrimitives,maxLQN.exponent,l2atoms,length(l2atoms)))
  end
  function LibInt2EngineNuclearAttraction(maxNumberPrimitives::Int,maxLQN::LQuantumNumber,geo::Geometry)
    l2atoms = Vector{Tuple{Int,Float64,Float64,Float64}}()
    for atom in geo.atoms
      push!(l2atoms,Tuple{Int,Float64,Float64,Float64}((atom.element.atomicNumber,atom.position.x,atom.position.y,atom.position.z)))
    end
    reinterpret(LibInt2Engine,ccall(dlsym(libint2jl,:_Z29createEngineNuclearAttractioniiPN7libint24AtomEi),Ptr{Nothing},(Cint,Cint,Ptr{Nothing},Cint),maxNumberPrimitives,maxLQN.exponent,l2atoms,length(l2atoms)))
  end
  function LibInt2EngineOverlap(maxNumberPrimitives::Int,maxLQN::LQuantumNumber)
    reinterpret(LibInt2Engine,ccall(dlsym(libint2jl,:_Z19createEngineOverlapii),Ptr{Nothing},(Cint,Cint),maxNumberPrimitives,maxLQN.exponent))
  end

  #  Destructor
  function destroy!(engine::LibInt2Engine)
    ccall(dlsym(libint2jl,:_Z13destroyEnginePN7libint26EngineE),Nothing,(LibInt2Engine,),engine)
  end


  ## library functionality
  function libInt2Initialize()
    ccall(dlsym(libint2jl,:_Z12libint2startv),Nothing,())
  end

  function libInt2Finalize()
    ccall(dlsym(libint2jl,:_Z11libint2stopv),Nothing,())
  end

  function computeFactorsNonAxial(shells::Vector{LibInt2Shell})
    naxfacs = Vector{Float64}()
    for sh in shells
      push!(naxfacs,computeFactorsNonAxial(LQuantumNumber(lqn(sh)))...)
    end
    return naxfacs
  end

  """
  If a LibInt2Engine is specified, it is used without assertion if it is of the correct type.
  """
  function IntegralsModule.computeMatrixBlockOverlap(engine::LibInt2Engine, μ::LibInt2Shell, ν::LibInt2Shell)
    buf = ccall(dlsym(libint2jl,:_Z13compute2cIntsPN7libint26EngineEPNS_5ShellES3_),Ptr{Cdouble},(LibInt2Engine,LibInt2Shell,LibInt2Shell), engine, μ,ν)
    block = reshape(unsafe_wrap(Array,buf,nbf(μ)*nbf(ν)),(nbf(ν),nbf(μ)))' # the inner loop within LibInt2 is over the second shell
    ShellModule.renormNonAxial!(block,[μ],[ν])
    return block
  end

  function IntegralsModule.computeMatrixBlockOverlap(μlib::LibInt2Shell,νlib::LibInt2Shell)
    totaldim, maxprims, maxlqn = computeDimensions([μlib,νlib])
    engine = LibInt2EngineOverlap(maxprims,maxlqn)

    result = copy(computeMatrixBlockOverlap(engine,μlib,νlib))

    destroy!(engine)
    return result
  end

  """
  If a LibInt2Engine is specified, it is used without assertion if it is of the correct type.
  """
  function computeMatrixBlockKinetic(engine::LibInt2Engine, μ::LibInt2Shell, ν::LibInt2Shell)
    buf = ccall(dlsym(libint2jl,:_Z13compute2cIntsPN7libint26EngineEPNS_5ShellES3_),Ptr{Cdouble},(LibInt2Engine,LibInt2Shell,LibInt2Shell), engine, μ,ν)
    block = reshape(unsafe_wrap(Array,buf,nbf(μ)*nbf(ν)),(nbf(ν),nbf(μ)))' # the inner loop within LibInt2 is over the second shell
    ShellModule.renormNonAxial!(block,[μ],[ν])
    return block
  end

  function computeMatrixBlockKinetic(μlib::LibInt2Shell,νlib::LibInt2Shell)
    totaldim, maxprims, maxlqn = computeDimensions([μlib,νlib])
    engine = LibInt2EngineKinetic(maxprims,maxlqn)

    result = copy(computeMatrixBlockKinetic(engine,μlib,νlib))

    destroy!(engine)
    return result
  end

  """
  If a LibInt2Engine is specified, it is used without assertion if it is of the correct type.
  """
  function computeMatrixBlockNuclearAttraction(engine::LibInt2Engine, μ::LibInt2Shell, ν::LibInt2Shell)
    buf = ccall(dlsym(libint2jl,:_Z13compute2cIntsPN7libint26EngineEPNS_5ShellES3_),Ptr{Cdouble},(LibInt2Engine,LibInt2Shell,LibInt2Shell), engine, μ,ν)
    block = reshape(unsafe_wrap(Array,buf,nbf(μ)*nbf(ν)),(nbf(ν),nbf(μ)))' # the inner loop within LibInt2 is over the second shell
    ShellModule.renormNonAxial!(block,[μ],[ν])
    return block
  end

  function computeMatrixBlockNuclearAttraction(μlib::LibInt2Shell,νlib::LibInt2Shell, atom::Atom)
    totaldim, maxprims, maxlqn = computeDimensions([μlib,νlib])
    engine = LibInt2EngineNuclearAttraction(maxprims,maxlqn,atom)

    result = copy(computeMatrixBlockNuclearAttraction(engine,μlib,νlib))

    destroy!(engine)
    return result
  end

  function computeMatrixBlockNuclearAttraction(μlib::LibInt2Shell,νlib::LibInt2Shell, geo::Geometry)
    totaldim, maxprims, maxlqn = computeDimensions([μlib,νlib])
    engine = LibInt2EngineNuclearAttraction(maxprims,maxlqn,geo)

    result = copy(computeMatrixBlockNuclearAttraction(engine,μlib,νlib))

    destroy!(engine)
    return result
  end

  """
  If a LibInt2Engine is specified, it is used without assertion if it is of the correct type.
  """
  function computeTensorBlockElectronRepulsionIntegrals(engine::LibInt2Engine, μ::LibInt2Shell, ν::LibInt2Shell, λ::LibInt2Shell, σ::LibInt2Shell)
    buf = ccall(dlsym(libint2jl,:_Z13compute4cIntsPN7libint26EngineEPNS_5ShellES3_S3_S3_),Ptr{Cdouble},(LibInt2Engine,LibInt2Shell,LibInt2Shell,LibInt2Shell,LibInt2Shell), engine, σ,λ,ν,μ)
    block = reshape(unsafe_wrap(Array,buf,nbf(μ)*nbf(ν)*nbf(λ)*nbf(σ)),(nbf(μ),nbf(ν),nbf(λ),nbf(σ))) # Apparently for ERIs the order is correct already (in contrast to the matrix functions
    ShellModule.renormNonAxial!(block,[μ],[ν],[λ],[σ])
    return block
  end

  function computeTensorBlockElectronRepulsionIntegrals(μlib::LibInt2Shell, νlib::LibInt2Shell, λlib::LibInt2Shell, σlib::LibInt2Shell)
    totaldim, maxprims, maxlqn = computeDimensions([μlib,νlib,λlib,σlib])
    engine = LibInt2EngineCoulomb(maxprims,maxlqn)

    result = copy(computeTensorBlockElectronRepulsionIntegrals(engine, μlib,νlib,λlib,σlib))

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








else # if the libint2 is not available we have to stub the relevant functions
  include("LibInt2Module.stubs")
end

end # Module





