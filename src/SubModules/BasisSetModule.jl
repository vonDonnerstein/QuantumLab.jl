module BasisSetModule
using StringParserPEG
export BasisSet, readBasisSetTX93
using ..BaseModule
using ..AtomModule
import Base.show
import Base.==
import Base.≈
using Printf

struct PrimitiveGaussianBasisFunctionDefinition
    exponent::Float64
    prefactor::Float64
end
show(io::IO,pgbf::PrimitiveGaussianBasisFunctionDefinition) = print(io,@sprintf("lin:%9.5f exp:%12.5f",pgbf.prefactor, pgbf.exponent))
≈(pgbfdef1::PrimitiveGaussianBasisFunctionDefinition,pgbfdef2::PrimitiveGaussianBasisFunctionDefinition;kwargs...) = (≈(pgbfdef1.exponent,pgbfdef2.exponent;kwargs...) && ≈(pgbfdef1.prefactor,pgbfdef2.prefactor;kwargs...))
≈(vec1::Vector{PrimitiveGaussianBasisFunctionDefinition},vec2::Vector{PrimitiveGaussianBasisFunctionDefinition};kwargs...) = all(.≈(vec1,vec2;kwargs...))

struct ContractedGaussianBasisFunctionDefinition
    lQuantumNumber::LQuantumNumber
    primitives::Array{PrimitiveGaussianBasisFunctionDefinition,1}
end
function show(io::IO,cgbf::ContractedGaussianBasisFunctionDefinition)
  str  = cgbf.lQuantumNumber.symbol*"─ $(cgbf.primitives[1])\n"
  for prim in cgbf.primitives[2:end-1]
    str *= "│  $(prim)\n"
  end
  str *= "╰─ $(cgbf.primitives[end])\n"
  print(io,str)
end
function show(io::IO,vec::Vector{ContractedGaussianBasisFunctionDefinition})
  for cgbf in vec
    show(io,cgbf)
  end
end
==(cgbfdef1::ContractedGaussianBasisFunctionDefinition,cgbfdef2::ContractedGaussianBasisFunctionDefinition) = (cgbfdef1.lQuantumNumber==cgbfdef2.lQuantumNumber && cgbfdef1.primitives == cgbfdef2.primitives)
==(vec1::Vector{ContractedGaussianBasisFunctionDefinition},vec2::Vector{ContractedGaussianBasisFunctionDefinition}) = all(vec1 .== vec2)
≈(cgbfdef1::ContractedGaussianBasisFunctionDefinition,cgbfdef2::ContractedGaussianBasisFunctionDefinition;kwargs...) = (cgbfdef1.lQuantumNumber==cgbfdef2.lQuantumNumber && all(.≈(cgbfdef1.primitives, cgbfdef2.primitives; kwargs...)))
≈(vec1::Vector{ContractedGaussianBasisFunctionDefinition},vec2::Vector{ContractedGaussianBasisFunctionDefinition};kwargs...) = all(.≈(vec1, vec2;kwargs...))

struct BasisSet
    definitions::Dict{Element,Array{ContractedGaussianBasisFunctionDefinition,1}}
end
function show(io::IO,basset::BasisSet)
  for (key,val) in basset.definitions
    println(io,"$key => ")
    show(val)
  end
end
≈(bassset1::BasisSet,basset2::BasisSet;kwargs...) = keys(basset1.definitions)==keys(basset2.definitions) && all(≈(basset1.defintions[el],basset2.definitons[el];kwargs...) for el in keys(basset1.definitions))

tx93 = Grammar("""
	       start   => (basset & -(*(comment))) {liftchild}
	       basset  => *((-(*(comment)) & elementgroup) {liftchild}) {"basset"}
	       comment => ?(r(^!.*)r | space) & eol

	       elementgroup  => (elementheader & CGBFs) {"entry"}
	       elementheader => (-('FOR' & space) & element & -(?(space) & eol)) {liftchild}
	       element       => r([A-Za-z]*)r {"element"}
	       CGBFs	       => +(CGBF) {"cgbfs"}

	       lqn => ('SPD' | 'SP' | 'L' | 'S' | 'P' | 'D' | 'F' | 'G' | 'H') {nodevalue}
	       CGBF    => (lqn & +(PGBF)) {"cgbf"}
	       PGBF    => (-(space) & float & +(-(space) & float) & -(?(space) & eol)) {"pgbf"}
	       """,standardrules)

fromTX93(node,children,::MatchRule{:pgbf}) = [children[1],[n.children[1] for n in children[2].children]]
function QLcgbf(lins,exps,lqn,index=1)
  ContractedGaussianBasisFunctionDefinition(
     LQuantumNumber(lqn),
     [PrimitiveGaussianBasisFunctionDefinition(exp,pre[index]) for (pre,exp) in zip(lins,exps)]
  )
end
function fromTX93(node,children,::MatchRule{:cgbf})
  lqn  = children[1]
  exps = [p[1] for p in children[2].children]
  lins = [p[2] for p in children[2].children]
  if lqn == "SPD"
    @assert length(lins[1]) == 3
    sshell = QLcgbf(lins,exps,"S",1)
    pshell = QLcgbf(lins,exps,"P",2)
    dshell = QLcgbf(lins,exps,"D",3)
    return [sshell,pshell,dshell]
  elseif lqn == "SP" || lqn == "L"
    @assert length(lins[1]) == 2
    sshell = QLcgbf(lins,exps,"S",1)
    pshell = QLcgbf(lins,exps,"P",2)
    return [sshell,pshell]
  elseif lqn in ["S","P","D","F","G","H"]
    @assert length(lins[1]) == 1
    shell = QLcgbf(lins,exps,lqn,1)
    return [shell]
  else
    @error("Unexpected sharedsp lqn ($lqn)")
  end
end
fromTX93(node,children,::MatchRule{:cgbfs}) = vcat(children...)
fromTX93(node,children,::MatchRule{:element}) = Element(node.value)
fromTX93(node,children,::MatchRule{:entry}) =(children[1] => children[2])
fromTX93(node,children,::MatchRule{:basset}) = BasisSet(Dict(children...))

function BasisSetTX93(str::AbstractString)
  (ast,pos,err) = parse(tx93,str)
  if err != nothing
    @error err
  end
  return transform(fromTX93,ast)
end

function readBasisSetTX93(filename::AbstractString)
  fd = open(filename)
  basSet = BasisSetTX93(read(fd,String))
  close(fd)
  return basSet
end

end # module
