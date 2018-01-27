module BasisSetModule
export BasisSet, readBasisSetTX93
using ..BaseModule
using ..AtomModule

immutable PrimitiveGaussianBasisFunctionDefinition
    exponent::Float64
    prefactor::Float64
end

immutable ContractedGaussianBasisFunctionDefinition
    lQuantumNumber::LQuantumNumber
    primitives::Array{PrimitiveGaussianBasisFunctionDefinition,1}
end

immutable BasisSet
    definitions::Dict{Element,Array{ContractedGaussianBasisFunctionDefinition,1}}
end


function readBasisSetTX93(filename::AbstractString)
  basSet = BasisSet(Dict())                                 # initialize return value

  elem = Element("H")                                       # functionglobal elem
  lqn = "S"                                         # functionglobal lqn
  contractedDefinition = ContractedGaussianBasisFunctionDefinition(LQuantumNumber("S"),[])  # functionglobal contractedDefinition

  fd = open(filename)
  for line in eachline(fd)
    # comments
    if ismatch(r"^(!|\s*$)",line) continue end              # skip

    # new element
    if ismatch(r"^FOR",line)
      if(contractedDefinition.primitives != [])                 #   and not the first
    append!(basSet.definitions[elem],[contractedDefinition])                    #   finish contractedDefinition
    contractedDefinition = ContractedGaussianBasisFunctionDefinition(LQuantumNumber(lqn),[])    #   start new contractedDefinition
      end
      elem = Element(match(r"^FOR\s*(\w*)",line).captures[1])           #   update elem
      basSet.definitions[elem] = ContractedGaussianBasisFunctionDefinition[]    #   initialize elem-definitions empty

    # definition line
    else
      regmatch = match(Regex(raw"(?P<lqn>\w*)?\s*(?P<exp>"*floatregex*raw")\s*(?P<lin>"*floatregex*raw")"),line)
      # new contracted
      if regmatch[:lqn] != ""                   # if new contracted
    lqn = regmatch[:lqn]
    if(contractedDefinition.primitives != [])               #    not a new element
      append!(basSet.definitions[elem],[contractedDefinition])                  #   finish contractedDefinition
    end
      contractedDefinition = ContractedGaussianBasisFunctionDefinition(LQuantumNumber(lqn),[])  #   start new contractedDefinition
      end
      # always (add primitive to current contracted)
      exponent = float(regmatch[:exp])
      linfactor = float(regmatch[:lin])
      append!(contractedDefinition.primitives,[PrimitiveGaussianBasisFunctionDefinition(exponent,linfactor)])
    end
  end
  if(contractedDefinition.primitives != [])                 #   if still contracteds left
    append!(basSet.definitions[elem],[contractedDefinition])        #   finish contractedDefinition
  end
  close(fd)
  return basSet
end

end # module
