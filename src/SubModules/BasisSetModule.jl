module BasisSetModule
export BasisSet, readBasisSetTX93, readBasisSetGaussian, readBasisSetGAMESSUS, readBasisSetNWChem, readBasisSetDalton
using ..BaseModule
using ..AtomModule
using StringParserPEG

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


# Supplementary Functions for the PEG BasisSet Parsers

# MergeDictionaries - used to construct one large dictionary from an array of dictionaries created for each element beforehand.
function MergeDictionaries(children) 
dicty = Dict() 

for i = 1:length(children) 
  dicty = merge(dicty,children[i]) 
end 

return dicty
end

# CreateDictionariesOfElement - used in the NWChem Parser in order to create a dictionary of contracted Gaussians for each element. These dictionaries are later merged by the MergeDictionaries function in one large dictionary.
function CreateDictionariesOfElements(children)
Elem = children[1][1]
contracteds = Vector()

for i=1:length(children)
 cont = BasisSetModule.ContractedGaussianBasisFunctionDefinition(children[i][2],children[i][3][1])
 push!(contracteds,cont)
end

dictyElem = Dict(zip([Elem],[contracteds]))

return dictyElem
end

# MergeContractedFunctions - used by the Gaussian Parser to create a vector of all the contracted Gaussians of one element. This is employed due to the Gaussian basis set file defining both S and P orbital functions in one line (SP).
function MergeContractedFunctions(children)
bases = Vector()
for i=1:length(children)
if typeof(children[i]) == BasisSetModule.ContractedGaussianBasisFunctionDefinition
  push!(bases,children[i])
elseif typeof(children[i]) == Array{BasisSetModule.ContractedGaussianBasisFunctionDefinition,1}
  for j=1:length(children[i])
  push!(bases,children[i][j]) 
  end
end
end

return bases
end

#SPOrbitals - Separates S and P orbital functions from the SP syntax used in the Gaussian Basis Set files.
function SPOrbitals(children)
OrbS = Vector()
OrbP = Vector()
S = Vector()
P = Vector()
for i=1:length(children)
  S1 = children[i][1]
  push!(S,S1)
end
for i=1:length(children)
  S2 = children[i][2]
  push!(P,S2)
end

OrbS = BasisSetModule.ContractedGaussianBasisFunctionDefinition(LQuantumNumber("S"),S)
OrbP = BasisSetModule.ContractedGaussianBasisFunctionDefinition(LQuantumNumber("P"),P)
#push!(Primitives,OrbS)
#push!(Primitives,OrbP)


Primitives = vcat(OrbS,OrbP)
return Primitives
end

#=
function MergeDictionariesDalton(children) 
  dicty = Dict() 

for i = 1:length(children)
  dicty = merge(dicty,children[i])
end 

return dicty
end
=#

#MergeContractedsDalton - creates a vector containing all contracted Gaussians of one element.
function MergeContractedsDalton(children)
bases = Vector()
for i=1:length(children)
  for j=1:length(children[i])
  push!(bases,children[i][j]) 
  end
end

return bases
end

# CreateContractedFunctions - this function first creates a vector of vectors that contains the primitives of each contracted function of one element. These are then used to construct the contracted Gaussians.
function CreateContractedFunctions(children)
Primitives = Vector()
A = zeros(length(children[2]),length(children[2][1]))
for i = 1:length(children[2]) 
  for j = 1:length(children[2][i])
  A[i,j] = children[2][i][j]
  end
end 
B = 0

a = size(A)
C = [Vector{Any}() for _ in 1:a[2]-1]

for j = 1:a[2]-1
  for i = 1:a[1]
    if A[i,j] == 0 
    continue
    else 
    B = BasisSetModule.PrimitiveGaussianBasisFunctionDefinition(A[i,end],A[i,j])
    push!(C[j],B)
    end
  end
end 



D = Vector(length(C))
for i=1:length(C)
  D[i]=BasisSetModule.ContractedGaussianBasisFunctionDefinition(children[1],C[i])
end
return D
end


############################################################################################################################################
# Main Functions

# Dalton 

function readBasisSetDalton(filename::AbstractString)
readBasisSetStringDalton = Grammar("""
start => (-(comments) & -(databegin) & bases & ?(end)) {"start"}
comments => *(*(comment) & newline) {"comments"}
comment => -(exclam & (r([ a-zA-Z0-9\-\(\)\*\:\.\/\,\+]*)r) & ?(newline)) 
bases => +(base) {"bases"}
base => (finalelement & -(finalelement2) & sets) {"base"}
convertedelement => (element) {"convertedelement"}
element => (r([a-zA-Z]+)r) {"element"}
sets => +(set) {"sets"}
set => (-(exclam) & -(space) & QuantumNumber & -(space) & -(funcs) & -(newline) & -(nofun) & primitives) {"set"}
primitives => +(primitive) {"primitives"} 
primitive => (expon & -(space) & +( (lin & -(onlyspace)) {"linandzeros"}) {"remaininglinear"} & -(space)) {"primitive"}
QuantumNumber => ((opl)) {"convopl"}
expon => (num) {"expon"}
lin => (num) {"lin"}
opl => ('s' {"S"} |'d' {"D"} | 'f' {"F"} | 'p' {"P"} | 'g' {"g"} | 'h' | 'i' | 'j' | 'k') {"opl"}
end => -('!' & r([^&]+)r) 

finalelement => (-(exclam) & -(space) & convertedelement & -(space) & -(orbitalcount) & -(space)) {"finalelement"}
finalelement2 => (-(exclam) & -(space) & convertedelement & -(space) & -(orbitalcount) & -(space)) {"finalelement2"}

onlyspace => r([ ]+)r
funcs => ('functions')
funnum => (num)
nofun => ('H' & space & num & space & num & space)
databegin => ('! Basis = ' & r([^&\\n]+)r & '\n')
NumberOfPrimitives => (num)
ScalingFactor => (num)
null => ('0')
num => (r([0-9\+\E\-]+)r & '.' & r([0-9\+\E\-]+)r | r([0-9\+\E\-]+)r |minus & r([0-9\+\E\-]+)r & '.' & r([0-9\+\E\-]+)r |minus & r([0-9\+\E\-]+)r) {"num"}
minus => ('-')
space => r([ \\t\\n\\r]*)r
sep => ('#BASIS SET: ') {"sep"}
orbitalcount => (r([a-z0-9\(\)\,]+)r & ' -> ' & r([^A-Z\\n]+)r & '\n') {"orbitalcount"}
exclam => ('!')
newline => ('\n')
doublenewline => -('\n\n')
""")


toresult(node,children,::MatchRule{:default}) =  node.value
toresult(node,children,::MatchRule{:expon}) =  parse(Float64,node.value) 
toresult(node,children,::MatchRule{:lin}) =  parse(Float64,node.value)
toresult(node,children,::MatchRule{:linandzeros}) = children[1]
toresult(node,children,::MatchRule{:remaininglinear}) = children 
toresult(node,children,::MatchRule{:element}) =  string(""", node.value, """)  
toresult(node,children,::MatchRule{:convertedelement}) =  Element(ElementSymbols[node.value]) 
toresult(node,children,::MatchRule{:finalelement}) =  children[1]
toresult(node,children,::MatchRule{:opl}) =  string(""", node.value, """)
toresult(node,children,::MatchRule{:convopl}) = LQuantumNumber(uppercase(children[1]))
toresult(node,children,::MatchRule{:bases})= MergeDictionaries(children)
toresult(node,children,::MatchRule{:base})= Dict(zip([children[1]],[children[2]]))
toresult(node,children,::MatchRule{:sets}) =  MergeContractedsDalton(children)
toresult(node,children,::MatchRule{:set}) = CreateContractedFunctions(children)
toresult(node,children,::MatchRule{:primitives}) = children
toresult(node,children,::MatchRule{:primitive}) = push!(children[2],children[1])
toresult(node,children,::MatchRule{:start}) = BasisSet(children[1])


fileString = read(filename, String)
(ast,pos,error) = parse(readBasisSetStringDalton, fileString)

if error == nothing
  geom = transform(toresult,ast)
  return geom
  else
  return error
end

end

function readBasisSetGAMESSUS(filename::AbstractString)
readBasisSetStringGAMESSUS = Grammar("""
start => (-(comments) & -(databegin) & bases & ?(end)) {"start"}
comments => *(*(comment) & newline) {"comments"}
comment => -(exclam & (r([ a-zA-Z0-9\-\(\)\*\:\.\/\,\+]*)r) & ?(newline)) 
bases => +(base) {"bases"}
base => (convertedelement & -(space) & sets) {"base"}
convertedelement => (element) {"convertedelement"}
element => (r([a-zA-Z]+)r) {"element"}
sets => +(set) {"sets"}
set => (QuantumNumber & -(space) & -(numberofprimitives) & -(space) & primitives) {"set"}
primitives => +(primitive) {"primitives"} 
primitive => (-(funnum) & -(space) & expon & -(space) & lin & -(space)) {"primitive"}
QuantumNumber => ((opl)) {"convopl"}
expon => (num) {"expon"}
lin => (num) {"lin"}
opl => ('S' {"S"} |'D' {"D"} | 'F' {"F"} | 'P' {"P"} | 'G' {"G"} | 'H' | 'I' | 'J' | 'K') {"opl"}
end => -('\$END' & r([^&]+)r) 

funnum => (num)
numberofprimitives => (num)
databegin => ('\$DATA\n')
linS => (num) {"linS"}
linP => (num) {"linP"}
NumberOfPrimitives => (num)
ScalingFactor => (num)
null => ('0')
num => (r([0-9\+\E\-]+)r & '.' & r([0-9\+\E\-]+)r | r([0-9\+\E\-]+)r |minus & r([0-9\+\E\-]+)r & '.' & r([0-9\+\E\-]+)r |minus & r([0-9\+\E\-]+)r) {"num"}
minus => ('-')
space => r([ \\t\\n\\r]*)r
sep => ('#BASIS SET: ') {"sep"}
orbitalcount => (r([a-z0-9\(\)\,]+)r & ' -> ' & r([^A-Z\\n]+)r & '\n') {"orbitalcount"}
exclam => ('!')
newline => ('\n')
doublenewline => -('\n\n')
""")

toresult(node,children,::MatchRule{:default}) =  node.value
toresult(node,children,::MatchRule{:expon}) =  parse(Float64,node.value) 
toresult(node,children,::MatchRule{:lin}) =  parse(Float64,node.value)
toresult(node,children,::MatchRule{:linS}) =  parse(Float64,node.value)
toresult(node,children,::MatchRule{:linP}) =  parse(Float64,node.value)
toresult(node,children,::MatchRule{:element}) =  string(""", node.value, """)  
toresult(node,children,::MatchRule{:convertedelement}) =  Element(ElementSymbols[node.value]) 
toresult(node,children,::MatchRule{:opl}) =  string(""", node.value, """)
toresult(node,children,::MatchRule{:convopl}) = LQuantumNumber(children[1])
toresult(node,children,::MatchRule{:bases})= MergeDictionaries(children)
toresult(node,children,::MatchRule{:base})= Dict(zip([children[1]],[children[2][1]]))
toresult(node,children,::MatchRule{:sets}) = [children]
toresult(node,children,::MatchRule{:set}) = BasisSetModule.ContractedGaussianBasisFunctionDefinition(children[1],children[2][1])
toresult(node,children,::MatchRule{:primitives}) = [children]
toresult(node,children,::MatchRule{:primitive}) = BasisSetModule.PrimitiveGaussianBasisFunctionDefinition(children[1],children[2])
toresult(node,children,::MatchRule{:start}) = BasisSet(children[1])

fileString = read(filename, String)
(ast,pos,error) = parse(readBasisSetStringGAMESSUS, fileString)

if error == nothing
  geom = transform(toresult,ast)
  return geom
  else
  return error
end

end
#TX93

function readBasisSetTX93(filename::AbstractString)

readBasisSetStringTX93 = Grammar("""
start => (-(comments) & bases & ?(end))  {"start"}
comments => *(*(comment) & newline)
comment => -(exclam & (r([ a-zA-Z0-9\-\(\)\*\:\.\/\,\+]*)r) & ?(newline)) 
bases => *(base) {"bases"}
base => (-(for) & -(space) & convertedelement & -(newline) & sets) {"base"}
convertedelement => (element) {"convertedelement"}
element => (r([a-zA-Z]+)r) {"element"}
sets => +(set) {"sets"}
set => (convopl & primitives) {"set"}
primitives => +(primitive) {"primitives"} 
primitive => (-(space) & expon & -(space) & lin & -(space)) {"primitive"}
convopl => (opl) {"convopl"}
expon => (num) {"expon"}
lin => (num) {"lin"}
opl => ('S' {"S"} |'D' {"D"} | 'F' {"F"} | 'P' {"P"} | 'G' {"G"} | 'H' | 'I' | 'J' | 'K') {"opl"}
end => -(-(exclam) & r([^&]+)r) 

num => (r([0-9\+\E\-]+)r & '.' & r([0-9\+\E\-]+)r | r([0-9\+\E\-]+)r |minus & r([0-9\+\E\-]+)r & '.' & r([0-9\+\E\-]+)r |minus & r([0-9\+\E\-]+)r) {"num"}
minus => ('-')
space => r([ \\t\\n\\r]*)r
for => ('FOR')
exclam => ('!')
newline => ('\n')
doublenewline => -('\n\n')
""")

toresult(node,children,::MatchRule{:default}) =  node.value
toresult(node,children,::MatchRule{:expon}) =  parse(Float64,node.value) 
toresult(node,children,::MatchRule{:lin}) =  parse(Float64,node.value)
toresult(node,children,::MatchRule{:element}) =  string(""", node.value, """)  
toresult(node,children,::MatchRule{:convertedelement}) =  Element(node.value) 
toresult(node,children,::MatchRule{:opl}) =  string(""", node.value, """)
toresult(node,children,::MatchRule{:convopl}) = LQuantumNumber(node.value)
toresult(node,children,::MatchRule{:bases})= MergeDictionaries(children)
toresult(node,children,::MatchRule{:base})= Dict(zip([children[1]],[children[2][1]]))
toresult(node,children,::MatchRule{:sets}) = [children]
toresult(node,children,::MatchRule{:set}) = BasisSetModule.ContractedGaussianBasisFunctionDefinition(children[1],children[2][1])
toresult(node,children,::MatchRule{:primitives}) = [children]
toresult(node,children,::MatchRule{:primitive}) = BasisSetModule.PrimitiveGaussianBasisFunctionDefinition(children[1],children[2])
toresult(node,children,::MatchRule{:start}) = BasisSet(children[1])


fileString = read(filename, String)
(ast,pos,error) = parse(readBasisSetStringTX93, fileString)

if error == nothing
  geom = transform(toresult,ast)
  return geom
  else
  return error
end

end


#Gaussian

function readBasisSetGaussian(filename::AbstractString)

readBasisSetStringGaussian94 = Grammar("""
start => (-(comments) & bases & ?(end)) {"start"}
comments => *(*(comment) & newline)
comment => -(exclam & (r([ a-zA-Z0-9\-\(\)\*\:\.\/\,\+]*)r) & ?(newline)) 
bases => *(base) {"bases"}
base => (-(sep) & -(newline) & convertedelement & -(space) & -(null) & -(space) & sets) {"base"}
convertedelement => (element) {"convertedelement"}
element => (r([a-zA-Z]+)r) {"element"}
sets => *(set | SPSet) {"sets"}
set => ((QuantumNumber &  primitives)) {"set"}
primitives => +(primitive) {"primitives"} 
primitive => (-(space) & expon & -(space) & lin & -(space)) {"primitive"}
QuantumNumber => ((opl) & -(space) & -(NumberOfPrimitives) & -(space) & -(ScalingFactor) & -(space)) {"convopl"}
expon => (num) {"expon"}
lin => (num) {"lin"}
opl => ('S' {"S"} |'D' {"D"} | 'F' {"F"} | 'P' {"P"} | 'G' {"G"} | 'H' | 'I' | 'J' | 'K') {"opl"}
end => -(-(sep) & r([^&]+)r) 

SPSet => (-(SP) & SPPrimitives) {"SPSet"}
SP => ((SPOrb) & -(space) & -(NumberOfPrimitives) & -(space) & -(ScalingFactor))
SPOrb => ('SP') {"SPOrb"}
SPPrimitives => *(SPPrimitive) {"SPPrimitives"}
SPPrimitive => (-(space) & expon & -(space) & linS & -(space) & linP & -(space)) {"SPPrimitive"}


linS => (num) {"linS"}
linP => (num) {"linP"}
NumberOfPrimitives => (num)
ScalingFactor => (num)
null => ('0')
num => (r([0-9\+\E\-]+)r & '.' & r([0-9\+\E\-]+)r | r([0-9\+\E\-]+)r |minus & r([0-9\+\E\-]+)r & '.' & r([0-9\+\E\-]+)r |minus & r([0-9\+\E\-]+)r) {"num"}
minus => ('-')
space => r([ \\t\\n\\r]*)r
sep => ('****') {"sep"}
exclam => ('!')
newline => ('\n')
doublenewline => -('\n\n')
""")

toresult(node,children,::MatchRule{:default}) =  node.value
toresult(node,children,::MatchRule{:expon}) =  parse(Float64,node.value) 
toresult(node,children,::MatchRule{:lin}) =  parse(Float64,node.value)
toresult(node,children,::MatchRule{:linS}) =  parse(Float64,node.value)
toresult(node,children,::MatchRule{:linP}) =  parse(Float64,node.value)
toresult(node,children,::MatchRule{:element}) =  string(""", node.value, """)  
toresult(node,children,::MatchRule{:convertedelement}) =  Element(node.value) 
toresult(node,children,::MatchRule{:opl}) =  string(""", node.value, """)
toresult(node,children,::MatchRule{:convopl}) = LQuantumNumber(children[1])
toresult(node,children,::MatchRule{:bases})= MergeDictionaries(children)
toresult(node,children,::MatchRule{:base})= Dict(zip([children[1]],[children[2]]))
toresult(node,children,::MatchRule{:sets}) = MergeContractedFunctions(children)
toresult(node,children,::MatchRule{:set}) = BasisSetModule.ContractedGaussianBasisFunctionDefinition(children[1],children[2][1])
toresult(node,children,::MatchRule{:primitives}) = [children] 
toresult(node,children,::MatchRule{:primitive}) = BasisSetModule.PrimitiveGaussianBasisFunctionDefinition(children[1],children[2])
toresult(node,children,::MatchRule{:SPPrimitive}) = [BasisSetModule.PrimitiveGaussianBasisFunctionDefinition(children[1],children[2]),BasisSetModule.PrimitiveGaussianBasisFunctionDefinition(children[1],children[3])]
toresult(node,children,::MatchRule{:SPPrimitives}) = SPOrbitals(children)
toresult(node,children,::MatchRule{:SPSet}) = children[1]
toresult(node,children,::MatchRule{:start}) = BasisSet(children[1])

fileString = read(filename, String)
(ast,pos,error) = parse(readBasisSetStringGaussian94, fileString)

if error == nothing
  geom = transform(toresult,ast)
  return geom
  else
  return error
end

end

#NWChem

function readBasisSetNWChem(filename::AbstractString)

readBasisSetStringNWChem = Grammar("""
start => (-(comments) & printline & bases & ?(end)) {"start"}
printline => -('BASIS "ao basis" PRINT\n') {"printline"}
comments => *(*(comment) & newline) {"comments"}
comment => -(exclam & (r([ a-zA-Z0-9\-\(\)\*\:\.\/\,\+]*)r) & ?(newline)) 
bases => +(base) {"bases"}
base => (-(sep) & -(orbitalcount) & sets) {"base"}
convertedelement => (element) {"convertedelement"}
element => (r([a-zA-Z]+)r) {"element"}
sets => +(set) {"sets"}
set => (convertedelement & -(space) & QuantumNumber & primitives) {"set"}
primitives => +(primitive) {"primitives"} 
primitive => (-(space) & expon & -(space) & lin & -(space)) {"primitive"}
QuantumNumber => ((opl)) {"convopl"}
expon => (num) {"expon"}
lin => (num) {"lin"}
opl => ('S' {"S"} |'D' {"D"} | 'F' {"F"} | 'P' {"P"} | 'G' {"G"} | 'H' | 'I' | 'J' | 'K') {"opl"}
end => ('END' & r([^&]+)r) 



linS => (num) {"linS"}
linP => (num) {"linP"}
NumberOfPrimitives => (num)
ScalingFactor => (num)
null => ('0')
num => (r([0-9\+\E\-]+)r & '.' & r([0-9\+\E\-]+)r | r([0-9\+\E\-]+)r |minus & r([0-9\+\E\-]+)r & '.' & r([0-9\+\E\-]+)r |minus & r([0-9\+\E\-]+)r) {"num"}
minus => ('-')
space => r([ \\t\\n\\r]*)r
sep => ('#BASIS SET: ') {"sep"}
orbitalcount => (r([a-z0-9\(\)\,]+)r & ' -> ' & r([^A-Z\\n]+)r & '\n') {"orbitalcount"}
exclam => ('#')
newline => ('\n')
doublenewline => -('\n\n')
""")

toresult(node,children,::MatchRule{:default}) =  node.value
toresult(node,children,::MatchRule{:expon}) =  parse(Float64,node.value) 
toresult(node,children,::MatchRule{:lin}) =  parse(Float64,node.value)
toresult(node,children,::MatchRule{:linS}) =  parse(Float64,node.value)
toresult(node,children,::MatchRule{:linP}) =  parse(Float64,node.value)
toresult(node,children,::MatchRule{:element}) =  string(""", node.value, """)  
toresult(node,children,::MatchRule{:convertedelement}) =  Element(node.value) 
toresult(node,children,::MatchRule{:opl}) =  string(""", node.value, """)
toresult(node,children,::MatchRule{:convopl}) = LQuantumNumber(children[1])
toresult(node,children,::MatchRule{:bases})= MergeDictionaries(children)
toresult(node,children,::MatchRule{:base})= children[1]
toresult(node,children,::MatchRule{:sets}) = CreateDictionariesOfElements(children)
toresult(node,children,::MatchRule{:set}) = children
toresult(node,children,::MatchRule{:primitives}) = [children] 
toresult(node,children,::MatchRule{:primitive}) = BasisSetModule.PrimitiveGaussianBasisFunctionDefinition(children[1],children[2])
toresult(node,children,::MatchRule{:start}) = BasisSet(children[1])

fileString = read(filename, String)
(ast,pos,error) = parse(readBasisSetStringNWChem, fileString)

if error == nothing
  geom = transform(toresult,ast)
  return geom
  else
  return error
end

end


end # module
