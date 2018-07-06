module BasisSetModule
export BasisSet, readTX93, readGaussian, readGAMESSUS, readNWChem, readDalton
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

# Element Dictionary : Full name to Symbol
 ElementSymbols = Dict{String,String}(
"HYDROGEN" =>                 "H",
"HELIUM" =>                  "He",
"LITHIUM" =>                 "Li",
"BERYLLIUM" =>               "Be",
"BORON" =>                    "B",
"CARBON" =>                   "C",
"NITROGEN" =>                 "N",
"OXYGEN" =>                   "O",
"FLUORINE" =>                 "F",
"NEON" =>                    "Ne",
"SODIUM" =>                  "Na",
"MAGNESIUM" =>               "Mg", 
"ALUMINUM" =>                "Al",
"ALUMINIUM" =>         "Al", 
"SILICON" =>                 "Si", 
"PHOSPHOROUS" =>               "P", 
"SULFUR" =>                   "S", 
"CHLORINE" =>                "Cl", 
"ARGON" =>                   "Ar", 
"POTASSIUM" =>                "K", 
"CALCIUM" =>                 "Ca", 
"SCANDIUM" =>                "Sc", 
"TITANIUM" =>                "Ti", 
"VANADIUM" =>                 "V", 
"CHROMIUM" =>                "Cr", 
"MANGANESE" =>               "Mn", 
"IRON" =>                    "Fe", 
"COBALT" =>                  "Co", 
"NICKEL" =>                  "Ni", 
"COPPER" =>                  "Cu", 
"ZINC" =>                    "Zn", 
"GALLIUM" =>                 "Ga", 
"GERMANIUM" =>               "Ge", 
"ARSENIC" =>                 "As", 
"SELENIUM" =>                "Se", 
"BROMINE" =>                 "Br", 
"KRYPTON" =>                 "Kr", 
"RUBIDIUM" =>                "Rb", 
"STRONTIUM" =>               "Sr", 
"YTTRIUM" =>                  "Y", 
"ZIRCONIUM" =>               "Zr", 
"NIOBIUM" =>                 "Nb", 
"MOLYBDENUM" =>              "Mo", 
"TECHNETIUM" =>              "Tc", 
"RUTHENIUM" =>               "Ru", 
"RHODIUM" =>                 "Rh", 
"PALLADIUM" =>               "Pd", 
"SILVER" =>                  "Ag", 
"CADMIUM" =>                 "Cd", 
"INDIUM" =>                  "In", 
"TIN" =>                     "Sn", 
"ANTIMONY" =>                "Sb", 
"TELLURIUM" =>               "Te", 
"IODINE" =>                   "I",
"XENON" =>                   "Xe", 
"CESIUM" =>                  "Cs", 
"BARIUM" =>                  "Ba", 
"LANTHANUM" =>               "La", 
"CERIUM" =>                  "Ce", 
"PRASEODYMIUM" =>            "Pr", 
"NEODYMIUM" =>               "Nd", 
"PROMETHIUM" =>              "Pm", 
"SAMARIUM" =>                "Sm", 
"EUROPIUM" =>                "Eu", 
"GADOLINIUM" =>              "Gd", 
"TERBIUM" =>                 "Tb", 
"DYSPROSIUM" =>              "Dy", 
"HOLMIUM" =>                 "Ho", 
"ERBIUM" =>                  "Er", 
"THULIUM" =>                 "Tm", 
"YTTERBIUM" =>               "Yb", 
"LUTETIUM" =>                "Lu", 
"HAFNIUM" =>                 "Hf", 
"TANTALUM" =>                "Ta", 
"TUNGSTEN" =>                 "W", 
"RHENIUM" =>                 "Re", 
"OSMIUM" =>                  "Os", 
"IRIDIUM" =>                 "Ir", 
"PLATINUM" =>                "Pt", 
"GOLD" =>                    "Au", 
"MERCURY" =>                 "Hg", 
"THALLIUM" =>                "Tl", 
"LEAD" =>                    "Pb", 
"BISMUTH" =>                 "Bi", 
"POLONIUM" =>                "Po", 
"ASTATINE" =>                "At", 
"RADON" =>                   "Rn", 
"FRANCIUM" =>                "Fr", 
"RADIUM" =>                  "Ra", 
"ACTINIUM" =>                "Ac", 
"THORIUM" =>                 "Th", 
"PROTACTINIUM" =>            "Pa", 
"URANIUM" =>                  "U",
"NEPTUNIUM" =>               "Np", 
"PLUTONIUM" =>               "Pu", 
"AMERICIUM" =>               "Am", 
"CURIUM" =>                  "Cm", 
"BERKELIUM" =>               "Bk", 
"CALIFORNIUM" =>             "Cf", 
"EINSTEINIUM" =>             "Es", 
"FERMIUM" =>                 "Fm", 
"MENDELEVIUM" =>             "Md", 
"NOBELIUM" =>                "No", 
"LAWRENCIUM" =>              "Lr", 
"RUTHERFORDIUM" =>           "Rf", 
"DUBNIUM" =>                 "Db", 
"SEABORGIUM" =>              "Sg", 
"BOHRIUM" =>                 "Bh", 
"HASSIUM" =>                 "Hs", 
"MEITNERIUM" =>              "Mt", 
"DARMSTADTIUM" =>            "Ds", 
"ROENTGENIUM" =>             "Rg", 
"COPERNICIUM" =>             "Cn", 
"NIHONIUM" =>                "Nh", 
"FLEROVIUM" =>               "Fl", 
"MOSCOVIUM" =>               "Mc", 
"LIVERMORIUM" =>             "Lv", 
"TENNESSINE" =>              "Ts", 
"OGANESSON" =>               "Og", 
)

# Supplementary Functions for the PEG BasisSet Parsers


function dict(children) 
  dicty = Dict() 
for i = 1:length(children) 
  dicty = merge(dicty,children[i]) 
end 

return dicty
end

function DictElements(children)
Elem = children[1][1]
contracteds = Vector()
for i=1:length(children)
 cont = BasisSetModule.ContractedGaussianBasisFunctionDefinition(children[i][2],children[i][3][1])
 push!(contracteds,cont)
end

dictyElem = Dict(zip([Elem],[contracteds]))

return dictyElem
end

function bases(children)
bases = Vector()
for i=1:length(children)
if typeof(children[i]) == QuantumLab.BasisSetModule.ContractedGaussianBasisFunctionDefinition
  push!(bases,children[i])
elseif typeof(children[i]) == Array{QuantumLab.BasisSetModule.ContractedGaussianBasisFunctionDefinition,1}
  for j=1:length(children[i])
  push!(bases,children[i][j]) 
  end
end
end

return bases
end

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

function dictDalton(children) 
  dicty = Dict() 
for i = 1:length(children)
  dicty = merge(dicty,children[i])
end 

return dicty
end


function basesDalton(children)
bases = Vector()
for i=1:length(children)
  for j=1:length(children[i])
  push!(bases,children[i][j]) 
  end
end

return bases
end


function CreateCorrectPrimitives(children)
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

function LinearFactors(children)

LinearFactor=0
for i=1:length(children)
  if children[i][1] == 0
  continue
  else 
  LinearFactor=children[i][1]
  end
end
return LinearFactor
end
############################################################################################################################################
# Main Functions

# Dalton 

function readDalton(filename::AbstractString)
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
opl => ('s' {"S"} |'d' {"D"} | 'f' {"F"} | 'p' {"P"} | 'g' {"g"}) {"opl"}
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
toresult(node,children,::MatchRule{:bases})= dict(children)
toresult(node,children,::MatchRule{:base})= Dict(zip([children[1]],[children[2]]))
toresult(node,children,::MatchRule{:sets}) =  basesDalton(children)
toresult(node,children,::MatchRule{:set}) = CreateCorrectPrimitives(children)
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

function readGAMESSUS(filename::AbstractString)
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
opl => ('S' {"S"} |'D' {"D"} | 'F' {"F"} | 'P' {"P"} | 'G' {"G"}) {"opl"}
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
toresult(node,children,::MatchRule{:bases})= dict(children)
toresult(node,children,::MatchRule{:base})= Dict(zip([children[1]],[children[2][1]]))
toresult(node,children,::MatchRule{:sets}) = [children]
toresult(node,children,::MatchRule{:set}) = BasisSetModule.ContractedGaussianBasisFunctionDefinition(children[1],children[2][1])
toresult(node,children,::MatchRule{:primitives}) = [children]
toresult(node,children,::MatchRule{:primitive}) = BasisSetModule.PrimitiveGaussianBasisFunctionDefinition(children[1],children[2])
#Dict(zip([children[1]],[BasisSetModule.ContractedGaussianBasisFunctionDefinition(children[2],children[3][1])]))



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

function readTX93(filename::AbstractString)

readBasisSetStringTX93 = Grammar("""
start => (-(comments) & bases) {"start"}
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
opl => ('S' {"S"} | 'P' {"P"} | 'D' {"D"} | 'F' {"F"}) {"opl"}



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
toresult(node,children,::MatchRule{:bases})= dict(children)
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

function readGaussian(filename::AbstractString)

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
opl => ('S' {"S"} |'D' {"D"} | 'F' {"F"} | 'P' {"P"}) {"opl"}
end => (-(sep) & -(space))

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
toresult(node,children,::MatchRule{:bases})= dict(children)
toresult(node,children,::MatchRule{:base})= Dict(zip([children[1]],[children[2]]))
toresult(node,children,::MatchRule{:sets}) = bases(children)
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

function readNWChem(filename::AbstractString)

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
opl => ('S' {"S"} |'D' {"D"} | 'F' {"F"} | 'P' {"P"} | 'G' {"G"}) {"opl"}
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
toresult(node,children,::MatchRule{:bases})= dict(children)
toresult(node,children,::MatchRule{:base})= children[1]
toresult(node,children,::MatchRule{:sets}) = DictElements(children)
toresult(node,children,::MatchRule{:set}) = children
toresult(node,children,::MatchRule{:primitives}) = [children] 
toresult(node,children,::MatchRule{:primitive}) = BasisSetModule.PrimitiveGaussianBasisFunctionDefinition(children[1],children[2])
#Dict(zip([children[1]],[BasisSetModule.ContractedGaussianBasisFunctionDefinition(children[2],children[3][1])]))



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
