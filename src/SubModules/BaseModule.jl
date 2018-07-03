module BaseModule
export Position, Ionization, LQuantumNumber, MQuantumNumber, MQuantumNumbers, distance, origin, floatregex, @T_str, doublefactorial, 𝐈, evaluateFunction, trlog
import Base.*, Base.+, Base./, Base.-, Base.isless, Base.convert

immutable Position
	x::Float64
	y::Float64
	z::Float64
end

const origin = Position(0.,0.,0.)

"""
    floatregex
is the regex-string that matches any floating point number with or without scientific notation

use it by combining it like
```
Regex(raw"5\\+" * floatregex)
```
to match e.g. '5+3.7'
"""
const floatregex = raw"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?" 

*(c::Real,p::Position) = Position(c*p.x,c*p.y,c*p.z)
*(p::Position,c::Real) = Position(c*p.x,c*p.y,c*p.z)
/(p::Position,c::Real) = Position(p.x/c,p.y/c,p.z/c)
+(p::Position,q::Position) = Position(p.x+q.x,p.y+q.y,p.z+q.z)
-(p::Position,q::Position) = Position(p.x-q.x,p.y-q.y,p.z-q.z)

convert(::Type{Vector{Float64}},p::Position) = [p.x,p.y,p.z]

function distance(p::Position,q::Position)
  pq = p - q
  return sqrt((pq.x)^2+(pq.y)^2+(pq.z)^2)
end

immutable Ionization
  charge::Int64
end

"""
MQuantumNumber describes a single orientation of an atomic orbital in cartesian space by denoting the exponent of x, y and z of the radial part
E.g.:   a d_xy orbital can be described by MQuantumNumber(1,1,0)
"""
immutable MQuantumNumber
  x::Int
  y::Int
  z::Int
end

lqnExponent = Dict{String,Int}(
  "S" => 	0,
  "P" => 	1,
  "D" => 	2,
  "F" =>	3,
  "G" =>	4,
  "H" =>	5,
  "I" =>	6,
  "K" =>	7,
)

exponentLqn = Dict{Int,String}(
  0 => 	"S",
  1 => 	"P",
  2 => 	"D",
  3 =>	"F",
  4 =>	"G",
  5 =>	"H",
  6 =>	"I",
  7 =>	"K",
)


"""
The LQuantumNumber desribes whether a function is spherical (s, lqn=0), polarized (p, lqn=1), etc.
As is generally known the total MQuantumNumber in pure coordinates lies between -lqn and +lqn.
This means that the sum of the three cartesian components can at most b. LQuantumNumber.
"""
immutable LQuantumNumber
	symbol::String
	exponent::Int

	LQuantumNumber(sym::AbstractString) = new(sym,lqnExponent[sym])
	LQuantumNumber(exp::Int) = new(exponentLqn[exp],exp)
end

"""
    numberMQNsCartesian(l::LQuantumNumber)
returns the number of (cartesian, i.e. 'xxy' etc.) `MQuantumNumber`s with the given `LQuantumNumber`. This is equivalent to
`length(MQuantumNumbers(lqn))`.
"""
function numberMQNsCartesian(lqn::LQuantumNumber)
  l = lqn.exponent
  div((l+1)^2+(l+1),2)
end

"""
    numberMQNsPure(l::LQuantumNumber)
returns the number of pure m-Quantum numbers with the given `LQuantumNumber`: -l <= m <= l.
"""
function numberMQNsPure(lqn::LQuantumNumber)
  2*lqn.exponent + 1
end

"""
MQuantumNumber**s** returns all MQuantumNumber objects corresponding to the LQuantumNumber object given to the constructor as an iterable collection.
Order is for example xx,xy,xz,yy,yz,zz for LQuantumNumber = 2
"""
immutable MQuantumNumbers
  lqn::LQuantumNumber
  mqnarray::Array{MQuantumNumber,1}
  count::Int

  function MQuantumNumbers(lqn::LQuantumNumber)
    maxmqn = lqn.exponent
    mqnarray=MQuantumNumber[]
		for x in maxmqn:-1:0
			for y in (maxmqn-x):-1:0
				z = maxmqn-x-y
				append!(mqnarray,[MQuantumNumber(x,y,z)])
			end
		end
		count = ((maxmqn+1)^2+(maxmqn+1))/2
		new(lqn,mqnarray,count)
	end
end
Base.start(mqns::MQuantumNumbers) = 1
Base.next(mqns::MQuantumNumbers,state) = (mqns.mqnarray[state],state+1)
Base.done(mqns::MQuantumNumbers,state) = state > mqns.count
Base.length(mqns::MQuantumNumbers) = mqns.count
#Base.eltype(::Type{MQuantumNumbers}) = Int

isless(lqn1::LQuantumNumber, lqn2::LQuantumNumber) = lqn1.exponent<lqn2.exponent

function evaluateFunction(x::Position, f::Function)
  return f(x)
end

doublefactorial(n::Int) = prod(n:-2:1)

type 𝐈type
  scalar::Float64
end
const 𝐈 = 𝐈type(1.)
+(𝐈t::𝐈type,m::Matrix) = 𝐈t.scalar*eye(m) + m
+(m::Matrix,𝐈t::𝐈type) = 𝐈t.scalar*eye(m) + m
*(𝐈t::𝐈type,scalar::Real) = 𝐈type(𝐈t.scalar*scalar)
*(scalar::Real,𝐈t::𝐈type) = 𝐈type(𝐈t.scalar*scalar)

"""
    trLog(m::Matrix)
trace(logm(m::Matrix)) efficiently via Cholesky-decomposition.
"""
function trlog(m::Matrix)
  L = chol(m)
  result = 0.
  for idx in 1:size(L)[1]
    result += 2*log(L[idx,idx])
  end
  return result
end


end # module
