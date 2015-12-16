module BaseModule
export Position, LQuantumNumber, MQuantumNumber, MQuantumNumbers, distance
import Base.*, Base.+, Base./, Base.-

immutable Position
	x::Float64
	y::Float64
	z::Float64
end

*(c::Real,p::Position) = Position(c*p.x,c*p.y,c*p.z)
*(p::Position,c::Real) = Position(c*p.x,c*p.y,c*p.z)
/(p::Position,c::Real) = Position(p.x/c,p.y/c,p.z/c)
+(p::Position,q::Position) = Position(p.x+q.x,p.y+q.y,p.z+q.z)
-(p::Position,q::Position) = Position(p.x-q.x,p.y-q.y,p.z-q.z)

function distance(p::Position,q::Position) 
  pq = p - q
  return sqrt((pq.x)^2+(pq.y)^2+(pq.z)^2)
end

immutable MQuantumNumber
  x::Int
  y::Int
  z::Int
end

lqnExponent = Dict{ASCIIString,Int}(
  "S" => 	0,
  "P" => 	1,
  "D" => 	2,
  "F" =>	3,
  "G" =>	4,
  "H" =>	5,
)
immutable LQuantumNumber
	symbol::ASCIIString
	exponent::Int

	LQuantumNumber(sym::AbstractString) = new(sym,lqnExponent[sym])
end

immutable MQuantumNumbers
  lqn::LQuantumNumber
  mqnarray::Array{MQuantumNumber,1}
  count::Int

  function MQuantumNumbers(lqn::LQuantumNumber)
    maxmqn = lqn.exponent
    mqnarray=MQuantumNumber[]
    for z in maxmqn:-1:0
      for y in (maxmqn-z):-1:0
	x = maxmqn-z-y
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
Base.length(mqns::MQuantumNumbers) = MQuantumNumbers.count
#Base.eltype(::Type{MQuantumNumbers}) = Int

end # module
