module GeometryModule
export Geometry, readGeometryXYZ, angstrom2bohr, angstrom2bohr!, computeEnergyInteratomicRepulsion

using ..BaseModule
using ..AtomModule

type Geometry
	atoms::Array{Atom,1}
end

function readGeometryXYZ(filename::AbstractString; unit::Symbol=:Bohr)
	geo = Geometry([])
	lines = open(readlines,filename)
	for line in lines[3:end]	# line 1 is n. of atoms, line 2 is comment
		columns = split(line)
		elem = Element(columns[1])
		pos = Position(map(parse,columns[2:4])...)	# float is not recommended for use on tuples
		append!(geo.atoms,[Atom(elem,pos)])
	end
	if (unit == :Bohr) angstrom2bohr!(geo) end	# Should be default as Bohr is an atomic unit
	return geo
end

angstrom2bohr(x::Real) = (return x/0.529177210)
angstrom2bohr(pos::Position) = (return Position(angstrom2bohr(pos.x),angstrom2bohr(pos.y),angstrom2bohr(pos.z)))
function angstrom2bohr!(geo::Geometry)
  for (atom in geo.atoms)
    atom.position = angstrom2bohr(atom.position)
  end
end

function computeEnergyInteratomicRepulsion(
  atom1::Atom,
  atom2::Atom)
  q1 = atom1.element.atomicNumber
  q2 = atom2.element.atomicNumber
  A = atom1.position
  B = atom2.position
  return q1*q2/(distance(A,B))
end

function computeEnergyInteratomicRepulsion(geo::Geometry)
  result = 0.::Float64
  for (A in 1:length(geo.atoms))
    for (B in 1:A-1)
      result += computeEnergyInteratomicRepulsion(geo.atoms[A],geo.atoms[B])
    end
  end
  return result
end

end #module
