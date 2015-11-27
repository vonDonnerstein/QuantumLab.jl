module GeometryModule
export Geometry, readGeometryXYZ, angstrom2bohr, angstrom2bohr!

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

angstrom2bohr(x::Float64) = (return x/0.529177210)
angstrom2bohr(pos::Position) = (return Position(angstrom2bohr(pos.x),angstrom2bohr(pos.y),angstrom2bohr(pos.z)))
function angstrom2bohr!(geo::Geometry)
  for (atom in geo.atoms)
    atom.position = angstrom2bohr(atom.position)
  end
end

end #module
