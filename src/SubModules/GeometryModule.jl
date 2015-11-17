module GeometryModule
export Geometry, readGeometryXYZ

using ..BaseModule
using ..AtomModule

type Geometry
	atoms::Array{Atom,1}
end

function readGeometryXYZ(filename::String)
	geo = Geometry([])
	lines = open(readlines,filename)
	for line in lines[3:end]	# line 1 is n. of atoms, line 2 is comment
		columns = split(line)
		elem = Element(columns[1])
		pos = Position(map(float64,tuple(columns[2:4]...))...)	# float is not recommended for use on tuples
		append!(geo.atoms,[Atom(elem,pos)])
	end
	return geo
end


end #module
