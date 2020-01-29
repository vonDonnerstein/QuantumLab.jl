module GeometryModule
export Geometry, readGeometryXYZ, bohr2angstrom, angstrom2bohr, angstrom2bohr!, computeEnergyInteratomicRepulsion, computeNumberElectrons, computeRotationalConstants, computeMatrixInertiaTensor

using ..BaseModule
using ..AtomModule
using LinearAlgebra

using Printf
using Markdown
import Base.show

struct Geometry
	atoms::Array{Atom,1}
end

function show(io::IO,mime::MIME"text/plain",geo::Geometry)
  show(io,mime,Markdown.parse("    \t      x                  y                  z")); print(io,"\n")
  for at in geo.atoms
    show(io,mime,Markdown.parse("**$(at.element.symbol)**")) 
    print(io,@sprintf("\t% 18.10f % 18.10f % 18.10f\n",at.position.x,at.position.y,at.position.z))
  end
end

"""
    readGeometryXYZ(filename::AbstractString; unit::Symbol=:Bohr)

reads an .xyz-Geometry file (acc. to the unofficial standard the coordinates 
of which are interpreted as Angstrom) and yields a Geometry object in `unit`.
"""
function readGeometryXYZ(filename::AbstractString; unit::Symbol=:Bohr)
	geo = Geometry([])
	lines = open(readlines,filename)
	for line in lines[3:end]	# line 1 is n. of atoms, line 2 is comment
		columns = split(line)
		elem = Element(columns[1])
		pos = Position(map(Meta.parse,columns[2:4])...)	# float is not recommended for use on tuples
		append!(geo.atoms,[Atom(elem,pos)])
	end
	if (unit == :Bohr) angstrom2bohr!(geo) end	# Should be default as Bohr is an atomic unit
	return geo
end

function Geometry(filename::AbstractString)
  if isfile(filename)
    return readGeometryXYZ(filename)
  end
  for ext in ("xyz", "XYZ")
    if isfile("$filename.$ext")
      return readGeometryXYZ("$filename.$ext")
    end
  end
  error("Couldn't find file: $filename")
end

angstrom2bohr(x::Real) = x/0.529177210
bohr2angstrom(x::Real) = x*0.529177210
angstrom2bohr(pos::Position) = (return Position(angstrom2bohr(pos.x),angstrom2bohr(pos.y),angstrom2bohr(pos.z)))
function angstrom2bohr!(geo::Geometry)
  for atom in geo.atoms
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
  result = 0.
  for A in 1:length(geo.atoms), B in 1:A-1
      result += computeEnergyInteratomicRepulsion(geo.atoms[A],geo.atoms[B])
  end
  return result
end

function computeNumberElectrons(geo::Geometry,charge::Int=0)
  sum(at.element.atomicNumber for at in geo.atoms) - charge
end

function computePositionCenterOfMass(geo)
  weight = AtomModule.atomicWeight
  sum(weight[at.element.atomicNumber] * at.position for at in geo.atoms)/sum(weight[at.element.atomicNumber] for at in geo.atoms)
end

function computeMatrixInertiaTensor(geo)
  weight = AtomModule.atomicWeight

  com = computePositionCenterOfMass(geo)
  X,Y,Z = com.x,com.y,com.z

  Ixx = sum([weight[at.element.atomicNumber] * ((at.position.y-Y)^2 + (at.position.z-Z)^2) for at in geo.atoms])
  Iyy = sum([weight[at.element.atomicNumber] * ((at.position.x-X)^2 + (at.position.z-Z)^2) for at in geo.atoms])
  Izz = sum([weight[at.element.atomicNumber] * ((at.position.x-X)^2 + (at.position.y-Y)^2) for at in geo.atoms])

  Ixy = Iyx = sum([weight[at.element.atomicNumber] * -((at.position.x-X)*(at.position.y-Y)) for at in geo.atoms])
  Ixz = Izx = sum([weight[at.element.atomicNumber] * -((at.position.x-X)*(at.position.z-Z)) for at in geo.atoms])
  Iyz = Izy = sum([weight[at.element.atomicNumber] * -((at.position.y-Y)*(at.position.z-Z)) for at in geo.atoms])

  return [Ixx Ixy Ixz; Iyx Iyy Iyz; Izx Izy Izz]
end

function computeRotationalConstants(geo;unit=:MHz)
  I = computeMatrixInertiaTensor(geo) * 1822.89 # amu Bohr^2 * mₑ/amu (kg/mol==amu)
  B = 1 ./ (2*eigvals(I)) # 1/(mₑ Bohr^2)  (* 1 Hartree*s)^2  => Hartree^2/Hartree => Hartree
  if unit==:hartree
    return B
  elseif unit==:cm1
    return B*2.1947e5 # Hartree -> cm-1
  elseif unit==:MHz
    return B*2.1947e5*2.99792458e4 # Hartree -> cm-1 -> MHz
  end
end

end #module
