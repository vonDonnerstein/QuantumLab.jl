"""
The *LaplaceModule* contains all functionality to work with the Laplace transformation
1/x ≈ ∑ₜweight[t] * exp(-x node[t]).
The nodes and weights are contained within the *LaplacePoints* type. Hackbusch and coworkers
have optimized and pretabulated these for different x-ranges according to minmax
(Takatsuka, Ten-No, Hackbusch, JCP, 129 (2008), 044112).
Their pretabulated files can be obtained with *downloadLaplacePointsHackbusch*. Use *readLaplacePointsHackbusch* to
generate a LaplacePoints object from that file. After downloading, one can also use *findLaplacePointsHackbuschPretableLarger*/*findLaplacePointsHackbuschPretableSmaller* to obtain the best guesses for a given range. The ideal Laplace
problem considers x∈ [1,R]. For the transformation to the x∈ [A,B]-type Laplace problem utility functions are provided (*transformRangeToIdealLaplace*,*transformWeightFromIdealLaplace*,*transformNodeFromIdealLaplace*,*transformLaplacePointFromIdealLaplace*).
"""
module LaplaceModule
export LaplacePoints, downloadLaplacePointsHackbusch, findLaplacePointsHackbuschPretableLarger, findLaplacePointsHackbuschPretableSmaller, transformRangeToIdealLaplace, transformNodeFromIdealLaplace, transformWeightFromIdealLaplace, transformLaplacePointFromIdealLaplace
using ..DocumentationModule
using ..BaseModule
using ZipFile
import Base.display
import Base.==

"""
LaplacePoints contain the collection of weights and corresponding nodes for the quadrature
according to
1/x = ∫ exp(-xt) ≈ ∑ₜweight[t] * exp(-x node[t])
"""
type LaplacePoints
  weights::Array{Float64,1}
  nodes::Array{Float64,1}
end

function display(lp::LaplacePoints)
  assert(length(lp.weights) == length(lp.nodes))
  println("  Exponential Factor         Weight")
  println("  ------------------         ------")
  for (node,weight) in zip(lp.nodes,lp.weights)
    @printf("  %.16f         %f\n",node,weight)
  end
end

"""
Evaluates to the name of the Laplace point file which corresponds to the given number of points and limit
in Hackbusch nomenclature. These files are contained in the zip-file obtained via downloadLaplacePointsHackbusch.
"""
function computeFilenameHackbusch(numberOfPoints::Integer,R::Real)
  k = numberOfPoints
  kstr = @sprintf("k%02d",k)
  Rstr = replace(replace(@sprintf(".%.0E",R),r"\+0*(\d+)",s"\1"),r"(\d+)E(\d+)",s"\1_\2")
  filename = "1_x$kstr$Rstr"
  return filename
end

"""
Extract the number of points (k) and the range (1,R) from the filename generated with computeFilenameHackbusch
"""
function computeFilenameHackbuschReverse(filename::AbstractString)
  m = match(r"1_xk(\d+).(\d+)_(\d+)",filename)
  k = parse(m.captures[1])
  R = parse("$(m.captures[2])E$(m.captures[3])")
  return (k,R)
end

"""
Download a zip-file of pretabulated Laplace points to target location (default:"hackbusch")
"""
function downloadLaplacePointsHackbusch(target::AbstractString="hackbusch")
  download("https://www.mis.mpg.de/scicomp/EXP_SUM/1_x/1_xData.zip",target)
end
@doc GenericCitation("https://www.mis.mpg.de/scicomp/EXP_SUM/1_x/") downloadLaplacePointsHackbusch

"""
Generate a LaplacePoints object by extracting from the pretable zip-file obtained with downloadLaplacePointsHackbusch.
"""
function readLaplacePointsHackbusch(numberOfPoints::Integer,R::Real,zipfile::AbstractString="hackbusch")
  k = numberOfPoints
  dir = ZipFile.Reader(zipfile)
  file = ZipFile.ReadableFile
  for f in dir.files
    if f.name==computeFilenameHackbusch(k,R)
      file = f
    end
  end

  result = LaplacePoints(Array{Float64,1}(),Array{Float64,1}())
  for line in 1:k
    regmatch = match(Regex(" *(?P<weight>$floatregex) .*omega.*"),readline(file))
    push!(result.weights,float(regmatch[:weight]))
  end
  for line in k+1:2k
    regmatch = match(Regex(" *(?P<node>$floatregex) .*alpha.*"),readline(file))
    push!(result.nodes,float(regmatch[:node]))
  end
  close(file)
  return result
end

"""
Returns the LaplacePoints for the smallest range larger than the requested one for which
pretabulated Laplace points can be found. If the requested range is larger than any pretabulated one,
we expect the largest pretabulated one to be Rₖ, so the function issues a warning and
returns the LaplacePoints for that range.
"""
function findLaplacePointsHackbuschPretableLarger(numberOfPoints::Integer,R::Float64,zipfile::AbstractString="hackbusch")
  dir = ZipFile.Reader(zipfile)
  Rtight = Inf
  Rmax = 0
  for f in dir.files
    if ismatch(r"1_xk(\d+).(\d+)_(\d+)",f.name)
      k,range = computeFilenameHackbuschReverse(f.name)
      if k == numberOfPoints
	Rmax = (Rmax > range) ? Rmax : range;
	Rtight = (R <= range < Rtight) ? range : Rtight;
      end
    end
  end

  if Rtight < Inf
    return (readLaplacePointsHackbusch(numberOfPoints,Rtight,zipfile),Rtight)
  else
    # Expect R > Rₖ, so just return largest available table for the given number of points
    warn("numberOfPoints not sufficiently large ?!")
    return (readLaplacePointsHackbusch(numberOfPoints,Rmax,zipfile),Rmax)
  end
end

"""
Returns the LaplacePoints for the largest range smaller than the requested one for which
pretabulated Laplace points can be found.
"""
function findLaplacePointsHackbuschPretableSmaller(numberOfPoints::Integer,R::Float64,zipfile::AbstractString="hackbusch")
  dir = ZipFile.Reader(zipfile)
  Rtight = 0.
  for f in dir.files
    if ismatch(r"1_xk(\d+).(\d+)_(\d+)",f.name)
      k,range = computeFilenameHackbuschReverse(f.name)
      if k == numberOfPoints
	Rtight = (Rtight < range <= R) ? range : Rtight;
      end
    end
  end

  if Rtight > 0
    return (readLaplacePointsHackbusch(numberOfPoints,Rtight,zipfile),Rtight)
  else
    error("No pretabulated value for $numberOfPoints points and a Range of [1., $R] or smaller.")
  end
end

"""
For the ideal Laplace problem the values x must be known to be in the range [1,R].
If 1/x can in general be within [A,B], then all nodes and weights must be scaled with 1/A,
so that x' := x/A lies within [1,B/A (=:R)].

Consider: 1/x = ∑ w exp(-α x)  =>  A/x = ∑ A w exp(-α x)
                               =>  A/x = ∑ (A w) exp(-(A α) x/A)
                   =>  1/x' = ∑ w' exp(-α' x')
"""
function transformRangeToIdealLaplace(A::Float64,B::Float64)
  return (1,B/A)
end

"""
cmp. help(transformRangeToIdealLaplace)
"""
function transformNodeFromIdealLaplace(α::Float64,A::Float64)
  return α/A
end

"""
cmp. help(transformRangeToIdealLaplace)
"""
function transformWeightFromIdealLaplace(w::Float64,A::Float64)
  return w/A
end

"""
cmp. help(transformRangeToIdealLaplace)
"""
function transformLaplacePointFromIdealLaplace(lp::LaplacePoints, A)
  LaplacePoints(map(x->transformWeightFromIdealLaplace(x,A),lp.weights), map(x->transformNodeFromIdealLaplace(x,A),lp.nodes))
end

"""
1/x ≈ ∑ᵢ wᵢ exp(-αᵢ x)
"""
function computeInverseByLaplaceApproximation(x,lp::LaplacePoints)
  sum=0.
  for (node,weight) in zip(lp.nodes,lp.weights)
    sum += weight*exp(-node*x)
  end
  return sum
end
@doc JournalCitation(["A. Takatsuka", "S. Ten-No", "W. Hackbusch"],"J. Chem. Phys.",129,044112,2008) computeInverseByLaplaceApproximation

"""
1/(x² + ω²) ≈ ∑ᵢ wᵢ sin(ω αᵢ)/ω exp(-αᵢ x)
"""
function computeRPADenominatorByDoubleLaplace(x::Float64,ω::Float64,lp::LaplacePoints)
  sum = 0.
  for (node,weight) in zip(lp.nodes,lp.weights)
    sum += weight * sin(ω*node)/ω * exp(-node*x)
  end
  return sum
end
@doc JournalCitation(["H. F. Schurkus","C. Ochsenfeld"],"J. Chem. Phys.",144,031101,2016) computeRPADenominatorByDoubleLaplace

==(x::LaplacePoints, y::LaplacePoints) = x.weights == y.weights && x.nodes == y.nodes

end # LaplaceModule
