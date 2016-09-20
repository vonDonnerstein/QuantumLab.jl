"""
The *LaplaceModule* contains all functionality to work with the Laplace transformation
1/x ≈ ∑ₜweight[t] * exp(-x node[t]).
The nodes and weights are contained within the *LaplacePoints* type. Hackbusch and coworkers
have optimized and pretabulated these for different x-ranges according to minmax
(Takatsuka, Ten-No, Hackbusch, JCP, 129 (2008), 044112).
Their pretabulated files can be obtained with *downloadLaplacePointsHackbusch*. Use *readLaplacePointsHackbusch* to
generate a LaplacePoints object from such a file. After downloading one can also use *findLaplacePointsHackbuschPretableLarger*/*findLaplacePointsHackbuschPretableSmaller* to obtain the best guesses for a given range. The ideal Laplace
problem considers x∈ [1,R]. For the transformation to the x∈ [A,B]-type Laplace problem utility functions are provided (*transformRangeToIdealLaplace*,*transformWeightFromIdealLaplace*,*transformNodeFromIdealLaplace*,*transformLaplacePointFromIdealLaplace*).
"""
module LaplaceModule
export LaplacePoints, downloadLaplacePointsHackbusch, findLaplacePointsHackbuschPretableLarger, findLaplacePointsHackbuschPretableSmaller, transformRangeToIdealLaplace, transformNodeFromIdealLaplace, transformWeightFromIdealLaplace, transformLaplacePointFromIdealLaplace
using ..DocumentationModule
using ..BaseModule
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
A dictionary mapping the number of Laplace Points to the ranges for which Hackbusch and
coworkers list pretabulated values on their webpage.
"""
const hackbuschpretables = Dict{Int64,Array{Int64,1}}(
  1=>[2:10;],
  2=>[2:10; 20:10:50],
  3=>[2:10; 20:10:100; 200],
  4=>[2:10; 20:10:100; 200:100:500],
  5=>[2:10; 20:10:100; 200:100:1000; 2000],
  6=>[2:10; 20:10:100; 200:100:1000; 2000; 3000],
  7=>[2:10; 20:10:100; 200:100:1000; 2000:1000:7000],
  8=>[      10:10:100; 200:100:500; 700; 1000:1000:5000; 7000; 1e4; 2e4],
  9=>[      10:10:100; 200:100:500; 700; 1000:1000:5000; 7000; 1e4:1e4:3e4],
 10=>[      10:10:100; 200:100:500; 700; 1000:1000:5000; 7000; 1e4:1e4:3e4; 5e4; 1e5],
 11=>[      10:10:100; 200:100:500; 700; 1000:1000:5000; 7000; 1e4:1e4:3e4; 5e4; 1e5; 2e5],
 12=>[      10:10:100; 200:100:500; 700; 1000:1000:5000; 7000; 1e4:1e4:3e4; 5e4; 1e5:1e5:3e5],
 13=>[      10:10:100; 200:100:500; 700; 1000:1000:5000; 7000; 1e4:1e4:3e4; 5e4; 1e5:1e5:4e5],
 14=>[      10:10:100; 200:100:500; 700; 1000:1000:5000; 7000; 1e4:1e4:3e4; 5e4; 1e5:1e5:5e5; 7e5],
 15=>[      10:10:100; 200:100:1000; 2000:1000:1e4; 2e4:1e4:5e4; 7e4; 1e5; 2e5; 5e5; 1e6; 2e6],
 16=>[      10:10:100; 200:100:1000; 2000:1000:1e4; 2e4:1e4:5e4; 7e4; 1e5; 2e5; 5e5; 1e6:1e6:3e6],
 17=>[      20:10:100; 200:100:1000; 2000:1000:1e4; 2e4:1e4:5e4; 7e4; 1e5; 2e5; 5e5; 1e6:1e6:4e6],
 18=>[      40:10:100; 200:100:1000; 2000:1000:1e4; 2e4:1e4:5e4; 7e4; 1e5; 2e5; 5e5; 1e6:1e6:5e6; 7e6],
 19=>[      50:10:100; 200:100:1000; 2000:1000:1e4; 2e4:1e4:5e4; 7e4; 1e5; 2e5; 5e5; 1e6:1e6:5e6; 7e6; 1e7],
 20=>[      60:10:100; 200:100:1000; 2000:1000:1e4; 2e4:1e4:5e4; 7e4; 1e5; 2e5; 5e5; 1e6:1e6:5e6; 7e6; 1e7; 2e7],
 21=>[      90:10:100; 200:100:1000; 2000:1000:1e4; 2e4:1e4:5e4; 7e4; 1e5; 2e5; 5e5; 1e6:1e6:5e6; 7e6; 1e7:1e7:3e7],
 22=>[             100:100:500; 700; 1000:1000:5000; 7000; 1e4:1e4:5e4; 7e4; 1e5:1e5:5e5; 7e5; 1e6:1e6:5e6; 7e6; 1e7:1e7:4e7],
 23=>[             200:100:500; 700; 1000:1000:5000; 7000; 1e4:1e4:5e4; 7e4; 1e5:1e5:5e5; 7e5; 1e6:1e6:5e6; 7e6; 1e7:1e7:5e7; 7e7],
 24=>[             200:100:500; 700; 1000:1000:5000; 7000; 1e4:1e4:5e4; 7e4; 1e5:1e5:5e5; 7e5; 1e6:1e6:5e6; 7e6; 1e7:1e7:5e7; 7e7; 1e8],
 25=>[             200:100:500; 700; 1000:1000:10000;       2e4; 5e4;        1e5; 2e5; 5e5;    1e6; 2e6; 5e6;    1e7; 5e7;         1e8; 2e8],
 26=>[             400:100:1000;     2000:1000:10000;2e4:1e4:1e5;      2e5:1e5:5e5; 7e5;       1e6:1e6:5e6; 7e6; 1e7:1e7:5e7; 7e7; 1e8:1e8:3e8],
 27=>[                                               1e4:1e4:1e5;      2e5:1e5:1e6;            2e6:1e6:1e7; 2e7:1e7:5e7; 7e7; 1e8:1e8:4e8],
 28=>[                                                                                                                        1e8:1e8:4e8; 7e8],
 29=>[7e8; 1e9],
 30=>[7e8; 1e9; 2e9],
 31=>[7e8; 1e9; 2e9],
 32=>[7e8; 1e9:1e9:3e9],
 33=>[7e8; 1e9:1e9:4e9],
 34=>[7e8; 1e9:1e9:5e9; 7e9],
 35=>[7e8; 1e9:1e9:3e9; 5e9; 1e10],
 36=>[1e10; 2e10],
 37=>[1e10; 2e10],
 38=>[1e10; 3e10],
 39=>[1e10; 2e10; 4e10],
 40=>[1e10; 2e10; 5e10],
 41=>[1e10; 2e10; 3e10; 7e10],
 42=>[1e10; 2e10; 5e10; 1e11],
 43=>[1e11; 2e11],
 44=>[1e11; 2e11],
 45=>[1e11; 3e11],
 46=>[1e11; 2e11; 4e11],
 47=>[1e11; 5e11],
 48=>[1e11; 2e11; 3e11; 7e11],
 49=>[1e11; 2e11; 3e11; 5e11; 1e12],
 50=>[2e8:1e8:5e8; 7e8; 1e9:1e9:7e9; 1e10:1e10:5e10; 7e10; 1e11:1e11:5e11; 7e11; 1e12; 2e12],
 51=>[1e12; 2e12],
 52=>[1e12; 3e12],
 53=>[4e11; 1e12:1e12:4e12]
)
@doc GenericCitation("http://www.mis.mpg.de/scicomp/EXP_SUM/1_x/") hackbuschpretables

"""
Evaluates to the name of the Laplace point file which corresponds to the given number of points and limit
in Hackbusch nomenclature. See also downloadLaplacePointsHackbusch for how to obtain these files.
"""
function computeFilenameHackbusch(numberOfPoints::Integer,R::Real)
  k = numberOfPoints
  kstr = @sprintf("k%02d",k)
  Rstr = replace(@sprintf("_%.0E",R),r"\+0*(\d+)",s"\1")
  filename = "1_x$kstr$Rstr"
  return filename
end

"""
Download all files with pretabulated Laplace points to the given folder (default:"hackbusch")
"""
function downloadLaplacePointsHackbusch(targetdir::AbstractString="hackbusch")
  mkdir(targetdir)
  for (k,Rlist) in hackbuschpretables
    for R in Rlist
      filename = computeFilenameHackbusch(k,R)
      println(filename)
      download("http://www.mis.mpg.de/scicomp/EXP_SUM/1_x/$filename","$targetdir/$filename")
    end
  end
end
@doc GenericCitation("http://www.mis.mpg.de/scicomp/EXP_SUM/1_x/") downloadLaplacePointsHackbusch

"""
Generate a LaplacePoints object by reading in the specified file.
"""
function readLaplacePointsHackbusch(filename::AbstractString)
  k = countlines(filename)/2
  fd = open(filename)
  result = LaplacePoints(Array{Float64,1}(),Array{Float64,1}())
  for line in 1:k
    regmatch = match(Regex(" *(?P<weight>$floatregex) .*omega.*"),readline(fd))
    push!(result.weights,float(regmatch[:weight]))
  end
  for line in k+1:2k
    regmatch = match(Regex(" *(?P<node>$floatregex) .*alpha.*"),readline(fd))
    push!(result.nodes,float(regmatch[:node]))
  end
  close(fd)
  return result
end

function readLaplacePointsHackbusch(numberOfPoints::Integer,R::Real,dir::AbstractString="hackbusch")
  filename = computeFilenameHackbusch(numberOfPoints,R)
  readLaplacePointsHackbusch("$dir/$filename")
end
@doc """
Alternative to specifying the file directly one can also specify the numberOfPoints, upper limit,
and directory name (default:"hackbusch") to load the corresponding pretable file..
""" readLaplacePointsHackbusch

"""
Returns the LaplacePoints for the smallest range larger than the requested one for which
pretabulated Laplace points can be found. If the requested range is larger than any pretabulated one,
we expect the largest pretabulated one to be Rₖ, so the function issues a warning and
returns the LaplacePoints for that range.
"""
function findLaplacePointsHackbuschPretableLarger(numberOfPoints::Integer,R::Float64,dir::AbstractString="hackbusch")
  for r in hackbuschpretables[numberOfPoints]
    if R <= r
      return (readLaplacePointsHackbusch(numberOfPoints,r,dir),r)
    end
  end
  # Expect R > Rₖ, so just return largest available table for the given number of points
  warn("numberOfPoints not sufficiently large ?!")
  r = hackbuschpretables[numberOfPoints][end]
  return (readLaplacePointsHackbusch(numberOfPoints,r,dir),r)
end

"""
Returns the LaplacePoints for the largest range smaller than the requested one for which
pretabulated Laplace points can be found.
"""
function findLaplacePointsHackbuschPretableSmaller(numberOfPoints::Integer,R::Float64,dir::AbstractString="hackbusch")
  for r in reverse(hackbuschpretables[numberOfPoints])
    if R >= r
      return (readLaplacePointsHackbusch(numberOfPoints,r,dir),r)
    end
  end
  error("No pretabulated value for $numberOfPoints points and a Range of [1., $R] or smaller.")
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
