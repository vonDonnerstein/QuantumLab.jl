module MatrixModule
import Base.show
import Base.*
import Base.+
import QuantumLab
export SparseMatrixBCSR, SparseMatrixBCSRSymmetric, purgeSparseMatrix!, convertSpMToMBCSR

function readMatrix(filename::String,dimensions::Tuple{Integer,Integer}=(0,0))
  floats = reinterpret(Float64,read(filename))

  if dimensions==(0,0) # try quadratic matrix
    dim = Int64(sqrt(length(floats)))
    if dim*dim != length(floats)
      error("Non-square matrix. Specify dimensions.")
    end
    dimensions = (dim,dim)
  end

  return reshape(floats,dimensions)
end

function ASCIIArt(matrix::Matrix)
  result = Matrix{Char}(size(matrix))
  for i in 1:length(matrix)
    if matrix[i] >= 10e-0
      result[i] = 'M'
    elseif matrix[i] >= 10e-1
      result[i] = 'N'
    elseif matrix[i] >= 10e-2
      result[i] = '@'
    elseif matrix[i] >= 10e-3
      result[i] = '%'
    elseif matrix[i] >= 10e-4
      result[i] = '#'
    elseif matrix[i] >= 10e-5
      result[i] = 't'
    elseif matrix[i] >= 10e-6
      result[i] = '+'
    elseif matrix[i] >= 10e-7
      result[i] = ';'
    elseif matrix[i] >= 10e-8
      result[i] = ':'
    elseif matrix[i] >= 10e-9
      result[i] = ','
    elseif matrix[i] >= 10e-10
      result[i] = '.'
    else
      result[i] = ' '
    end
  end
  return result
end

function writePixmap(filename::String,matrix::Matrix)
  # taking the specification from wikipedia
  dimensions = size(matrix)
  colors = ['M' => "#ff0000",
            'N' => "#ff3300",
            '@' => "#ff3333",
            '%' => "#ff6633",
            '#' => "#ff6666",
            't' => "#ff9966",
            '+' => "#ff9999",
            ';' => "#ffbb99",
            ':' => "#ffbbbb",
            ',' => "#ffeebb",
            '.' => "#ffeeee",
            ' ' => "#ffffff",
            '0' => "#bbbbbb", # zero by matrix format (blocking)
	   ]
  asciiart = ASCIIArt(matrix)
  open(filename,"w") do fn
    write(fn,"/* XPM */\n")
    write(fn,"static char* $(split(filename,'.')[1])[] = {\n")
    write(fn,"\"$(dimensions[2]) $(dimensions[1]) $(length(colors)) 1\",\n")
    for (sym,color) in colors
      write(fn,"\"$sym c $color\",\n")
    end
    for line in 1:dimensions[1]
      write(fn,"\"$(join(asciiart[line,:]))\",\n")
    end
  end
end

#SPARSE MATRIX CONSTRUCTORS

#@vonDonnerstein: Hyperlink not working don't know why
"""
    SparseMatrixBCSR(val,col,rowptr,rowBlockingPattern,colBlockingPattern)
Type definition for sparse matrix in BCSR format, for reference see: [ResearchGate]
(https://www.researchgate.net/figure/Compressed-Row-Storage-CRS-and-Block-Compressed-Row-Storage-BCRS-examples-for-the_fig1_29617855)

# Example
```julia-repl
julia> SparseMatrixBCSR([[1. 2.; 3. 4.], [5. 6.; 7. 8.]], [1,2], [0,1,2], [(1,2),(3,4)], [(1,2),(3,4)])
4×4 Array{Any,2}:
 1.0  2.0   *    * 
 3.0  4.0   *    * 
  *    *   5.0  6.0
  *    *   7.0  8.0
```
Note that the blocking pattern has to fulfill the following requirements: \n
	1) the first tuple has to have 1 as the first entry, e.g. [(1,...] \n
	2) the last tuple has to have the matrix dimension as last entry, e.g. [...,dim)] \n
	3) the last integer+1 of every tuple has to be the first integer of the next tuple, e.g. [...,(4,7),(8,10),...]
"""
type SparseMatrixBCSR
	val::Array{Array{Float64,2},1}
	col::Array{Int,1}
	rowptr::Array{Int,1}
	rowpattern::Array{Tuple{Int64,Int64},1}
	colpattern::Array{Tuple{Int64,Int64},1}
end

"""
    SparseMatrixBCSR(M,rowBlockingPattern,colBlockingPattern)
returns a matrix M in BCSR format with the given row- and column blocking pattern

# Example
```julia-repl
julia> M = [1.0 2.0 0.0 0.0 9.0 10.0; 3.0 0.0 4.0 0.0 11.0 12.0; 5.0 6.0 0.0 0.0 0.0 0.0; 0.0 7.0 8.0 0.0 0.0 0.0]
4×6 Array{Float64,2}:
 1.0  2.0  0.0  0.0   9.0  10.0
 3.0  0.0  4.0  0.0  11.0  12.0
 5.0  6.0  0.0  0.0   0.0   0.0
 0.0  7.0  8.0  0.0   0.0   0.0

julia> rowBlockingPattern = [(1,2),(3,4)]
2-element Array{Tuple{Int64,Int64},1}:
 (1, 2)
 (3, 4)

julia> colBlockingPattern = [(1,3),(4,6)]
2-element Array{Tuple{Int64,Int64},1}:
 (1, 3)
 (4, 6)

julia> SparseMatrixBCSR(M,rowBlockingPattern,colBlockingPattern)
4×6 Array{Any,2}:
 1.0  2.0  0.0  0.0   9.0  10.0
 3.0  0.0  4.0  0.0  11.0  12.0
 5.0  6.0  0.0   *     *     * 
 0.0  7.0  8.0   *     *     * 
```
"""
function SparseMatrixBCSR(M::Array{Float64,2},rowpattern::Array{Tuple{Int64,Int64},1},colpattern::Array{Tuple{Int64,Int64},1})
	(checkBlockingPattern(M,rowpattern,1) && checkBlockingPattern(M,colpattern,2)) || error("The blocking pattern does not match the matrix or is flawed,
	please see blocking pattern requirements: help?> SparseMatrixBCSR")
    val,col,rowptr = computeSparseMatrixBCSRFields(M,rowpattern,colpattern)
	SparseMatrixBCSR(val,col,rowptr,rowpattern,colpattern)
end

"""
    SparseMatrixBCSRSymmetric(val,col,rowptr,blockingPattern)
Type definition for symmetric sparse matrix in BCSR format, special case of ```SparseMatrixBCSR```, further help:
```julia-repl 
help?> SparseMatrixBCSR
```
"""
type SparseMatrixBCSRSymmetric
	val::Array{Union{Array{Float64,2},LowerTriangular{Float64,Array{Float64,2}}},1}
	col::Array{Int,1}
	rowptr::Array{Int,1}
	pattern::Array{Tuple{Int64,Int64},1}
end

"""
    SparseMatrixBCSRSymmetric(M,blockingPattern)
returns a matrix M in symmetric BCSR format with the given blocking pattern

# Example
```julia-repl
julia> M = [1.0 2.0 0.0 6.0; 2.0 3.0 0.0 7.0; 0.0 0.0 0.0 0.0; 6.0 7.0 0.0 0.0]
 4×4 Array{Float64,2}:
 1.0  2.0  0.0  6.0
 2.0  3.0  0.0  7.0
 0.0  0.0  0.0  0.0
 6.0  7.0  0.0  0.0

julia> blockingPattern = [(1,2),(3,4)]
2-element Array{Tuple{Int64,Int64},1}:
 (1, 2)
 (3, 4)

julia> SparseMatrixBCSRSymmetric(M,blockingPattern)
4×4 Array{Any,2}:
 1.0  2.0  0.0  6.0
 0.0  3.0  0.0  7.0
 0.0  0.0   *    * 
 6.0  7.0   *    * 
```
"""
function SparseMatrixBCSRSymmetric(M::Array{Float64,2},pattern::Array{Tuple{Int64,Int64},1})
	checkBlockingPattern(M,pattern,1) || error("The blocking pattern does not match the matrix or is flawed,
	please see blocking pattern requirements: ```julia-repl help?> SparseMatrixBCSR ```")
	val,col,rowptr = computeSparseMatrixBCSRSymmetricFields(M,pattern,true)
	SparseMatrixBCSRSymmetric(val,col,rowptr,pattern)
end

"""
	show(io,SpM)
displays a sparse matrix in BCSR format with * when norm(block) < 1e-10
"""
function show(io::IO,SpM::Union{SparseMatrixBCSR,SparseMatrixBCSRSymmetric})
	M = fillSpMWithChar(SpM,*)
	display(M)
end

#HELPERS (+ functions for advanced users)

"""
    computeBlockingPattern(blockingPattern)
returns the matrix blocking pattern itself (for fallback reasons)
"""
function computeBlockingPattern(pattern::Vector{Tuple{Int64,Int64}})
    return pattern
end

"""
    checkBlockingPattern(M,blockingPattern,dim)
checks blocking pattern for (1) first element = 1, (2) last element = dimension of matrix, (3) gaplessness
"""
function checkBlockingPattern(M::Array{Float64,2},pattern::Array{Tuple{Int64,Int64},1},dim::Int)
	check::Bool = true
	if pattern[1][1] != 1 return false end
	if pattern[end][2] != size(M,dim) return false end
	s::Int32 = length(pattern)-1
	
	for i = 1:s
		if pattern[i][2]+1 != pattern[i+1][1] return false end
	end

	return check
end

"""
    computeSparseMatrixBCSRFields(M,rowBlockingPattern,colBlockingPattern)
return val, col and rowptr for a matrix M with a given blocking pattern
"""
function computeSparseMatrixBCSRFields(M::Array{Float64,2},rowpattern::Array{Tuple{Int64,Int64},1},colpattern::Array{Tuple{Int64,Int64},1})
	s1::Int32						= length(rowpattern)
	s2::Int32						= length(colpattern)
	nnzb::Int64						= 0
	val::Array{Array{Float64,2},1}	= []
	col::Array{Int64,1}				= []
	rowptr::Array{Float64,1}		= [0]

	for i = 1:s1
		for j = 1:s2
			block = computeBlock(M,rowpattern,colpattern,i,j)
			if norm(block) != 0.
				push!(val, block)
				push!(col,j)
				nnzb += 1
			end
		end
		push!(rowptr,nnzb)
	end

	return val,col,rowptr
end

"""
    computeBCSRSpMFields(M,blockingPattern,Bool)
return val, col and rowptr for a matrix M with a given blocking pattern, when Bool = true then diagonal elements will be stored as LowerTriangular
"""
function computeSparseMatrixBCSRSymmetricFields(M::Array{Float64,2},pattern::Array{Tuple{Int64,Int64},1},tri::Bool)
	s::Int64						= length(pattern)
	nnzb::Int64						= 0
	val::Array{Array{Float64,2},1}	= []
	col::Array{Int64,1}				= []
	rowptr::Array{Float64,1}		= [0]

	for i = 1:s
		for j = 1:s
			if j <= i
				block = computeBlock(M,pattern,pattern,i,j)
				if tri == true && i == j
					block = LowerTriangular(block)
				end
				if norm(block) != 0.
					push!(val,block)
					push!(col,j)
					nnzb += 1
				end
			end
		end
		push!(rowptr,nnzb)
	end
		
	return val,col,rowptr
end

"""
    computeDifferenceRowptr(Integer,rowptr)
returns the running indices for running over all blocks in a row; b-a = # of blocks in a row
"""
function computeDifferenceRowptr(i::Int,rowptr)
    a = rowptr[i]+1
    b = rowptr[i+1]
	return a,b
end

"""
    computeBlock(M,rowBlockingPattern,colBlockingPattern,i,j)
returns the indices for the block to be stored in when converting SpM to M

# Example
```julia-repl
julia>  M = [1.0 2.0 0.0 6.0; 2.0 3.0 0.0 7.0; 0.0 0.0 0.0 0.0; 6.0 7.0 0.0 0.0]
4×4 Array{Float64,2}:
 1.0  2.0  0.0  6.0
 2.0  3.0  0.0  7.0
 0.0  0.0  0.0  0.0
 6.0  7.0  0.0  0.0

julia> computeBlock(M,[(1,2),(3,4)],[(1,2),(3,4)],1,1)
2×2 Array{Float64,2}:
 1.0  2.0
 2.0  3.0

julia> computeBlock(M,[(1,2),(3,4)],[(1,2),(3,4)],1,2)
 2×2 Array{Float64,2}:
 0.0  6.0
 0.0  7.0
```
Note that i and j are the row/column indices of the returned block
"""
function computeBlock(M::Array{Float64,2},rowpattern::Array{Tuple{Int64,Int64},1},colpattern::Array{Tuple{Int64,Int64},1},i::Int,j::Int)
	block::Array{Float64,2} = M[rowpattern[i][1]:rowpattern[i][2],colpattern[j][1]:colpattern[j][2]]
	return block
end 

"""
    computeSegmentForBlock(SpM,i,j)
returns the matrix indices M[a:b,c:d] computed from the block matrix indices i (row) and j (columnrow) and j (column))

# Example
```julia-repl
julia> M = [1.0 2.0 0.0 6.0; 2.0 3.0 0.0 7.0; 0.0 0.0 0.0 0.0; 6.0 7.0 0.0 0.0]
4×4 Array{Float64,2}:
 1.0  2.0  0.0  6.0
 2.0  3.0  0.0  7.0
 0.0  0.0  0.0  0.0
 6.0  7.0  0.0  0.0

julia> SpM = SparseMatrixBCSR(M,[(1,2),(3,4)],[(1,2),(3,4)])
4×4 Array{Any,2}:
 1.0  2.0  0.0  6.0
 2.0  3.0  0.0  7.0
 0.0  0.0   *    * 
 6.0  7.0   *    * 

julia> computeSegmentForBlock(SpM,1,1)
(1, 2, 1, 2)
```
"""
function computeSegmentForBlock(SpM::Union{SparseMatrixBCSR,SparseMatrixBCSRSymmetric},i,j)
	if isa(SpM,SparseMatrixBCSR)
		rowpattern = SpM.rowpattern
		colpattern = SpM.colpattern
	elseif isa(SpM,SparseMatrixBCSRSymmetric)
		rowpattern = SpM.pattern
		colpattern = SpM.pattern
	end
	a = rowpattern[i][1]
	b = rowpattern[i][2]
	c = colpattern[j][1]
	d = colpattern[j][2]

	return a,b,c,d
end
	
"""
    fillSpMWithChar(SpM,Char)
returns the corresponding dense matrix M for a sparse matrix SpM and every 'zero block' gets filled with the specified Char
"""
function fillSpMWithChar(SpM::Union{SparseMatrixBCSR,SparseMatrixBCSRSymmetric},a::Any)
	if isa(SpM,SparseMatrixBCSR)
		dim1::Int32		= SpM.rowpattern[end][2]
		dim2::Int32		= SpM.colpattern[end][2]
	elseif isa(SpM,SparseMatrixBCSRSymmetric)
		dim1			= SpM.pattern[end][2]
		dim2			= SpM.pattern[end][2]
	end
		
	M::Array{Any,2}		= fill(a,dim1,dim2)
	len::Int32			= length(SpM.rowptr)-1
	δ::Array{Int,1}		= [SpM.rowptr[i+1]-SpM.rowptr[i] for i in 1:len]

	j = 1
	for i = 1:len
		z = 0
		while z < δ[i]
			a,b,c,d = computeSegmentForBlock(SpM,i,SpM.col[j])
			el = SpM.val[j]
			if isa(SpM,SparseMatrixBCSR)
				M[a:b,c:d] = el
			elseif isa(SpM,SparseMatrixBCSRSymmetric)
				if a == c && b == d && typeof(SpM.val[j]) == LowerTriangular{Float64,Array{Float64,2}} 
					el = Symmetric(el,:L)
				end
				M[a:b,c:d] = el
				M[c:d,a:b] = el'
			end
			j += 1
			z += 1
		end
	end
		
	return M
end

"""
    computeRowFromRowptr(SpM)
returns the attribute row (analog to SparseMatrixBCSR field col) for the specified sparse matrix
"""
function computeRowFromRowptr(SpM::SparseMatrixBCSR)
	len::Int32          = length(SpM.rowptr)-1
	δ::Array{Int,1}     = [SpM.rowptr[i+1]-SpM.rowptr[i] for i in 1:len]
	row::Array{Int,1}	= []
	
	for i = 1:length(δ), j = 1:δ[i]
		push!(row,i)
	end
	
	return row
end

"""
    computeRowptrFromRow(row,SpM)
returns the SparseMatrixBCSR field rowptr computed from row for the specified sparse matrix
"""
function computeRowptrFromRow(row::Array{Int64,1},SpM::SparseMatrixBCSR)
	rowptr::Array{Int64,1}	= zeros(length(SpM.rowpattern)+1)
	δ::Array{Int64,1}		= computeNumberOfDistinctIntegersInVector(row,SpM)
	
	for i = 1:length(rowptr)
		rowptr[i] = δ[i]+sum(δ[1:i-1])
	end
	
	return rowptr
end

"""
    computeNumberOfDistinctIntegersInVector(Vector,SpM)
returns the number of distinct integers in the specified vector
"""
function computeNumberOfDistinctIntegersInVector(vec::Array{Int64,1},SpM::SparseMatrixBCSR)
	res::Array{Int64,1} = zeros(SpM.rowpattern[end][2]+1)
	
	for i in vec
		res[i+1] += 1
	end
	
	return res
end

#@vonDonnerstein: I want something like 7/2 = 3 or 5/2 = 2 for len, like in JAVA when you divide Integers, this: Int64(round(length(pattern)/2,0)) is very clunky
"""
    computeTupleBlockingPattern(blockingPattern)
returns the blocking pattern in format [(x,y),...] from a blocking pattern in the format [x,y,...]
"""
function computeTupleBlockingPattern(pattern::Array{Int,1})
	len								= Int64(round(length(pattern)/2,0))
	a::Array{Tuple{Int64,Int64},1}	= []
	[push!(a,(pattern[2*(i-1)+1],pattern[2*(i-1)+2])) for i in 1:len]
	return a
end

"""
    convertSpMToMBCSR(SpM)
returns the dense matrix representation M of a sparse matrix SpM
"""
function convertSpMToMBCSR(SpM::Union{SparseMatrixBCSR,SparseMatrixBCSRSymmetric})
	M = fillSpMWithChar(SpM,0.)
	return M
end

"""
    convertVToBV(Vector,SpM)
returns the specified vector with the blocking pattern of a sparse matrix SpM applied
"""
function convertVToBV(vec::Array{Float64,1},SpM::SparseMatrixBCSR)
	len::Int64							= length(SpM.colpattern)
	blockvec::Array{Array{Float64,1}}	= []
	
	for i = 1:len
		push!(blockvec,vec[SpM.colpattern[i][1]:SpM.colpattern[i][2]])
	end	
	
	return blockvec
end

#MATRIX OPERATIONS

"""
    purgeSparseMatrix!(SpM)
eliminates all blocks with norm(block) < Θ from sparse matrix
"""
function purgeSparseMatrix!(SpM::SparseMatrixBCSR,Θ::Float64)
    len                     = length(SpM.rowptr)-1
    del::Array{Int64}		= []
    delrow::Array{Int64}	= [] 
    
    for i = 1:len
        a,b = computeDifferenceRowptr(i,SpM.rowptr)
        for j = a:b
		if norm(SpM.val[j]) < Θ
                push!(del,j)
                push!(delrow,i)
            end
        end
    end

    deleteat!(SpM.val,del)
    deleteat!(SpM.col,del)
    for i = 1:length(delrow)
        SpM.rowptr[(delrow[i]+1):end] -= 1
    end

    return SpM
end

function *(SpM::SparseMatrixBCSR,vec::Array{Float64,1})
	if SpM.colpattern[end][2] != length(vec) error("Dimensions do not match (ง •̀_•́)ง ผ(•̀_•́ผ) \n
	the matrix has a row length of $(SpM.colpattern[end][2]) while the vectors is of length $(length(vec)).") end
	len						= length(SpM.rowptr)-1
	res::Array{Float64,1}	= zeros(length(vec)) 
	blockVec				= convertVToBV(vec,SpM)
	
	for i = 1:len
		a,b = computeDifferenceRowptr(i,SpM.rowptr)
		for j = a:b
			res[SpM.colpattern[i][1]:SpM.colpattern[i][2]] += SpM.val[j] * blockVec[SpM.col[j]]
		end
	end
	
	return res
end

#TODO: findlast für Array von Tuples vom Typ findlast(pattern,(Int,3)), wuerde computeTupleBlockingPattern sparen
function *(SpM::SparseMatrixBCSR,M::Array{Float64,2})
	if SpM.colpattern[end][2] != size(M,1) error("Dimensions do not match ԅ(≖‿≖ԅ) \n
	the sparse matrix has a row length of $(SpM.colpattern[end][2]) while the dense matrix has a column length of $(size(M,1)).") end
	s		= size(M,2)
	pattern = collect(Iterators.flatten(SpM.rowpattern))
	pattern = pattern[1:findlast(pattern,s)]
	b		= computeTupleBlockingPattern(pattern)
	res		= zeros(SpM.rowpattern[end][2],size(M,2))	
	row		= computeRowFromRowptr(SpM)	
	
	for i = 1:length(SpM.val)
		for j = 1:length(b)
			res[SpM.colpattern[row[i]][1]:SpM.colpattern[row[i]][2],b[j][1]:b[j][2]] += SpM.val[i]*M[SpM.colpattern[SpM.col[i]][1]:SpM.colpattern[SpM.col[i]][2],b[j][1]:b[j][2]]
		end
	end
	
	return res
end

function *(SpM1::SparseMatrixBCSR,SpM2::SparseMatrixBCSR)
	if SpM1.colpattern != SpM2.rowpattern error("The blocking patterns of the two matrices do not match (ﾉ◕ヮ◕)ﾉ*:・ﾟ✧") end
	row1 = computeRowFromRowptr(SpM1)
	row2 = computeRowFromRowptr(SpM2)
	
	dict::Dict{Tuple{Int64,Int64},Array{Float64,2}} = Dict{Tuple{Int64,Int64},Array{Float64,2}}() 
	
	for i = 1:length(SpM1.val), j = 1:length(SpM2.val)
			if SpM1.col[i] == row2[j]
				if haskey(dict,(row1[i],SpM2.col[j])) 
					dict[(row1[i],SpM2.col[j])] += SpM1.val[i]*SpM2.val[j]
				else 
					dict[(row1[i],SpM2.col[j])] = SpM1.val[i]*SpM2.val[j]
				end
			end
	end
	
	dict2 = sort(collect(dict), by = x->(x[1],x[2]))
	rowptr = computeRowptrFromRow([x[1][1] for x in dict2],SpM1)
	col = [x[1][2] for x in dict2]
	val = [x[2] for x in dict2]
	
	return SparseMatrixBCSR(val,col,rowptr,SpM1.colpattern,SpM2.rowpattern)
end

end #end of module
