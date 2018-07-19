module MatrixModule
import Base.show
import Base.*
import Base.+
import QuantumLab
export BCSRSpM

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



"""
    BCSRSpM(val,col,rowptr,rowpattern,colpattern)
Type definition for sparse matrix in BCSR format
"""
type BCSRSpM
	val::Array{Array{Float64,2},1}
	col::Array{Int,1}
	rowptr::Array{Int,1}
	rowpattern::Array{Tuple{Int64,Int64},1}
	colpattern::Array{Tuple{Int64,Int64},1}
    BCSRSpM(val,col,rowptr,rowpattern,colpattern) = new(val,col,rowptr,rowpattern,colpattern)
end

"""
    BCSRSpM(M,rowpattern,colpattern)
returns a matrix M in BCSR format with the given blocking row- and colpattern
"""
function BCSRSpM(M::Array{Float64,2},rowpattern::Array{Tuple{Int64,Int64},1},colpattern::Array{Tuple{Int64,Int64},1})
	(checkPattern(M,rowpattern,1) && checkPattern(M,colpattern,2)) || error("Hey idiot your pattern doesn't match the matrix ᕦ(ò_óˇ)ᕤ")
    val,col,rowptr = computeBCSRSpMFields(M,rowpattern,colpattern)
	BCSRSpM(val,col,rowptr,rowpattern,colpattern)
end

"""
    symBCSRSpM(val,col,rowptr,pattern)
Type definition for symmetric sparse matrix in BCSR format
"""
type symBCSRSpM
	val::Array{Union{Array{Float64,2},LowerTriangular{Float64,Array{Float64,2}}},1}
	col::Array{Int,1}
	rowptr::Array{Tuple{Int64,Int64},1}
	pattern::Array{Tuple{Int64,Int64},1}
	symBCSRSpM(val,col,rowptr,pattern) = new(val,col,rowptr,pattern)
end

"""
    symBCSRSpM(M,pattern)
returns a matrix M in symmetric BCSR format with the given blocking pattern
"""
function symBCSRSpM(M::Array{Float64,2},pattern::Array{Tuple{Int64,Int64},1})
	checkPattern(M,pattern,1) || error("In case of error flip table ╯(°□°）╯︵ ┻━┻")
	val,col,rowptr = computeSymBCSRSpMFields(M,pattern,true)
	symBCSRSpM(val,col,rowptr,pattern)
end

"""
	show(io,SpM)
displays a sparse matrix in BCSR format with * when norm(block) < 1e-10
"""
function show(io::IO,SpM::Union{BCSRSpM,symBCSRSpM})
	M = fillSpMWithChar(SpM,*)
	display(M)
end

"""
    checkPattern(M,pattern,dim)
checks blocking pattern for (1) first element = 1, (2) last element = dimension of matrix, (3) gaplessness
"""
function checkPattern(M::Array{Float64,2},pattern::Array{Tuple{Int64,Int64},1},dim::Int)
	check::Bool = true
	if pattern[1][1] != 1 return false end
	if pattern[end][2] != size(M,dim) return false end
	s::Int32 = length(pattern)-1
	
	for i = 1:s
		if pattern[i][2]+1 != pattern[i+1][1] return false end
	end

	return check
end

"Define multiplication for BCSRSpM sparse matrix with vector"
*(SpM::BCSRSpM,vec::Array{Float64,1}) = multiplySpMV(SpM,vec)

"Define multiplication for BCSRSpM sparse matrix with dense matrix"
*(SpM::BCSRSpM,M::Array{Float64,2}) = multiplySpMM(SpM,M)

"Define multiplication for BCSRSpM sparse matrix with BCSRSpM sparse matrix"
*(SpM1::BCSRSpM,SpM2::BCSRSpM) = multiplySpMSpM2(SpM1,SpM2)

"Define addition for tuples of type (Int64,Int64,Array{Float64,2})"
+(A::Tuple{Int64,Int64,Array{Float64,2}},B::Tuple{Int64,Int64,Array{Float64,2}}) = (A[1],A[2],A[3]+B[3])

"""
    computePattern(pattern)
returns the matrix blocking pattern itself (for fallback reasons)
"""
function computePattern(pattern::Vector{Tuple{Int64,Int64}})
    return pattern
end

"""
    purgeSparseMatrix!(SpM)
eliminates all blocks with norm(block) < 1e-10 from sparse matrix
"""
function purgeSparseMatrix!(SpM::BCSRSpM)
    len						= length(SpM.rowptr)-1
    del::Array{Int64}		= []
    delrow::Array{Int64}	= [] 
    
    for i = 1:len
        a,b = computeDifferenceRowptr(i,SpM.rowptr)
        for j = a:b
            if norm(SpM.val[j]) < 1e-10
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

#TODO for loop kuerzer schreiben
"""
    symmetrizeLowerTriangular(lt)
returns symmetric of lower triangular matrix
"""
function symmetrizeLowerTriangular(lt::LowerTriangular{Float64,Array{Float64,2}})
	M::Array{Float64,2} = Array{Float64,2}(lt)
	n					= size(M,1)
	N					= zeros(n,n)
	
	for i = 1:n, j = 1:n
		N[i,j] = M[i,j]
		if i == j N[i,j] = 0. end
	end

	return M+N'
end

"""
    computeBlock(M,rowpattern,colpattern,i,j)
returns the indices for the block to be stored in when converting SpM to M
"""
function computeBlock(M::Array{Float64,2},rowpattern::Array{Tuple{Int64,Int64},1},colpattern::Array{Tuple{Int64,Int64},1},i::Int,j::Int)
	block::Array{Float64,2} = M[rowpattern[i][1]:rowpattern[i][2],colpattern[j][1]:colpattern[j][2]]
	return block
end 

"""
    computeBCSRSpMFields(M,rowpattern,colpattern)
return val, col and rowptr for a matrix M with a given blocking pattern
"""
function computeBCSRSpMFields(M::Array{Float64,2},rowpattern::Array{Tuple{Int64,Int64},1},colpattern::Array{Tuple{Int64,Int64},1})
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
    computeBCSRSpMFields(M,pattern,Bool)
return val, col and rowptr for a matrix M with a given blocking pattern, when Bool = true then diagonal elements will be stored as LowerTriangular
"""
function computeSymBCSRSpMFields(M::Array{Float64,2},pattern::Array{Tuple{Int64,Int64},1},tri::Bool)
	s::Int64						= length(pattern)
	nnzb::Int64						= 0
	val::Array{Array{Float64,2},1}	= []
	col::Array{Int64,1}				= []
	rowptr::Array{Float64,1}		= [0]
	block::Array{Float64,2}			= []

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
    convertSpMToMBCSR(SpM)
returns the dense matrix representation M of a sparse matrix SpM
"""
function convertSpMToMBCSR(SpM::Union{BCSRSpM,symBCSRSpM})
	M = fillSpMWithChar(SpM,0.)
	return M
end

#TODO: macht aehnliches wie computeBlock -> vereinigen
"""
    computeSegmentForBlock(SpM,i,j)
returns indices a,b,c,d where a given block from SpM is stored
"""
function computeSegmentForBlock(SpM::Union{BCSRSpM,symBCSRSpM},i,j)
	if isa(SpM,BCSRSpM)
		rowpattern = SpM.rowpattern
		colpattern = SpM.colpattern
	elseif isa(SpM,symBCSRSpM)
		rowpattern = SpM.pattern
		colpattern = SpM.pattern
	end
	a = rowpattern[i][1]
	b = rowpattern[i][2]
	c = colpattern[j][1]
	d = colpattern[j][2]

	return a,b,c,d
end
	
#TODO: gefaellt mir noch nicht, andere Lösung für Anzeigen der Matrizen finden
"""
    fillSpMWithChar(SpM,Char)
returns the corresponding dense matrix M for a sparse matrix SpM and every 'zero block' gets filled with the specified Char
"""
function fillSpMWithChar(SpM::Union{BCSRSpM,symBCSRSpM},a::Any)
	if isa(SpM,BCSRSpM)
		dim1::Int32		= SpM.rowpattern[end][2]
		dim2::Int32		= SpM.colpattern[end][2]
	elseif isa(SpM,symBCSRSpM)
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
			if isa(SpM,BCSRSpM)
				M[a:b,c:d] = el
			elseif isa(SpM,symBCSRSpM)
				if a == c && b == d && typeof(SpM.val[j]) == LowerTriangular{Float64,Array{Float64,2}} 
					el = symmetrizeLowerTriangular(el)
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
    computeDifferenceRowptr(Integer,rowptr)
returns the running indices to run over all blocks in a row; b-a = # of blocks in a row
"""
function computeDifferenceRowptr(i::Int,rowptr)
    a = rowptr[i]+1
    b = rowptr[i+1]
	return a,b
end

#"Function applies pattern of sparse matrix to a vector"
"""
    convertVToBV(Vector,SpM)
returns the specified vector with the blocking pattern of a sparse matrix SpM applied
"""
function convertVToBV(vec::Array{Float64,1},SpM::BCSRSpM)
	len::Int64							= length(SpM.colpattern)
	blockvec::Array{Array{Float64,1}}	= []
	
	for i = 1:len
		push!(blockvec,vec[SpM.colpattern[i][1]:SpM.colpattern[i][2]])
	end	
	
	return blockvec
end

"""
    multiplySpMV(SpM,Vector)
returns the product Vector2 = SpM*Vector
"""
function multiplySpMV(SpM::BCSRSpM,vec::Array{Float64,1})
	if SpM.colpattern[end][2] != length(vec) error("Dimensions do not match (ง •̀_•́)ง ผ(•̀_•́ผ)") end
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

"""
    computeRowFromRowptr(SpM)
returns the attribute row (analog to BCSRSpM field col) for the specified sparse matrix
"""
function computeRowFromRowptr(SpM::BCSRSpM)
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
returns the BCSRSpM field rowptr computed from row for the specified sparse matrix
"""
function computeRowptrFromRow(row::Array{Int64,1},SpM::BCSRSpM)
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
function computeNumberOfDistinctIntegersInVector(vec::Array{Int64,1},SpM::BCSRSpM)
	res::Array{Int64,1} = zeros(SpM.rowpattern[end][2]+1)
	
	for i in vec
		res[i+1] += 1
	end
	
	return res
end

#TODO: Int64(round(length(pattern)/2,0)) kuerzer schreiben
"""
    computeTuplePatternFromPattern(pattern)
returns the blocking pattern in format [(x,y),...] from a blocking pattern in the format [x,y,...]
"""
function computeTuplePatternFromPattern(pattern::Array{Int,1})
	len								= Int64(round(length(pattern)/2,0))
	a::Array{Tuple{Int64,Int64},1}	= []
	[push!(a,(pattern[2*(i-1)+1],pattern[2*(i-1)+2])) for i in 1:len]
	return a
end

#TODO: findlast für Array von Tuples vom Typ findlast(pattern,(Int,3)), wuerde computeTuplePatternFromPattern sparen
"""
    multiplySpMM(SpM,M)
returns the product M2 = SpM*M
"""
function multiplySpMM(SpM::BCSRSpM,M::Array{Float64,2})
	if SpM.colpattern[end][2] != size(M,1) error("Dimensions do not match ԅ(≖‿≖ԅ)") end
	s		= size(M,2)
	pattern = collect(Iterators.flatten(SpM.rowpattern))
	pattern = pattern[1:findlast(pattern,s)]
	b		= computeTuplePatternFromPattern(pattern)
	res		= zeros(SpM.rowpattern[end][2],size(M,2))	
	row		= computeRowFromRowptr(SpM)	
	
	for i = 1:length(SpM.val)
		for j = 1:length(b)
			res[SpM.colpattern[row[i]][1]:SpM.colpattern[row[i]][2],b[j][1]:b[j][2]] += SpM.val[i]*M[SpM.colpattern[SpM.col[i]][1]:SpM.colpattern[SpM.col[i]][2],b[j][1]:b[j][2]]
		end
	end
	
	return res
end

"""
    multiplySpMSpM(SpM1,SpM2)
returns the product M = SpM1*SpM2
"""
function multiplySpMSpM(SpM1::BCSRSpM,SpM2::BCSRSpM)
	if SpM1.colpattern != SpM2.rowpattern error("Patterns do not match (ﾉ◕ヮ◕)ﾉ*:・ﾟ✧") end
	res		= zeros(SpM1.rowpattern[end][2],SpM2.colpattern[end][2])
	row1	= computeRowFromRowptr(SpM1)
	row2	= computeRowFromRowptr(SpM2)

	for i = 1:length(SpM1.val), j = 1:length(SpM2.val)
			if SpM1.col[i] == row2[j]
				res[SpM1.colpattern[row1[i]][1]:SpM1.colpattern[row1[i]][2],SpM2.rowpattern[SpM2.col[j]][1]:SpM2.rowpattern[SpM2.col[j]][2]] += SpM1.val[i]*SpM2.val[j]
			end
	end
	
	return res
end

"""
    multiplySpMSpM2(SpM1,SpM2)
returns the product SpM3 = SpM1*SpM2
"""
function multiplySpMSpM2(SpM1::BCSRSpM,SpM2::BCSRSpM)
	if SpM1.colpattern != SpM2.rowpattern error("Patterns do not match (ﾉ◕ヮ◕)ﾉ*:・ﾟ✧") end
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
	
	return BCSRSpM(val,col,rowptr,SpM1.colpattern,SpM2.rowpattern)
end

#TODO: push! geht nicht und col muss gedreht werden
function computeMatrixTransposeOfSpM(SpM::BCSRSpM)
	row = computeRowFromRowptr(SpM)
	SpMT::BCSRSpM = BCSRSpM([],[],[],SpM.rowpattern,SpM.colpattern)
	for i = 1:length(SpM.val)
		push!(SpMT.val,SpM.val[i]')
		push!(SpMT.col,row[i])
		row[i]			= SpM.col[i]
	end
	SpMT.rowptr = computeRowptrFromRow(row,SpM2)
	return SpMT
end

#@Henry: laesst sich so nicht in MatrixModule aufrufen weil ShellModule benutzt
#wird, gehoert aber mMn eher hier hin als in ShellModule
#"Function does a benchmark for SpMSpM vs SpMM vs MM"
#"""
#    benchMarkMatrixMultiplication()
#returns a dictionary with computation times for M*M, SpM*M and SpM*SpM for different amylose molecules (Glc)n
#"""
#function benchMarkMatrixMultiplication()
#	
#	testMolecule	= ["amylose2.xyz","amylose4.xyz","amylose8.xyz","amylose16.xyz","amylose32.xyz","amylose64.xyz"]
#	basisSet		= "STO-3G"
#	results::Dict{String,Tuple{Float64,Float64,Float64}} = Dict{String,Tuple{Float64,Float64,Float64}}()
#	
#	for i = 1:6
#		basis			= computeBasisShells(BasisSet(basisSet),Geometry(testMolecule[i]))
#		S				= computeMatrixOverlap(basis)
#		SpM				= BCSRSpM(S,basis,true)
#		purgeSparseMatrix!(SpM)
#		timingMM		= []
#		timingSpMM		= []
#		timingSpMSpM	= []
#		for j = 1:3
#			push!(timingMM,@elapsed S*S)
#			push!(timingSpMM,@elapsed SpM*S)
#			push!(timingSpMSpM,@elapsed SpM*SpM)
#		end
#		results[testMolecule[i]] = (mean(timingMM),mean(timingSpMM),mean(timingSpMSpM))
#	end
#	
#	return sort(collect(results), by = x->x[2])
#end
end
