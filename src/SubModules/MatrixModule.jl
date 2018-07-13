module MatrixModule
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
end
