using Colors
using Vega

function heatmap(mat::Array{Float64,2})
  Nrow, Ncol = size(mat)
  xvals = Array(Float64,length(mat))
  yvals = Array(Float64,length(mat))
  colorvals = Array(Float64,length(mat))
  index = 1
  for x in 1:Ncol
    for y in 1:Nrow
      xvals[index] = x
      yvals[index] = y
      colorvals[index] = mat[y,x]
      index+=1
    end
  end
  Vega.heatmap(x=xvals, y=yvals, color=colorvals)
#  plot(x=xvals,y=yvals,color=colorvals, Geom.rectbin, 
#  	Scale.x_discrete, Scale.y_discrete, 
#	Coord.cartesian(yflip=true,fixed=true),
#	Guide.xlabel(""),Guide.ylabel(""),
#	Scale.ContinuousColorScale(Scale.lab_gradient(colorant"white",colorant"red"),minvalue=0.))
end

function writeMatrix(filename::String,matrix::Matrix)
  stream = open(filename,"w+")
  write(stream,mat)
  close(stream)
end

function readMatrix{T<:Real}(filename::String,format::Type{T},dims::Tuple{Vararg{Int64}})
  stream = open(filename,"r")
  matrix = read(stream,format,dims)
  close(stream)
  return matrix
end
