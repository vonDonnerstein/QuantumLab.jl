"""
Although the BSE exposes a special SOAP-based API for programmatic access to its database, I had too much trouble getting the desired results from it. Therefore the solution in this module is based around regexing of the
human readable HTML output.
"""
module BasisSetExchangeModule
using ..BasisSetModule

export downloadBasisSetExchange, computeBasisSetsAvailable
import ..BasisSetModule.BasisSet
using HTTP
using Printf
using ZipFile

bseZip    = homedir()*"/pqs.zip"
bseFormat = "pqs"

function Base.getindex(dir::ZipFile.Reader,fn::String)
  for f in dir.files
    if f.name == fn
      return f
    end
  end
  @error "$fn not found."
end

function downloadBasisSetExchange(;bseFormat=bseFormat, bseZip=bseZip)
  response = HTTP.request("GET","http://www.basissetexchange.org/download/current/$bseFormat/zip")
  write(bseZip,response.body)
  @assert isfile(bseZip)
  println("Downloaded $bseFormat formatted basis sets to $bseZip.")
end

function computeBasisSetsAvailable(;bseZip=bseZip,bseFormat=bseFormat)
  list = Dict{String,String}()
  zip = ZipFile.Reader(bseZip)
  for f in zip.files
    if occursin(Regex("\\.$bseFormat\$"),f.name)
      shortname = match(Regex(".*/([^\\.]*).*$bseFormat"), f.name).captures[1]
      push!(list, shortname => f.name)
    end
  end
  return list
end

function computeBasisSetExchangeEntry(name::AbstractString; bseZip=bseZip)
  dict = computeBasisSetsAvailable(bseZip=bseZip)
  if haskey(dict,name)
    return dict[name]
  else
    result = String[]
    for (shortname,entry) in dict
      if (occursin(Regex(name,"i"),shortname))
        push!(result,entry)
      end
    end
  end
  if length(result) == 0
    error("No matching basis set entry found (searching for $name in $bseZip)")
  elseif length(result)>1
    display(result)
    error("More than one matching entry found.")
  end
  return result[1]
end

"""
    BasisSet(name)

creates a BasisSet object from the specified file (possibly with a corresponding file extension).
If no such file can be found the argument is interpreted as the name of a basis set and that basis
set obtained from Basis Set Exchange.
"""
function BasisSet(name::AbstractString; bseZip=bseZip)

  if isfile(name)
    return readBasisSetTX93(name)
  end

  for ext in ("tx93", "TX93", "Tx93")
    if isfile("$name.$ext")
      return readBasisSetTX93("$name.$ext")
    end
  end

  @info("BasisSet file $name not found on disk. Obtaining from BasisSetExchange...")
  if !isfile(bseZip)
    downloadBasisSetExchange(bseZip=bseZip)
  end
  location   = computeBasisSetExchangeEntry(name)
  str = ""
  open(bseZip) do bse # trying to get around ZipFile Issue #14
    str      = read(ZipFile.Reader(bse)[location], String)
  end
  return BasisSetModule.BasisSetTX93(str)
end

end # module
