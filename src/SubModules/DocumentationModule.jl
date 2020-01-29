module DocumentationModule
export Citation, JournalCitation, BookCitation, GenericCitation, summarize
import Base.show
import Base.convert
import Base.promote_rule
using Markdown
using InteractiveUtils

"""Abstract class supertyping all objects that are entries of a bibliography (for documenting functions). Each subtype has its corresponding fields and defines its own show function."""
abstract type Citation end

"""
Documentation is a generic type to which all types which are to be used with `@doc` can be automatically converted.
At the moment this means: Any type of Citation and Markdown
"""
struct Documentation
  doc::Union{Citation,Markdown.MD}
end

"""Citation class for articles in scientific jounals."""
mutable struct JournalCitation <: Citation
  authors::Array{AbstractString,1}
  journalabbr::AbstractString
  vol::UInt
  page::UInt
  year::UInt
end

function show(io::IO,mime::MIME"text/plain",cite::JournalCitation)
  show(io,mime,Markdown.parse(""" - $(join(cite.authors,", ")), $(cite.journalabbr) **$(cite.vol)**, *$(cite.page)* ($(cite.year))"""))
  print(io,"\n")
end

"""Citation class for complete Books. Refer to BookSectionCitation for more specific parts of a book."""
mutable struct BookCitation <: Citation
  authors::Array{AbstractString,1}
  title::AbstractString
  ISBN::AbstractString
end

function show(io::IO,mime::MIME"text/plain",cite::BookCitation)
  show(io,mime,Markdown.parse(""" - $(join(cite.authors,", ")): *$(cite.title)* (ISBN $(cite.ISBN))"""))
  print(io,"\n")
end

"""Citation class for anything that can not be described by any of the other subclasses of Citation. Make sure none of the others fits the purpose before using this one."""
struct GenericCitation <: Citation
  string::AbstractString
end

function show(io::IO,mime::MIME"text/plain",cite::GenericCitation)
  show(io,mime,Markdown.parse(""" - $(cite.string)"""))
  print(io,"\n")
end

function show(io::IO,mime::MIME"text/plain",docuArray::Array{Documentation,1})
  local cites = Array{Citation,1}()
  show(io,mime,Markdown.parse("""
  # Documentation
  """))
  print(io,"\n")
  for document in docuArray
    docu = document.doc
    if typeof(docu)==Markdown.MD
      show(io,mime,docu)
    elseif isa(docu,Citation)
      push!(cites,docu)
    end
  end

  print(io,"\n\n")

  show(io,mime,Markdown.parse("""
  # Citations:
  """))
  print(io,"\n")
  for citekey in cites
    show(io,mime,citekey)
  end
end

function show(io::IO,mime::MIME"text/plain",docuArray::Array{Citation,1})
  show(io,mime,Markdown.parse("""
  # Citations:
  """))
  print(io,"\n")
  for citekey in docuArray
    show(io,mime,citekey)
  end
end

function show(io::IO,mime::MIME"text/plain",docuArray::Array{T,1}) where T<:Citation
  show(io,mime,Markdown.parse("""
  # Citations:
  """))
  print(io,"\n")
  for citekey in docuArray
    show(io,mime,citekey)
  end
end

"""
Taken from base/docs/Docs.jl, this function shows a summary of the fields, type and supertypes of a DataType.
In Docs.jl the function `typesummary` is only meant to be called when no other documentation can be found for it.
If `typesummary` doesn't print "No documentation found" anymore in the future, then this function becomes obsolete.
"""
function summarize(f::DataType)
    parts = [
    """
    **Summary:**
    ```julia
    $(f.abstract ? "abstract" : f.mutable ? "type" : "immutable") $f <: $(supertype(f))
    ```
    """
    ]
    if !f.abstract && !isempty(fieldnames(f))
        pad    = maximum([length(string(f)) for f in fieldnames(f)])
        fields = ["$(rpad(f, pad)) :: $(t)" for (f, t) in zip(fieldnames(f), f.types)]
        push!(parts,
        """
        **Fields:**
        ```julia
        $(join(fields, "\n"))
        ```
        """)
    end
    if !isempty(subtypes(f))
        push!(parts,
        """
        **Subtypes:**
        ```julia
        $(join(subtypes(f), "\n"))
        ```
        """)
    end
    Markdown.parse(join(parts, "\n"))
end

# promotion and conversion
Base.promote_rule(::Type{String},::Type{Documentation})         = Documentation
Base.promote_rule(::Type{String},::Type{Citation})         	= Documentation
Base.promote_rule(::Type{String},::Type{GenericCitation}) 	= Documentation
Base.promote_rule(::Type{String},::Type{BookCitation}) 		= Documentation
Base.promote_rule(::Type{String},::Type{JournalCitation}) 	= Documentation
Base.promote_rule(::Type{Markdown.MD},::Type{Documentation}) 	= Documentation
Base.promote_rule(::Type{Markdown.MD},::Type{Citation})         = Documentation
Base.promote_rule(::Type{Markdown.MD},::Type{GenericCitation}) 	= Documentation
Base.promote_rule(::Type{Markdown.MD},::Type{BookCitation}) 	= Documentation
Base.promote_rule(::Type{Markdown.MD},::Type{JournalCitation}) 	= Documentation
Base.promote_rule(::Type{Documentation},::Type{GenericCitation}) 	= Documentation
Base.promote_rule(::Type{Documentation},::Type{BookCitation}) 		= Documentation
Base.promote_rule(::Type{Documentation},::Type{JournalCitation}) 	= Documentation
Base.convert(::Type{Documentation},str::String) = Documentation(Markdown.parse(str))
Base.convert(::Type{Documentation},md::Markdown.MD) = Documentation(md)
Base.convert(::Type{Documentation},cite::Citation) = Documentation(cite)

macro ignore_warnings(expr)
  quote
    STDERRorig = stderr
    drop = redirect_stderr()
    result = $expr
    redirect_stderr(STDERRorig)
    return result
  end
end

end # module
