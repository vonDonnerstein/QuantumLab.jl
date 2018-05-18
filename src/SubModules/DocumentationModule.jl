module DocumentationModule
export Citation, JournalCitation, BookCitation, GenericCitation, summarize, @add_doc
import Base.display
import Base.convert
import Base.promote_rule

"""Abstract class supertyping all objects that are entries of a bibliography (for documenting functions). Each subtype has its corresponding fields and defines its own display function."""
abstract type Citation end

"""
Documentation is a generic type to which all types which are to be used with `@doc` can be automatically converted.
At the moment this means: Any type of Citation and Markdown
"""
type Documentation
  doc::Union{Citation,Base.Markdown.MD}
end

"""Citation class for articles in scientific jounals."""
type JournalCitation <: Citation
  authors::Array{AbstractString,1}
  journalabbr::AbstractString
  vol::UInt
  page::UInt
  year::UInt
end

function display(cite::JournalCitation)
  display(Base.Markdown.parse(""" - $(join(cite.authors,", ")), $(cite.journalabbr) **$(cite.vol)**, *$(cite.page)* ($(cite.year))"""))
end

"""Citation class for complete Books. Refer to BookSectionCitation for more specific parts of a book."""
type BookCitation <: Citation
  authors::Array{AbstractString,1}
  title::AbstractString
  ISBN::AbstractString
end

function display(cite::BookCitation)
  display(Base.Markdown.parse(""" - $(join(cite.authors,", ")): *$(cite.title)* (ISBN $(cite.ISBN))"""))
end

"""Citation class for anything that can not be described by any of the other subclasses of Citation. Make sure none of the others fits the purpose before using this one."""
type GenericCitation <: Citation
  string::AbstractString
end

function display(cite::GenericCitation)
  display(Base.Markdown.parse(""" - $(cite.string)"""))
end

function display(docuArray::Array{Documentation,1})
  local cites = Array{Citation,1}()
  display(Base.Markdown.parse("""
  # Documentation
  """))
  for document in docuArray
    docu = document.doc
    if typeof(docu)==Base.Markdown.MD
      display(docu)
    elseif isa(docu,Citation)
      push!(cites,docu)
    end
  end

  print('\n')

  display(Base.Markdown.parse("""
  # Citations:
  """))
  for citekey in cites
    display(citekey)
  end
end

function display(docuArray::Array{Citation,1})
  display(Base.Markdown.parse("""
  # Citations:
  """))
  for citekey in docuArray
    display(citekey)
  end
end

function display{T<:Citation}(docuArray::Array{T,1})
  display(Base.Markdown.parse("""
  # Citations:
  """))
  for citekey in docuArray
    display(citekey)
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
    if !isempty(fieldnames(f))
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
Base.promote_rule(::Type{Base.Markdown.MD},::Type{Documentation}) 	= Documentation
Base.promote_rule(::Type{Base.Markdown.MD},::Type{Citation}) 		= Documentation
Base.promote_rule(::Type{Base.Markdown.MD},::Type{GenericCitation}) 	= Documentation
Base.promote_rule(::Type{Base.Markdown.MD},::Type{BookCitation}) 	= Documentation
Base.promote_rule(::Type{Base.Markdown.MD},::Type{JournalCitation}) 	= Documentation
Base.promote_rule(::Type{Documentation},::Type{GenericCitation}) 	= Documentation
Base.promote_rule(::Type{Documentation},::Type{BookCitation}) 		= Documentation
Base.promote_rule(::Type{Documentation},::Type{JournalCitation}) 	= Documentation
Base.convert(::Type{Documentation},md::Base.Markdown.MD) = Documentation(md)
Base.convert(::Type{Documentation},cite::Citation) = Documentation(cite)

macro ignore_warnings(expr)
  quote
    STDERRorig = STDERR
    drop = redirect_stderr()
    result = $expr
    redirect_stderr(STDERRorig)
    return result
  end
end

macro add_doc(newentry,obj)
  _curmod = current_module()
  _obj = :($_curmod.$obj)
  _olddoc = eval(:(@doc $_obj))
  return quote
    if isa($_olddoc, Vector)
      @ignore_warnings @doc [$newentry,$_olddoc...] $_obj  # making sure we promote to the correct type
    else
      @ignore_warnings @doc [$newentry,$_olddoc] $_obj     # making sure we promote to the correct type
    end
  end
end

end # module
