#######################################################################################
"""
The following code is taken from the MIT licensed https://github.com/tkf/DisplayAs.jl

Copyright (c) 2019 Takafumi Arakaki <aka.tkf@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy of this
software and associated documentation files (the "Software"), to deal in the Software
without restriction, including without limitation the rights to use, copy, modify, merge,
publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons
to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.

------------------------------------------------------------------------------------
In Catalyst this code is considered **internal**, and solely used for forcing a given image
type in documentation due to limitations in Documenter.jl and Plots.jl. It should never be
used for any other purpose or exported.
"""

# _showables = [
#     (:Text, "text/plain")
#     (:MD, "text/markdown")
#     (:HTML, "text/html")
#     (:JSON, "application/json")
#     (:SVG, "image/svg+xml")
#     (:PNG, "image/png")
#     (:GIF, "image/gif")
#     (:PDF, "application/pdf")
#     (:EPS, "application/eps")
#     (:JPEG, "image/jpeg")
#     (:PS, "application/postscript")
#     (:LaTeX, "text/latex")
#     (:CSV, "text/csv")
#     (:TSV, "text/tab-separated-values")
# ]

# currently we only use PNG conversion, but others can be added from above as needed
_showables = [(:PNG, "image/png")]

struct Showable{mime <: MIME}
    content::Any
    options::NamedTuple
end

function Showable{T}(content; options...) where {T <: MIME}
    @nospecialize
    return Showable{T}(content, (; options...))
end

for (_, mime) in _showables
    MIMEType = typeof(MIME(mime))
    @eval Base.show(io::IO, ::$MIMEType, s::Showable{>:$MIMEType}; options...) = show(
        io, $MIMEType(), s.content; s.options..., options...)
end

macro mime_str(s)
    :(Showable{MIME{Symbol($s)}})
end

for (name, mime) in _showables
    @eval const $name = @mime_str $mime
end

#######################################################################################
