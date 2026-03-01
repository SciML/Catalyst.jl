# SymLatexWrapper metadata type — controls LaTeX rendering of multi-character symbol names.
# Not exported by Symbolics, so we access it via option_to_metadata_type and wrap
# all usage behind accessor functions.
const _SymLatexWrapper = Symbolics.option_to_metadata_type(Val{:latexwrapper}())

"""
    _has_latex_wrapper(s)

Return `true` if the symbolic variable `s` has a custom LaTeX wrapper set.
"""
_has_latex_wrapper(s) = hasmetadata(unwrap(s), _SymLatexWrapper)

"""
    _get_latex_wrapper(s)

Return the LaTeX wrapper function for the symbolic variable `s`.
"""
_get_latex_wrapper(s) = getmetadata(unwrap(s), _SymLatexWrapper)

"""
    _set_latex_wrapper(s, f)

Set the LaTeX wrapper function for the symbolic variable `s` to `f`.
"""
_set_latex_wrapper(s, f) = setmetadata(s, _SymLatexWrapper, f)

"""
    _add_latex_wrapper(s)

Add `SymLatexWrapper` metadata (set to `string`) so that multi-character names render
as plain math italic instead of `\\mathtt{}` in LaTeX. No-op if already set.
"""
function _add_latex_wrapper(s)
    _has_latex_wrapper(s) && return s
    _set_latex_wrapper(s, string)
end
