# Structural equality helper for `LatentGaussianModel` comparison in
# tests. Core LGM does not define `Base.:(==)` on `LatentGaussianModel`
# / components / mappings (default falls back to `===`), so a recursive
# field-walking helper is needed to assert that the macro expansion
# produces a model identical to the hand-written form.
#
# When LGM core gains proper `Base.:(==)` overloads (likely PR-3 or a
# follow-up Tier-1 PR), this helper can be removed in favour of `==`
# directly.

function _struct_isequal(a, b)
    typeof(a) === typeof(b) || return false
    if a isa AbstractArray
        return a == b
    elseif a isa Tuple
        length(a) == length(b) || return false
        for (p, q) in zip(a, b)
            _struct_isequal(p, q) || return false
        end
        return true
    elseif a isa Number || a isa AbstractString || a isa Symbol || a isa Bool
        return a == b
    elseif a isa Function || a isa Module
        return a === b
    elseif isstructtype(typeof(a))
        for f in fieldnames(typeof(a))
            _struct_isequal(getfield(a, f), getfield(b, f)) || return false
        end
        return true
    else
        return a == b
    end
end
