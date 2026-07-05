"""
    LGMIntegralsExt

Scaffold weakdep extension. SciML-style quadrature backend selection
(see `plans/dependencies.md`) lands with the posterior-expectation
milestone; this module exists so that any package in the tree may load
`Integrals` without Pkg trying to precompile a missing source file.
"""
module LGMIntegralsExt

end # module
