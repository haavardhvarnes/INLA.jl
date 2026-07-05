"""
    LGMHCubatureExt

Scaffold weakdep extension. Adaptive cubature for posterior expectations
(see `plans/dependencies.md`) lands with the posterior-expectation
milestone; this module exists so that any package in the tree may load
`HCubature` without Pkg trying to precompile a missing source file.
"""
module LGMHCubatureExt

end # module
