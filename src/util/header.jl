using Printf: @sprintf
using Dates: now

"""
Default header for writing wannier90 files.
"""
default_header() = @sprintf "Created by WannierIO.jl %s" string(now())
