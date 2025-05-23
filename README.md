# SDSHFringeMitigation

<!-- [![Build Status](https://github.com/jasper9000/SDSHFringeMitigation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jasper9000/SDSHFringeMitigation.jl/actions/workflows/CI.yml?query=branch%3Amain) -->

The code and data accompanying the publication

> J. Riebesehl, D. C. Nak and D. Zibar, "Interference Fringe Mitigation in Short-Delay Self-Heterodyne Laser Phase Noise Measurements" (2025)


# Installation steps

1. Install the Julia programming language. Preferably, install the latest version via [juliaup](https://github.com/JuliaLang/juliaup).


2. Download this repository. Launch a terminal and navigate to the root of directory of *SDSHFringeMitigation*.

3. Launch julia in the terminal in this directory.

4. Install Pluto:

    `julia> import Pkg; Pkg.add("Pluto")`

6. Launch the notebook of choice using Pluto:
   
    `julia> using Pluto`

    `julia> Pluto.run(notebook="examples/example_simulation.jl")`

