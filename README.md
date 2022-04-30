# ZeissMicroscopyFormat

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaIO.github.io/ZeissMicroscopyFormat.jl/stable)
[![Build Status](https://github.com/JuliaIO/ZeissMicroscopyFormat.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaIO/ZeissMicroscopyFormat.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaIO/ZeissMicroscopyFormat.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaIO/ZeissMicroscopyFormat.jl)

This package provides support for the [CZI image format for microscopy](https://www.zeiss.com/microscopy/int/products/microscope-software/zen/czi.html) in the Julia programming language. It is registered with the [FileIO](https://github.com/JuliaIO/FileIO.jl) package.

# Installation

Within Julia, perform the following:

```julia
julia> using Pkg

julia> Pkg.add(["ZeissMicroscopyFormat", "FileIO"])
```

This would add the package to whatever environment is currently active. Generally it is recommended to use [project-specific environments](https://pkgdocs.julialang.org/v1/environments/). More information about package management can be found in the [Pkg documentation](https://pkgdocs.julialang.org/v1/).

# Usage

In an environment with these packages installed, usage is as simple as

```julia
julia> using FileIO

julia> img = load("my_image_file.czi")
```

The format will be automatically detected by FileIO, which will load this package to import the specified file.

A general tool for viewing multidimensional images is [ImageView](https://github.com/JuliaImages/ImageView.jl).
