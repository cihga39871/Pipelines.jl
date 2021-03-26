# Pipelines.jl

*A lightweight Julia package for computational pipelines.*

Building reusable programs and pipelines is easier than you have ever thought.

| **Documentation**                                                               |
|:-------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://cihga39871.github.io/Pipelines.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://cihga39871.github.io/Pipelines.jl/dev) |

## Package Features

- Easy to build both simple and complex tasks.

- Supports external command lines and pure Julia functions.

- Supports **resuming** interrupted tasks, **skipping** finished tasks.

- Supports dependency check.

## Installation

Pipelines.jl can be installed using the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run

```julia
pkg> add Pipelines
# If it fails, use
pkg> add https://github.com/cihga39871/Pipelines.jl.git
```

To use the package, type

```julia
using Pipelines
```

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **documentation of the most recently tagged version.**
- [**DEVEL**][docs-dev-url] &mdash; *documentation of the in-development version.*
