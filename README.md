# Pipelines.jl

*A lightweight Julia package for computational pipelines.*

Building reusable pipelines and workflows is easier than you have ever thought.

| **Documentation**                                                               |
|:-------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://cihga39871.github.io/Pipelines.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://cihga39871.github.io/Pipelines.jl/dev) |

## Package Features

- Easy to build both simple and complex tasks.

- Supports external command lines and pure Julia functions.

- Supports **resuming** interrupted tasks, **retrying** failed tasks, and **skipping** finished tasks.

- Supports dependency check.

- Supports inputs, outputs validation, and so on.

- Supports program queuing and workload management with [JobSchedulers.jl](https://github.com/cihga39871/JobSchedulers.jl)

## Installation

Pipelines.jl can be installed using the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run

```julia
pkg> add Pipelines
```

To use the package, type

```julia
using Pipelines
```

## Documentation


- [**STABLE**](https://cihga39871.github.io/Pipelines.jl/stable) &mdash; **documentation of the most recently tagged version.**
- [**DEVEL**](https://cihga39871.github.io/Pipelines.jl/dev) &mdash; *documentation of the in-development version.*
