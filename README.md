[![Build Status](https://travis-ci.com/mewertd2/KVectors.jl.svg?branch=master)](https://travis-ci.com/mewertd2/KVectors.jl)
[![codecov.io](https://codecov.io/github/mewertd2/KVectors.jl/coverage.svg?branch=master)](https://codecov.io/github/mewertd2/KVectors.jl?branch=master)

# KVectors

The `KVectors` [Julia](http://julialang.org) package defines the `KVector` Type
to represent linear combinations of [k-blades](https://en.wikipedia.org/wiki/Blade_(geometry)).

All Blades in a given KVector have the same grade.

This is the main object in [Grassmann Algebra](https://en.wikipedia.org/wiki/Exterior_algebra).

KVectors essentially extends the algebras and Types defined in [Blades](https://github.com/mewertd2/Blades.jl) with the `+` operator.  This allows for a vector space of Blades.  

Most operators that act on Blades can act on KVectors.  Notable exceptions are the inner product `⋅` and geometric product `*`.
These can not, in general, operate on KVectors as they can result in a mixed grade vector.  The KVectors algebra is not closed under such operators.  For mixed grades you need [Multivectors](https://github.com/mewertd2/Multivectors.jl)

See the documentation of [Blades](https://github.com/mewertd2/Blades.jl) for more information.

## Project Information

### Contributing

Please read [CONTRIBUTING.md](./CONTRIBUTING.md) for details.

### Authors

* **Michael Alexander Ewert** - Developer - [Digital Domain](https://digitaldomain.com)

### License

This project is licensed under a modified Apache 2.0 license - see the [LICENSE](./LICENSE) file for details
