# Doubles.jl

[![Build Status](https://travis-ci.org/musm/Doubles.jl.svg?branch=master)](https://travis-ci.org/musm/Doubles.jl)

[![Coverage Status](https://coveralls.io/repos/musm/Doubles.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/musm/Doubles.jl?branch=master)

[![codecov.io](http://codecov.io/github/musm/Doubles.jl/coverage.svg?branch=master)](http://codecov.io/github/musm/Doubles.jl?branch=master)


This package implements Double arithmetic in Julia. In particular this package focuses on an efficient implementation and high performance.

Double arithmetic is based on Dekker's technique from 1971. The basic idea is to store an unevaluated sum of a number in terms of its high and low bits.

This package is inspired by [`DoubleDouble.jl`](https://github.com/simonbyrne/DoubleDouble.jl). See notes on differences on why this package was created.

# Usage
The only thing that you need to do to use this package is `x = Double(2.0)` and multiplication, division, sqrt, abs, addition, subtraction, negation should all work.

# Notes
Why is this package efficient compared to  other double  arithmetic packages? In other packages operations between Doubles and Float types may promote the Float type to a Double with a zero in the low bit field. This is wasteful and leads to an unnecessary number of calculations. In this package we are smart and promote float types to a wrapper type (ie. `Single`) which has specialized methods for all the implemented arithmetic operations to avoid the unnecessary calculations if the low bit field was non zero. To make this work without having to define methods for float types directly (only on the wrapper type `Single`), requires a clever hack on julia's promotion and conversion machinery.

In other words the user can simply write `1.0 * Double(1.0)` and this will properly call the method that operates on `Single  * Single`. Similarly, `1.0 * Double(1.0,1e-14)` will also properly get promoted and call the specialized method `Single * Double` .