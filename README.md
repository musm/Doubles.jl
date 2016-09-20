# Doubles.jl

[![Build Status](https://travis-ci.org/musm/Doubles.jl.svg?branch=master)](https://travis-ci.org/musm/Doubles.jl)

[![Coverage Status](https://coveralls.io/repos/musm/Doubles.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/musm/Doubles.jl?branch=master)

[![codecov.io](http://codecov.io/github/musm/Doubles.jl/coverage.svg?branch=master)](http://codecov.io/github/musm/Doubles.jl?branch=master)


This package implements Double arithmetic in Julia. In particular this package focuses on an effecient implementation and high performance.

Double arithmetic is based on Dekker's technique from 1971. The basic idea is to store an unevaluated sum of a number in terms of its high and low bits.
