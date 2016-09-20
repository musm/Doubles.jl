# Doubles.jl

[![Build Status](https://travis-ci.org/musm/Doubles.jl.svg?branch=master)](https://travis-ci.org/musm/Doubles.jl)

[![Coverage Status](https://coveralls.io/repos/musm/Doubles.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/musm/Doubles.jl?branch=master)

[![codecov.io](http://codecov.io/github/musm/Doubles.jl/coverage.svg?branch=master)](http://codecov.io/github/musm/Doubles.jl?branch=master)


This package implements Double and TwoFold arithmetic in Julia. In particular' this package's main goal is efficiency and performance.

Double arithmetic is based on Dekker's technique from 1971. The basic idea is to store an unevaluated sum of a number in terms of its high and low bits.

Twofold arithmetic is similar to Double arithmetic and is also based on Dekker's ideas. Twofolds share much of the same advantages of Dekker arithmetic, but have a different emphasis; the goal of TwoFold arithmetic is to monitor error, thus the normalization step is omitted in twofold arithmetic since it is irrelevant to the goal of tracking errors since it eliminates useful information about inaccuracies in calculations.

The main advantage of double style arithmetic and their most important feature is speed and accuracy. While they cannot fully simulate quad precision (the range is restricted to the corresponding float type they are based upon), they are accurate up to roughly 32 decimal places and are fast, compared to quad precision which is typically 100x slower than binary64 type. The arithmetic in this package is typically no worse than 10x of the corresponding float type. Since they store two float types, they are limited to 2x the speed of the corresponding float type.
