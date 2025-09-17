# Halton.m
This function generates Halton sequences for multidimensional SMLE. It follows Train (2003, Chapter 9.3.3-9.3.5).
Some rudimentary checks are included. Nevertheless be careful with the primes, dimensions, and settings that you choose.

The function returns the results as (N x dimensions x draws)-arrays.

```yaml
Function:
[H, Z] = halton(N, dimensions, draws, PROP, VAL, ...)

Returns:
H - Halton draws
Z - corresponding values from a standard normal distribution

Required inputs:
N - number of observational units for which Halton draws are needed 
dimensions - dimensions of the integral
draws - number of draws

Optional inputs:
'prime' - vector of primes, default are primes in ascending order.
'burn' - specify the number of points at the beginning of the sequence to be dropped/burnt. Default are 50 (the first point, zero, will always be dropped).
'leap' - specify the points to miss out between returned points. Default is no leap.
'random' - Set to 1 if you want randomized Halton draws a la Bhat (2003). It is recommended to set a SEED before running halton.m with the 'random' option.
'scramble' - Set to 1 if you want a scrambled Halton set. Recommended for high-dimensional sets. (For details see http://mathworks.com/help/stats/qrandset.scramble.html)

References:
See http://www.jyu.fi/ersa2003/cdrom/papers/406.pdf on correlations of Halton sequences for larger primes.
Train, K. 2003. Discrete Choice Methods with Simulation. Cambridge: Cambridge University Press.
```
