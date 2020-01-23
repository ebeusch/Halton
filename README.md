# Halton.m
This function generates Halton sequences for multidimensional SMLE. It follows Train (2003, Chapter 9.3.3-9.3.5).
Some rudimentary checks are included. Nevertheless be careful with the primes, dimensions, and settings that you choose.

The function returns the results as (N x dimensions x draws)-arrays.

```yaml
Function:
[H, S] = halton(N, dimensions, draws, PROP, VAL, ...)

Returns:
H - Halton draws
S - corresponding values from a standard normal distribution

Required inputs:
N - number of observational units for which Halton draws are needed 
dimensions - dimensions of the integral
draws - number of draws

Optional inputs:
'prnum' - vector of primes, default are primes in ascending order.
'burn' - specify the number of points at the beginning of the sequence to be dropped/burnt. Default are 50 (the first point, zero, will always be dropped).
'leap' - specify the points to miss out between returned points. Default is no leap.
'random' - Set to 1 if you want randomized Halton draws a la Bhat (2003). It is recommended to set a SEED before running halton.m with the 'random' option.
'scramble' - Set to 1 if you want a scrambled Halton set. Recommended for high-dimensional sets. (For details see http://mathworks.com/help/stats/qrandset.scramble.html)

References:
See http://www.jyu.fi/ersa2003/cdrom/papers/406.pdf on correlations of Halton sequences for larger primes.
Train, K. 2003. Discrete Choice Methods with Simulation. Cambridge: Cambridge University Press.
```

### Matlab Code
```matlab
function [H, S] = halton(N, dimensions, draws, varargin)

%% Optional input arguments
opts = inputParser;
    opt1 = 'prnum';
    val1 = 0;
        addParameter(opts,opt1,val1,@isnumeric);
    opt2 = 'burn';
    val2 = 50;
        addParameter(opts,opt2,val2,@isnumeric);
    opt3 = 'leap';
    val3 = 0;
        addParameter(opts,opt3,val3,@isnumeric);
    opt4 = 'random';
    val4 = 0;
        addParameter(opts,opt4,val4,@isnumeric);
    opt5 = 'scramble';
    val5 = 0;
        addParameter(opts,opt5,val5,@isnumeric);
parse(opts, varargin{:});
    prnum = opts.Results.prnum;
    burn = opts.Results.burn;
    leap = opts.Results.leap;
    randhalt = opts.Results.random;
    toscramble = opts.Results.scramble;

%% Defining the dimensions to be used
if prnum == 0
     p = dimensions;  
else
    if size(prnum,1) > size(prnum,2)
        prnum = prnum';
    end
    %% some checks for dimensions and primes
    if dimensions ~= size(prnum,2)
        error('Dimensions do not match number of primes supplied')
    end
    if min(isprime(prnum)) == 0
        error('Non-prime was supplied')
    end
    %% set dimensions for Matlab's haltonset fct based on max prime
    p = size(primes(max(prnum)+1),2);
end
    if dimensions == 1 && (prnum==0 || prnum==2)
        error('For dimension==1 a prime > 2 has to be supplied')
    end

%% Some warnings
if dimensions >= 7 && toscramble == 0
    warning('With high dimensional Halton draws scrambling is recommended')
end
if (burn - max(prnum) <=10) || (size(prnum,2)==1 && dimensions>=13)
    warning('The default burn setting might be too short for your primes')
end

%% Halton sequences
if leap == 0
    h = haltonset(p,'Skip',burn+1); % sequence starts with zero; no leaps
else
    h = haltonset(p,'Skip',burn+1,'Leap',leap);
end
 % Scramble
if toscramble == 1
    h = scramble(h,'RR2');
end
hp = net(h,N*draws);
    if prnum ~= 0
        [~, primi] = ismember(prnum, primes(max(prnum)+1));
        hp = hp(: , primi);
    end
    
%% Randomization a la Bhat (see Train (2003, page 264))
if randhalt == 1
    mu = rand(1,dimensions);
        mu = repmat(mu,N*draws,1);
    hp = hp + mu;
        c = hp > 1;
        hp = hp - c;
end

%% Reshaping
H = reshape(hp,draws,N,dimensions); % NxRxJ 
H = permute(H,[2,3,1]); % NxJxR

%% Standard normal correspondence
if nargout >= 2
    S = norminv(H);
end
```
