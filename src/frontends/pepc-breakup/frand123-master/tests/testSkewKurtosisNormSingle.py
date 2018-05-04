import sys
import numpy as np
import scipy.stats as stats
import argparse

# parse arguments
parser = argparse.ArgumentParser(description='Carry out statistical tests for skew and kurtosis' )
parser.add_argument( '--threefry', action = 'store_true', help = 'use values for Threefry RNG' )
parser.add_argument( '--ars', action = 'store_true', help = 'use values for ARS RNG' )
parser.add_argument( '--mu', help = 'expectation. default: 0', default = 0, type = float )
parser.add_argument( '--sigma', help = 'variance. default: 1', default = 1, type = float )
args = parser.parse_args()

# check at least one RNG is set
if( not ( args.threefry or args.ars ) ):
   print( 'No information on RNG used provided.' )
   sys.exit( 1 )

# check only one RNG is set
if( args.threefry and args.ars ):
   print( 'Please choose only one RNG' )
   sys.exit( 1 )

# Do we use threefry?
if( args.threefry ):
   pLimSkew = 0.5
   pLimKurtosis = 0.5
   pLimNormal = 0.78
   RNG = 'Threefry'
# use ARS
else:
   pLimSkew = 0.6
   pLimKurtosis = 0.8
   pLimNormal = 0.8
   RNG = 'ARS'

# read data from file
f = open( 'tests/rand_norm_single.out', 'rb' )
rnd_var = np.fromfile( f, 'single' )

# run skewtest
res = stats.skewtest( rnd_var )

if( res.pvalue < pLimSkew ):
   print( 'Skew test for %s with mu = %e and sigma = %e failed.\np-value: %e, expected p-value: > %e' % ( RNG, args.mu, args.sigma, res.pvalue, pLimSkew ) )
   sys.exit( 1 )

# run kurtosistest
res = stats.kurtosistest( rnd_var )

if( res.pvalue < pLimKurtosis ):
   print( 'Kurtosis test for %s with mu = %e and sigma = %e failed.\np-value: %e, expected p-value> > %e' % ( RNG, args.mu, args.sigma, res.pvalue, pLimKurtosis ) )
   sys.exit( 1 )

# run normaltest
res = stats.normaltest( rnd_var )

if( res.pvalue < pLimNormal ):
   print( 'Normal test for %s with mu = %e and sigma = %e failed.\np-value: %e, expected p-value> > %e' % ( RNG, args.mu, args.sigma, res.pvalue, pLimNormal ) )
   sys.exit( 1 )
