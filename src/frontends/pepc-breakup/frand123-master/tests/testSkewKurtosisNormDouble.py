import sys
import numpy as np
import scipy.stats as stats
import argparse

# parse arguments
parser = argparse.ArgumentParser(description='Carry out statistical tests for skew and kurtosis' )
parser.add_argument( '--polar', action = 'store_true', help = 'use values for polar version of Box-Muller transformation.' )
parser.add_argument( '--wichura', action = 'store_true', help = 'use values for inverse transformation sampling by Wichura.' )
parser.add_argument( '--mu', help = 'expectation. default: 0', default = 0, type = float )
parser.add_argument( '--sigma', help = 'variance. default: 1', default = 1, type = float )
args = parser.parse_args()

# check at least one is set
if( not ( args.polar or args.wichura ) ):
   print( 'No information on transformation used provided.' )
   sys.exit( 1 )

# Do we use polar Box-Muller's transformation?
if( args.polar ):
   # check that no other is set
   if( args.wichura ):
      print( 'Provide one transformation only.' )
      sys.exit( 1 )
   # set values
   transformation = 'Polar Box-Muller'
   pLimSkew = 0.27
   pLimKurtosis = 0.6
   pLimNormal = 0.5

# Do we use Wichura's transformation?
if( args.wichura ):
   # set values
   transformation = 'Wichura'
   pLimSkew = 0.35
   pLimKurtosis = 0.6
   pLimNormal = 0.58

# read data from file
f = open( 'tests/rand_norm_double.out', 'rb' )
rnd_var = np.fromfile( f, 'double' )

# run skewtest
res = stats.skewtest( rnd_var )

if( res.pvalue < pLimSkew ):
   print( 'Skew test for %s with mu = %e and sigma = %e failed.\np-value: %e, expected p-value: > %e' % ( transformation, args.mu, args.sigma, res.pvalue, pLimSkew ) )
   sys.exit( 1 )

# run kurtosistest
res = stats.kurtosistest( rnd_var )

if( res.pvalue < pLimKurtosis ):
   print( 'Kurtosis test for %s with mu = %e and sigma = %e failed.\np-value: %e, expected p-value> > %e' % ( transformation, args.mu, args.sigma, res.pvalue, pLimKurtosis ) )
   sys.exit( 1 )

# run normaltest
res = stats.normaltest( rnd_var )

if( res.pvalue < pLimNormal ):
   print( 'Normal test for %s with mu = %e and sigma = %e failed.\np-value: %e, expected p-value> > %e' % ( transformation, args.mu, args.sigma, res.pvalue, pLimNormal ) )
   sys.exit( 1 )
