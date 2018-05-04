function testCentralMomentsSingle()
   % compare central moments of the random numbers against analytical central moments

   % read data
   f = fopen( 'tests/rand_single.out' );
   a = fread( f, Inf, 'single' );
   fclose( f );

   % comute the mean
   mu = mean( a );

   % shift a by mean (building block of central moments)
   a = a - mu;

   % initialize b to 1 / length( a ) to reduce computational cost
   b = ones( size( a ) ) / length( a );

   % iterate over the central moments
   for i = 1:75
      % calculate b = a^i
      b = b .* a;

      % compute analytic central moment (2^(-i-1)((-1)^i+1)/(i+1) is the analytical value for the i-th central moment)
      analyticalCentralMoment = ( 2^(-i-1) * ( (-1)^i + 1 ) ) / ( i + 1 );

      % compute numic central moment
      numericCentralMoment = sum( b );

      % compare even and odd central moments as odd central moments should be zero
      if( mod( i, 2 ) == 0 )
         % compute difference between the analytical and numeric central moment
         diffMoment = analyticalCentralMoment - numericCentralMoment;

         % compute scaled error
         scaledError = diffMoment / analyticalCentralMoment;

         if( scaledError > 1e-3 )
            disp( sprintf( 'Failure in testCentralMomentsSingle' ) )
            disp( sprintf( '%d-th central moment: scaledError = %e', i, scaledError ) )
            exit( 1 )
         end
      else
         if( abs( numericCentralMoment > 1e-6 ) )
            disp( sprintf( 'Failure in testCentralMomentsSingle' ) )
            disp( sprintf( '%d-th central moments too far off 0, numerical central moment: %e'. i, numericCentralMoment ) )
            exit( 1 )
         end
      end
   end
   disp( 'testCentralMomentsSingle passed' )
   exit( 0 )
end
