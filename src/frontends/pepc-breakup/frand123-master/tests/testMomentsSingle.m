function testMomentsSingle()
   % compare moments of the random numbers against analytical moments
   % view deviation of more than 1e-4 as error

   % read data
   f = fopen( 'tests/rand_single.out' );
   a = fread( f, Inf, 'single' );
   fclose( f );

   % initialize b to 1 / length( a ) to reduce computational cost
   b = ones( size( a ) ) / length( a );

   % iterate over the moments
   for i = 1:75
      % calculate b = a^i
      b = b .* a;

      % compute analytical moment
      analyticalMoment = 1 / ( i + 1 );

      % compare difference (1/(i+1) is the analytical value for the i-th moment
      relativeError = ( analyticalMoment - sum( b ) ) / analyticalMoment;

      if( abs( relativeError ) > 2e-3 )
         disp( sprintf( 'Failure un testMomentsSingle' ) )
         disp( sprintf( 'relative error for %d-th moment is larger than 1e-3.', i ) )
         exit( 1 )
      end
   end
   disp( 'testMomentsSingle passed' )
   exit( 0 )
end
