function testRandSingle()
   pkg load statistics
   f = fopen( './tests/rand_single.out' );
   a = fread( f, Inf, 'single' );
   fclose( f );
   if( runstest( a, 'Alpha', 0.00001, 'Method', 'exact' ) == 0 )
      disp('Test passed' )
      exit(0)
   else
      disp('Test tests/testRandSingle.x failed' )
      exit(1)
   end
end
