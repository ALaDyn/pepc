module rng_wrapper
   use iso_c_binding

   interface

      integer(c_int) function gen_norm_double_rng(ctr_seed, key_seed, output) bind(c, name="gen_norm_double_rng")
         import :: c_int, c_double
         integer(c_int), dimension(4), intent(in) :: key_seed
         integer(c_int), dimension(4), intent(inout) :: ctr_seed
         real(c_double), dimension(8), intent(inout) :: output
      end function gen_norm_double_rng

   end interface

end module rng_wrapper
