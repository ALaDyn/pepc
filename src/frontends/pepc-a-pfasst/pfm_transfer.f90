!> Transfer (interpolate, restrict) routines.                                                                           
module pfm_transfer
  use module_debug
  use iso_c_binding
  use pfm_encap
  implicit none
  private
  save
  
  public interpolate
  public restrict

contains


   subroutine interpolate(qFp, qGp, levelF, ctxF, levelG, ctxG)
      type(c_ptr), intent(in), value :: qFp, qGp, ctxF, ctxG
      integer,     intent(in)        :: levelF, levelG

      type(app_data_t), pointer :: qF, qG

      call pepc_status('|--> interpolate()')

      call c_f_pointer(qFp,qF)
      call c_f_pointer(qGp,qG)

      ! TODO: Interpolate coarse (qG) to fine (qF) nodes (or copy if nothing has to be done)

   end subroutine interpolate


  subroutine restrict(qFp, qGp, levelF, ctxF, levelG, ctxG)
      type(c_ptr), intent(in), value :: qFp, qGp, ctxF, ctxG
      integer,     intent(in)        :: levelF, levelG

      type(app_data_t), pointer :: qF, qG

      call pepc_status('|--> restrict()')

      call c_f_pointer(qFp,qF)
      call c_f_pointer(qGp,qG)

      ! TODO: restriction by point injection form qF to qG (or nothing if nothing has to be done)

   end subroutine restrict


end module pfm_transfer

