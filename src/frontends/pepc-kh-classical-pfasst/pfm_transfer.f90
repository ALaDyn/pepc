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

   !> Interpolate coarse (qGp) to fine (qFp) nodes
   subroutine interpolate(qFp, qGp, levelF, ctxF, levelG, ctxG, t)
      type(c_ptr), intent(in), value :: qFp, qGp, ctxF, ctxG
      integer,     intent(in)        :: levelF, levelG
      real(pfdp),  intent(in) :: t

      type(app_data_t), pointer :: qF, qG

      call pepc_status('|--> interpolate()')

      call c_f_pointer(qFp,qF)
      call c_f_pointer(qGp,qG)

      ! FIXME: we just copy here, i.e. no coarsening in space for now
      qF%x(:,:) = qG%x
      qF%v(:,:) = qG%v

   end subroutine interpolate


   !> Restriction from qFp  to qGp
   subroutine restrict(qFp, qGp, levelF, ctxF, levelG, ctxG, t)
      type(c_ptr), intent(in), value :: qFp, qGp, ctxF, ctxG
      integer,     intent(in)        :: levelF, levelG
      real(pfdp),  intent(in) :: t

      type(app_data_t), pointer :: qF, qG

      call pepc_status('|--> restrict()')

      call c_f_pointer(qFp,qF)
      call c_f_pointer(qGp,qG)

      ! FIXME: we just copy here, i.e. no coarsening in space for now
      qG%x(:,:) = qF%x
      qG%v(:,:) = qF%v

   end subroutine restrict


end module pfm_transfer

