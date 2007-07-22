! ==============================================
!
!                SPLIT_DUMP
!
!  Read particle distribution from dump file and
!  split it into Ncpu sub-directories
!
! ==============================================

program split_dump


  implicit none

  integer :: i, ipe, idummy, ierr, me

  character(30) :: cinfile, cdump, cdumpfile, cpu_infile, cpu_dumpfile
  character(9) :: ct
  character(11) :: cme
  integer :: ner, nir, np_beamr, npartr, iconf, iens, start_step
  real, allocatable :: x(:),y(:),z(:),ux(:),uy(:),uz(:),q(:),m(:), &
       ax(:),ay(:),az(:),phi(:)
  integer, allocatable ::  pid(:),label(:)
  real :: epsr, thetar, xlr, ylr, zlr, boxr, qe, qi, mass_e, mass_i
  real :: vplas, aplas, qplas, trun, tlaser, Zion, a_ii, eps
  real :: omegar, lambdar, plasma_center(1:3), focus(1:3)
  real :: axdum, aydum, azdum,phidum, bdum
  integer :: ioffset, i1, i2, npp_partial, npp_total, ipass, me_read, j, nrest, nadd
  integer :: ncpu, nmax, npp, timestamp


  write(6,*) 'Give # cpus to split data set by:'
  read(5,*) ncpu
  write(6,*) 'Give timestamp:'
  read(5,*) timestamp
  ! dump file 

  ! get destination filename suffix from dump counter
  do i=0,4
     cdump(6-i:6-i) =  achar(mod(timestamp/10**i,10) + 48)  
  end do
  cdump(1:1) = achar(timestamp/10**5 + 48)

  cinfile="dumps/parts_info."//cdump(1:6)
  open(80,file=cinfile)

  !  Dump file to read from
  cdumpfile="dumps/parts_dump."//cdump(1:6)
  open (81,file=cdumpfile) 

  ! Root reads info block in run directory to check that run parameters correct

  read(80,'(7(9x,i8/),10(9x,f12.5/),9(9x,1pe12.5/),2(9x,3f12.5/))')  &    ! info block - skip variable names
       start_step, npartr, &
       ner, nir, np_beamr, iconf, iens, &
       xlr, ylr, zlr, boxr, &
       epsr, thetar, &
       tlaser, trun, omegar, lambdar, &
       qe, qi, mass_e, mass_i, Zion, a_ii, Vplas, Aplas, Qplas, &
       plasma_center(1:3), focus(1:3)
  close(80)

  ! echo info block
  write(*,'(7(a12,i12/),9(a12,f12.5/),9(a12,1pe12.5/),2(a12,3(1pe12.5)/))')  &   
       'Start time: ', start_step, & 
       '# particles: ', npartr, &
       '# electrons: ', ner, &
       '# ions: ', nir, & 
       '# beam: ', np_beamr, &
       'Geom  : ', iconf, & 
       'Scheme: ', iens, &
       'Box_x: ',xlr,&
       'Box_y: ', ylr, &
       'Box_z: ', zlr, &
       'eps: ', eps,  &
       'theta: ', thetar, &
       'tlaser: ',tlaser, &
       'trun: ',trun, &
       'omega:',omegar, &
       'lambda:', lambdar, &
       'qe:',qe, &
       'qi:',qi, &
       'mass_e:',mass_e, &
       'mass_i:',mass_i, &
       'Zion:',Zion, &
       'a_ii:', a_ii, &
       'Vplas: ',Vplas, &
       'Aplas: ',Aplas, &
       'Qplas: ',Qplas, &
       'center: ',plasma_center(1:3), &
       'focus: ',focus(1:3)



  write (*,*) 'Splitting sets by 1:',ncpu

  npp_total = ner+nir  
  npp = npp_total/ncpu

  nrest = mod(npp_total,ncpu)  ! Remainder

  nmax = npp+nrest
  allocate (x(nmax),y(nmax),z(nmax),ux(nmax),uy(nmax),uz(nmax),ax(nmax),ay(nmax),az(nmax), &
       q(nmax), m(nmax), pid(nmax), label(nmax), phi(nmax) )



  ! Filename (directory) prefix
  do me=0,ncpu-1
     cme = "data/pe" &   
          // achar(me/1000+48) &
          // achar(mod(me/100,10)+48) &
          // achar(mod(me/10,10)+48) &
          // achar(mod(me,10)+48)  ! Convert 4-digit PE number into character string


     cpu_infile=cme//"/parts_info."//cdump(1:6)
     cpu_dumpfile=cme//"/parts_dump."//cdump(1:6)

     open(60,file=cpu_infile)
     open(61,file=cpu_dumpfile)


     if (mod(me,ncpu)==ncpu-1) then
        nadd = nrest  ! Put remainder on last PE
     else
        nadd = 0
     endif

     write(*,*) 'CPU ',me,': Copying ',npp+nadd,' particles out of ',npp_total,' to ',cme



     !  Read next npp particles 

     do i=1,npp+nadd
        read(81,*) x(i),y(i),z(i),ux(i),uy(i),uz(i),q(i),m(i), &
             ax(i),ay(i),az(i),phi(i), pid(i),label(i)
     end do

     ! Info block
     write(60,'(7(a9,i8/),10(a9,f12.5/),9(a9,1pe12.5/),2(a9,3(1pe12.5)/))')  &   
          'current_step=',start_step, 'npp=',npp, &
          'ne=',ner, 'ni=',nir, 'npbeam=',np_beamr, 'geometry=', iconf, &
          'scheme=',iens, &
          'xl=',xlr, 'yl=',ylr, 'zl=',zlr, 'boxsize=',zlr, &
          'eps=', epsr, 'theta=',thetar,' tlaser= ',tlaser,' trun= ',trun, &
          'omega=',omegar,'lambda=',lambdar,'  qe=',qe,'  qi=',qi, &
          'mass_e=',mass_e,'mass_i=',mass_i,'Zion=',Zion,'a_ii=',a_ii, &
          'Vplas=',Vplas,'Aplas=',Aplas,'Qplas=',Qplas, &
          'center=',plasma_center(1:3),'focus=',focus(1:3)

     do i=1,npp+nadd
        write(61,'((12(1pe14.5),2i9))')  &
             x(i),y(i),z(i),ux(i),uy(i),uz(i),q(i),m(i), &
             ax(i),ay(i),az(i),phi(i), pid(i),label(i)
     end do

     close(60)
     close(61)

  end do




end program split_dump



