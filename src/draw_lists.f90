
! ======================
!
!   DRAW_LISTS
!
!   Draw selected particle lists
!
! ======================

subroutine draw_lists


  use treevars

  implicit none

  integer*8 :: key_list(nppm)

  integer, dimension(nppm) ::  addr_list, level_list, node_list, owner_list       ! list data

  character(15) :: cfile
  character(1) :: csnap
  character(3) :: cme,cip
  character(7), parameter :: colors(0:9) = (/"orange ", "cyan   ", "magenta", "blue   ", "green  ", &
       "red    ","yellow ","grey20 ","violet ","brown  "/) 


  integer :: i, ip, j, ilev, isnap, ix, iy, nbits
  real :: s, xt, yt
  logical :: write_keys=.false.
  integer :: key2addr        ! Mapping function to get hash table address from key


!  do ip=1,npp,npp/2
  ip=1
  !  Header file 
  cme = achar(me/100+48) // achar(me/10+48) // achar(mod(me,10)+48)  ! Convert 3-digit PE number into character string
  cip = achar(ip/100+48) // achar(ip/10+48) // achar(mod(ip,10)+48)  ! Convert 3-digit particle number into character string
  cfile="list_"//cme//cip//".gle"
  write (ipefile,'(a)') cfile
 open(60,file=cfile)

  !  initialise graphics filter

!  write (60,'(a,2(/a),2(/a,2f13.4))') &
!       'size 18 18' &
!       , 'set font rm' &
!       , 'set lwidth 0.001 lstyle 1' &
!       , 'begin scale ',15./boxsize,15./boxsize &
!       , 'begin translate ',.01*boxsize,.01*yl

!  write (60,'(a/a,2f13.4)') 'amove 0 0','box ',boxsize,yl


  write (60,'(a)') 'set color black hei .3'
  write (60,'(a)')  'psize=0.4' 
  write (60,'(a,2f13.4)') 'amove ',x(ip),y(ip)
  write (60,'(a,f10.5,2a)') 'marker otimes '
  write (60,'(a,2f13.4)') 'rmove ',.1,.1

  write (60,'(a,i6)') 'text ',pelabel(ip)         ! write out label


  !  Now do selected particle lists

  nlist = nterm(ip)

  ! get keys of twig nodes from hash table
  key_list(1:nlist) = intlist(1:nlist,ip)
  addr_list(1:nlist) = (/ ( key2addr( key_list(i) ),i=1,nlist) /)   !  Table address

  ! get levels of twigs
  level_list(1:nlist) = log( 1.*key_list(1:nlist) )/log(8.)
  node_list(1:nlist) =   htable( addr_list(1:nlist))%node   !  nodes
  owner_list(1:nlist) =   htable( addr_list(1:nlist))%owner      ! owner


  write (60,'(a)') 'set lwidth .01 hei .02'

  do j=1,nlist
    write (60,'(a,a7)') 'set color ',colors( mod(owner_list(j),10) )
     ! get box coords from keys
     nbits = level_list(j)
     ix = SUM( (/ (2**i*ibits( key_list(j),idim*i,1 ), i=0,nbits-1) /) )
     iy = SUM( (/ (2**i*ibits( key_list(j),idim*i+1,1 ), i=0,nbits-1) /) )

     s = boxsize/2**(level_list(j))          !  box length
     xt=ix*s + xmin
     yt=iy*s + ymin
     write (60,'(a,2f13.4)') 'amove ',xt,yt
     write (60,'(a,2f13.4)') 'box ',s,s

     write (60,'(a,2f13.4)') 'amove ',xcoc( node_list(j) ),ycoc( node_list(j) )        ! Centre of charge of twig node
     write (60,'(a,f10.5,2a)') 'circle ',(abs_charge( node_list(j)))**(.33), &
     		'*psize fill ',colors( mod(owner_list(j),10) )
    ! write (60,'(a,f10.5,2a)') 'circle ',psize,' fill ',colors( mod(owner_list(j),10) )
     !
     write (60,'(a)') 'set color black'
     write (60,'(a,2f13.4)') 'amove ',xt+s/2,yt+s/2
     if (write_keys) write (60,'(a,i6)') 'text ',key_list(j)  ! write out key

  end do

!  write (60,'(2(/a))') 'end scale ','end translate'

  close(60)
! end do
end subroutine draw_lists
