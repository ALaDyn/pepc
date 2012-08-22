!!!!!!!!!!!!!!!!!!!!
!! Diagnose module
!!
!! Enthaelt Methoden fuer die Diagnose des Plasmas
!!!!!!!!!!!!!!!!!!!!

!erster test

MODULE diagnostics


    use variables

    implicit none

    CONTAINS

!===============================================================================

    subroutine write_wallpotential(p)

        implicit none 

        type(t_particle), allocatable, intent(in) :: p(:)
        integer                                   :: n,ip

        n=size(p)

        DO ip=1,n
            IF (p(ip)%label<0) THEN                  !only wall particles            
                write(my_rank,'(i8,4es8.2)')p(ip)%label,p(ip)%x,p(ip)%results%pot
            END IF
        END DO 

    end subroutine 


!===============================================================================

    subroutine write_velocities(p)

        implicit none 

        type(t_particle), allocatable, intent(in) :: p(:)
        integer                                   :: n,ip

        n=size(p)

        open(20, file='test.ion', action='write')
        open(21, file='test.electron', action='write')
        DO ip=1,n
            IF (p(ip)%label>0) THEN                  !only particles  
                IF (p(ip)%data%q>0) THEN             !ions
                    write(20,*)p(ip)%label,p(ip)%x,p(ip)%data%v,p(ip)%data%q,p(ip)%data%m
                ELSE
                    write(20,*)p(ip)%label,p(ip)%x,p(ip)%data%v,p(ip)%data%q,p(ip)%data%m
                END IF     
            END IF
        END DO 
        close(20)
        close(21)

    end subroutine 




END MODULE
