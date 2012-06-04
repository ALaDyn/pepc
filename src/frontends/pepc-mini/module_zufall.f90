!!!!!!!!!!!!!!!!!!!!
!! zufallszahlen module
!!
!! Enthaelt Methoden zur Erzeugung von verschieden verteilten Zufallszahlen
!!!!!!!!!!!!!!!!!!!!

MODULE zufall

    implicit none

    CONTAINS

!======================================================================================
!Füllt "list" mit gauss-verteilten zufallszahlen 
    SUBROUTINE random_gauss_list(list,mu,sigma)
        implicit none
  
        real*8, intent(inout) :: list(:)
        real*8  :: v(2), pi, r, p, mu, sigma
        integer :: n, i
    
        pi = 2.0_8*acos(0.0_8)
        n  = size(list)

        DO i=1, n, 2
            call random_number(v)
            r = sqrt(-2.0_8 * log(v(1)))
            p = 2.0_8*pi*v(2)
            list(i)                = r * sin(p)*sigma+mu
            if((i+1)<=n) list(i+1) = r * cos(p)*sigma+mu
        END DO
  
  END SUBROUTINE

!======================================================================================
!Füllt "list" mit Zufallszahlen die gemäß f(v)=m/T * v * exp(-m/(2T)*v**2) verteilt sind, vtherm=sqrt(T/m)
    SUBROUTINE random_gaussian_flux_list(list,vtherm)
        implicit none
  
        real*8, intent(inout) :: list(:)
        real*8  :: vtherm,y
        integer :: n,i
    
        n  = size(list)

        DO i=1, n
            call random_number(y)
            list(i)=vtherm*sqrt(-2.*log(1-y))          
        END DO
  
    END SUBROUTINE


!======================================================================================
!Erzeugt Zufallszahlen die gemäß f(v)=m/T * v * exp(-m/(2T)*v**2) verteilt sind, vtherm=sqrt(T/m)
    SUBROUTINE random_gaussian_flux(x,vtherm)
        implicit none
  
        real*8, intent(inout) :: x
        real*8  :: vtherm,y
                
        call random_number(y)
        x=vtherm*sqrt(-2.*log(1-y))          
          
    END SUBROUTINE

END MODULE
