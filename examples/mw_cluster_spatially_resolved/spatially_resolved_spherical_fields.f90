module FFTW3
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
end module FFTW3

module progress_bar
  implicit none
  
  logical, private, parameter :: showprogress            = .true.
  integer, private, parameter :: barlen                  = 100
  character, private, parameter :: progress_char         = '='
  character, private, parameter :: progress_char_special = '+'

  contains
  
  subroutine progress(j,k,special)  
    implicit none  
    integer, intent(in):: j,k  
    logical, intent(in), optional :: special
    character(len=18) :: barfmt = "(a1,I3,a,x,???a,$)"
    
    integer :: value
    integer, save :: lastvalue = -1
    
    character, save :: bar(0:barlen+1)
    character :: pchar
    
    if (.not. showprogress) return
    value = (barlen*j)/k
    if (value .eq. lastvalue) return
    lastvalue = value
    
    pchar = progress_char
    
    if (present(special)) then
      if (special) then
        pchar = progress_char_special
      endif
    endif
    
    if (value==0) then
      bar(0)        = "|"
      bar(1:barlen) = " "
      bar(barlen+1) = "|"
    else
      bar(value)        = pchar
    endif
    

    write(barfmt(12:14), '(I3.3)') barlen+2
    
    ! print the progress bar.  
    write(6, barfmt) char(13), (j*100)/k, "%", bar
    
  end subroutine progress  

end module



module module_data
  implicit none

  integer :: Na, Nt, Nt_max
  
  integer, parameter :: NUMCOMPONENTS = 4
  
  integer :: NR, NTheta, NPhi
  real*8 :: maxR, rcl
  
  ! indices: spatial component (iR,iT,iP), time/frequency index (t/w)
  real*8, allocatable, dimension(:,:,:,:,:) :: observable
  real*8, allocatable, dimension(:) :: abscissa ! will hold time/frequency values
  real*8, allocatable, dimension(:,:,:,:) :: grid
  
  contains
  
    subroutine free_data()
      implicit none
      
      if (allocated(observable))    deallocate(observable)
      if (allocated(abscissa))      deallocate(abscissa)
      if (allocated(grid))          deallocate(grid)
      
    end subroutine free_data
    
    
    subroutine load_grid(filename)
      use progress_bar
      implicit none
      character(*), intent(in) :: filename
      
      integer :: idx, iR, iTheta, iPhi
      real*8 :: r
      
      write(*,'("[STATUS] ", "load_grid")')

      open(87,file=trim(filename),status='old',position='rewind',action='read', FORM='unformatted')
      
      read(87) Nt_max, maxR, NR, NTheta, NPhi

      Na = (NR+1) * (NTheta+1) * (NPhi+1)

      write(*, '("# maxR     = ", g15.5)') maxR
      write(*, '("# NR       = ",   i15)') NR
      write(*, '("# NTheta   = ",   i15)') NTheta
      write(*, '("# NPhi     = ",   i15)') NPhi
      
      allocate(grid(0:NR, 0:NTheta, 0:NPhi, 1:3))
      
      idx = 0
      do iR = 0,NR
        do iTheta = 0,NTheta
          do iPhi = 0,NPhi
	    idx = idx + 1
            call progress(idx, Na)
	    read (87) grid(iR, iTheta, iPhi, 1:3)
          end do
        end do
      end do
 
      write(*,*) ! progress bar line break
      
      close(87)

    end subroutine
    
    
    logical function load_data(filename, starttime_fs)
      use progress_bar
      implicit none
      character(*), intent(in) :: filename
      real*8, intent(in) :: starttime_fs
      integer :: itime
      real*8 :: time_fs
      integer :: ios, iR, iTheta, iPhi
      integer :: datindex, lineindex
      real*8, allocatable :: dummy(:,:,:,:)

      write(*,'("[STATUS] ", "load_data")')

      open(87,file=trim(filename),status='old',position='rewind',action='read', FORM='unformatted')
      
      read(87) Nt_max, maxR, rcl, NR, NTheta, NPhi

      Na = (NR+1) * (NTheta+1) * (NPhi+1)

      write(*, '("# Nt_max   = ",   i15)') Nt_max
      write(*, '("# maxR     = ", g15.5)') maxR
      write(*, '("# rcl      = ", g15.5)') rcl
      write(*, '("# NR       = ",   i15)') NR
      write(*, '("# NTheta   = ",   i15)') NTheta
      write(*, '("# NPhi     = ",   i15)') NPhi
      write(*, '("# Na       = ",   i15)') Na
      
      allocate(observable(1:NUMCOMPONENTS, 0:NR, 0:NTheta, 0:NPhi, Nt_max))
      allocate(dummy(1:NUMCOMPONENTS, 0:NR, 0:NTheta, 0:NPhi))
      allocate(abscissa(1:Nt_max))
      
      datindex  = 0
      lineindex = 0
      
      do
        read(87, iostat=ios) itime, time_fs
	
	lineindex = lineindex + 1
	
	if (ios .ne. 0) exit
	
        do iR = 0,NR
          do iTheta = 0,NTheta
            do iPhi = 0,NPhi
              read (87) dummy(1:NUMCOMPONENTS, iR, iTheta, iPhi)
            end do
          end do
        end do
	
	if (time_fs >= starttime_fs) then
	  datindex = datindex + 1
            observable(1:NUMCOMPONENTS, 0:NR, 0:NTheta, 0:NPhi, datindex) = dummy(1:NUMCOMPONENTS, 0:NR, 0:NTheta, 0:NPhi)
	    abscissa(datindex) = time_fs
	endif
	
	call progress(lineindex, Nt_max, time_fs >= starttime_fs)

      end do
      
      write(*,*) ! progress bar line break

      close(87)
      
      Nt = datindex ! store number of actually available timesteps
      write(*, '("# Nt       = ",   i15)') Nt

      load_data = (Nt > 2)

      if (.not. load_data) then
        write(*,'(a,/,a)') 'Did not find a sufficient number of timesteps.', 'Reduce starttime_fs or check your data'
      endif
      
    end function load_data
    
    subroutine observable_spectrum_to_file(filename)
      use progress_bar
      implicit none
      character(*), intent(in) :: filename

      integer :: i, iR, iTheta, iPhi
 
      write(*,'("[STATUS] ", "observable_spectrum_to_file: ",a)') filename
 
      open(24,file=trim(filename),status='unknown',position='rewind',action='write')

      do i=1,Nt
        call progress(i, Nt)

        write(24,'(g15.5)', advance='no') abscissa(i)

        do iR = 0,NR
          do iTheta = 0,NTheta
            do iPhi = 0,NPhi
	      write(24,'(3(x,g15.5))', advance='no') observable(4, iR, iTheta, iPhi, i) ! TODO: for now we only dump the potential
            end do	
	  end do
        end do
	
	write(24,*) ! line break in output file
      end do

      close(24)
      
      write(*,*) ! progress bar line break
      
    end subroutine
    
    
    subroutine observable_write_to_vtk(fileprefix)
      use progress_bar
      implicit none
      character(*), intent(in) :: fileprefix
      
      character*256 :: filename
      integer :: i
      integer :: iR, iTheta, iPhi
      real*8 :: theta, phi, r
      real*8, parameter :: pi = 3.14159265358979323846_8
      integer, parameter :: VTK_HEXAHEDRON = 12
      
      write(*,'("[STATUS] ", "observable_write_to_vtk: ",a)') fileprefix

      do i = 1,Nt
        call progress(i, Nt)

        write(filename,'(a,"_",I7.7,".vtk")') fileprefix, i
	
        open(24,file=trim(filename),status='unknown',position='rewind',action='write')
	
	! write VTK header
	write(24,'("# vtk DataFile Version 2.0")')
	write(24,'("spherical fields")')
	write(24,'("ASCII")')
	! Grid definition
	write(24,'("DATASET UNSTRUCTURED_GRID")')
	! grid point coordinates
	write(24,'("POINTS ", I0, " double")') (NR+1)*(NTheta+1)*(NPhi+1)
        do iR = 0,NR
          do iTheta = 0,NTheta
            do iPhi = 0,NPhi
	      write(24,'(3(x,g25.12))') grid(iR, iTheta, iPhi, 1:3)
            end do	
	  end do
        end do
	  write(24,*)
	! grid cell types
	write(24,'("CELL_TYPES ", I0)') NR*NTheta*NPhi
        do iR = 0,NR-1
          do iTheta = 0,NTheta-1
            do iPhi = 0,NPhi-1
	      write(24,'(x,I0)', advance='no') VTK_HEXAHEDRON
            end do	
	  end do
        end do
	  write(24,*)
	! grid cell definitions
	write(24,'("CELLS ", I0,x,I0)') NR*NTheta*NPhi, (8+1)*NR*NTheta*NPhi
        do iR = 0,NR-1
          do iTheta = 0,NTheta-1
            do iPhi = 0,NPhi-1
	      write(24,'(9(x,I0))') 8, &
	                            idx(iR  , mod(iPhi+1, NPhi), iTheta+1) , &
	                            idx(iR  ,     iPhi         , iTheta+1) , &
	                            idx(iR+1,     iPhi         , iTheta+1) , &
	                            idx(iR+1, mod(iPhi+1, NPhi), iTheta+1) , &
	                            idx(iR  , mod(iPhi+1, NPhi), iTheta  ) , &
	                            idx(iR  ,     iPhi         , iTheta  ) , &
	                            idx(iR+1,     iPhi         , iTheta  ) , &
	                            idx(iR+1, mod(iPhi+1, NPhi), iTheta  )
            end do	
	  end do
        end do
	  write(24,*)
	! point data
	write(24,'("POINT_DATA ", I0)') (NR+1)*(NTheta+1)*(NPhi+1)
	! electrical field
	write(24,'("SCALARS e_field double 3")')
	write(24,'("LOOKUP_TABLE default")')

        do iR = 0,NR
          do iTheta = 0,NTheta
            do iPhi = 0,NPhi
	      write(24,'(3(x,g25.12))') observable(1:3, iR, iTheta, iPhi, i)
            end do	
	  end do
        end do
	
	write(24,*)

	! electrical potential
	write(24,'("SCALARS potential double 1")')
	write(24,'("LOOKUP_TABLE default")')

        do iR = 0,NR
          do iTheta = 0,NTheta
            do iPhi = 0,NPhi
	      write(24,'(g25.12)') observable(4, iR, iTheta, iPhi, i)
            end do	
	  end do
        end do
	
        close(24)
      end do

      write(*,*) ! progress bar line break
      
    contains
    
      integer function idx(iR, iPhi, iTheta)
        implicit none
	integer, intent(in) :: iR, iTheta, iPhi
	idx = iR*(NTheta+1)*(NPhi+1) + iTheta * (NPhi+1) + iPhi
      end function
      
    end subroutine
    
  
end module


module math
  implicit none

    real*8, parameter :: pi = 3.14159265358979323846_8
    
    integer, parameter :: MAXN = 16
    integer, parameter :: MAXL = 12
    
    ! from besselstuff.nb Mathematica-Worksheet
    real*8, parameter :: SphericalBesselJZero(0:MAXL, 1:MAXN) = reshape([ &
      3.14159265359_8, 4.49340945791_8, 5.76345919689_8, 6.98793200050_8, 8.18256145259_8, 9.35581211105_8, 10.5128354081_8, 11.6570321925_8, 12.7907817121_8, 13.9158226104_8, 15.0334693039_8, 16.1447429423_8, 17.2504547841_8, &
      6.28318530718_8, 7.72525183694_8, 9.09501133064_8, 10.4171185475_8, 11.7049071546_8, 12.9665301728_8, 14.2073924588_8, 15.4312892103_8, 16.6410028815_8, 17.8386431992_8, 19.0258535361_8, 20.2039426328_8, 21.3739721812_8, &
      9.42477796077_8, 10.9041216594_8, 12.3229409706_8, 13.6980231532_8, 15.0396647076_8, 16.3547096394_8, 17.6479748702_8, 18.9229991985_8, 20.1824707649_8, 21.4284869721_8, 22.6627206581_8, 23.8865307560_8, 25.1010385210_8, &
      12.5663706144_8, 14.0661939128_8, 15.5146030109_8, 16.9236212852_8, 18.3012559595_8, 19.6531521018_8, 20.9834630689_8, 22.2953480191_8, 23.5912748180_8, 24.8732139239_8, 26.1427676434_8, 27.4012592589_8, 28.6497962617_8, &
      15.7079632679_8, 17.2207552719_8, 18.6890363554_8, 20.1218061745_8, 21.5254177334_8, 22.9045506479_8, 24.2627680424_8, 25.6028559538_8, 26.9270407788_8, 28.2371343600_8, 29.5346341078_8, 30.8207940865_8, 32.0966767254_8, &
      18.8495559215_8, 20.3713029593_8, 21.8538742227_8, 23.3042469889_8, 24.7275655478_8, 26.1277501372_8, 27.5078683649_8, 28.8703733470_8, 30.2172627094_8, 31.5501883818_8, 32.8705345977_8, 34.1794746665_8, 35.4780131752_8, &
      21.9911485751_8, 23.5194524987_8, 25.0128032023_8, 26.4767636645_8, 27.9155761994_8, 29.3325625786_8, 30.7303807316_8, 32.1111962397_8, 33.4768008195_8, 34.8286965377_8, 36.1681571359_8, 37.4962736358_8, 38.8139888812_8, &
      25.1327412287_8, 26.6660542588_8, 28.1678297080_8, 29.6426045403_8, 31.0939332141_8, 32.5246612886_8, 33.9371083026_8, 35.3331941827_8, 36.7145291272_8, 38.0824790873_8, 39.4382144800_8, 40.7827470981_8, 42.1169585325_8, &
      28.2743338823_8, 29.8115987909_8, 31.3201417074_8, 32.8037323852_8, 34.2653900861_8, 35.7075769531_8, 37.1323317249_8, 38.5413648517_8, 39.9361278109_8, 41.3178646902_8, 42.6876512847_8, 44.0464252109_8, 45.3950094399_8, &
      31.4159265359_8, 32.9563890398_8, 34.4704883313_8, 35.9614058047_8, 37.4317367682_8, 38.8836309555_8, 40.3188925092_8, 41.7390528671_8, 43.1454250176_8, 44.5391446334_8, 45.9212017638_8, 47.2924656053_8, 48.6537041118_8, &
      34.5575191895_8, 36.1006222444_8, 37.6193657536_8, 39.1164701903_8, 40.5941896534_8, 42.0544164128_8, 43.4987571413_8, 44.9285896767_8, 46.3451060653_8, 47.7493457344_8, 49.1422214247_8, 50.5245397256_8, 51.8970175245_8, &
      37.6991118431_8, 39.2444323612_8, 40.7671158214_8, 42.2695149778_8, 43.7536054311_8, 45.2210650159_8, 46.6733329250_8, 48.1116545550_8, 49.5371160745_8, 50.9506714535_8, 52.3531638708_8, 53.7453428658_8, 55.1278782244_8, &
      40.8407044967_8, 42.3879135681_8, 43.9139818114_8, 45.4209639723_8, 46.9106054901_8, 48.3844038606_8, 49.8436551888_8, 51.2894900803_8, 52.7229017097_8, 54.1447680563_8, 55.5558697207_8, 56.9569043527_8, 58.3484984443_8, &
      43.9822971503_8, 45.5311340140_8, 47.0601416128_8, 48.5711298516_8, 50.0656518347_8, 51.5450520426_8, 53.0105034817_8, 54.4630367429_8, 55.9035630242_8, 57.3328925814_8, 58.7517496714_8, 60.1607847680_8, 61.5605846363_8, &
      47.1238898038_8, 48.6741442320_8, 50.2057283367_8, 51.7202484304_8, 53.2190952897_8, 54.7034825077_8, 56.1744764965_8, 57.6330202625_8, 59.0799524622_8, 60.5160228380_8, 61.9419048385_8, 63.3582060271_8, 64.7654767358_8, &
      50.2654824574_8, 51.8169824873_8, 53.3508435853_8, 54.8685009575_8, 56.3712071531_8, 57.8600629728_8, 59.3360419635_8, 60.8000100548_8, 62.2527414665_8, 63.6949317159_8, 65.1272083475_8, 66.5501398542_8, 67.9642431481_8  &
      ], shape(SphericalBesselJZero))
      
    real*8, parameter :: BesselJZero(0:MAXL, 1:MAXN) = reshape([ &
      2.40482555770_8, 3.83170597021_8, 5.13562230184_8, 6.38016189592_8, 7.58834243451_8, 8.77148381596_8, 9.93610952423_8, 11.0863700192_8, 12.2250922637_8, 13.3543004779_8, 14.4755006861_8, 15.5898478845_8, 16.6982499339_8, &
      5.52007811029_8, 7.01558667043_8, 8.41724414069_8, 9.76102312968_8, 11.0647094885_8, 12.3386041975_8, 13.5892901706_8, 14.8212687270_8, 16.0377741909_8, 17.2412203825_8, 18.4334636670_8, 19.6159669040_8, 20.7899063601_8, &
      8.65372791291_8, 10.1734681351_8, 11.6198411721_8, 13.0152007217_8, 14.3725366716_8, 15.7001740797_8, 17.0038196678_8, 18.2875828325_8, 19.5545364310_8, 20.8070477893_8, 22.0469853647_8, 23.2758537263_8, 24.4948850439_8, &
      11.7915344390_8, 13.3236919363_8, 14.7959517824_8, 16.2234661603_8, 17.6159660498_8, 18.9801338752_8, 20.3207892136_8, 21.6415410199_8, 22.9451731319_8, 24.2338852578_8, 25.5094505542_8, 26.7733225455_8, 28.0267099500_8, &
      14.9309177085_8, 16.4706300509_8, 17.9598194950_8, 19.4094152264_8, 20.8269329570_8, 22.2177998966_8, 23.5860844356_8, 24.9349278877_8, 26.2668146412_8, 27.5837489636_8, 28.8873750635_8, 30.1790611788_8, 31.4599600353_8, &
      18.0710639679_8, 19.6158585105_8, 21.1169970530_8, 22.5827295931_8, 24.0190195248_8, 25.4303411542_8, 26.8201519834_8, 28.1911884595_8, 29.5456596710_8, 30.8853789677_8, 32.2118561997_8, 33.5263640756_8, 34.8299869903_8, &
      21.2116366299_8, 22.7600843806_8, 24.2701123136_8, 25.7481666993_8, 27.1990877660_8, 28.6266183073_8, 30.0337223866_8, 31.4227941923_8, 32.7958000373_8, 34.1543779239_8, 35.4999092054_8, 36.8335713419_8, 38.1563775047_8, &
      24.3524715308_8, 25.9036720876_8, 27.4205735500_8, 28.9083507809_8, 30.3710076671_8, 31.8117167240_8, 33.2330417628_8, 34.6370893521_8, 36.0256150639_8, 37.4000999772_8, 38.7618070179_8, 40.1118232710_8, 41.4510923079_8, &
      27.4934791320_8, 29.0468285349_8, 30.5692044955_8, 32.0648524071_8, 33.5371377118_8, 34.9887812946_8, 36.4220196683_8, 37.8387173829_8, 39.2404479952_8, 40.6285537190_8, 42.0041902367_8, 43.3683609475_8, 44.7219435432_8, &
      30.6346064684_8, 32.1896799110_8, 33.7165195092_8, 35.2186707386_8, 36.6990011287_8, 38.1598685620_8, 39.6032394161_8, 41.0307736916_8, 42.4438877433_8, 43.8438014203_8, 45.2315741035_8, 46.6081326763_8, 47.9742935313_8, &
      33.7758202136_8, 35.3323075501_8, 36.8628565113_8, 38.3704724348_8, 39.8576273022_8, 41.3263832540_8, 42.7784816132_8, 44.2154085053_8, 45.6384441822_8, 47.0487007377_8, 48.4471513873_8, 49.8346535104_8, 51.2119670041_8, &
      36.9170983537_8, 38.4747662348_8, 40.0084467335_8, 41.5207196704_8, 43.0137377234_8, 44.4893191232_8, 45.9490159980_8, 47.3941657556_8, 48.8259303816_8, 50.2453269553_8, 51.6532516682_8, 53.0504989591_8, 54.4377769283_8, &
      40.0584257646_8, 41.6170942128_8, 43.1534537784_8, 44.6697431166_8, 46.1678535129_8, 47.6493998067_8, 49.1157737248_8, 50.5681846798_8, 52.0076914567_8, 53.4352271570_8, 54.8516190760_8, 56.2576047151_8, 57.6538448119_8, &
      43.1997917132_8, 44.7593189977_8, 46.2979966772_8, 47.8177856915_8, 49.3203606864_8, 50.8071652030_8, 52.2794539036_8, 53.7383253720_8, 55.1847479393_8, 56.6195802665_8, 58.0435879282_8, 59.4574569084_8, 60.8618046825_8, &
      46.3411883717_8, 47.9014608872_8, 49.4421641104_8, 50.9650299062_8, 52.4715513985_8, 53.9630265584_8, 55.4405920689_8, 56.9052499920_8, 58.3578890253_8, 59.7993016310_8, 61.2301979773_8, 62.6512173882_8, 64.0629378249_8, &
      49.4826098974_8, 51.0435351836_8, 52.5860235068_8, 54.1116155698_8, 55.6216509098_8, 57.1173027815_8, 58.5996056312_8, 60.0694769983_8, 61.5277351668_8, 62.9751135342_8, 64.4122724129_8, 65.8398088044_8, 67.2582645563_8  &      
      ], shape(BesselJZero))
  
  contains
    
    real*8 function factorial(n)
        implicit none
        integer, intent(in) :: n

        if (n<0) then
          write(*,*) "Tried to calculate factorial of negative argument."
	  stop
        end if

        if (n>170) then
          write(*,*) "Tried to calculate factorial with n>170. This would lead to numeric overflow."
	  stop
        end if

        select case (n)
          case ( 0)
            factorial =                            1._8
          case ( 1)
            factorial =                            1._8
          case ( 2)
            factorial =                            2._8
          case ( 3)
            factorial =                            6._8
          case ( 4)
            factorial =                           24._8
          case ( 5)
            factorial =                          120._8
          case ( 6)
            factorial =                          720._8
          case ( 7)
            factorial =                         5040._8
          case ( 8)
            factorial =                        40320._8
          case ( 9)
            factorial =                       362880._8
          case (10)
            factorial =                      3628800._8
          case (11)
            factorial =                     39916800._8
          case (12)
            factorial =                    479001600._8
          case (13)
            factorial =                   6227020800._8
          case (14)
            factorial =                  87178291200._8
          case (15)
            factorial =                1307674368000._8
          case (16)
            factorial =               20922789888000._8
          case (17)
            factorial =              355687428096000._8
          case (18)
            factorial =             6402373705728000._8
          case (19)
            factorial =           121645100408832000._8
          case (20)
            factorial =          2432902008176640000._8
          case (21)
            factorial =         51090942171709440000._8
          case (22)
            factorial =       1124000727777607680000._8
          case (23)
            factorial =      25852016738884976640000._8
          case (24)
            factorial =     620448401733239439360000._8
          case default
            factorial = gamma(real(n+1, 8))
        end select
    end function factorial   
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Computes the associated Legendre polynomial \f$P_{lm}(x)\f$.
    !> Here m and l are integers satisfying  \f$0 \leq m \leq l\f$,
    !> while x lies in the range \f$-1 \leq x \leq 1\f$.
    !>
    !> Code fragment for \f$P_l^m(x)\f$ taken from
    !>
    !> Numerical Recipes in Fortran 77: The Art of Scientific Computing
    !>              (ISBN 0-521-43064-X)
    !> pg. 246ff
    !>
    !> and modified to give \f$P_{lm}(x)\f$:
    !> \f$ P_{lm}(x) = (-1)^m P_l^m (x) \f$, see
    !>
    !> Abramowitz and Stegun: Handbook of Mathematical Functions
    !> Section 8. Legendre Functions (pg. 332)
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real*8 function LegendreP(l,m,x)
      implicit none
      integer, intent(in) :: l, m
      real*8 ::x

      integer :: i,ll
      real*8 :: fact,pll,pmm,pmmp1,somx2

      pll = 0._8

      if ( (m < 0) .or. (m > l) .or. (abs(x) > 1) .or.  (l<0)) then
        write(*,*) 'Invalid arguments for LegendreP(',l,m,x,')'
	stop
      endif

      pmm = 1._8     ! Compute P_m^m

      if (m > 0) then
        somx2 = sqrt((1._8-x)*(1._8+x))
        fact  = 1._8

        do i = 1,m
           pmm  = -pmm * fact * somx2
           fact = fact+2._8
        enddo
      endif

      if (l == m) then
        LegendreP = pmm
      else
        pmmp1 = x*(2._8*m+1._8)*pmm  ! Compute P_m+1^m

        if (l == m+1) then
          LegendreP = pmmp1
        else                  ! Compute P_l^m , l > m + 1
          do ll = m+2,l
            pll   = (x*(2._8*ll-1._8)*pmmp1-(ll+m-1._8)*pmm)/real(ll-m,kind=8)
            pmm   = pmmp1
            pmmp1 = pll
          enddo

          LegendreP = pll
        endif
      endif

      LegendreP = (-1._8)**m * LegendreP

    end function LegendreP
    
    
    ! see notes sheet [D] and besselstuff.nb
    real*8 function SphericalBesselJ(l,x)
      implicit none
      real*8, intent(in) :: x
      integer, intent(in) :: l
      
      integer :: p
      real*8 :: cx, sxx
      
      integer, parameter :: maxp = 6
      
      integer*8, parameter :: Ccoeffs(0:maxp, 0:12) = reshape([ &
                              0,           0,           0,        0,      0,  0, 0, &
                             -1,           0,           0,        0,      0,  0, 0, &
                             -3,           0,           0,        0,      0,  0, 0, &
                            -15,           1,           0,        0,      0,  0, 0, &
                           -105,          10,           0,        0,      0,  0, 0, &
                           -945,         105,          -1,        0,      0,  0, 0, &
                         -10395,        1260,         -21,        0,      0,  0, 0, &
                        -135135,       17325,        -378,        1,      0,  0, 0, &
                       -2027025,      270270,        6930,       36,      0,  0, 0, &
                      -34459425,     4729725,     -135135,      990,     -1,  0, 0, &
                     -654729075,    91891800,    -2837835,    25740,    -55,  0, 0, &
		   -13749310575,  1964187225,   -64324260,   675675,  -2145,  1, 0, &
  	          -316234143225, 45831035250, -1571349780, 18378360, -75075, 78, 0  & 
            ], shape(Ccoeffs))
      
      integer*8, parameter :: Scoeffs(0:maxp, 0:12) = reshape([ &
                              1,             0,           0,          0,       0,     0, 0, &
                              1,             0,           0,          0,       0,     0, 0, &
                              3,            -1,           0,          0,       0,     0, 0, &
                             15,            -6,           0,          0,       0,     0, 0, &
                            105,           -45,           1,          0,       0,     0, 0, &
                            945,          -420,          15,          0,       0,     0, 0, &
                          10395,         -4725,         210,         -1,       0,     0, 0, &
                         135135,        -62370,        3150,        -28,       0,     0, 0, &
                        2027025,       -945945,       51975,       -630,       1,     0, 0, &
                       34459425,     -16216200,      945945,     -13860,      45,     0, 0, &
                      654729075,    -310134825,    18918900,    -315315,    1485,    -1, 0, &
                    13749310575,   -6547290750,   413513100,   -7567560,   45045,   -66, 0, &
                   316234143225, -151242416325,  9820936125, -192972780, 1351350, -3003, 1  &
            ], shape(Scoeffs))

      if (x==0._8) then
        if (l==0) then
	  sphericalBesselJ = 1._8
	else
	  sphericalBesselJ = 0._8
	endif
      else
	cx   = cos(x)
	sxx  = sin(x)/x
      
        sphericalBesselJ = 0.
	
	do p=0,maxp
	  sphericalBesselJ = sphericalBesselJ + x**(2*p-l) * ( real(Ccoeffs(p,l),8) * cx +  real(Scoeffs(p,l),8) * sxx )
	end do
	
      endif
      
    end function
    
    
    real*8 function Rfunc_Raitza(n,l,a,r)
      implicit none
      integer, intent(in) :: n, l
      real*8, intent(in) :: r, a
      
      real*8 :: xnl, NNnl
      
      xnl          = SphericalBesselJZero(l,n)
      NNnl         = (a*a*a)/2.*SphericalBesselJ(l+1,xnl)
      Rfunc_Raitza = NNnl * SphericalBesselJ(l, xnl/a * r)
      
    end function
    
    
    real*8 function Rfunc_Wang(n,l,a,r)
      implicit none
      integer, intent(in) :: n, l
      real*8, intent(in) :: r, a
      
      real*8 :: xnl, NNnl
      
      xnl        = BesselJZero(l,n)
      NNnl       = a*a/2. * ( bessel_jn(l+1, xnl) )**2
      Rfunc_Wang = bessel_jn(l, xnl/a * r) / sqrt( NNnl)
      
    end function    
    
end module


module spherical_fourier
  implicit none
  
    logical, parameter :: use_Raitza_definition = .false. ! if false, we use Wang's definition
  
    real*8, allocatable, dimension(:,:,:) :: Rtilda
    real*8, allocatable, dimension(:,:)   :: MM
    real*8, allocatable, dimension(:,:,:) :: P
    real*8, allocatable, dimension(:,:)   :: E
  
  contains
  
    subroutine spherical_fourier_decomposition(freqindex, Scnlm_w)
      use module_data
      implicit none
      integer, intent(in) :: freqindex
      real*8, intent(out) :: Scnlm_w(1:NUMCOMPONENTS, 1:NR/2, 0:NTheta/2, 0:NTheta/2)
      integer :: n,l,m
      integer :: iR, iTheta, iPhi, ic
      
      real*8 :: tmp
      
      Scnlm_w(:, :, :, :) = 0._8
      
        do n = 1,NR/2
          do l = 0,NTheta/2
	    do m=0,NTheta/2
	  
	      if (m <= l) then
                do iR = 0,NR
                  do iTheta = 0,NTheta
                    do iPhi = 0,NPhi
		    
		      tmp = Rtilda(n, l, iR) * MM(l, m) * P(l, m, iTheta) * E(m, iPhi)
		    
                      do ic = 1,NUMCOMPONENTS
		        Scnlm_w(ic, n, l, m) = Scnlm_w(ic, n, l, m) + observable(ic, iR, iTheta, iPhi, freqindex) * tmp
                      end do

		    end do
		  end do
	        end do
	      
	      endif

	    end do
	  end do
        end do
      
      write(*,*) freqindex, abscissa(freqindex), Scnlm_w(1,:,:,:)
      
    end subroutine
    

    subroutine spherical_fourier_decomposition_init()
      use module_data
      use math
      implicit none
      integer :: n, l, m, iR, iTheta, iPhi
      real*8 :: r, r2, rv(3), costheta, Phi
      
      write(*,'("[STATUS] ", "spherical_fourier_decomposition_init")')
      
      allocate(Rtilda(1:NR/2, 0:NTheta/2, 0:NR))
      allocate(MM(0:Ntheta/2, 0:Ntheta/2))
      allocate(P(0:NTheta/2,  0:Ntheta/2, 0:NTheta))
      allocate(E(0:Ntheta/2,  0:NPhi))
      
      ! initialize Rtilda
      do iR=0,NR
	rv = grid(iR, NTheta, NPhi, 1:3)
        r2 = dot_product(rv, rv)
	r  = sqrt(r2)

        do n=1,NR/2
	  do l=0,NTheta/2
            if (use_Raitza_definition) then
	      Rtilda(n, l, iR) = Rfunc_Raitza(n,l,maxR,r) * r2
	    else
	      Rtilda(n, l, iR) = Rfunc_Wang(n,l,maxR,r) * r2
	    endif
	  end do
	end do
      end do
      
      ! initialize M
      do l=0,NTheta/2
        do m=0,NTheta/2
	  if (m <= l) then
    	    MM(l, m) = sqrt( (2._8*l+1._8)/(4._8*pi) * factorial(l-m)/factorial(l+m) )
	  else
	    MM(l, m) = 0._8
	  endif
	end do
      end do
      
      ! initialize P
      do iTheta=0,NTheta
	rv = grid(NR, iTheta, NPhi, 1:3)
        r2 = dot_product(rv, rv)
	r  = sqrt(r2)
	
	if (r > 0._8) then
  	  costheta = rv(3)/r
	else
	  costheta = 1._8
	endif
	
	do l=0,NTheta/2
	  do m=0,NTheta/2
	    if (m <= l) then
	      P(l, m, iTheta) = LegendreP(l,m,costheta) * sqrt(1-costheta*costheta)
	    else
	      P(l, m, iTheta) = 0._8
	    endif
	  end do
	end do
      end do

      ! initialize E
      do iPhi=0,NPhi
	rv = grid(NR, NTheta, iPhi, 1:3)

        if ((rv(1) == 0._8) .and. (rv(2) == 0._8)) then
          Phi = 0._8
        else
         Phi = atan2(rv(2), rv(1))
        end if
      
        do m=0,NTheta/2
	  E(m, iPhi) = cos(m*Phi)
	end do
      end do        
      
      
      
    end subroutine
    
    
    subroutine spherical_fourier_decomposition_free()
      implicit none
      
      if (allocated(Rtilda)) deallocate(Rtilda)
      if (allocated(MM))     deallocate(MM)
      if (allocated(P))      deallocate(P)
      if (allocated(E))      deallocate(E)
    
    end subroutine
    
    
    subroutine spherical_fourier_decomposition_for_all(component, filename)
      use progress_bar
      use module_data
      implicit none
      real*8 :: Scnlm_w(1:NUMCOMPONENTS, 1:NR/2, 0:NTheta/2, 0:NTheta/2)
      integer, intent(in) :: component
      character(*), intent(in) :: filename
      character*11 :: formatstring = '????(g15.5)'
      
      integer :: i
      
      call spherical_fourier_decomposition_init()


      write(*,'("[STATUS] ", "spherical_fourier_decomposition_for_all")')
      
      write(formatstring(1:4),'(I4.4)') NR/2*(NTheta/2+1)*(NTheta/2+1) + 1


      open(24,file=trim(filename),status='unknown',position='rewind',action='write')

      do i=1,Nt/2 ! second half of spectrum is boring anyway
        call progress(i, Nt/2)
        call spherical_fourier_decomposition(i, Scnlm_w)
        write(24,formatstring, advance='no') abscissa(i), Scnlm_w(component, :, :, :)
      end do
      
      close(24)

      call spherical_fourier_decomposition_free()
      
      write(*,*) ! progress bar line break

    end subroutine 

end module



module computations
  use module_data
  use progress_bar
  implicit none
  
  contains
  
  subroutine observable_fourier_forward()
    use, intrinsic :: iso_c_binding
    use FFTW3
    implicit none
    
    type(C_PTR) :: fftw_plan
    complex(C_DOUBLE_COMPLEX), dimension(Nt) :: tmparray
    integer :: c, iR, iTheta, iPhi
    integer :: nstep, nstep_max
    
    write(*,'("[STATUS] ", "observable_fourier_forward")')

    fftw_plan = fftw_plan_dft_1d(size(tmparray), tmparray, tmparray, FFTW_BACKWARD, FFTW_ESTIMATE)

    nstep = 0
    nstep_max = (NR+1) * (NTheta+1) * (NPhi+1)
    
    do iR = 0,NR
      do iTheta = 0,NTheta
        do iPhi = 0,NPhi
	  nstep = nstep + 1
	  call progress(nstep, nstep_max)
	  
          do c=1,NUMCOMPONENTS
            tmparray(1:Nt) = observable(c,iR, iTheta, iPhi, 1:Nt)
     
            ! perform FFT
            call fftw_execute_dft(fftw_plan, tmparray, tmparray)
      
            observable(c,iR, iTheta, iPhi, 1:Nt) = abs(tmparray(1:Nt))
	  end do
	end do
      end do
    end do
    
    call fftw_destroy_plan(fftw_plan)
    
    write(*,*) ! progress bar line break

  end subroutine


  subroutine abscissa_convert_to_frequencies()
    implicit none
    real*8 :: delta_t, wmax, delta_w
    integer :: i
    real*8, parameter :: pi = 3.14159265358979323846_8
    
    write(*,'("[STATUS] ", "abscissa_convert_to_frequencies")')

    delta_t = abscissa(2)-abscissa(1)
    wmax = 2*pi*1/delta_t
    delta_w = wmax/Nt
    
    abscissa(1:Nt)        = [((i-1)*delta_w,i=1,Nt)]
    abscissa(Nt+1:Nt_max) = 0.
    
  end subroutine

  
    
  subroutine analysis_workflow(filename_in, filename_grid, starttime_fs)
    use spherical_fourier
    implicit none
    character(*), intent(in) :: filename_in, filename_grid
    real*8, intent(in) :: starttime_fs
    
    if (load_data(filename_in, starttime_fs)) then
      call load_grid(filename_grid)
      !call observable_write_to_vtk('fields/spherical_fields')
      call observable_fourier_forward()
      call abscissa_convert_to_frequencies()
      !call observable_write_to_vtk('fields/spherical_fields_fft')
      call observable_spectrum_to_file('./spectrum.dat') ! TODO: rename this file
      call spherical_fourier_decomposition_for_all(4,'./sperical_fourier_coeffs_phi.dat') ! TODO: rename this file
    endif

    call free_data()

  end subroutine
  
end module


program spatially_resolved_spherical_fields
  use module_data
  use computations
  implicit none
  integer :: argc, i
  character*256 :: dirname_in
  character*256 :: filename_in, filename_grid
  
  argc = command_argument_count()
  
  call system('mkdir -p fields')

  if (argc < 1) then
    write(*,*) 'Call with file to be processed as first argument. Exiting.'
    stop
  endif

  do i=1,argc
    call get_command_argument(i, dirname_in)
    
    filename_in   = trim(dirname_in)//'/field_spherical.dat'
    filename_grid = trim(dirname_in)//'/field_spherical.dat_grid.dat'
    
    write(*,'(2/"Reading data from directory ",a)') trim(dirname_in)

    call analysis_workflow(trim(filename_in), trim(filename_grid), 100.0_8)
  end do

   

end program
