module module_cmdline
  implicit none

  logical :: do_resume, input_specified, do_convert
  character(len=255) :: resume_file,input_file, convert_file
  integer :: resume_step, convert_step

  contains

!==============================================================
    subroutine read_args()
      implicit none

      character(len=255) :: arg,arg_val

      logical :: file_exists, arg_is_val

      integer :: i,j

      resume_file=""
      input_file=""

      do_convert=.false.
      do_resume=.false.
      file_exists=.false.
      arg_is_val=.false.
      input_specified=.false.

      if (command_argument_count()==0) then
        write (*,*)
        write (*,*) "********************************************"
        write (*,*) "***** No Command line arguments given. *****"
        write (*,*) "****** At least input (-i) is needed. ******"
        write (*,*) "********************************************"
        write (*,*)
        stop
      end if

      do i = 1, command_argument_count()
        if (arg_is_val) then
          arg_is_val=.false.
          cycle
        end if

        call get_command_argument(i, arg)
        select case (arg)
          case ('-h', '--help')
            call print_help()
            stop
          case ('-c', '--convert')
            do_convert = .true.
            call get_command_argument(i+1, arg_val)
            inquire(FILE=arg_val, EXIST=file_exists)
            if (file_exists) then
                convert_file = trim(arg_val)
                do j=1,255
                  if (convert_file(j:j+3)==".mpi") then
                    read(convert_file(j-6:j),"(i6)")convert_step
                    exit
                  end if
                end do
                file_exists=.false.
                arg_is_val=.true.                 !set so that next cycle will be skipped
            else
              write (*,*)
              write (*,*) "********************************************"
              write (*,*) "*** Chosen file to convert does not exist ***"
              write (*,*) "********************************************"
              write (*,*)
              stop
            end if
          case ('-r', '--resume')
            do_resume = .true.
            call get_command_argument(i+1, arg_val)
            inquire(FILE=arg_val, EXIST=file_exists)
            if (file_exists) then
              resume_file=trim(arg_val)
              do j=1,255
                if (resume_file(j:j+3)==".mpi") then
                  read(resume_file(j-6:j),"(i6)")resume_step
                  exit
                end if
              end do
              file_exists=.false.
              arg_is_val=.true.                 !set so that next cycle will be skipped
            else
              write (*,*)
              write (*,*) "********************************************"
              write (*,*) "*** Chosen file to resume does not exist ***"
              write (*,*) "********************************************"
              write (*,*)
              stop
            end if
          case ('-i', '--input')
            call get_command_argument(i+1, arg_val)
            inquire(FILE=arg_val, EXIST=file_exists)
            if (file_exists) then
              input_file=trim(arg_val)
              input_specified=.true.
              file_exists=.false.
              arg_is_val=.true.                 !set so that next cycle will be skipped
            else
              write (*,*)
              write (*,*) "*******************************************"
              write (*,*) "***** Chosen inputfile does not exist *****"
              write (*,*) "*******************************************"
              write (*,*)
              stop
            end if
          case default
            write (*,*)
            write (*,*) "*******************************************"
            write (*,'(a,a)') 'Unrecognized command-line option: ', arg
            write (*,*) "*******************************************"
            write (*,*)
            call print_help()
            stop
        end select
      end do

      if ((input_specified .eqv. .false.) .and. (do_convert .eqv. .false.)) then
        write (*,*)
        write (*,*) "********************************************"
        write (*,*) "** Inputfile (-i) is needed to run pepc-f **"
        write (*,*) "********************************************"
        write (*,*)
        stop
      end if


    end subroutine read_args

!==============================================================

    subroutine print_help()
      print '(a)', 'usage: pepc-f -i inputfile [OPTIONS]'
      print '(a)', ''
      print '(a)', 'pepc-f needs at least an input file. Other arguments are optional'
      print '(a)', ''
      print '(a)', 'pepc-f options:'
      print '(a)', ''
      print '(a)', '  -c(--convert) convertfile convert checkpoint to npy format'
      print '(a)', '  -i(--input) inputfile     start run with specified inputfile'
      print '(a)', '  -r(--resume) resumefile   resume run from specified checkpoint'
      print '(a)', '                            the .mpi file has to be selected'
      print '(a)', '  -h, --help                print usage information and exit'
    end subroutine print_help

!==============================================================
end module module_cmdline
