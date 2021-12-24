module arguments 
    use kind_m
    use mod_mpg
    implicit none
    private
    public :: get_option
contains
    !------------------------------------------------------------------
    subroutine get_option(mpg, fn_in)
        type (mpeg_parameters), intent(   out) :: mpg
        character (len = *)   , intent(in out) :: fn_in
        integer   (kind = 4) :: narg, iarg, istatus
        character (len = 40) :: buffer
        character (len =  6) :: fmt
        buffer = ''
        iarg = 0
        narg = command_argument_count()
        if (narg == 0) then 
           call print_option()
           stop
        end if 
        do
            iarg = iarg + 1
            if (iarg > narg) call option_error( trim(buffer) )
            call get_command_argument(iarg, buffer)
            if (buffer(1:1) /= '-') exit  
            select case(trim(buffer))
            case ('-b') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(i', istatus, ')' 
                read(buffer, fmt) mpg%ibit_rate
                if (mpg%ibit_rate < 1 .or. mpg%ibit_rate > 14) call option_error( trim(buffer) )  
            case ('-crc') 
                mpg%icrc       = 0 ! crc16 on 
            case ('-cut') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(i', istatus, ')' 
                read(buffer, fmt) icut
                if (icut < 0 .or. icut > 32) call option_error( trim(buffer) )
            case ('-athmin') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(f', istatus, '.0)' 
                read(buffer, fmt) ath_min
            case ('-athmax') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(f', istatus, '.0)' 
                read(buffer, fmt) ath_max
            case ('-switch') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(f', istatus, '.0)' 
                read(buffer, fmt) switch
            case ('-xsm') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(f', istatus, '.0)' 
                read(buffer, fmt) xsm
            case ('-xms') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(f', istatus, '.0)' 
                read(buffer, fmt) xms
            case ('-offset') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(f', istatus, '.0)' 
                read(buffer, fmt) offset_l
            case ('-pow') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(f', istatus, '.0)' 
                read(buffer, fmt) pow
            case ('-r0') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(f', istatus, '.0)' 
                read(buffer, fmt) r0
            case ('-r1') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(f', istatus, '.0)' 
                read(buffer, fmt) r1
            case ('-skip') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(f', istatus, '.0)' 
                read(buffer, fmt) skip
                if (skip < 1.0_kd) then
                    write(*, *) 'input out of range: skip >= 1.0'
                    call option_error( trim(buffer) )
                end if           
            case ('-cuth')
                cut_factor = 0.0_kd
            case ('-longfft')
                q_short_fft = .false.
            case ('-l') 
                mblock_type_param = 0
            case ('-s') 
                q_sm = .false.
                mblock_type_param = 20
            case ('-m') 
                q_sm = .false.
                mblock_type_param = 21
            case ('-sm') 
                q_sm = .true.
            case ('-c') 
                mpg%icopyright = 1 ! copyrigt on
            case ('-o') 
                mpg%ioriginal  = 1 ! original on
            case ('-ms')
                qms_stereo = .true.
            case ('-ns')
                qms_stereo = .false.
            case ('-nomask')
                q_mask = .false.
            case ('-inner')
                q_inner = .true.
            case ('-helper')
                q_help = .true.
            case ('-noalias')
                q_alias = .false.
            case ('-noinfo')
                q_info = .false.
            case ('-about')
                call print_about()
            case default
                call option_error(trim(buffer))
            end select
        end do
        fn_in = trim(buffer)
    end subroutine get_option
    !---------------------------------------------------------------
    subroutine option_error(text)
        character (len = *), intent(in) :: text
        call print_option()
        write(*, *)
        write(*, '(3a)') ' >>>>>> error near option >>>>>>', text, '<<<<<<'
        stop
    end subroutine option_error
    !---------------------------------------------------------------
    subroutine print_option()
        write(*, *) 'usage : uzura3 -options file_name '
        write(*, *) '      : file_name.wav -> file_name.mp3'
        write(*, *) 'option:-b  1..14    bitrate for cbr mode                 (default 9 : 128kbps)'
        write(*, *) '       -crc         crc16 error protection on            (default off)'
        write(*, *) '       -c           copyright flag on                    (default off)'
        write(*, *) '       -o           original  flag on                    (default off)'
        write(*, *) '       -cut 1..32   band cut-off : place after -b option (default 32: 22.05khz)'
        write(*, *) '       -l           long-block-only                      (default off)'
        write(*, *) '       -s           short-mode for short-block           (default off)'
        write(*, *) '       -m           mixed-mode for short-block           (default off)'
        write(*, *) '       -sm          short & mixed-mode for short-block   (default on )'
        write(*, *) '       -xsm xx      short / mixed switching parameter    (default 2.0)'
        write(*, *) '       -switch  xx  long/short switching parameter       (default 2.0)'
        write(*, *) '       -skip    xx  speeds up outer loop                 (default 10.0)'
        write(*, *) '       -ms/-ns      stereo mode (ms/normal)              (default ms)'
        write(*, *) '       -xms xx      ms/ns stereo switching parameter     (default 0.5)'
        write(*, *) '       -nomask      masking off                          (default on)'
        write(*, *) '       -athmin   xx minimum of ath  [ db ]               (default -119.0)'
        write(*, *) '       -athmax   xx ceiling of ath  [ db ]               (default 0.0)'
        write(*, *) '       -offset   xx offset for spreading function [ db ] (default -38.0)'
        write(*, *) '       -pow      xx non-linear sum parameter             (default 0.3)'
        write(*, *) '       -longfft     use long fft for short/mixed blocks  (default off)'
        write(*, *) '       -inner       inner loop only                      (default off)'
        write(*, *) '       -noalias     anti-alias for mixed-block off       (default on)'
        write(*, *) '       -debug       print debug info                     (default on) '
        write(*, *) '       -about       about uzura3 '
        write(*, *) '        '
        write(*, *) 'example cbr 128kbps crc on : uzura3 -crc   file_name '
        write(*, *) 'example short-block normal stereo  : uzura3 -s -ns file_name '
    end subroutine print_option
    !---------------------------------------------------------------
    subroutine print_about()
        implicit none
        write(*, *) ' ******* UZURA3 is an MPEG-I/layer-III encoder written in Fortran95. ********'
        write(*, *) ' **** this program comes with absolutely no warranty. (c) h.o. 2000-2004 **** '
        write(*, *)
        write(*, *)
        write(*, *) ' http://members.tripod.co.jp/kitaurawa/index.html   (Japanese page)'
        write(*, *) ' http://members.tripod.co.jp/kitaurawa/index_e.html (English  page)'
        stop 
    end subroutine print_about
    !---------------------------------------------------------------
end module arguments