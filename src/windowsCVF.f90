module arguments 
    ! module for command line option  ! compaq (dec) visual fortran for intel (windows) 
!    use dflib
!    use ifport ! intel visual fortran
    use mod_mpg
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
            case ('-switch') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(f', istatus, '.0)' 
                read(buffer, fmt) switch
            case ('-ath_min') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(f', istatus, '.0)' 
                read(buffer, fmt) ath_min
            case ('-ath_max') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(f', istatus, '.0)' 
                read(buffer, fmt) ath_max
            case ('-crc') 
                mpg%icrc       = 0 ! crc16 on 
            case ('-cut') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(i', istatus, ')' 
                read(buffer, fmt) icut
                if (icut < 0 .or. icut > 32) call option_error( trim(buffer) )
            case ('-xms') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(f', istatus, '.0)' 
                read(buffer, fmt) xms
            case ('-xsm') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(f', istatus, '.0)' 
                read(buffer, fmt) xsm
            case ('-offset') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(f', istatus, '.0)' 
                read(buffer, fmt) offset
            case ('-tempo') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(f', istatus, '.0)' 
                read(buffer, fmt) tempo
            case ('-factor') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(f', istatus, '.0)' 
                read(buffer, fmt) factor
                if (factor < 0.0d0 .or. factor > 1.0d0) then 
                    write(*, *) 'input out of range: 0 <= factor <= 1'
                    call option_error( trim(buffer) )
                end if
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
            case ('-pm') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(f', istatus, '.0)' 
                read(buffer, fmt) pm_factor
            case ('-skip') 
                iarg = iarg + 1
                if ( iarg >= narg ) call option_error( trim(buffer) )
                call get_command_argument(iarg, buffer, istatus)
                write(fmt, '(a, i1, a)') '(f', istatus, '.0)' 
                read(buffer, fmt) skip
                if (skip < 1.0d0) then
                    write(*, *) 'input out of range: skip >= 1.0'
                    call option_error( trim(buffer) )
                end if
            case ('-cuth') 
                cut_factor = 0.0d0 
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
            case ('-v') ! vbr mode
                q_vbr = .true. 
                qms_stereo = .false.
                icut = 32
            case ('-rio500') 
                q_rio500 = .true.
            case ('-c') 
                mpg%icopyright = 1 ! copyrigt on
            case ('-o') 
                mpg%ioriginal  = 1 ! original on
            case ('-ms')
                qms_stereo = .true.
            case ('-ns')
                qms_stereo = .false.
                pm_factor = 1.1d0 ! for vbr mode
            case ('-nomask')
                q_mask = .false.
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
        write(*, '(2a)') ' >>>>>> error near >>>>>>', text
        stop
    end subroutine option_error
    !---------------------------------------------------------------
    subroutine print_option()
        write(*, *) 'usage : uzura3 -options file_name '
        write(*, *) '      : file_name.wav -> file_name.mp3'
        write(*, *) 'option:-b  1..14   bitrate for cbr mode                 (default 9 : 128kbps)'
        write(*, *) '       -crc        crc16 error protection on            (default off)'
        write(*, *) '       -c          copyright flag on                    (default off)'
        write(*, *) '       -o          original  flag on                    (default off)'
        write(*, *) '       -cut 1..32  band cut-off : place after -b option (default 26: 17.9khz)'
        write(*, *) '       -cuth       cut band 21 (l) or 12 (s/m)          (default off) '
        write(*, *) '       -v          vbr mode  (ns, icut = 32)            (default off)'
        write(*, *) '       -rio500     avoid rio500 vbr skip bug            (default off)'
        write(*, *) '       -l          long-block-only                      (default off)'
        write(*, *) '       -s          short-mode for short-block           (default off)'
        write(*, *) '       -m          mixed-mode for short-block           (default off)'
        write(*, *) '       -sm         short & mixed-mode for short-block   (default on )'
        write(*, *) '       -xsm xx     short / mixed switching parameter    (default 1.5)'
        write(*, *) '       -switch xx  long/short switching parameter       (default 2.0)'
        write(*, *) '       -skip   xx  speeds up outer loop                 (default 1.3)'
        write(*, *) '       -ms/-ns     stereo mode (ms/normal)              (default ms)'
        write(*, *) '       -xms xx     ms/ns      switching parameter       (default 0.5)'
        write(*, *) '       -nomask     masking off                          (default on)'
        write(*, *) '       -ath_min xx minimum of ath  [ db ]               (default -125)'
        write(*, *) '       -ath_max xx ceiling of ath  [ db ]               (default  0.0)'
        write(*, *) '       -offset  xx offset for mask [ db ]               (default 40.0)'
        write(*, *) '       -tempo   xx temporal masking factor              (default 0.85)'
        write(*, *) '       -factor  xx bit distribution among (gr, ch)      (default 0.4)'
        write(*, *) '       -noalias    anti-alias for mixed-block off       (default on)'
        write(*, *) '       -debug      print debug info                     (default on) '
        write(*, *) '       -about      about uzura3 '
        write(*, *) '        '
        write(*, *) 'example cbr 128kbps crc on : uzura3 -crc   file_name '
        write(*, *) 'example vbr normal stereo  : uzura3 -v -ns file_name '
    end subroutine print_option
    !---------------------------------------------------------------
    subroutine print_about()
        implicit none
        write(*, *) ' ******* uzura3 is an mpeg-i/layer-iii encoder written in fortran90. ********'
        write(*, *) ' **** this program comes with absolutely no warranty. (c) h.o. 2000-2004 **** '
        write(*, *)
        write(*, *)
        write(*, *) ' http://members.tripod.co.jp/kitaurawa/index.html   (japanese page)'
        write(*, *) ' http://members.tripod.co.jp/kitaurawa/index_e.html (english  page)'
        stop 
    end subroutine print_about
    !---------------------------------------------------------------
end module arguments