program uzura3
    use arguments
    use wav_io
    use bit_io
    use mod_mpg
    use mod_polyphase
    use mod_psycho
    use mod_mdct
    use mod_layer3
    use mod_huffman
    use mod_encode
    implicit none
    type (riff_chunk) :: riff
    type (mpeg_parameters) :: mpg
    real (kind = 8), allocatable :: pcm(:, :), subband(:, :, :), r_mdct(:, :, :, :)
    integer        , allocatable :: i_mdct(:, :, :), mblock_type(:, :), mblock_prev(:)
    integer :: iframe_length, ianc
    integer :: nchannel, itotal_frames, i !, iframe
    character (len = 40) :: file_name, fn_in, fn_out
    !
    mpg%mtype           =  3 ! 0:mpeg2.5, 1:---, 2:mpeg2, 3:mpeg1
    mpg%layer           =  1 ! layer { 1 = iii, 2 = ii, 3 = i }
    mpg%ibit_rate       =  9 ! 128kbps
    mpg%isample_rate    =  0
    mpg%ipadding        =  0
    mpg%icrc            =  1 ! crc16  0 enabled / 1 disabled
    mpg%iprivate        =  0
    mpg%mode            =  1 ! joint stereo
    mpg%mode_extension  =  2 ! intensity_stereo off / ms_stereo on
    mpg%icopyright      =  0
    mpg%ioriginal       =  0
    mpg%iemphasis       =  0 ! 
    call init_huffman()
    call init_huffman_quadruples()
    call get_option(mpg, file_name)
    fn_in  = trim(file_name) // '.wav'
    fn_out = trim(file_name) // '.mp3'
    call pr_info(mpg)
    call open_wav_file(10, fn_in)
    call check_riff_chunk(riff)
    call play_time(riff)
    call sample_rate(riff, mpg)
    call init_scalefactor_bands(mpg%isample_rate)
    nchannel = riff%fmt%ichannels
    itotal_frames = riff%dat%ichunk_size / (mpeg_frame_size(mpg%layer) * nchannel * 2)   ! iso 2.4.2.3 padding bit (p22)
    if ( mod(riff%dat%ichunk_size, mpeg_frame_size(mpg%layer) * nchannel * 2) > 480 ) itotal_frames = itotal_frames + 1
    allocate( pcm(2304, 2), subband(32, 36, nchannel) )
    allocate( r_mdct(32, 18, 2, nchannel), i_mdct(576, 2, nchannel), mblock_type(2, nchannel), mblock_prev(2) )
    call open_mpg_file(9, fn_out)
    mblock_type = mblock_type_param
    pcm = 0.0d0
    call read_pcm_1frame(pcm)
    do iframe = 1, itotal_frames  ! main loop
        call read_pcm_1frame(pcm)
        call psycho(pcm, mpg, mblock_type)
        call polyphase_filter36(pcm(1:1632, :), subband)
        call sub_mdct(subband, r_mdct, mblock_type, q_alias)
        call alloc_bits(mblock_type, r_mdct, i_mdct, mpg, iframe_length, ianc)
        call encode_all(mpg, i_mdct, nchannel, ianc)
        call write_bits_1frame(iframe_length)
        if ( mod(iframe, 200) == 1 ) call update_status(riff, iframe, itotal_frames) 
    end do
    call update_status(riff, iframe - 1, itotal_frames) 
    write(*, *) 'total frames', iframe
    call close_wav_file()
    call close_mpg_file()
    write(*, '(a)') ' normal end. '
    call print_debug_info()
    stop
contains
!------------------------------------------------------------------
    subroutine print_debug_info()
        if (q_info) then
            print '(a)', ' ======== parameters ========================================================='
            print '(4(a, f8.2))', ' ath_min  =', ath_min,  ' ath_max  =', ath_max, ' offset =', offset_l
            print '(3(a, f8.2))', ' switch   =', switch,   ' xsm      =', xsm    , ' xms    =', xms                
            print '(3(a, f8.2))', ' masking: p-norm =', pow
            if (q_short_fft) print '(a)', ' short fft used for short/mixed blocks'
            print '(a)', ' ======== info ==============================================================='
            print '(5(a, i7))', ' block type:long', long, ':short', nshort, ':mixed', mix, ':type1', m1, ':type3', m3
            print '(a, 2i7, a, 2i6, a, 3i7)', ' ms/ns select   ', ms, ns, ' [', ns1, ns2, '] :long/short switch', nn1, nn2
            print '(a)', '..............................................................................'
            print '(a, i7)', ' average scale factor (long 0-20)', n_sc_l
            print '(11f6.2)', tot_sc_l( 0: 9)  / real( max(n_sc_l, 1), kind = 8)  
            print '(11f6.2)', tot_sc_l(10:20)  / real( max(n_sc_l, 1), kind = 8) 
            print '(a, i7)', ' average scale factor (short 0-11)', n_sc_s
            print '(12f6.2)', ( tot_sc_s(:, 1) + tot_sc_s(:, 2) + tot_sc_s(:, 3) ) / real( max(3 * n_sc_s, 1), kind = 8) 
            print '(a)', ' selected huffman table (table 0-31) ' 
            print '( 4(1x, 10i7/) )', ntable 
            print '( (1x, a, 10i7/) )', 'selected count1 table a, b', ntab_ab 
            print '(a)', '..............................................................................'
            print '(3(a,  i6))', ' pre-emphasis ', n_emph, ' :scalefactor_scale ', n_scale, ' :subblock_gain ', n_sub_gain 
            print '(a)', '..............................................................................' 
            print '(a, 14i5  )', ' bit:', (i, i = 1, 14)
            print '(a, 14f5.1)', ' (%):', real(nbits * 100) / real(sum(nbits))
            print '(a, f7.2, a, es10.4)', ' average bit rate (kbps)', sum(real(nbits * mpeg_bit_rates(1:14,mpg%layer))) &
                                         / real(sum(nbits)), '      :maximum distortion ', distortion_max
            print '(a)', ' ============================================================================='
        end if
    end subroutine print_debug_info
!------------------------------------------------------------------
    subroutine sample_rate(riff, mpg)
        type (riff_chunk     ), intent(in ) :: riff
        type (mpeg_parameters), intent(out) :: mpg
        select case (riff%fmt%isamples_per_sec)
        case (44100)
            mpg%isample_rate = 0
        case (48000)
            mpg%isample_rate = 1
        case (32000)
            mpg%isample_rate = 2
        case default
            write(*, *) 'sampling rate ', riff%fmt%isamples_per_sec, ' is not supported in MPEG-I.'
            stop 'abnormal end'
        end select
    end subroutine sample_rate
!------------------------------------------------------------------
    subroutine play_time(riff)
        type (riff_chunk), intent(in) :: riff
        integer :: itot_time, ihour, imin, isec
        itot_time = itotal_play_time(riff)
        ihour =          itot_time / 3600
        imin  =      mod(itot_time, 3600) / 60
        isec  = mod( mod(itot_time, 3600) , 60 )
        write(*, '(a, i3, a, i2, a, i2)') ' playtime ', ihour, ':', imin, ':', isec
        write(*, *)
    end subroutine play_time
!------------------------------------------------------------------
    function itotal_play_time(riff) result(ires) ! return total play time in second
        type (riff_chunk), intent(in) :: riff
        integer :: ires
        ires = riff%dat%ichunk_size / riff%fmt%ibytes_per_sec
    end function itotal_play_time    
!------------------------------------------------------------------
    subroutine pr_info(mpg)
        type (mpeg_parameters), intent(in) :: mpg
        write(*, *) 'Uzura3 (mp3 encoder/Fortran95) ver.0.6 (c) H.O. psychoacoustic model Baumgarte et al.'
        if (mpg%icrc == 0) write(*, *) 'crc16 error protection enabled'
    end subroutine pr_info
!------------------------------------------------------------------
    subroutine update_status(riff, iframe, itot_frames)
        type (riff_chunk), intent(in) :: riff
        integer, intent(in) :: iframe, itot_frames
        integer :: it(8), ielapsed, iel_min, iel_sec
        integer, save :: istart
        logical :: qfirst = .true.
        real    :: percent
        character (len = 10) :: time, date, zone
        call date_and_time(date, time, zone, it)
        if (qfirst) then
           istart   = it(5) * 3600 + it(6) * 60 + it(7)
           qfirst   = .false.
        end if
        ielapsed = it(5) * 3600 + it(6) * 60 + it(7) - istart
        iel_min  =     ielapsed / 60
        iel_sec  = mod(ielapsed , 60)
        percent = real(100 * iframe) / real(itot_frames)
        write(*, '(a, f6.2, a, i4, a, i4, 2(a, i2), 3(a, i2.2), a, i4.2, a, i2.2, a)')  &
              '+processed...', percent, '% [', nint(percent * itotal_play_time(riff) / 100), 'secs] ', &
              it(1), '/', it(2), '/', it(3), ' ', it(5), ':', it(6), ':', it(7), &
              ' time elapsed ', iel_min, 'min ', iel_sec, 'sec'
    end subroutine update_status
!----------------------------------------------------------------------------
end program uzura3
