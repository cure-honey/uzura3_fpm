module mod_psycho
    use kind_m
    use mod_mpg
    use mod_fft
    implicit none
    private
    public :: psycho, calc_mask  ! subroutine
    real (kind = kd), save :: pi, pi2 
    real (kind = kd), save :: ath_l(576), ath_s(192, 3), ath0_l(100), ath0_s(100), dbath_l(100), dbath_s(100)
    real (kind = kd), save :: sf_l(100, 100),  sf_s(100, 100)
    real (kind = kd), save :: afft_l(576, 2, 2), afft_s(192, 3, 2, 2)
    real (kind = kd), save :: dbfft_l(100, 2, 2), dbfft_s(100, 3, 2, 2)
    real (kind = kd), save :: freq_l(576), freq_s(192), bark_l(576), bark_s(192), bw_l(576), bw_s(192)
    real (kind = kd), save :: weight_l(576), weight_s(192)
    complex (kind = kd), save :: cfft_l(1152, 2, 2), cfft_s(384, 3, 2, 2)
    integer, save          :: ibark_l(576), ibark_s(192), ifb_l(100, 0:2), ifb_s(100, 0:2), bw0_l(100), bw0_s(192)
contains
    !----------------------------------------------------------------
    subroutine init_absolute_threshold(isample_rate)
        integer, intent(in) :: isample_rate
        real(kind = kd):: freq, temp !, ath(576)
        integer :: i 
        pi  = 4.0_kd * atan(1.0_kd)
        pi2 = 2.0_kd * pi
        !
        do i = 1, 576
            freq = real(isample_rate, kind = kd) / 2.0_kd / 1000.0_kd * (real(i - 1, kind = kd) + 0.0_kd) / 576.0_kd
        !    temp =  3.64_kd  * freq ** (-0.8_kd) & 
        !         -  6.50_kd  * exp(-0.6_kd * (freq -  3.3_kd)**2.0_kd) &
        !         +  0.001_kd * freq ** 4.0_kd 
        !         + ath_min
            temp = 3.64_kd * freq ** (-0.8_kd) &                      ! alternative ATH function 
                 - 6.50_kd * exp(-0.6_kd * (freq - 3.3_kd)**2.0_kd) & ! reference: lame ath-type 3  
            !      + 5.1_kd 
            !     + 0.001_kd * freq ** 3.80_kd + 4.92_kd  ! 12.61    
            !  temp = 0.001_kd * freq ** 3.81_kd + 4.60_kd  ! 12.93     
            !  temp = 0.001_kd * freq ** 3.82_kd + 4.27_kd  ! 13.26     
            !  temp = 0.001_kd * freq ** 3.83_kd + 3.94_kd  ! 13.59     
            !  temp = 0.001_kd * freq ** 3.84_kd + 3.61_kd  ! 13.93     
                  + 0.001_kd * freq ** 3.85_kd + 3.27_kd  ! 14.28     
            !      + 0.001_kd * freq ** 4.00_kd - 2.86_kd  ! 20.70     
            temp = min(temp + ath_min, ath_max)
            ath_l(i) = 10.0_kd**(temp / 20.0_kd)   
        end do
        !
        do i = 1, 192
            freq = real(isample_rate, kind = kd) / 2.0_kd / 1000.0_kd * (real(i - 1, kind = kd) + 0.0_kd) / 192.0_kd
            temp = 3.64_kd * freq ** (-0.8_kd) &                      ! alternative ATH function 
                 - 6.50_kd * exp(-0.6_kd * (freq - 3.3_kd)**2.0_kd) & ! reference: lame ath-type 3  
            !     + 0.001_kd * freq ** 3.80_kd + 4.92_kd  ! 12.61    
                 + 0.001_kd * freq ** 3.85_kd + 3.27_kd  ! 14.28    
            !     + 0.001_kd * freq ** 4.00_kd - 2.86_kd  ! 20.70    
            temp = min(temp + ath_min, ath_max)
            ath_s(i, :) = 10.0_kd**(temp / 20.0_kd)   
        end do
    end subroutine init_absolute_threshold
!--------------------------------------------------------------------------------------------------------------------
    function switch_q(wx, igranule, ichannel) result(ires) ! attack detection 
        real (kind = kd), intent(in) :: wx(:, :, :)
        integer, intent(in) :: igranule, ichannel
        integer :: ires
        real (kind = kd), save :: sum0a(2, 2), sum1a(2, 2) = 0.0_kd, sum0b(2, 2), sum1b(2, 2) = 0.0_kd
        sum0a(igranule, ichannel) = sum1a(igranule, ichannel)
        sum1a(igranule, ichannel) = sum( abs(wx(1:36, igranule, ichannel)) )
        sum0b(igranule, ichannel) = sum1b(igranule, ichannel)
        sum1b(igranule, ichannel) = sum( abs(wx(37:, igranule, ichannel)) )
        ! 
        if ( sum0a(igranule, ichannel) < 1.25_kd .or. sum0b(igranule, ichannel) < 0.70_kd ) then  ! when silent, use long
            ires = 0 
        else if (sum1a(igranule, ichannel) > xsm * switch * sum0a(igranule, ichannel) ) then  ! band 1:2 attack 
            ires = mblock_type_param ! 20
            nn1 = nn1 + 1 ! debug info
        else if ( sum1b(igranule, ichannel) >   switch * sum0b(igranule, ichannel) ) then    ! band 3:32 attack  
            ires = mblock_type_param ! 20
            if ( q_sm ) ires = 21 ! mixed 
            nn2 = nn2 + 1 ! debug info
        else    ! long block         
            ires = 0
        end if
        if (mblock_type_param == 0) ires = 0 ! force long-only mode
    end function switch_q
!--------------------------------------------------------------------------------------------------------------------
    subroutine attack(wx, mblock_type) !..... ! iso figure c.7 (p.95) window switching state diagram
        real (kind = kd), intent(in ) :: wx(:, :, :)
        integer         , intent(out) :: mblock_type(:, :)
        integer :: igranule, ichannel, iattack
        integer, save :: mblock_prev(2, 2) = 0
        do ichannel = 1, size(wx, dim = 3)    
            do igranule = 1, 2
                iattack = switch_q(wx, igranule, ichannel)
                select case (iattack)
                case (0) ! no-attack
                    select case (mblock_prev(igranule, ichannel))
                    case ( 0, 30, 31) ! long
                        mblock_type(igranule, ichannel) =  0
                    case (10) 
                        mblock_type(igranule, ichannel) = 20
                    case (11) 
                        mblock_type(igranule, ichannel) = 21
                    case (20) ! short
                        mblock_type(igranule, ichannel) = 30
                    case (21) ! mixed
                        mblock_type(igranule, ichannel) = 31
                    case default
                        write(*, *) 'error: psycho : unexpected block type: case 0-x',  igranule, &
                                     mblock_prev(igranule, ichannel), mblock_type(igranule, ichannel) 
                        stop
                    end select
                case (20) ! attack-short
                    select case (mblock_prev(igranule, ichannel))
                    case ( 0, 30, 31)  
                        mblock_type(igranule, ichannel) = 10
                    case (10) 
                        mblock_type(igranule, ichannel) = 20
                    case (11) 
                        mblock_type(igranule, ichannel) = 21
                    case (20) 
                        mblock_type(igranule, ichannel) = 20
                    case (21) 
                        mblock_type(igranule, ichannel) = 21 
                    case default
                       write(*, *) 'error: psycho : unexpected block type: case 20-x',  igranule, &
                                    mblock_prev(igranule, ichannel), mblock_type(igranule, ichannel) 
                       stop
                    end select
                case (21) ! attack-mixed
                    select case (mblock_prev(igranule, ichannel))
                    case ( 0, 30, 31)  
                        mblock_type(igranule, ichannel) = 11
                    case (10) 
                        mblock_type(igranule, ichannel) = 20
                    case (11) 
                        mblock_type(igranule, ichannel) = 21
                    case (20) 
                        mblock_type(igranule, ichannel) = 20 
                    case (21) 
                        mblock_type(igranule, ichannel) = 21
                    case default
                        write(*, *) 'error: psycho : unexpected block type: case 21-x', igranule, &
                                     mblock_prev(igranule, ichannel), mblock_type(igranule, ichannel) 
                        stop
                    end select
                case default
                    write(*, *) 'error: psycho : unexpected block type'
                    stop
                end select
                mblock_prev(igranule, ichannel) = mblock_type(igranule, ichannel)
                !....... debug info ...............................................
                select case ( mblock_type(igranule, ichannel) )
                case ( 0) 
                    long = long + 1
                case (20)
                    nshort = nshort + 1
                case (21) 
                    mix = mix + 1
                case (10, 11)
                    m1 = m1 + 1
                case (30, 31)
                    m3 = m3 + 1
                end select
            end do
        end do
    end subroutine attack    
!-----------------------------------------------------------------------------------------------------------------------
    subroutine fft_long(nchannel, pcm, afft_l, fft576)
        integer            , intent(in ) :: nchannel
        real (kind = kd)   , intent(in ) :: pcm(:, :)
        real (kind = kd)   , intent(out) :: afft_l(:, :, :)
        complex (kind = kd), intent(out) :: fft576(1152, 2, 2)
        integer :: ichannel, igranule, m1, m2
        do igranule = 1, 2
            do ichannel = 1, nchannel
                m1 = 1 + 480 * (igranule - 1)
                m2 = m1 + 1152 - 1
                fft576(:, igranule, ichannel) = cmplx(pcm(m1:m2, ichannel), 0.0_kd, kind = kd) ! put 1152 real data -> get 576 complex fft
                call fft23han(7, 2, indx576, omega576, fft576(:, igranule, ichannel), han576 ) ! 2^7 * 3^2 = 1152
            !    call fft23(7, 2, indx576, omega576, fft576(:, igranule, ichannel) )
                afft_l(:, igranule, ichannel) = abs(fft576(1:576, igranule, ichannel))
            end do
        end do
    end subroutine fft_long
!------------------------------------------------------------------------------------------------
    subroutine fft_short(nchannel, pcm, afft_s, fft192)
        integer            , intent(in ) :: nchannel
        real (kind = kd)   , intent(in ) :: pcm(:, :)
        real (kind = kd)   , intent(out) :: afft_s(:, :, :, :)
        complex (kind = kd), intent(out) :: fft192(384, 3, 2, 2)
        integer :: ichannel, igranule, iwin, m1, m2
        do igranule = 1, 2
            do ichannel = 1, nchannel
               do iwin = 1, 3
                    m1 = 1 + 480 * (igranule - 1) + 384 * (iwin - 1)
                    m2 = m1 + 384 - 1
                    fft192(:, iwin, igranule, ichannel) = cmplx(pcm(m1:m2, ichannel), 0.0_kd, kind = kd) ! put 384 real data -> get 192 complex fft
                    call fft23han(7, 1, indx192, omega192, fft192(:, iwin, igranule, ichannel), han192 ) ! 2^7 * 3^1 = 384
                !    call fft23(7, 1, indx192, omega192, fft192(:, iwin, igranule, ichannel) )
                    afft_s(:, iwin, igranule, ichannel) = abs(fft192(1:192, iwin, igranule, ichannel))
                end do
            end do
        end do
    end subroutine fft_short
!------------------------------------------------------------------------------------------------
    subroutine calc_wx(nchannel, wx)
        integer        , intent(in ) :: nchannel
        real(kind = kd), intent(out) :: wx(:, :, :)
        integer :: igranule, ichannel
        forall (igranule = 1:2, ichannel = 1:nchannel) wx(:, igranule, ichannel) = afft_l(:, igranule, ichannel) * weight_l
    end subroutine calc_wx
!------------------------------------------------------------------------------------------------
    subroutine psycho(pcm, mpg, mblock_type)
        type (mpeg_parameters), intent(in out) :: mpg
        integer               , intent(   out) :: mblock_type(:, :)
        real (kind = kd)      , intent(in    ) :: pcm(:, :)
        logical, save :: qfirst = .true.
        integer       :: nchannel
        real (kind = kd) :: wx(576, 2, 2)
        !..... initialization .........................................................................................
        if (qfirst) then
            qfirst = .false.               
            call init_absolute_threshold( mpeg_sample_rates(mpg%isample_rate) )
            call init_mask( mpeg_sample_rates(mpg%isample_rate) )
            call init_fft()
        end if
        !..... fft 576 / 192 ................................................................................................
        nchannel = size(mblock_type, 2) 
        call fft_long (nchannel, pcm, afft_l, cfft_l)
        call fft_short(nchannel, pcm, afft_s, cfft_s)
        !..... weighted intensity .....................................................................................
        call calc_wx(nchannel, wx)
        !..... attack detection .......................................................................................
        call attack(wx, mblock_type)
        !  psycho calaulates fft & decides long/short/mixed block.  
        !  mid-side/normal stereo decision (layer3.f90)
        !
        ! debug info
        nbits(mpg%ibit_rate) = nbits(mpg%ibit_rate) + 1
    end subroutine psycho
!------------------------------------------------------------------------------------------------
    subroutine init_mask(nsample_rate)
        integer, intent(in) :: nsample_rate
        integer :: i
        real (kind  = 8) :: f0, f1

        do i = 1, 576
            freq_l(i) = real(nsample_rate, kind = kd) / 2.0_kd * (real(i - 1, kind = kd) + 0.0_kd) / 576.0_kd ! khz
            bark_l(i) = bark(freq_l(i) / 1000.0_kd) 
            ibark_l(i) = int(bark_l(i) + 0.1_kd) + 1 
            f0 = real(nsample_rate, kind = kd) / 2000.0_kd * real(i - 1, kind = kd) / 576.0_kd
            f1 = real(nsample_rate, kind = kd) / 2000.0_kd * real(i    , kind = kd) / 576.0_kd
            bw_l(i) = bark(f1) - bark(f0)
            weight_l(i) = bw_l(i) / (f1 - f0) 
        end do
        do i = 1, 192
            freq_s(i) = real(nsample_rate, kind = kd) / 2.0_kd * (real(i - 1, kind = kd) + 0.0_kd) / 192.0_kd ! khz
            bark_s(i) = bark(freq_s(i) / 1000.0_kd)  
            ibark_s(i) = int(bark_s(i) + 0.1_kd) + 1
            f0 = real(nsample_rate, kind = kd) / 2000.0_kd * real(i - 1, kind = kd) / 192.0_kd
            f1 = real(nsample_rate, kind = kd) / 2000.0_kd * real(i    , kind = kd) / 192.0_kd
            bw_s(i) = bark(f1) - bark(f0)
            weight_s(i) = bw_s(i) / (f1 - f0) 
        end do
!
        do i = 576, 1, -1
            ifb_l(max(1, ceiling(4 * bark_l(i))), 1) = i
        end do
        do i = 1, 576
            ifb_l(max(1, ceiling(4 * bark_l(i))), 2) = i 
        end do
        bw_l = -huge(0.0_kd)        
        do i = 1, 100
            ifb_l(i, 0) = ifb_l(i, 2) - ifb_l(i, 1) + 1  
            if (ifb_l(i, 1) /= 0) bw0_l(i) = ifb_l(i, 0) * (freq_l(2) - freq_l(1))
            if (ifb_l(i, 1) /= 0) ath0_l(i) = ath_l( (ifb_l(i, 2) + ifb_l(i, 1)) / 2)
            if (ifb_l(i, 1) /= 0) dbath_l(i) = 20.0_kd * log10(ath0_l(i))
        end do    
!        
        do i = 192, 1, -1
            ifb_s(max(1, ceiling(4 * bark_s(i))), 1) = i
        end do
        do i = 1, 192
            ifb_s(max(1, ceiling(4 * bark_s(i))), 2) = i 
        end do
        bw_s = -huge(0.0_kd)        
        do i = 1, 100
            ifb_s(i, 0) = ifb_s(i, 2) - ifb_s(i, 1) + 1  
            if (ifb_s(i, 1) /= 0) bw0_s(i) = ifb_s(i, 0) * (freq_s(2) - freq_s(1))
            if (ifb_s(i, 1) /= 0) ath0_s(i) = ath_s( (ifb_s(i, 2) + ifb_s(i, 1)) / 2, 1)
            if (ifb_s(i, 1) /= 0) dbath_s(i) = 20.0_kd * log10(ath0_s(i))
        end do    
    end subroutine init_mask
!------------------------------------------------------------------------------------------------
    subroutine calc_mask(igranule, ichannel, mblock_type, xmask)
        integer         , intent(in ) :: igranule, ichannel, mblock_type
        real (kind = kd), intent(out) :: xmask(:, :)
        real (kind = kd) :: x0_l(100), x0_s(100, 3)
        real (kind = kd) :: y_l(576), y_s(192, 3)
        real (kind = kd) :: yn_l(576), yn_s(192, 3)
        real (kind = kd) :: tmp_l(100), tmp_s(100, 3)
        real (kind = kd) :: v_l(100), v_s(100, 3) 
        integer :: i, j, iwin, iband
        integer, parameter :: il = 1, iu = 1 ! average over bark: i-il:i+iu
        != reference =======================================================================================================
        ! F. Baumgarte, C. Ferekidis and H. Fuchs, "A Nonlinear Psychoacoustic Model Applied to the ISO MPEG Layer 3 Coder"
        !===================================================================================================================
        ! sum up to 1/4 bark
        v_l = 1.0e-32_kd 
        forall(i = 1:100, ifb_l(i, 1) /= 0) v_l(i) = sum(afft_l(ifb_l(i, 1):ifb_l(i, 2), igranule, ichannel)**2) / ifb_l(i, 0) !* 4
        ! non-linear sum
        dbfft_l(:, igranule, ichannel) = 10.0_kd * log10(2.0_kd * 1152.0_kd * v_l / 2.0_kd) ! 1/N FFT.f90 ! AFFT 1152 -> 576 
        forall (i = 1:100, j = 1:100) sf_l(i, j) = spreading_function( 0.25_kd * (i - j), dbfft_l(j, igranule, ichannel) )   
        do i = 1, 100   !  pcm^2 : MDCT^2 ~ 300 : 1 ; pcm^2 : FFT^2 ~ sqrt(1152/384) : 1
            tmp_l = sf_l(i, :) + dbfft_l(:, igranule, ichannel) - dbath_l + ath_min + offset_l &![dB]  49.2dB
                                                                                    + 20.0_kd * log10( 576.0_kd / 2.0_kd)            
            x0_l(i) = sum( 10.0_kd**(tmp_l * pow / 10.0_kd) )**( 0.5_kd / pow) * ath0_l(i)
        end do
        forall (i = 2:99, ifb_l(i, 1) /= 0) & ! average over nearest neighbors ~ 1bark
               tmp_l(i) = minval(x0_l(i - il:i + iu), mask = ifb_l(i - il:i + iu, 1) /= 0)
        forall (i = 1:100, ifb_l(i, 1) /= 0) yn_l(ifb_l(i, 1):ifb_l(i, 2)) = tmp_l(i)
        y_l = max( ath_l, yn_l )
        ! help for parameter decision  
        if (q_help .and. mod(iframe, 200) == 0) print *, count(ath_l < yn_l), ' ath_l < yn_l'
        ! short block : uses psychoacoustic analysis for long blocks
        forall (i = 1:192) yn_s(i, :) = minval( yn_l(3 * (i - 1) + 1:3 * (i - 1) + 3) ) 
        y_s = max( ath_s, yn_s ) 
        ! short block (if q_short_fft == .true.) : independent psychoacoustic analysis for short blocks
        if (q_short_fft .and. (mblock_type == 20 .or. mblock_type == 21)) then
            v_s = 1.0e-32_kd 
            forall (i = 1:100, iwin = 1:3, ifb_s(i, 1) /= 0) &
                  v_s(i, iwin) = sum(afft_s(ifb_s(i, 1):ifb_s(i, 2), iwin, igranule, ichannel)**2) / ifb_s(i, 0) 
            dbfft_s(:, :, igranule, ichannel) = 10.0_kd * log10(2.0_kd * 384.0_kd * v_s * 3 * 3 / 2.0_kd) !normalization factor 9 ?? ! 1/N FFT.f90 
            do iwin = 1, 3    
                forall (i = 1:100, j = 1:100) sf_s(i, j) = &
                       spreading_function( 0.25_kd * (i - j), dbfft_s(j, iwin, igranule, ichannel))  
                do i = 1, 100                                                               
                    tmp_s(:, iwin) = sf_s(i, :) + dbfft_s(:, iwin, igranule, ichannel) & 
                                   - dbath_s + ath_min  + offset_l + 20.0_kd * log10( 576.0_kd / 2.0_kd) ! 49.2dB shift polyphase+MDCT : FFT
                    x0_s(i, iwin) = sum( 10.0_kd**(tmp_s(:, iwin) * pow / 10.0_kd) )**( 0.5_kd / pow) * ath0_s(i) 
                end do
            end do   
            forall (i = 2:99,  iwin = 1:3, ifb_s(i, 1) /= 0) &  ! average over nearest neighbors ~ 1bark
                    tmp_s(i, iwin) = minval(x0_s(i - il:i + iu, iwin), mask = ifb_s(i - il:i + iu, 1) /= 0) 
            forall (i = 1:100, iwin = 1:3, ifb_s(i, 1) /= 0) yn_s(ifb_s(i, 1):ifb_s(i, 2), iwin) = tmp_s(i, iwin)
            y_s = max( ath_s, yn_s ) 
        end if
        ! reorder 32x18 --> 576 (long), 32x3x6 --> 576 (short), 2x18 + 30x3x6 -> 576 (mixed)
        select case (mblock_type)
        case (0, 10, 11, 30, 31) ! long block  ! transition block
            forall (iband = 1:32, i = 1:18) xmask(iband, i) = y_l(18 * (iband - 1) + i)
        case (20) ! short block
            forall (iband = 1:32, iwin = 1:3, i = 1:6) xmask(iband, 6 * (iwin - 1) + i) = y_s(6 * (iband - 1) + i, iwin)
        case (21) ! mixed block   
            forall (iband = 1:2, i = 1:18) xmask (iband, i) = y_l(18 * (iband - 1) + i)
            forall (iband = 3:32, iwin = 1:3, i = 1:6) xmask(iband, 6 * (iwin - 1) + i) = y_s(6 * (iband - 1) + i, iwin)
        case default
            print *, 'mblock_type', mblock_type
            stop 'calc_mask: should not come here'
        end select
    end subroutine calc_mask
 !------------------------------------------------------------------------------------------------
    pure elemental function bark(f) result(res)
        real (kind = kd), intent(in) :: f
        real (kind = kd) :: res
        res = 13.0_kd * atan(0.76_kd * f) + 3.5_kd * atan( (f / 7.5_kd)**2.0_kd )
    end function bark
 !------------------------------------------------------------------------------------------------
    pure elemental function spreading_function(z, x) result(res)
        real (kind = kd), intent(in) :: z, x
        real (kind = kd) :: res
        if ( z > 0.0_kd ) then
            res = -(22.0_kd - 0.2_kd * (x + ath_min)) !93.3_kd)) ! 10.0_kd * log10( 1.0_kd / 2.0_kd**(15*2 + 1) ))) ! -93.31929866
        else
            res =  30.0_kd * z 
        end if
    end function spreading_function
 !------------------------------------------------------------------------------------------------
end module mod_psycho