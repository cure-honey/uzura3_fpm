module mod_psycho
    use kind_m
    use mod_mpg
    use mod_fft
    implicit none
    private
    public :: psycho, calc_mask, alpha  ! subroutine
    real (kind = kd), parameter :: pi = 4 * atan(1.0_kd), pi2 = 2 * pi
    real (kind = kd) :: ath_l(576) , ath_s(192, 3)
    real (kind = kd) :: sf_l(576, 576), sf_s(192, 192) 
    real (kind = kd) :: afft_l(576, 2, 2), afft_s(192, 3, 2, 2)
    real (kind = kd) :: phi_l(576, 2, 2), phi_s(192, 3, 2, 2)
    real (kind = kd) :: freq_l(576), freq_s(192), bark_l(576), bark_s(192), bw_l(576), bw_s(192)
    real (kind = kd) :: weight_l(576), weight_s(192)
    integer          :: ibark_l(576), ibark_s(192), ifb_l(25, 0:2), ifb_s(25, 0:2)
    real (kind = kd), parameter :: alpha = 0.27_kd ! non-linear factor
    !
contains
    !----------------------------------------------------------------
    subroutine init_absolute_threshold(isample_rate)
        integer, intent(in ) :: isample_rate
        real(kind = kd):: freq, temp !, ath(576)
        integer :: i, k
        do i = 1, 576
            freq = real(isample_rate, kind = kd) / 2.0_kd / 1000.0_kd * (real(i, kind = kd) - 0.0_kd) / 576.0_kd
        !  temp =  3.64_kd  * freq ** (-0.8_kd) & 
        !       -  6.50_kd  * exp(-0.6_kd * (freq -  3.3_kd)**2.0_kd) &
        !       +  0.001_kd * freq ** 4.0_kd &      
        !       + ath_min
            temp = 3.64_kd * freq ** (-0.8_kd) &                     ! alternative ath function 
                 - 6.50_kd * exp(-0.6_kd * (freq - 3.3_kd)**2.0_kd) &  ! reference: lame ath-type 3  
                 + 5.1_kd 
        !
            if      (freq >  5.0_kd .and. freq <=  5.5_kd) then
                temp =  5.1_kd + 4.0_kd * (freq -  5.0_kd) 
            else if (freq >  5.5_kd .and. freq <=  8.0_kd) then 
                temp =  7.1_kd + 2.76_kd * (freq -  5.5_kd) 
            else if (freq >  8.0_kd .and. freq <= 10.0_kd) then
                temp = 14.0_kd + 1.00_kd * (freq -  8.0_kd) 
            else if (freq > 10.0_kd .and. freq <= 11.5_kd) then
                temp = 16.0_kd + 0.33_kd * (freq - 10.0_kd) 
            else if (freq > 11.5_kd .and. freq <= 12.0_kd) then
                temp = 16.5_kd + 2.0_kd * (freq - 11.5_kd) 
            else if (freq > 12.0_kd) then
                temp = 0.001_kd * freq ** 3.81_kd + 4.60_kd  ! 12.93     
            !
            !  temp = 0.001_kd * freq ** 3.80_kd + 4.92_kd  ! 12.61    
            !  temp = 0.001_kd * freq ** 3.81_kd + 4.60_kd  ! 12.93     
            !  temp = 0.001_kd * freq ** 3.82_kd + 4.27_kd  ! 13.26     
            !  temp = 0.001_kd * freq ** 3.83_kd + 3.94_kd  ! 13.59     
            !  temp = 0.001_kd * freq ** 3.84_kd + 3.61_kd  ! 13.93     
            !  temp = 0.001_kd * freq ** 3.85_kd + 3.27_kd  ! 14.28     
            !  temp = 0.001_kd * freq ** 4.00_kd - 2.86_kd  ! 20.70     
            end if
            temp = min(temp + ath_min, ath_max)  
            ath_l(i) = 10.0_kd**(temp / 20.0_kd) 
        end do
        !
        do i = 1, 192
            k = 3 * (i - 1) + 1
            ath_s(i, :) = minval( ath_l(k:k + 2) ) 
        end do
    end subroutine init_absolute_threshold
    !------------------------------------------------------------------------------------------------
    function switch_q(wx) result(ires) ! attack detection (uzura original)
        real (kind = kd), intent(in) :: wx(:, :)
        integer :: ires
        real (kind = kd), save :: sum0a, sum1a, sum0b, sum1b
        sum0a = sum1a
!        sum1a = sum( wx(1:36, :) )
        sum1a = sum( wx(1:36, :)**alpha )**(1 / alpha)
        sum0b = sum1b
!        sum1b = sum( wx(37: , :) ) 
        sum1b = sum( wx(37: , :)**alpha )**(1 / alpha)
        if  ( sum1a > switch * sum0a .or. sum1b > switch * sum0b ) then
            ires = mblock_type_param
            if (q_sm .and. sum1a < xsm * sum0a .and. sum0a < xsm * sum1a ) ires = 21 ! mixed 
        else 
            ires = 0
        end if
        if (mblock_type_param == 0) ires = 0 ! force long-only mode

        !debug info
        if ( sum1a >= switch * sum0a ) nn1 = nn1 + 1
        if ( sum1b >= switch * sum0b ) nn2 = nn2 + 1
    end function switch_q
!--------------------------------------------------------------------------------------------------------------------
    subroutine attack(wx, mblock_type) !..... ! iso figure c.7 (p.95) window switching state diagram
        real (kind = kd), intent(in ) :: wx(:, :, :)
        integer        , intent(out) :: mblock_type(:, :)
        integer :: i, iattack
        integer, save :: mblock_prev = 0
        do i = 1, 2
            iattack = switch_q(wx(:, i, :))
            select case (iattack)
            case (0) ! no-attack
                select case (mblock_prev)
                case ( 0, 30, 31) ! long
                    mblock_type(i, :) =  0
                case (10) 
                    mblock_type(i, :) = 20
                case (11) 
                    mblock_type(i, :) = 21
                case (20) ! short
                    mblock_type(i, :) = 30
                case (21) ! mixed
                    mblock_type(i, :) = 31
                case default
                    write(*, *) 'error: psycho : unexpected block type: case 0-x', i, mblock_prev, mblock_type(i, 1) 
                    stop
                end select
            case (20) ! attack-short
                select case (mblock_prev)
                case ( 0, 30, 31)  
                    mblock_type(i, :) = 10
                case (10) 
                    mblock_type(i, :) = 20
                case (11) 
                    mblock_type(i, :) = 21
                case (20) 
                    mblock_type(i, :) = 20
                case (21) 
                   mblock_type(i, :) = 21 
                case default
                   write(*, *) 'error: psycho : unexpected block type: case 20-x', i, mblock_prev, mblock_type(i, 1) 
                   stop
                end select
            case (21) ! attack-mixed
                select case (mblock_prev)
                case ( 0, 30, 31)  
                    mblock_type(i, :) = 11
                case (10) 
                    mblock_type(i, :) = 20
                case (11) 
                    mblock_type(i, :) = 21
                case (20) 
                    mblock_type(i, :) = 20 
                case (21) 
                    mblock_type(i, :) = 21
                case default
                    write(*, *) 'error: psycho : unexpected block type: case 21-x', i, mblock_prev, mblock_type(i, 1) 
                    stop
                end select
            case default
                write(*, *) 'error: psycho : unexpected block type'
                stop
            end select
           !mblock_type = 20 ! for debug
            mblock_prev = mblock_type(i, 1)
            !....... debug info ...............................................
            select case ( mblock_type(i, 1) )
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
    end subroutine attack
!-----------------------------------------------------------------------------------------------------------------------
    ! not working correctly ; mdct based version is better  
    subroutine mid_side(mpg, wx, mblock_type) ! iso  c.2.4.3.4.9.2,  g.2 ms_stereo and intensity stereo coding layer iii
        real (kind = kd)       , intent(in    ) :: wx(:, :, :)
        type (mpeg_parameters), intent(in out) :: mpg
        integer        , intent(in) :: mblock_type(:, :)
        integer         :: igranule, nchannel, n0, n1
        real (kind = kd) :: tmp1, tmp2
        integer, save :: mode_old = 0, mode_ext_old = 0
        logical         :: qms
        nchannel = size(wx, 3)
        qms = qms_stereo
        select case (mpg%isample_rate) ! threshold ~ 7khz (empirical value)
        case (0) ! 44.1khz
            n0 = 183
        case (1) ! 48.0khz
            n0 = 168
        case (2) ! 32.0khz
            n0 = 252
        case default
            stop ' sample_rate error : subroutine mid_side '
        end select
        n1 = n0 + 1
        tmp1 = 0.0_kd
        tmp2 = 0.0_kd
        do igranule = 1, 2
        ! pop musics often have different l-r behavior in low and high frequency 
            tmp1 = tmp1 + sum( abs( wx(1:n0, igranule, 1) - wx(1:n0, igranule, 2) ) ) 
            tmp2 = tmp2 + sum(    ( wx(1:n0, igranule, 1) + wx(1:n0, igranule, 2) ) ) 
        end do
        if ( tmp1 > xms * tmp2 ) then 
            qms = .false.
            ns1 = ns1 + 1
            ns  = ns  + 1
        end if
        n1 = n0 + 1
        tmp1 = 0.0_kd
        tmp2 = 0.0_kd
        do igranule = 1, 2
        ! pop musics often have different l-r behavior in low and high frequency 
            tmp1 = tmp1 + sum( abs( wx(n1:, igranule, 1) - wx(n1:, igranule, 2) ) ) 
            tmp2 = tmp2 + sum(    ( wx(n1:, igranule, 1) + wx(n1:, igranule, 2) ) ) 
        end do
        if ( abs(tmp1) > xms * abs(tmp2) ) then 
            qms = .false.
            ns1 = ns1 + 1
            ns  = ns  + 1
        end if
        !
        if (mblock_type(1, 1) /= 0 .and. mblock_type(1, 1) /= 20 .and. mblock_type(1, 1) /= 21) then  
            mpg%mode = 0 !mode_old 
            mpg%mode_extension = 0 ! mode_ext_old
        else
            if (qms) then 
                 mpg%mode            =  1 ! joint stereo
                 mpg%mode_extension  =  2 ! intensity_stereo off / ms_stereo on
                 ms = ms + 1 
            else
                 mpg%mode            =  0 ! normal stereo
                 mpg%mode_extension  =  0 ! intensity_stereo off / ms_stereo off 
            end if
        end if 
        mode_old     = mpg%mode 
        mode_ext_old = mpg%mode_extension 
    end subroutine mid_side
!------------------------------------------------------------------------------------------------
    subroutine fft_long(nchannel, pcm, afft_l, phi_l)
        integer       , intent(in ) :: nchannel
        real (kind = kd), intent(in ) :: pcm(:, :)
        real (kind = kd), intent(out) :: afft_l(:, :, :), phi_l(:, :, :)
        complex (kind = kd) :: fft576(1152, 2, 2)
        integer :: ichannel, igranule, m1,m2
        do igranule = 1, 2
            do ichannel = 1, nchannel
                m1 = 1 + 480 * (igranule - 1)
                m2 = m1 + 1152 - 1
                fft576(:, igranule, ichannel) = cmplx(pcm(m1:m2, ichannel), 0.0_kd, kind = kd) ! put 1152 real data -> get 576 complex fft
             !  call fft23han(7, 2, indx576, omega576, fft576(:, igranule, ichannel), han576 ) ! 2^7 * 3^2 = 1152
                call fft23(7, 2, indx576, omega576, fft576(:, igranule, ichannel) )
                afft_l (:, igranule, ichannel) = abs(fft576(1:576, igranule, ichannel))
                afft_l (:, igranule, ichannel) = afft_l (:, igranule, ichannel) 
                phi_l(:, igranule, ichannel) = atan2(aimag(fft576(1:576, igranule, ichannel)), &
                                                      real(fft576(1:576, igranule, ichannel)))
            end do
        end do
        phi_l = phi_l + pi ! atan -pi~pi -> 0~2pi
    end subroutine fft_long
!------------------------------------------------------------------------------------------------
    subroutine fft_short(nchannel, pcm, afft_s, phi_s)
        integer       , intent(in ) :: nchannel
        real (kind = kd), intent(in ) :: pcm(:, :)
        real (kind = kd), intent(out) :: afft_s(:, :, :, :), phi_s(:, :, :, :)
        complex (kind = kd) :: fft192(384, 3, 2, 2)
        integer :: ichannel, igranule, iwin, m1, m2
        do igranule = 1, 2
            do ichannel = 1, nchannel
               do iwin = 1, 3
                   m1 = 1 + 480 * (igranule - 1) + 384 * (iwin - 1)
                   m2 = m1 + 384 - 1
                   fft192(:, iwin, igranule, ichannel) = cmplx(pcm(m1:m2, ichannel), 0.0_kd, kind = kd) ! put 384 real data -> get 192 complex fft
               !  call fft23han(7, 1, indx192, omega192, fft192(:, iwin, igranule, ichannel), han192 ) ! 2^7 * 3^1 = 384
                    call fft23(7, 1, indx192, omega192, fft192(:, iwin, igranule, ichannel) )
                    afft_s (:, iwin, igranule, ichannel) = abs(fft192(1:192, iwin, igranule, ichannel))
                    afft_s (:, iwin, igranule, ichannel) = afft_s (:, iwin, igranule, ichannel) 
                    phi_s(:, iwin, igranule, ichannel) = atan2(aimag(fft192(1:192, iwin, igranule, ichannel)), &
                                                                real(fft192(1:192, iwin, igranule, ichannel)))
                end do
            end do
        end do
        phi_s = phi_s + pi ! atan2 -pi~pi -> 0~2pi
    end subroutine fft_short
!------------------------------------------------------------------------------------------------
    subroutine calc_wx(nchannel, wx)
        integer        , intent(in ) :: nchannel
        real(kind = kd), intent(out) :: wx(:, :, :)
        integer :: igranule, ichannel
        do igranule = 1, 2
            do ichannel = 1, nchannel
!                wx(:, igranule, ichannel) = ( afft_l(:, igranule, ichannel) * weight_l )**2.0_kd 
                wx(:, igranule, ichannel) = ( afft_l(:, igranule, ichannel) * weight_l )**alpha 
            end do
        end do
    end subroutine calc_wx
!------------------------------------------------------------------------------------------------
    subroutine psycho(pcm, mpg, mblock_type)
        type (mpeg_parameters), intent(in out) :: mpg
        integer               , intent(   out) :: mblock_type(:, :)
        real (kind = kd)      , intent(in    ) :: pcm(:, :)
        logical, save :: qfirst = .true.
        integer        :: igranule, ichannel, nchannel
        real (kind = kd) :: wx(576, 2, 2), pm, pm0(2, 2)
        !..... initialization .........................................................................................
        if (qfirst) then
            qfirst = .false.               
            !!!- to avoid vbr bug of portable mp3 player diamond multimedia rio500. first frame bitrate must be less than average bitrate
            if (q_vbr .and. q_rio500) mpg%ibit_rate = 8 ! force first frame 112kbps for rio500 (firmware 1.15)  
            call init_absolute_threshold( mpeg_sample_rates(mpg%isample_rate) )
            call init_mask( mpeg_sample_rates(mpg%isample_rate) )
            call init_fft()
        end if
        !..... fft 576/192 ............................................................................................
        nchannel = size(mblock_type, 2) 
        call fft_long (nchannel, pcm, afft_l, phi_l)
        call fft_short(nchannel, pcm, afft_s, phi_s)
        !..... weighted intensity .....................................................................................
        call calc_wx(nchannel, wx)
        !..... attack detection .......................................................................................
        call attack(wx, mblock_type)
        !..... ms/ns selection ........................................................................................
        !call mid_side(mpg, wx, mblock_type) ! not working correctly : mdct based model works better
        !..... vbr  ...................................................................................................
        if (q_vbr) then ! simple implimentation 
            do igranule = 1, 2
                do ichannel = 1, nchannel
!                    pm0(igranule, ichannel) = sum(wx(:, igranule, ichannel) * bark_l)**2 !psychoacoustic moment (uzura original) 
                    pm0(igranule, ichannel) = sum(wx(:, igranule, ichannel) * bark_l)**(1 / alpha) !psychoacoustic moment (uzura original) 
                end do 
            end do
            pm = sum(pm0) / ( sum(wx) + 1.0d-9 )
            if (pm < 0.1_kd) then
                mpg%ibit_rate =  1 
            else if (pm             < 1.0_kd) then
                mpg%ibit_rate =  9  
            else if (pm * pm_factor < 2.5_kd) then
                mpg%ibit_rate = 10  
            else if (pm * pm_factor < 3.5_kd) then
                mpg%ibit_rate = 11  
            else if (pm * pm_factor < 5.0_kd) then
                mpg%ibit_rate = 12 
            else if (pm * pm_factor < 7.0_kd) then 
                mpg%ibit_rate = 13 
            else 
                mpg%ibit_rate = 14
            end if
        end if
        nbits(mpg%ibit_rate) = nbits(mpg%ibit_rate) + 1
    end subroutine psycho
!------------------------------------------------------------------------------------------------
    subroutine init_mask(nsample_rate)
        integer, intent(in) :: nsample_rate
        integer :: i, j
        real (kind  = 8) :: f0, f1
        do i = 1, 576
            freq_l(i) = real(nsample_rate, kind = kd) / 2.0_kd * (real(i, kind = kd) - 0.5_kd) / 576.0_kd ! khz
            bark_l(i) = bark(freq_l(i) / 1000.0_kd) 
            ibark_l(i) = int(bark_l(i) + 0.1_kd) + 1 
            f0 = real(nsample_rate, kind = kd) / 2000.0_kd * real(i - 1, kind = kd) / 576.0_kd
            f1 = real(nsample_rate, kind = kd) / 2000.0_kd * real(i    , kind = kd) / 576.0_kd
            bw_l(i) = bark(f1) - bark(f0)
            weight_l(i) = ( bark(f1) - bark(f0) ) / (f1 - f0) 
        end do
        do i = 1, 192
            freq_s(i) = real(nsample_rate, kind = kd) / 2.0_kd * (real(i, kind = kd) - 0.5_kd) / 192.0_kd ! khz
            bark_s(i) = bark(freq_s(i) / 1000.0_kd)  
            ibark_s(i) = int(bark_s(i)) + 1
            f0 = real(nsample_rate, kind = kd) / 2000.0_kd * real(i - 1, kind = kd) / 192.0_kd
            f1 = real(nsample_rate, kind = kd) / 2000.0_kd * real(i    , kind = kd) / 192.0_kd
            bw_s(i) = bark(f1) - bark(f0)
            weight_s(i) = ( bark(f1) - bark(f0) ) / (f1 - f0) 
        end do
        !
        ifb_l(1, 1) = 1
        do i = 1, 25
            do j = 1, 576
                if (ibark_l(j) == i - 1) ifb_l(i, 1) = j + 1 
                if (ibark_l(j) == i    ) ifb_l(i, 2) = j  
            end do
            ifb_l(i, 0) = ifb_l(i, 2) - ifb_l(i, 1) + 1
        end do
        ifb_s(1, 1) = 1
        do i = 1, 25
            do j = 1, 192
                if (ibark_s(j) == i - 1) ifb_s(i, 1) = j + 1 
                if (ibark_s(j) == i    ) ifb_s(i, 2) = j  
            end do
            ifb_s(i, 0) = ifb_s(i, 2) - ifb_s(i, 1) + 1
        end do
        !
        do i = 1, 576
            do j = 1, 576
               sf_l(i, j) = 10.0_kd ** ( spreading_function( bark_l(i) - bark_l(j) ) / 20.0_kd ) 
            end do
        end do 
        do i = 1, 192
            do j = 1, 192
                sf_s(i, j) = 10.0_kd ** ( spreading_function( bark_s(i) - bark_s(j) ) / 20.0_kd ) 
            end do
        end do 
    end subroutine init_mask
!------------------------------------------------------------------------------------------------
    subroutine calc_mask(igranule, ichannel, mblock_type, xmask, xnoise)
        integer       , intent(in    ) :: igranule, ichannel, mblock_type
        real (kind = kd), intent(   out) :: xmask(:, :), xnoise(:, :)
        real (kind = kd) :: x0_l(576), x0_s(192, 3), y0_l(576), y0_s(192, 3)
        real (kind = kd) :: d2phi_l(576), tone_l(576), fk_l(576), fl_l(576), tn_l(576)
        real (kind = kd) :: d2phi_s(192), tone_s(192), fk_s(192), fl_s(192), tn_s(192)
        real (kind = kd) :: yn(25)
        real (kind = kd), save :: x1_l(576) = 0.0_kd, x1_s(192, 3) = 0.0_kd
        real (kind = kd), save :: phi1_l(576, 2) = 0.0_kd, phi2_l(576, 2) = 0.0_kd
        real (kind = kd), save :: phi1_s(192, 2) = 0.0_kd, phi2_s(192, 2) = 0.0_kd
        real (kind = kd), save :: p0_l(576, 2) = 0.0_kd, p1_l(576, 2) = 0.0_kd, p2_l(576, 2) = 0.0_kd
        real (kind = kd), save :: p0_s(192, 2) = 0.0_kd, p1_s(192, 2) = 0.0_kd, p2_s(192, 2) = 0.0_kd
        real (kind = kd), save :: p3_s(192, 2) = 0.0_kd, p4_s(192, 2) = 0.0_kd, p5_s(192, 2) = 0.0_kd
        integer :: icritical_band, iwin
        !---------------------------------------------------------------------------------------------
        ! masking / allowed noise : reference bosse lincoln "an experimental high fidelity perceptual audio coder project in mus420 win 97"
        !---------------------------------------------------------------------------------------------
        d2phi_l = phi_l(:, igranule, ichannel) + phi2_l(:, ichannel) - 2.0_kd * phi1_l(:, ichannel) 
        p0_l(:, ichannel) = mod( abs(d2phi_l), pi2 ) / pi2  
        tone_l = 1.0_kd - max(p0_l(:, ichannel), p1_l(:, ichannel), p2_l(:, ichannel) )  
        ! mask for long block
        fk_l =  0.3_kd * tone_l +  0.5_kd * (1 - tone_l) 
        fl_l = 34.0_kd * tone_l + 20.0_kd * (1 - tone_l) 
        tn_l = 10.0_kd**( - ( fk_l * bark_l + fl_l + offset ) / 20.0_kd )  
!        x0_l = matmul(sf_l, afft_l(:, igranule, ichannel) * tn_l * weight_l) + ath_l
        x0_l = matmul(sf_l, afft_l(:, igranule, ichannel) * tn_l * weight_l ** alpha)**(1/alpha) + ath_l
        x1_l = tempo * x1_l + (1.0_kd - tempo) * x0_l ! temporal masking
        x0_l = max(x0_l, x1_l)
        ! allowed noise for long block : average mask over a critical band
        do icritical_band = 1, 25
!            yn(icritical_band) = sum( x0_l(ifb_l(icritical_band, 1):ifb_l(icritical_band, 2)) ) &
!                               / real(ifb_l(icritical_band, 0), kind = kd)
            yn(icritical_band) = sum( x0_l(ifb_l(icritical_band, 1):ifb_l(icritical_band, 2))**alpha)**(1/alpha) &
                                  / real(ifb_l(icritical_band, 0), kind = kd) 
        end do
        y0_l = max( ath_l, yn(ibark_l) ) 
        ! save old data
        p1_l(:, ichannel) = p0_l(:, ichannel)
        p2_l(:, ichannel) = p1_l(:, ichannel)
        phi1_l(:, ichannel) = phi_l(:, igranule, ichannel)
        phi2_l(:, ichannel) = phi1_l(:, ichannel)
        !
        do iwin = 1, 3
            d2phi_s = phi_s(:, iwin, igranule, ichannel) + phi2_s(:, ichannel) - 2.0_kd * phi1_s(:, ichannel)
            p0_s(:, ichannel) = mod( abs(d2phi_s), pi2 ) / pi2
            tone_s = 1.0_kd - max( p0_s(:, ichannel), p1_s(:, ichannel), p2_s(:, ichannel), & 
                                   p3_s(:, ichannel), p4_s(:, ichannel), p5_s(:, ichannel)  )
           ! mask for short block
            fk_s =  0.3_kd * tone_s +  0.5_kd * (1 - tone_s)  
            fl_s = 34.0_kd * tone_s + 20.0_kd * (1 - tone_s) 
            tn_s = 10.0_kd**( - ( fk_s * bark_s + fl_s + offset ) / 20.0_kd )
!            x0_s(:, iwin) = matmul(sf_s, afft_s(:, iwin, igranule, ichannel) * tn_s * weight_s) + ath_s(:, iwin)
            x0_s(:, iwin) = matmul(sf_s, afft_s(:, iwin, igranule, ichannel) * tn_s * weight_s **alpha)**(1/alpha) + ath_s(:, iwin)
           ! save old data
            p1_s(:, ichannel) = p0_s(:, ichannel)
            p2_s(:, ichannel) = p1_s(:, ichannel)
            p3_s(:, ichannel) = p2_s(:, ichannel)
            p4_s(:, ichannel) = p3_s(:, ichannel)
            p5_s(:, ichannel) = p4_s(:, ichannel)
            phi1_s(:, ichannel) = phi_s(:, iwin, igranule, ichannel)
            phi2_s(:, ichannel) = phi1_s(:, ichannel)
        end do
        x1_s(:, 1) = tempo * x1_s(:, 3) + (1.0_kd - tempo) * x0_s(:, 1) ! temporal masking   
        x1_s(:, 2) = tempo * x1_s(:, 1) + (1.0_kd - tempo) * x0_s(:, 2) ! 
        x1_s(:, 3) = tempo * x1_s(:, 2) + (1.0_kd - tempo) * x0_s(:, 3) !  
        x0_s  = max(x0_s, x1_s)
        ! allowed noise for short block  : average mask over a critical band
        do iwin = 1, 3
            do icritical_band = 1, 25
                yn(icritical_band) = sum( x0_s(ifb_s(icritical_band, 1):ifb_s(icritical_band, 2), iwin)**alpha ) **(1 / alpha) &
                                     / real(ifb_s(icritical_band, 0), kind = kd) 
            end do
            y0_s(:, iwin) = max( ath_s(:, iwin), yn(ibark_s(:)) ) 
        end do
        ! order to r_mdct style
        select case (mblock_type)
        case (0, 10, 11, 30, 31)
            call deorder_l(x0_l, xmask )
            call deorder_l(y0_l, xnoise)
        case (20)
            call deorder_s(x0_s, xmask )
            call deorder_s(y0_s, xnoise)
        case (21) 
            call deorder_m(x0_l, x0_s, xmask )
            call deorder_m(y0_l, y0_s , xnoise)
        case default
            stop 'subroutine: calc_noise'
        end select
    end subroutine calc_mask
!------------------------------------------------------------------------------------------------
    subroutine deorder_l(zl, xth)
        real (kind = kd), intent(in ) :: zl(:)
        real (kind = kd), intent(out) :: xth(:, :)
        integer :: i, iband
        forall(iband = 1:32, i = 1:18) xth(iband, i) = zl(18 * (iband - 1) + i)
    end subroutine deorder_l
!------------------------------------------------------------------------------------------------
    subroutine deorder_s(zs, xth)
        real (kind = kd), intent(in ) :: zs(:, :)
        real (kind = kd), intent(out) :: xth(:, :)
        integer :: i, iwin, iband
        forall(iband = 1:32, iwin = 1:3, i = 1:6) xth(iband, 6 * (iwin - 1) + i) = zs(6 * (iband - 1) + i, iwin) 
    end subroutine deorder_s
!------------------------------------------------------------------------------------------------
    subroutine deorder_m(zl, zs, xth)
        real (kind = kd), intent(in ) :: zl(:), zs(:, :)
        real (kind = kd), intent(out) :: xth(:, :)
        integer :: i, iwin, iband
        forall(iband = 1:2, i = 1:18) xth(iband, i) = zl( 18 * (iband - 1) + i)
        forall(iband = 3:32, iwin = 1:3, i = 1:6) xth(iband, 6 * (iwin - 1) + i) = zs(6 * (iband - 1) + i, iwin) 
    end subroutine deorder_m
 !------------------------------------------------------------------------------------------------
    function bark(f) result(res)
        real (kind = kd), intent(in) :: f
        real (kind = kd) :: res
        res = 13.0_kd * atan(0.76_kd * f) + 3.5_kd * atan( (f / 7.5_kd)**2.0_kd )
    end function bark
 !------------------------------------------------------------------------------------------------
    function spreading_function0(z) result(res)
        real (kind = kd), intent(in) :: z
        real (kind = kd) :: res
        res = 15.81_kd + 7.5_kd * (z + 0.474_kd) - 17.5_kd * sqrt(1.0_kd + (z + 0.474_kd)**2.0_kd) 
    end function spreading_function0
 !------------------------------------------------------------------------------------------------
    function spreading_function(z) result(res)
        real (kind = kd), intent(in) :: z
        real (kind = kd) :: res
       ! if (z < 0.0_kd) then
       !     res =  25.0_kd * z
       ! else
       !     res = -10.0_kd * z
       !end if 
        if ( z > 0.0_kd ) then
            res = -25.0_kd * z 
        else
            res =  75.0_kd * z 
        end if
    end function spreading_function
!------------------------------------------------------------------------------------------------
    function spreading_function2(z, y) result(res)
        real (kind = kd), intent(in) :: z, y
        real (kind = kd) :: res
        res = (15.81_kd - y) + 7.5_kd * (z + 0.474_kd) - (17.5_kd - y) * sqrt(1.0_kd + (z + 0.474_kd)**2.0_kd) 
    end function spreading_function2
!------------------------------------------------------------------------------------------------
end module mod_psycho