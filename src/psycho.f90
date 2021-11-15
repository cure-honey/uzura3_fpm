module mod_psycho
    use kind_m
    use mod_mpg
    use mod_fft
    implicit none
    private
    public :: psycho, calc_mask  ! subroutine
    real (kind = kd), save :: pi, pi2 
    real (kind = kd), save :: ath_l(576) , ath_s(192, 3)
    real (kind = kd), save :: sf_l(576, 576), sf_s(192, 192)
    real (kind = kd), save :: dbsf_l(576, 576), dbsf_s(192, 192) 
    real (kind = kd), save :: afft_l(576, 2, 2), afft_s(192, 3, 2, 2)
    real (kind = kd), save :: dbfft_l(576, 2, 2), dbfft_s(192, 3, 2, 2)
    real (kind = kd), save :: arg_l(576, 2, 2), arg_s(192, 3, 2, 2)
    real (kind = kd), save :: freq_l(576), freq_s(192), bark_l(576), bark_s(192), bw_l(576), bw_s(192)
    real (kind = kd), save :: weight_l(576), weight_s(192)
    complex (kind = kd), save :: cfft_l(1152, 2, 2),  cfft_s(384, 3, 2, 2)
    integer, save          :: ibark_l(576), ibark_s(192), ifb_l(25, 0:2), ifb_s(25, 0:2)
    !
contains
    !----------------------------------------------------------------
    subroutine init_absolute_threshold(isample_rate)
        integer, intent(in ) :: isample_rate
        real(kind = kd):: freq, temp !, ath(576)
        integer :: i, k
        pi  = 4.0_kd * atan(1.0_kd)
        pi2 = 2.0_kd * pi
        !
        do i = 1, 576
            freq = real(isample_rate, kind = kd) / 2.0_kd / 1000.0_kd * (real(i - 1, kind = kd) + 0.5_kd) / 576.0_kd
        !    temp =  3.64_kd  * freq ** (-0.8_kd) & 
        !         -  6.50_kd  * exp(-0.6_kd * (freq -  3.3_kd)**2.0_kd) &
        !         +  0.001_kd * freq ** 4.0_kd 
        !         + ath_min
            temp = 3.64_kd * freq ** (-0.8_kd) &                    ! alternative ath function 
                 - 6.50_kd * exp(-0.6_kd * (freq - 3.3_kd)**2.0_kd) !& ! reference: lame ath-type 3  
           !      + 5.1_kd 
            temp = 0.001_kd * freq ** 3.80_kd + 4.92_kd  ! 12.61    
            !  temp = 0.001_kd * freq ** 3.81_kd + 4.60_kd  ! 12.93     
            !  temp = 0.001_kd * freq ** 3.82_kd + 4.27_kd  ! 13.26     
            !  temp = 0.001_kd * freq ** 3.83_kd + 3.94_kd  ! 13.59     
            !  temp = 0.001_kd * freq ** 3.84_kd + 3.61_kd  ! 13.93     
            !  temp = 0.001_kd * freq ** 3.85_kd + 3.27_kd  ! 14.28     
        !      temp = 0.001_kd * freq ** 4.00_kd - 2.86_kd  ! 20.70     
            temp = min(temp + ath_min, ath_max)  
            ath_l(i) = 10.0_kd**(temp / 10.0_kd) ! intensity  
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
        real (kind = kd), save :: sum0a, sum1a = 0, sum0b, sum1b = 0
        sum0a = sum1a
        sum1a = sum( abs(wx(1:36, :)) )
        sum0b = sum1b
        sum1b = sum( abs(wx(37: , :)) )
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
    ! not used; not working correctly ; mdct based version is better  
    subroutine mid_side_fft(mpg, mblock_type) ! iso  c.2.4.3.4.9.2,  g.2 ms_stereo and intensity stereo coding layer iii
        type (mpeg_parameters), intent(in out) :: mpg
        integer               , intent(in    ) :: mblock_type(:, :)
        integer          :: igranule, nchannel, n0, n1
        real (kind = kd) :: cos_t, x1, x2
        integer, save :: mode_old = 0, mode_ext_old = 0
        logical       :: qms
        nchannel = size(mblock_type, 2)
        qms = qms_stereo
        !
        x1 = sqrt( sum(conjg(cfft_l(:, :, 1)) * cfft_l(:, :, 1)) )
        x2 = sqrt( sum(conjg(cfft_l(:, :, 2)) * cfft_l(:, :, 2)) )
        cos_t =    sum(conjg(cfft_l(:, :, 1)) * cfft_l(:, :, 2)) / (x1 * x2 + tiny(0.0_kd))
        if (abs(cos_t) < 0.76_kd .or. min(x1, x2) / max(x1, x2) < 0.5_kd) qms = .false.
        ns1 = ns1 + 1
        ns  = ns  + 1
!
!
        if (qms) then
            print *, qms, ':', cos_t, x1, x2, min(x1, x2) / max(x1, x2)
            do igranule = 1, 2
                x1 = sqrt( dot_product(cfft_l(:, igranule, 1), cfft_l(:, igranule, 1)) ) ! L
                x2 = sqrt( dot_product(cfft_l(:, igranule, 2), cfft_l(:, igranule, 2)) ) ! R
                cos_t = dot_product(cfft_l(:, igranule, 1), cfft_l(:, igranule, 2)) / (x1 * x2 + tiny(0.0_kd))
                if (abs(cos_t) < 0.76_kd .or. min(x1, x2) / max(x1, x2) < 0.5_kd) qms = .false.
            !    ns1 = ns1 + 1
            !    ns  = ns  + 1
                print *, qms, ':', cos_t, x1, x2, min(x1, x2) / max(x1, x2)
            end do 
            print  *
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
    end subroutine mid_side_fft
!------------------------------------------------------------------------------------------------
    subroutine fft_long(nchannel, pcm, afft_l, phi_l, fft576)
        integer         , intent(in ) :: nchannel
        real (kind = kd), intent(in ) :: pcm(:, :)
        real (kind = kd), intent(out) :: afft_l(:, :, :), phi_l(:, :, :)
        complex (kind = kd), intent(out) :: fft576(1152, 2, 2)
        integer :: ichannel, igranule, m1,m2
        do igranule = 1, 2
            do ichannel = 1, nchannel
                m1 = 1 + 480 * (igranule - 1)
                m2 = m1 + 1152 - 1
                fft576(:, igranule, ichannel) = cmplx(pcm(m1:m2, ichannel), 0.0_kd, kind = kd) ! put 1152 real data -> get 576 complex fft
             !  call fft23han(7, 2, indx576, omega576, fft576(:, igranule, ichannel), han576 ) ! 2^7 * 3^2 = 1152
                call fft23(7, 2, indx576, omega576, fft576(:, igranule, ichannel) )
                afft_l(:, igranule, ichannel) = abs(fft576(1:576, igranule, ichannel))
                arg_l (:, igranule, ichannel) = atan2(aimag(fft576(1:576, igranule, ichannel)), &
                                                       real(fft576(1:576, igranule, ichannel)))
            end do
        end do
        phi_l = phi_l + pi ! atan -pi~pi -> 0~2pi
    end subroutine fft_long
!------------------------------------------------------------------------------------------------
    subroutine fft_short(nchannel, pcm, afft_s, phi_s, fft192)
        integer       , intent(in ) :: nchannel
        real (kind = kd), intent(in ) :: pcm(:, :)
        real (kind = kd), intent(out) :: afft_s(:, :, :, :), phi_s(:, :, :, :)
        complex (kind = kd), intent(out) :: fft192(384, 3, 2, 2)
        integer :: ichannel, igranule, iwin, m1, m2
        do igranule = 1, 2
            do ichannel = 1, nchannel
               do iwin = 1, 3
                    m1 = 1 + 480 * (igranule - 1) + 384 * (iwin - 1)
                    m2 = m1 + 384 - 1
                    fft192(:, iwin, igranule, ichannel) = cmplx(pcm(m1:m2, ichannel), 0.0_kd, kind = kd) ! put 384 real data -> get 192 complex fft
                 !  call fft23han(7, 1, indx192, omega192, fft192(:, iwin, igranule, ichannel), han192 ) ! 2^7 * 3^1 = 384
                    call fft23(7, 1, indx192, omega192, fft192(:, iwin, igranule, ichannel) )
                    afft_s(:, iwin, igranule, ichannel) = abs(fft192(1:192, iwin, igranule, ichannel))
                    arg_s (:, iwin, igranule, ichannel) = atan2(aimag(fft192(1:192, iwin, igranule, ichannel)), &
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
                wx(:, igranule, ichannel) = afft_l(:, igranule, ichannel) * weight_l
            end do
        end do
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
        !..... fft 576/192 ............................................................................................
        nchannel = size(mblock_type, 2) 
        call fft_long (nchannel, pcm, afft_l, arg_l, cfft_l)
        call fft_short(nchannel, pcm, afft_s, arg_s, cfft_s)
        !..... weighted intensity .....................................................................................
        call calc_wx(nchannel, wx)
        !..... attack detection .......................................................................................
        call attack(wx, mblock_type)
        !..... ms/ns selection ........................................................................................
        !call mid_side_fft(mpg, mblock_type) ! not working correctly : mdct based model works better
        !
        ! debug info
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
            weight_l(i) = bw_l(i) / (f1 - f0) 
        end do
        do i = 1, 192
            freq_s(i) = real(nsample_rate, kind = kd) / 2.0_kd * (real(i, kind = kd) - 0.5_kd) / 192.0_kd ! khz
            bark_s(i) = bark(freq_s(i) / 1000.0_kd)  
            ibark_s(i) = int(bark_s(i)) + 1
            f0 = real(nsample_rate, kind = kd) / 2000.0_kd * real(i - 1, kind = kd) / 192.0_kd
            f1 = real(nsample_rate, kind = kd) / 2000.0_kd * real(i    , kind = kd) / 192.0_kd
            bw_s(i) = bark(f1) - bark(f0)
            weight_s(i) = bw_s(i) / (f1 - f0) 
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
                dbsf_l(i, j) = spreading_function( bark_l(i) - bark_l(j) )  
                sf_l(i, j) = 10.0_kd ** ( dbsf_l(i, j) / 20.0_kd ) 
            end do
        end do 
        do i = 1, 192
            do j = 1, 192
                dbsf_s(i, j) = spreading_function( bark_s(i) - bark_s(j) )
                sf_s(i, j) = 10.0_kd ** ( dbsf_s(i, j) / 20.0_kd ) 
            end do
        end do 
    end subroutine init_mask
!------------------------------------------------------------------------------------------------
    subroutine calc_mask(igranule, ichannel, mblock_type, xmask, xnoise)
        integer         , intent(in    ) :: igranule, ichannel, mblock_type
        real (kind = kd), intent(   out) :: xmask(:, :), xnoise(:, :)
        real (kind = kd) :: x0_l(576), x_l(576), y_l(576)
        real (kind = kd) :: tone_l(576), fk_l(576), fl_l(576), tn_l(576), sf(576, 576)
        real (kind = kd) :: yn(25)
        real (kind = kd), save :: x1_l(576) = 0.0_kd
        real (kind = kd), save :: phi1_l(576, 2) = 0.0_kd, phi2_l(576, 2) = 0.0_kd, afft_l0(576, 2) = 0.0_kd
        real (kind = kd), save :: arg1_l(576, 2) = 0.0_kd, arg2_l(576, 2) = 0.0_kd, arg3_l(576, 2) = 0.0_kd
        real (kind = kd), save :: err1_l(576, 2) = 0.0_kd, err2_l(576, 2) = 0.0_kd, errm_l(576, 2) = 0.0_kd
        integer :: icritical_band, iwin, iband, i, j
        real (kind = kd) :: tmp1
        !---------------------------------------------------------------------------------------------
        ! masking / allowed noise : reference Bosse Lincoln "an experimental high fidelity perceptual audio coder project in mus420 win 97"
        !---------------------------------------------------------------------------------------------
        !
        ! tonality estimation : angles (err1,2) & radius (errm)
        phi1_l(:, ichannel) = 2 * arg1_l(:, ichannel) - arg2_l(:, ichannel)   
        phi2_l(:, ichannel) = 2 * arg2_l(:, ichannel) - arg3_l(:, ichannel)   
        err1_l(:, ichannel) = mod(abs(arg_l (:, igranule, ichannel) - phi1_l(:, ichannel) ), pi) / pi 
        err2_l(:, ichannel) = mod(abs(arg1_l(:, ichannel)           - phi2_l(:, ichannel) ), pi) / pi 
        forall (i = 1:576)
            errm_l(i, ichannel) = abs( afft_l0(i, ichannel) - afft_l(i, igranule, ichannel) ) &
                          / (abs( max( afft_l0(i, ichannel),  afft_l(i, igranule, ichannel) ) ) + tiny(0.0_kd)) 
        end forall
        tone_l = 1.0_kd - max(err1_l(:, ichannel), err2_l(:, ichannel), errm_l(:, ichannel) )  
        !
        ! mask for long block
        fk_l =  0.3_kd * tone_l +  0.5_kd * (1 - tone_l) 
        fl_l = 34.0_kd * tone_l + 20.0_kd * (1 - tone_l) 
        if (.not. q_pnorm) then
            ! 1-norm 1 
            forall(i = 1:576, j = 1:576) sf(i, j) = sf_l(i, j) * 10.0_kd **(-fk_l(i) * bark_l(j) / 20.0_kd)
            tn_l = 10.0_kd**( - ( fl_l - offset_l - weight_l) / 20.0_kd ) ! [dB]
            x0_l = matmul(sf, afft_l(:, igranule, ichannel)**2) * tn_l 
        else               
            dbfft_l = 20.0_kd * log10(afft_l)
            do i = 1, 576
                tn_l = 2 * dbfft_l(:, igranule, ichannel)  + dbsf_l(i, :) - fk_l(i) * bark_l  - fl_l(i) + weight_l ![dB]
                ! 1-norm 2
!                x0_l(i) = sum( 10.0_kd**(tn_l - offset_l/ 20.0_kd) )  
                ! 2/3-norm
!                tmp1 = sum( 10.0_kd**(tn_l / 30.0_kd) )
!                x0_l(i) = sqrt(tmp1)*tmp1*fctr  
                ! 1/2-norm
!                tmp1 = sum( 10.0_kd**(tn_l / 40.0_kd) )
!                x0_l(i) = tmp1*tmp1*fctr
                ! p-norm  
                 x0_l(i) = sum( 10.0_kd**(tn_l * pow / 20.0_kd) )**( 1.0_kd / pow) * 10.0_kd**(offset_l / 20.0_kd)  
            end do
        end if
        !
        ! temporal masking        
        x1_l = tempo_l * x1_l + (1.0_kd - tempo_l) * x0_l 
        x0_l = max(x0_l, x1_l)
        !
        ! allowed noise : average mask over a critical band
        do icritical_band = 1, 25
            yn(icritical_band) = sum( x0_l(ifb_l(icritical_band, 1):ifb_l(icritical_band, 2)) ) &
                                    / real(ifb_l(icritical_band, 0), kind = kd)
        end do
        x_l = max( ath_l, x0_l )
        y_l = max( ath_l, yn(ibark_l) ) 
        !
        ! reorder
        select case (mblock_type)
        case (0) ! long block
            forall (iband = 1:32, i = 1:18) xmask (iband, i) = x_l(18 * (iband - 1) + i)
            forall (iband = 1:32, i = 1:18) xnoise(iband, i) = y_l(18 * (iband - 1) + i)
        case (10, 11, 30, 31) ! transition block
            forall (iband = 1:32, i = 1:18) xmask (iband, i) = x_l(18 * (iband - 1) + i) 
            forall (iband = 1:32, i = 1:18) xnoise(iband, i) = y_l(18 * (iband - 1) + i) 
        case (20) ! short block
            forall (iband = 1:32, iwin = 1:3, i = 1:6) 
                xmask (iband, 6 * (iwin - 1) + i) = minval(x_l(18 * (iband - 1) + 6 * (i - 1) + 1  &
                                                              :18 * (iband - 1) + 6 * (i - 1) + 3  ) )
                xnoise(iband, 6 * (iwin - 1) + i) = minval(y_l(18 * (iband - 1) + 6 * (i - 1) + 1  &
                                                              :18 * (iband - 1) + 6 * (i - 1) + 3  ) )
            end forall 
        case (21) ! mixed block   
            forall (iband = 1: 2, i = 1:18) xmask (iband, i) = x_l(18 * (iband - 1) + i) 
            forall (iband = 1: 2, i = 1:18) xnoise(iband, i) = y_l(18 * (iband - 1) + i)
            forall (iband = 3:32, iwin = 1:3, i = 1:6) 
                xmask (iband, 6 * (iwin - 1) + i) = minval(x_l(18 * (iband - 1) + 6 * (i - 1) + 1  &
                                                              :18 * (iband - 1) + 6 * (i - 1) + 3  ) )
                xnoise(iband, 6 * (iwin - 1) + i) = minval(y_l(18 * (iband - 1) + 6 * (i - 1) + 1  &
                                                              :18 * (iband - 1) + 6 * (i - 1) + 3  ) )
            end forall 
        case default
            print *, 'mblock_type', mblock_type
            stop 'calc_mask: should not come here'
        end select
        !
        ! save old data
        arg3_l(:, ichannel) = arg2_l(:, ichannel)
        arg2_l(:, ichannel) = arg1_l(:, ichannel)
        arg1_l(:, ichannel) = arg_l (:, igranule, ichannel)
        afft_l0(:, ichannel) = afft_l(:, igranule, ichannel) 
    end subroutine calc_mask
 !------------------------------------------------------------------------------------------------
    function bark(f) result(res)
        real (kind = kd), intent(in) :: f
        real (kind = kd) :: res
        res = 13.0_kd * atan(0.76_kd * f) + 3.5_kd * atan( (f / 7.5_kd)**2.0_kd )
    end function bark
 !------------------------------------------------------------------------------------------------
    pure elemental function spreading_function(z) result(res)
        real (kind = kd), intent(in) :: z
        real (kind = kd) :: res
        res = 15.81_kd + 7.5_kd * (z + 0.474_kd) - 17.5_kd * sqrt(1.0_kd + (z + 0.474_kd)**2.0_kd) 
    end function spreading_function
 !------------------------------------------------------------------------------------------------
    pure elemental function spreading_function0(z) result(res)
        real (kind = kd), intent(in) :: z
        real (kind = kd) :: res
        if ( z > 0.0_kd ) then
            res = -25.0_kd * z 
        else
            res =  75.0_kd * z 
        end if
    end function spreading_function0
 !------------------------------------------------------------------------------------------------
    pure elemental function spreading_function2(z, y) result(res)
        real (kind = kd), intent(in) :: z, y
        real (kind = kd) :: res
        res = (15.81_kd - y) + 7.5_kd * (z + 0.474_kd) - (17.5_kd - y) * sqrt(1.0_kd + (z + 0.474_kd)**2.0_kd) 
    end function spreading_function2
!------------------------------------------------------------------------------------------------
end module mod_psycho