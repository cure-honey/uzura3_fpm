module mod_layer3
    use kind_m
    use mod_mpg
    use mod_psycho
    use mod_inner_loop
    implicit none
    private
    public  :: alloc_bits, init_scalefactor_bands ! subroutine 
    public  :: side_info, scfct                   ! variable
    public  :: scale_factor                       ! type 
    integer, parameter :: npretab(0:20, 0:1) = &  ! iso table b.6 layer iii preemphasis (pretab)
             reshape((/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 3, 2/), shape(npretab))
    !                  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
    type :: scale_factor
        integer :: long(0:20)
        integer :: ishort(0:12, 3)
    end type scale_factor
    !
    type (scale_factor), save :: scfct(2, 2)
    type (side_info_  ), save :: side_info
    integer, save :: iscfband_l(0:20, 3, 0:2), iscfband_s(0:11, 3, 0:2)
contains
    !------------------------------------------------------------------------------------------------
    subroutine init_scalefactor_bands(n) ! iso table b.8 layer iii scalefactor bands
        integer, intent(in) :: n
        integer :: i, j
        iscfband_l = 0  
        iscfband_s = 0
        iscfband_l(:, 1, 2) = (/4, 4, 4, 4, 4, 4, 6, 6, 8, 10, 12, 16, 20, 24, 30, 38, 46, 56, 68, 84, 102/) !32.0khz
        iscfband_l(:, 1, 0) = (/4, 4, 4, 4, 4, 4, 6, 6, 8,  8, 10, 12, 16, 20, 24, 28, 34, 42, 50, 54,  76/) !44.1khz
        iscfband_l(:, 1, 1) = (/4, 4, 4, 4, 4, 4, 6, 6, 6,  8, 10, 12, 16, 18, 22, 28, 34, 40, 46, 54,  54/) !48.0khz
        iscfband_s(:, 1, 2) = (/4, 4, 4, 4, 6, 8, 12, 16, 20, 26, 34, 42/) !32.0khz
        iscfband_s(:, 1, 0) = (/4, 4, 4, 4, 6, 8, 10, 12, 14, 18, 22, 30/) !44.1khz
        iscfband_s(:, 1, 1) = (/4, 4, 4, 4, 6, 6, 10, 12, 14, 16, 20, 26/) !48.0khz
        do j = 0, 2
            do i = 1, 20
                iscfband_l(i    , 2, j) = iscfband_l(i - 1, 2, j) + iscfband_l(i - 1, 1, j)
                iscfband_l(i - 1, 3, j) = iscfband_l(i    , 2, j) - 1
            end do 
            iscfband_l(20, 3, j) = iscfband_l(20, 1, j) + iscfband_l(20, 2, j) - 1
            do i = 1, 11
                iscfband_s(i    , 2, j) = iscfband_s(i - 1, 2, j) + iscfband_s(i - 1, 1, j)
                iscfband_s(i - 1, 3, j) = iscfband_s(i    , 2, j) - 1
            end do 
            iscfband_s(11, 3, j) = iscfband_s(11, 1, j) + iscfband_s(11, 2, j) - 1
        end do
        iscalefactorband_l = iscfband_l(:, :, n)
        iscalefactorband_s = iscfband_s(:, :, n)
    end subroutine init_scalefactor_bands
    !------------------------------------------------------------------------------------------------
    subroutine init_scale_factor()
        integer :: igranule, ichannel
        do ichannel = 1, 2
            do igranule = 1, 2
                scfct(igranule, ichannel)%long   = 0
                scfct(igranule, ichannel)%ishort = 0
            end do
        end do
    end subroutine init_scale_factor
!------------------------------------------------------------------------------------------------
    subroutine init_side_info(mblock_type)
        integer, intent(in) :: mblock_type(:, :)
        integer :: igranule, ichannel
        side_info%main_data_begin = 0                                     !  9 bits
        side_info%iprivate_bits   = 0                                     !  5 / 3 bits
        side_info%iscfsi(4, 2)    = 0                                     !  1 bit 
        do igranule = 1, 2
            do ichannel = 1, size(mblock_type, 2)
                side_info%sub(igranule, ichannel)%ipart2_3_length        = 0    ! 12 bits
                side_info%sub(igranule, ichannel)%ibig_values            = 0    !  9 bits
                side_info%sub(igranule, ichannel)%iglobal_gain           = 0    !  8 bits
                side_info%sub(igranule, ichannel)%iscalefac_compress     = 0    !  4 bits   
                select case (mblock_type(igranule, ichannel))
                case ( 0) ! long-block
                    side_info%sub(igranule, ichannel)%iwindow_switching_flag = 0    !  1 bit    
                    side_info%sub(igranule, ichannel)%iblock_type            = 0    !  2 bits   
                    side_info%sub(igranule, ichannel)%mixied_block_flag      = 0    !  1 bit    
                case (10) ! start-block for short
                    side_info%sub(igranule, ichannel)%iwindow_switching_flag = 1    !  1 bit    
                    side_info%sub(igranule, ichannel)%iblock_type            = 1    !  2 bits   
                    side_info%sub(igranule, ichannel)%mixied_block_flag      = 0    !  1 bit    
                case (11) ! start-block for mixed
                    side_info%sub(igranule, ichannel)%iwindow_switching_flag = 1    !  1 bit    
                    side_info%sub(igranule, ichannel)%iblock_type            = 1    !  2 bits   
                    side_info%sub(igranule, ichannel)%mixied_block_flag      = 1    !  1 bit    
                case (30) ! stop-block for short
                    side_info%sub(igranule, ichannel)%iwindow_switching_flag = 1    !  1 bit    
                    side_info%sub(igranule, ichannel)%iblock_type            = 3    !  2 bits   
                    side_info%sub(igranule, ichannel)%mixied_block_flag      = 0    !  1 bit    
                case (31) ! stop-block for short
                    side_info%sub(igranule, ichannel)%iwindow_switching_flag = 1    !  1 bit    
                    side_info%sub(igranule, ichannel)%iblock_type            = 3    !  2 bits   
                    side_info%sub(igranule, ichannel)%mixied_block_flag      = 1    !  1 bit    
                case (20) ! short-block
                    side_info%sub(igranule, ichannel)%iwindow_switching_flag = 1    !  1 bit    
                    side_info%sub(igranule, ichannel)%iblock_type            = 2    !  2 bits   
                    side_info%sub(igranule, ichannel)%mixied_block_flag      = 0    !  1 bit    
                case (21) ! mixed-block
                    side_info%sub(igranule, ichannel)%iwindow_switching_flag = 1    !  1 bit    
                    side_info%sub(igranule, ichannel)%iblock_type            = 2    !  2 bits   
                    side_info%sub(igranule, ichannel)%mixied_block_flag      = 1    !  1 bit    
                case default
                    stop ' error : mblock_type '
                end select
                side_info%sub(igranule, ichannel)%itable_select(:)       = 0    !  5 bits   
                side_info%sub(igranule, ichannel)%isubblock_gain(:)      = 0    !  3 bits   
                side_info%sub(igranule, ichannel)%iregion0_count         = 0    !  4 bits
                side_info%sub(igranule, ichannel)%iregion1_count         = 0    !  3 bits
                side_info%sub(igranule, ichannel)%ipreflag               = 0    !  1 bit    
                side_info%sub(igranule, ichannel)%iscalefac_scale        = 0    !  1 bit   
                side_info%sub(igranule, ichannel)%icount1table_select    = 0    !  1 bit
                !.....
                side_info%sub(igranule, ichannel)%icount1                = 0    ! local use
            end do
        end do
    end subroutine init_side_info
!------------------------------------------------------------------------------------------------
    subroutine mid_side_select_mdct(mpg, mblock_type, r_mdct) ! iso  c.2.4.3.4.9.2,  g.2 ms_stereo and intensity stereo coding layer III
        type (mpeg_parameters), intent(in out) :: mpg
        integer               , intent(in    ) :: mblock_type(:, :)
        real (kind = kd)      , intent(in    ) :: r_mdct(:, :, :, :)
        integer          :: igranule, nchannel
        real (kind = kd) :: tmp1, tmp2, tmp3, tmp4, pi
        logical          :: qms
        integer :: n0
        pi = 4.0_kd * atan(1.0_kd)
        nchannel = size(r_mdct, 4)
        !
        if ( any(mblock_type(:, 1) /= mblock_type(:, 2)) & !  different block for LR channels 
        .or. any(mblock_type == 10) .or. any(mblock_type == 11) .or. any(mblock_type == 30) .or. any(mblock_type == 31) ) then 
            qms = .false.
            ! debug info
            ns1 = ns1 + 1
            ns  = ns  + 1
        else    
            qms = qms_stereo
            select case (mpg%isample_rate) ! threshold ~ 7khz (empirical value)
            case (0) ! 44.1khz
                n0 = 10
            case (1) ! 48.0khz
                n0 = 9
            case (2) ! 32.0khz
                n0 = 14
            case default
                stop ' sample_rate error : subroutine mid_side '
            end select
            !
            ! empirical 
            do igranule = 1, 2
                tmp1 = sum( abs( r_mdct(1:n0, :, igranule, 1) - r_mdct(1:n0, :, igranule, 2) ) ) 
                tmp2 = sum( abs( r_mdct(1:n0, :, igranule, 1) + r_mdct(1:n0, :, igranule, 2) ) )
                tmp3 = sum( abs( r_mdct(1:n0, :, igranule, 1)) )
                tmp4 = sum( abs( r_mdct(1:n0, :, igranule, 2)) )
            end do    

            if (  min(tmp1, tmp2) > xms * max(tmp1, tmp2) .or. &
            .not. min(tmp3, tmp4) > xms * max(tmp3, tmp4)        ) then ! normal stereo
                qms = .false.
            ! debug info
                ns2 = ns2 + 1
                ns  = ns  + 1
            end if 
            !
            do igranule = 1, 2
                tmp1 = sum( abs( r_mdct(:, :, igranule, 1) - r_mdct(:, :, igranule, 2) ) ) 
                tmp2 = sum( abs( r_mdct(:, :, igranule, 1) + r_mdct(:, :, igranule, 2) ) )
                tmp3 = sum( abs( r_mdct(:, :, igranule, 1)) )
                tmp4 = sum( abs( r_mdct(:, :, igranule, 2)) )
            end do    
            if (  min(tmp1, tmp2) > 0.95_kd * xms * max(tmp1, tmp2) .or. &
            .not. min(tmp3, tmp4) > 0.95_kd * xms * max(tmp3, tmp4)        ) then ! normal stereo
            ! debug info
                if (qms) then 
                    ns  = ns  + 1
                    ns2 = ns2 + 1
                end if    
            ! end debug info   
                qms = .false.
            end if
        end if     
        ! debug info
        if (qms) then 
            mpg%mode            =  1 ! joint stereo
            mpg%mode_extension  =  2 ! intensity_stereo off / ms_stereo on
            ms = ms + 1
        else
            mpg%mode            =  0 ! normal stereo
            mpg%mode_extension  =  0 ! intensity_stereo off / ms_stereo off
        end if
    end subroutine mid_side_select_mdct
!------------------------------------------------------------------------------------------------
    subroutine change_to_mid_side(mpg, r_mdct, ms_mdct) ! iso  c.2.4.3.4.9.2,  g.2 ms_stereo and intensity stereo coding layer III
        real (kind = kd)      , intent(in    ) ::  r_mdct(:, :, :, :)
        real (kind = kd)      , intent(   out) :: ms_mdct(:, :, :, :)
        type (mpeg_parameters), intent(in out) :: mpg
        real (kind = kd) :: tmp_ms(32, 18, 2, 2)
        integer :: nchannel
        nchannel = size(r_mdct, 4)
        if (mpg%mode ==  1) then 
            if (nchannel == 2) then
                tmp_ms(:, :, :, 1) = ( r_mdct(:, :, :, 1) + r_mdct(:, :, :, 2) ) / sqrt(2.0_kd)
                tmp_ms(:, :, :, 2) = ( r_mdct(:, :, :, 1) - r_mdct(:, :, :, 2) ) / sqrt(2.0_kd)
            else
                tmp_ms(:, :, :, 1) = r_mdct(:, :, :, 1) / sqrt(2.0_kd)
                tmp_ms(:, :, :, 2) = r_mdct(:, :, :, 1) / sqrt(2.0_kd)
            end if
            ms_mdct = tmp_ms
        else
            ms_mdct = r_mdct
        end if
    end subroutine change_to_mid_side
!------------------------------------------------------------------------------------------------
    subroutine calc_1norm(r_mdct, x_1norm)
        real (kind = kd) , intent(in ) :: r_mdct (:, :, :, :)
        real (kind = kd) , intent(out) :: x_1norm(:, :)
        integer :: igranule, ichannel
        do igranule = 1, 2
            do ichannel = 1, size(r_mdct, 4)
                x_1norm(igranule, ichannel) = sum( abs(r_mdct(:, :, igranule, ichannel)) ) 
!                x_1norm(igranule, ichannel) = sum( log(abs(r_mdct(:, :, igranule, ichannel))), &
!                                                mask = abs(r_mdct(:, :, igranule, ichannel)) > 0.0_kd ) 
            end do
        end do
!        x_1norm = exp(x_1norm / sum(x_1norm))
    end subroutine calc_1norm
!------------------------------------------------------------------------------------------------
    subroutine get_maxbits(mpg, max_bits) ! after iso 2.4.3.1 & 2.4.2.3 padding (p.22)
        type (mpeg_parameters), intent(in out) :: mpg
        integer               , intent(   out) :: max_bits
        integer                                :: idiff, islot_size
        integer, save                          :: irest = 0
        islot_size = 144000 * mpeg_bit_rates(mpg%ibit_rate, mpg%layer) / mpeg_sample_rates(mpg%isample_rate) 
        idiff  = mod(144000 * mpeg_bit_rates(mpg%ibit_rate, mpg%layer),  mpeg_sample_rates(mpg%isample_rate))
        irest  = irest - idiff
        if (irest < 0) then     
            mpg%ipadding = 1
            irest = irest + mpeg_sample_rates(mpg%isample_rate)
        else
            mpg%ipadding = 0
        end if
        max_bits = ( islot_size + mpg%ipadding ) * 8
    end subroutine get_maxbits 
!-------------------------------------------------------------------------------------------------
    subroutine mean_bits(mpg, max_bits, x_norm1, mbits, iused_bits)
        type (mpeg_parameters), intent(in ) :: mpg
        integer               , intent(in ) :: max_bits
        real (kind = kd)      , intent(in ) :: x_norm1(:, :)
        integer               , intent(out) :: mbits(:, :), iused_bits
        integer                             :: nbits, nchannel, limbits, ipos(2)
        nchannel = size(mbits, 2)
        if (nchannel == 1) then
            iused_bits =  32 + 136 ! monoral
        else
            iused_bits =  32 + 256 ! stereo 
        end if
        if (mpg%icrc == 0) iused_bits = iused_bits + 16
        nbits = max_bits - iused_bits
!  distributes bits between 2 granules * n channels according to sum of ideal quantization steps
        mbits = int(nbits * x_norm1 / sum(x_norm1 + epsilon(0.0_kd)))
        where (mbits < 20) mbits = 0
!  where (mbits /= 0) mbits = mbits + (nbits - sum(mbits)) / max(count(mbits /= 0), 1)
        ipos = minloc(mbits, mask = mbits /= 0)
        if (all(ipos /= 0)) mbits(ipos(1), ipos(2)) = mbits(ipos(1), ipos(2)) + nbits - sum(mbits)
!  avoids distributing too many bits at high bitrate   
        limbits = 4000  
        do while (maxval(mbits) > limbits) 
            where (mbits > limbits) mbits = limbits 
            where (mbits < limbits) mbits = (nbits - limbits * count(mbits == limbits)) / count(mbits /= limbits)
        end do
    end subroutine mean_bits
!------------------------------------------------------------------------------------------------
    subroutine alloc_bits(mblock_type, r_mdct, i_mdct, mpg, max_bits, ianc) ! iso c.1.5.4 
        integer               , intent(in    ) :: mblock_type(:, :)
        real (kind = kd)      , intent(in    ) :: r_mdct(:, :, :, :)
        integer               , intent(   out) :: i_mdct(:, :, :), max_bits, ianc
        type (mpeg_parameters), intent(in out) :: mpg
        integer :: ibit, ibit2, ichannel, nchannel, igranule, iused_bits, nused_bits, & 
                   mbits(size(r_mdct, 3), size(r_mdct, 4))
        real (kind = kd) ::  wk_mdct(size(r_mdct, 1), size(r_mdct, 2), size(r_mdct, 3), size(r_mdct, 4))
        real (kind = kd) ::  x_mdct (size(r_mdct, 1), size(r_mdct, 2), size(r_mdct, 3), size(r_mdct, 4))
        real (kind = kd) ::  x_mask (size(r_mdct, 1), size(r_mdct, 2), size(r_mdct, 3), size(r_mdct, 4))
        real (kind = kd) ::  z_mask (size(r_mdct, 1), size(r_mdct, 2))
        real (kind = kd) ::  x_1norm(size(r_mdct, 3), size(r_mdct, 4)), distortion
        logical :: qmask(2, size(r_mdct, 4))
        integer :: i, ipos(2)
        i_mdct = 0
        nchannel = size(r_mdct, 4) 
        call init_scale_factor()
        call init_side_info(mblock_type)
        do igranule = 1, 2
            do ichannel = 1, nchannel   
                call calc_mask(igranule, ichannel, mblock_type(igranule, ichannel), x_mask(:, :, igranule, ichannel))
            end do 
        end do
!--- LR <--> MS decision     
        call mid_side_select_mdct(mpg, mblock_type, r_mdct)
        x_mdct = r_mdct       
!--- cut masked signals       
        if (q_mask) where(abs(x_mdct) < x_mask) x_mdct = 0.0_kd 
        ! help for parameter decision
        if (q_help .and. mod(iframe, 200) == 0) print *, count(abs(x_mdct) < x_mask), 'abs(x_mdct) < x_mask'
!--- bit allocatiion among 2 granules * nchannel : critical 
        call change_to_mid_side(mpg, x_mdct / max(x_mask, epsilon(0.0_kd)), wk_mdct)    !
        call calc_1norm(wk_mdct, x_1norm)
        call get_maxbits(mpg, max_bits)
        call mean_bits  (mpg, max_bits, x_1norm, mbits, nused_bits)  
!--- LR => MS / LR => LR 
        call change_to_mid_side(mpg, x_mdct, wk_mdct) 
        ibit2 = 0 ! remain bits
        qmask = .true.
        do i = 1, 2 * nchannel
! maxbits to minbits  ! do igranule = 1, 2 ! do ichannel = 1, nchannel
            ipos = maxloc(mbits, mask = qmask)
            igranule = ipos(1)
            ichannel = ipos(2)
            qmask(igranule, ichannel) = .false.
            ! remain bits re-distribution
            ipos = minloc(mbits, mask = qmask)
            if (i /= 2 * nchannel) mbits(ipos(1), ipos(2)) = mbits(ipos(1), ipos(2)) + ibit2
            ibit2 = 0
            !
            ibit = mbits(igranule, ichannel) + ibit2 
            if ( mpg%mode == 1 ) then ! noise for MS: 1/sqrt(2)*min(noiseR, noiseL)
                z_mask(:, :) = sqrt(0.5_kd) * min(x_mask(:, :, igranule, 1), x_mask(:, :, igranule, 2)) 
            else
                z_mask(:, :) = x_mask(:, :, igranule, ichannel) 
            end if
            if ( x_1norm(igranule, ichannel) /= 0.0_kd ) then
                !--- scale factor loop   
                call outer_loop(ibit, mblock_type(igranule, ichannel), &
                                   wk_mdct (:, :, igranule, ichannel), z_mask(:, :), &
                                    i_mdct (:   , igranule, ichannel), side_info%sub(igranule, ichannel), &
                                      scfct(      igranule, ichannel), iused_bits, distortion)
            else 
                iused_bits = 0
                i_mdct(:, igranule, ichannel) = 0
            end if
            nused_bits = nused_bits + iused_bits
            side_info%sub(igranule, ichannel)%ipart2_3_length = iused_bits  
            ibit2 = ibit - iused_bits
            ! debug info
            ntable (side_info%sub(igranule, ichannel)%itable_select(1)) &
                     = ntable (side_info%sub(igranule, ichannel)%itable_select(1)) + 1
            ntable (side_info%sub(igranule, ichannel)%itable_select(2)) & 
                     = ntable (side_info%sub(igranule, ichannel)%itable_select(2)) + 1
            ntable (side_info%sub(igranule, ichannel)%itable_select(3)) &
                     = ntable (side_info%sub(igranule, ichannel)%itable_select(3)) + 1
            ntab_ab(side_info%sub(igranule, ichannel)%icount1table_select) &
                     = ntab_ab(side_info%sub(igranule, ichannel)%icount1table_select) + 1
        end do
        ianc = max_bits - nused_bits
    end subroutine alloc_bits
!------------------------------------------------------------------------------------------------    
    subroutine outer_loop(iallowed_bits, iblock_type, r_mdct, x_noise, iwk, side, sc_fac, iused_bits, &
                          distortion) ! iso c.1.5.4.3
        integer             , intent(in    ) :: iallowed_bits, iblock_type
        integer             , intent(   out) :: iused_bits
        real (kind = kd)    , intent(in    ) :: r_mdct(:, :), x_noise(:, :)
        real (kind = kd)    , intent(   out) :: distortion
        integer             , intent(   out) :: iwk(:)
        type (side_info_sub), intent(in out) :: side
        type (scale_factor) , intent(in out) :: sc_fac
        integer :: iover_l(0:20), iover_s(0:11, 3) 
        integer :: iscfac_bits, ihuff_bits, ibits_best, iscale, iwk_best(size(iwk))
        real (kind = kd)    :: wk(size(iwk)), th(size(iwk)), sc(size(iwk)), swk(size(iwk))
        real (kind = kd)    :: distortion_min
        type (side_info_sub):: side_best
        type (scale_factor) :: sc_fac_best
        logical             :: qexit, qfirst
        distortion_min = 1.0d10
        ibits_best = 0
        iwk_best = 0
        !
        scale0:do iscale = 0, 1                     ! scalefactor_scale loop                          ! iso c.1.5.4.3
            side%iscalefac_scale = iscale           ! 0: sqrt(2) or 1: 2
            side%ipreflag = 0                       ! start with preemphasis off                      ! iso c.1.5.4.3.4
            side%isubblock_gain = 0
            do                                      ! subblock_gain loop (only for short/mixed block) ! iso 2.4.3.4.7.1
                qfirst = .true.
                sc_fac%long = 0                        ! clear scale factor 
                sc_fac%ishort = 0                      ! clear scale factor
                do                                     ! scalefactor loop
                    call select_compress(iblock_type, sc_fac     , side%iscalefac_compress)
                    call calc_scfac_bit (iblock_type, iscfac_bits, side%iscalefac_compress) 
                    call reorder(iblock_type, r_mdct , wk, icut)
                    call reorder(iblock_type, x_noise, th,   32)

                    call calc_scale(iblock_type, iscale, sc_fac, npretab(:, side%ipreflag), side%isubblock_gain, sc)
                    swk = wk * sc ! scale mdct

                    if ( iallowed_bits - iscfac_bits < 0) exit scale0 ! no more bits 
                    ! quantization Huffman coding
                    call inner_loop(iallowed_bits - iscfac_bits, iblock_type, swk, iwk, side, ihuff_bits) ! iso c.1.5.4.3.2
                    ! calculate distortion ; scale factor change direction +1  
                    call calc_distortion(iblock_type, side%iglobal_gain, wk, iwk, th, sc, iover_l, iover_s, distortion) ! iso c.1.5.4.3.3
                    !
                    if ( distortion <= distortion_min ) then      ! save best parameters so far           ! iso c.1.5.4.3.1
                        distortion_min = distortion
                        sc_fac_best    = sc_fac                   ! scale factors
                        side_best      = side                     ! side informations
                        iwk_best       = iwk                      ! 576 quantized data
                        ibits_best     = iscfac_bits + ihuff_bits ! required bits for scale_factors & huffman codes
                    else
                        if ( distortion > skip * distortion_min ) exit ! short cut to avoid meaningless search
                    end if

                    if (q_inner) exit scale0 ! no scale factor fit

                    if ( sum(iover_l) + sum(iover_s) == 0 ) exit scale0 ! if converged return ! iso c.1.5.4.3.6
                    !
                    if ( qfirst .and. any(iover_l(12:20) > 0) ) then   ! iso c.1.5.4.3.4, Table B.6 11..20: 1 1 1 1 2 2 3 3 3 2
                        side%ipreflag = 1                              ! restart inner_loop with preemphasis on  
                        qfirst = .false.                               ! 
                        cycle
                    else
                        qfirst = .false.
                    end if
                    call increase_scale_factor(iblock_type, iover_l, iover_s, sc_fac, qexit)         ! iso c.1.5.4.3.5   
                    if (qexit) exit                           ! when scale_factor reached maximum value defined by iso 
                end do
                call calc_subblock_gain(iblock_type, iover_s, side%isubblock_gain, qexit) 
                if (qexit) exit                           ! if no subblock_gain increased exit       ! iso c.1.5.4.3.6  
            end do 
            exit
        end do scale0
        ! retrieve best parameters and return
        distortion = distortion_min                  
        sc_fac     = sc_fac_best
        side       = side_best 
        iwk        = iwk_best
        iused_bits = ibits_best
        ! debug info
        distortion_max = max(distortion_max, distortion)
        if ( side%iscalefac_scale == 1 ) n_scale = n_scale + 1
        if ( side%ipreflag        == 1 ) n_emph  = n_emph  + 1
        if ( any(sc_fac%long     /= 0) ) n_sc_l  = n_sc_l  + 1
        if ( any(sc_fac%ishort   /= 0) ) n_sc_s  = n_sc_s  + 1
        if ( sum(side%isubblock_gain) /= 0 ) n_sub_gain = n_sub_gain + 1
        tot_sc_l = tot_sc_l + real(sc_fac%long(0:20)       , kind = kd)
        tot_sc_s = tot_sc_s + real(sc_fac%ishort(0:11, 1:3), kind = kd)
    end subroutine outer_loop
!------------------------------------------------------------------------------------------------
    subroutine select_compress(iblock_type, sc_fac, icompress)  ! iso 2.4.3.4.5
        integer            , intent(in ) :: iblock_type
        type (scale_factor), intent(in ) :: sc_fac
        integer            , intent(out) :: icompress
        icompress = 0
        select case (iblock_type) 
        case (0, 10, 11, 30, 31) ! long block
            call select_compress_long (sc_fac, icompress)
        case (20) ! short block
            call select_compress_short(sc_fac, icompress)
        case (21) ! mixed block
            call select_compress_mixed(sc_fac, icompress)
        case default
            stop ' error : subroutine select_compress '
        end select
    end subroutine select_compress
!------------------------------------------------------------------------------------------------
    subroutine select_compress_long(sc_fac, icompress)     ! iso 2.4.3.4.5, 2.4.2.7 scalefac_compress
        type (scale_factor), intent(in ) :: sc_fac
        integer            , intent(out) :: icompress
        integer, parameter :: len_scale_compress(0:15, 2) = &  ! 2.4.2.7 scalefac_compress (iso p.26)
                 reshape( (/ 0, 0, 0, 0, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, &
                             0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 3  /), (/16, 2/) )
        integer :: i, n1, n2, iscband 
        n1 = 0
        n2 = 0
        do iscband =  0, 10
            n1 = max(n1, iget_len(sc_fac%long(iscband)) )  
        end do
        do iscband = 11, 20
            n2 = max(n2, iget_len(sc_fac%long(iscband)) )  
        end do
        do i = 0, 15
            if ( n1 == len_scale_compress(i, 1) .and. n2 <= len_scale_compress(i, 2) ) then
                icompress = i
                exit
            end if
        end do
        if (n2 >= 4) stop ' error : select_compress_long ' 
    end subroutine select_compress_long
!------------------------------------------------------------------------------------------------
    subroutine select_compress_short(sc_fac, icompress)     ! iso 2.4.3.4.5, 2.4.2.7 scalefac_compress
        type (scale_factor), intent(in ) :: sc_fac
        integer            , intent(out) :: icompress
        integer :: i, iwin, n1, n2, iscband 
        integer, parameter :: len_scale_compress(0:15, 2) = &   ! 2.4.2.7 scalefac_compress (iso p.26)
                 reshape( (/ 0, 0, 0, 0, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, &
                             0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 3  /), (/16, 2/) )
        n1 = 0
        n2 = 0
        do iscband =  0, 5
            do iwin = 1, 3
                n1 = max( n1, iget_len(sc_fac%ishort(iscband, iwin)) )  
            end do
        end do
        do iscband = 6, 11
            do iwin = 1, 3
                n2 = max( n2, iget_len(sc_fac%ishort(iscband, iwin)) )  
            end do
        end do
        do i = 0, 15
            if ( n1 == len_scale_compress(i, 1) .and. n2 <= len_scale_compress(i, 2) ) then
                icompress = i
                exit
            end if
        end do
        if (n2 >= 4) stop ' error : select_compress_short ' 
    end subroutine select_compress_short
!------------------------------------------------------------------------------------------------
    subroutine select_compress_mixed(sc_fac, icompress)      ! iso 2.4.3.4.5, 2.4.2.7 scalefac_compress
        type (scale_factor), intent(in ) :: sc_fac
        integer            , intent(out) :: icompress
        integer :: i, iwin, n1, n2, iscband 
        integer, parameter :: len_scale_compress(0:15, 2) = &    ! 2.4.2.7 scalefac_compress (iso p.26)   
                 reshape( (/ 0, 0, 0, 0, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, &
                             0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 3  /), (/16, 2/) )
        n1 = 0
        n2 = 0
        do iscband =  0, 7
            n1 = max(n1, iget_len(sc_fac%long(iscband)) )  
        end do
        do iscband =  3, 5
            do iwin = 1, 3
                n1 = max(n1, iget_len(sc_fac%ishort(iscband, iwin)) )  
            end do
        end do
        do iscband = 6, 11
            do iwin = 1, 3
                n2 = max(n2, iget_len(sc_fac%ishort(iscband, iwin)) )  
            end do
        end do
        icompress = 0
        do i = 0, 15
            if ( n1 == len_scale_compress(i, 1) .and. n2 <= len_scale_compress(i, 2) ) then
                icompress = i
                exit
            end if
        end do
        if (n2 >= 4) stop ' error : select_compress_mixed ' 
    end subroutine select_compress_mixed
!------------------------------------------------------------------------------------------------
    function iget_len(k) result(ires)
        integer, intent(in) :: k
        integer :: ires
        select case (k)
        case (0) 
            ires = 0
        case (1)
            ires = 1
        case (2:3)
            ires = 2
        case (4:7)
            ires = 3
        case (8:15)
            ires = 4
        case default
            write(*, *) 'input value', k
            stop ' error : function iget_len ' 
        end select
    end function iget_len
!------------------------------------------------------------------------------------------------
    subroutine increase_scale_factor(iblock_type, iover_l, iover_s, sc_fac, qexit) ! c.1.5.4.3.6
        integer            , intent(in    ) :: iblock_type, iover_l(0:20), iover_s(0:11, 1:3)
        type (scale_factor), intent(in out) :: sc_fac
        logical            , intent(   out) :: qexit
        integer, parameter :: max_4 = 15, max_3 = 7 ! limit scale factor ! max_4 = 2**4 - 1, max_3 = 2**3 - 1
        qexit = .false.
        select case (iblock_type)
        case (0, 10, 11, 30, 31) ! long-block
            where     (sc_fac%long( 0:10)        <  max_4) sc_fac%long( 0:10) = sc_fac%long( 0:10) + iover_l( 0:10)
            where     (sc_fac%long(11:20)        <  max_3) sc_fac%long(11:20) = sc_fac%long(11:20) + iover_l(11:20)
            if (maxval(sc_fac%long( 0:10))       >= max_4) qexit  = .true.
            if (maxval(sc_fac%long(11:20))       >= max_3) qexit  = .true.
        case (20) ! short-block
            where     (sc_fac%ishort(0: 5, 1:3)  <  max_4) sc_fac%ishort(0: 5, 1:3) = sc_fac%ishort(0: 5, 1:3) + iover_s(0: 5, 1:3)
            where     (sc_fac%ishort(6:11, 1:3)  <  max_3) sc_fac%ishort(6:11, 1:3) = sc_fac%ishort(6:11, 1:3) + iover_s(6:11, 1:3)
            if (maxval(sc_fac%ishort(0: 5, 1:3)) >= max_4) qexit = .true.
            if (maxval(sc_fac%ishort(6:11, 1:3)) >= max_3) qexit = .true.
        case (21) ! mixed-block
            where     (sc_fac%long(0:7)          <  max_4) sc_fac%long(0:7) = sc_fac%long(0:7) + iover_l(0:7)
            where     (sc_fac%ishort(3: 5, 1:3)  <  max_4) sc_fac%ishort(3: 5, 1:3) = sc_fac%ishort(3: 5, 1:3) + iover_s(3: 5, 1:3)
            where     (sc_fac%ishort(6:11, 1:3)  <  max_3) sc_fac%ishort(6:11, 1:3) = sc_fac%ishort(6:11, 1:3) + iover_s(6:11, 1:3)
            if (maxval(sc_fac%long(0:7))         >= max_4) qexit = .true.
            if (maxval(sc_fac%ishort(3: 5, 1:3)) >= max_4) qexit = .true.
            if (maxval(sc_fac%ishort(6:11, 1:3)) >= max_3) qexit = .true.
        case default
            stop ' error : increase_scale_factor'
        end select 
    end subroutine increase_scale_factor
!------------------------------------------------------------------------------------------------
    subroutine calc_scfac_bit(iblock_type, nbits, icompress)  ! iso 2.4.3.4.5
        integer, intent(out) :: nbits
        integer, intent(in ) :: iblock_type, icompress
        select case (iblock_type)
        case (0, 10, 11, 30, 31)
            nbits = 11 * len_scale_compress(icompress, 1) + 10 * len_scale_compress(icompress, 2)  
        case (20)
            nbits = 18 * len_scale_compress(icompress, 1) + 18 * len_scale_compress(icompress, 2)
        case (21)
            nbits = 17 * len_scale_compress(icompress, 1) + 18 * len_scale_compress(icompress, 2)
        case default
            stop 'error : subroutine calc_scfac_bit '
        end select
    end subroutine calc_scfac_bit 
!------------------------------------------------------------------------------------------------
    subroutine reorder(iblock_type, r_mdct, wk, icut)  ! iso 2.4.3.4.8
        integer         , intent(in ) :: iblock_type, icut
        real (kind = kd), intent(in ) :: r_mdct(:, :)
        real (kind = kd), intent(out) :: wk(:) 
        integer :: k
        k = 1
        wk = 0.0_kd
        select case (iblock_type)
        case (0) ! long-block
            call reorder_long (k, 1, icut,        r_mdct, wk)
        case (10, 11, 30, 31) ! start / stop
            call reorder_long (k, 1, icut,        r_mdct, wk)
        case (20) ! short-block
            call reorder_short(k, 1, icut, 0, 11, r_mdct, wk)
        case (21) ! mixed-block
            call reorder_long (k, 1,    2,        r_mdct, wk)
            call reorder_short(k, 3, icut, 3, 11, r_mdct, wk) 
        case default
            stop ' error : subroutine reorder '
        end select 
    end subroutine reorder
!..................................................................................................
    subroutine reorder_long(k, n0, n1, r_mdct, wk)
        integer         , intent(in out) :: k
        integer         , intent(in    ) :: n0, n1
        real (kind = kd), intent(in    ) :: r_mdct(:, :)
        real (kind = kd), intent(   out) :: wk(:) 
        integer :: iband
        wk = 0.0_kd
        do iband = n0, n1
            wk(k:k + 18 - 1) = r_mdct(iband, 1:18)
            k = k + 18
        end do
    end subroutine reorder_long
!..................................................................................................
    subroutine reorder_short(k, n0, n1, iscfb0, iscfb1, r_mdct, wk)
        integer         , intent(in out) :: k
        integer         , intent(in    ) :: n0, n1, iscfb0, iscfb1
        real (kind = kd), intent(in    ) :: r_mdct(:, :)
        real (kind = kd), intent(   out) :: wk(:) 
        real (kind = kd)                 :: wk0(size(wk) / 3, 3) 
        integer :: iwin, k0, m, n, iband, iscfb
        k0 = 1
        wk0 = 0.0_kd
        do iband = n0, n1
            do iwin = 1, 3
                m = 6 * (iwin - 1) + 1
                wk0(k0:k0 + 6 - 1, iwin) = r_mdct(iband, m:m + 6 - 1)
            end do
            k0 = k0 + 6
        end do
        ! reorder short-block
        k0 = 1
        do iscfb = iscfb0, iscfb1
            n = iscalefactorband_s(iscfb, 1)
            do iwin = 1, 3
               wk(k:k + n - 1) = wk0(k0:k0 + n - 1, iwin) 
               k = k + n
            end do
            k0 = k0 + n
        end do
        n = (576 - k + 1) / 3
        do iwin = 1, 3
            wk(k:k + n - 1) = wk0(k0:k0 + n - 1, iwin)
            k = k + n
        end do
    end subroutine reorder_short
!------------------------------------------------------------------------------------------
    subroutine calc_subblock_gain(iblock_type, iover_s, isubblock_gain, qexit) ! iso 2.4.2.7 subblock_gain, 2.4.3.4.7.1
        integer, intent(in    ) :: iblock_type, iover_s(0:11, 3)
        integer, intent(in out) :: isubblock_gain(:) 
        logical, intent(   out) :: qexit
        select case (iblock_type)
        case (0, 10, 11, 30, 31)
            qexit = .true.
        case (20) ! short-block
            call sub_subblock_gain(iover_s, isubblock_gain, qexit)
        case (21) ! mixed-block
            call sub_subblock_gain(iover_s, isubblock_gain, qexit) 
        case default
            stop ' error : subroutine reorder '
        end select 
    end subroutine calc_subblock_gain
!------------------------------------------------------------------------------------------
    subroutine sub_subblock_gain(iover_s, isubblock_gain, qexit)
        integer, intent(in    ) :: iover_s(0:11, 3)
        integer, intent(in out) :: isubblock_gain(:)
        logical, intent(   out) :: qexit
        integer :: iwin
        qexit = .true.
        do iwin = 1, 3
            if ( sum(iover_s(:, iwin)) /= 0 .and. isubblock_gain(iwin) < 7) then
                isubblock_gain(iwin) = isubblock_gain(iwin) + 1
                qexit = .false.
            end if
        end do
    end subroutine sub_subblock_gain
!------------------------------------------------------------------------------------------
    subroutine calc_scale(iblock_type, iscalefac_scale, sc_fac, ipretab, isubblock_gain, scale) ! iso 2.4.3.4.7.1
        integer             , intent(in ) :: iblock_type, iscalefac_scale, ipretab(0:20), isubblock_gain(:)
        type (scale_factor) , intent(in ) :: sc_fac
        real (kind = kd)    , intent(out) :: scale(:)
        integer :: k
        scale = 1.0_kd
        select case(iblock_type)
        case (0, 10, 11, 30, 31)
            k = 1
            call calc_scale_long (k, 0, 20, iscalefac_scale, sc_fac, ipretab       , scale) 
        case (20)
            k = 1
            call calc_scale_short(k, 0, 11, iscalefac_scale, sc_fac, isubblock_gain, scale) 
        case (21)
            k = 1
            call calc_scale_long (k, 0,  7, iscalefac_scale, sc_fac, ipretab,        scale) 
            call calc_scale_short(k, 3, 11, iscalefac_scale, sc_fac, isubblock_gain, scale) 
        case default
            stop ' error : subroutine rescale '
        end select
    end subroutine calc_scale
!..................................................................................................
    subroutine calc_scale_long(k, n0, n1, iscalefac_scale, sc_fac, ipretab, wk)
        integer             , intent(in out) :: k 
        integer             , intent(in    ) :: n0, n1, iscalefac_scale, ipretab(0:20)
        type (scale_factor) , intent(in    ) :: sc_fac
        real (kind = kd)    , intent(   out) :: wk(:)
        integer :: iscband, i
        do iscband = n0, n1
            do i = 1, iscalefactorband_l(iscband, 1)
                wk(k) = sqrt( real( 2**( (1 + iscalefac_scale) * &
                                       (sc_fac%long(iscband) + ipretab(iscband)) ), kind = kd )  )
                k = k + 1
            end do
        end do
    end subroutine calc_scale_long
!..................................................................................................
    subroutine calc_scale_short(k, n0, n1, iscalefac_scale, sc_fac, isubblock_gain, wk)
        integer             , intent(in out) :: k 
        integer             , intent(in    ) :: n0, n1, iscalefac_scale, isubblock_gain(:)
        type (scale_factor) , intent(in    ) :: sc_fac
        real (kind = kd)    , intent(   out) :: wk(:)
        integer :: iscband, i, iwin
        do iscband = n0, n1
            do iwin = 1, 3
                do i = 1, iscalefactorband_s(iscband, 1)
                    wk(k) = sqrt( real(2**((1 + iscalefac_scale) * sc_fac%ishort(iscband, iwin)), kind = kd ) ) & 
                          * real(2**(2 * isubblock_gain(iwin)), kind = kd) 
                    k = k + 1
                end do
            end do
        end do
    end subroutine calc_scale_short
!------------------------------------------------------------------------------------------
    subroutine calc_distortion(iblock_type, iglobal_gain, x, ix, th, sc, iover_l, iover_s, distortion) ! iso c.1.5.4.3.3
        real (kind = kd), intent(in ) :: x(:), th(:), sc(:)
        integer         , intent(in ) :: iblock_type, iglobal_gain, ix(:)
        integer         , intent(out) :: iover_l(0:), iover_s(0:, :) 
        real (kind = kd), intent(out) :: distortion 
        iover_l = 0
        iover_s = 0
        select case (iblock_type)
        case (0, 10, 11, 30, 31)
            call calc_dist_long (iglobal_gain, x, ix, th, sc, iover_l, distortion)  
        case (20)
            call calc_dist_short(iglobal_gain, x, ix, th, sc, iover_s, distortion)  
        case (21)
            call calc_dist_mixed(iglobal_gain, x, ix, th, sc, iover_l, iover_s, distortion)  
       case default
            stop ' error : calc_distortion '
        end select
    end subroutine calc_distortion
!..................................................................................................
    subroutine calc_dist_long(iglobal_gain, x, ix, th, sc, iover_l, distortion)
        real (kind = kd), intent(in ) :: x(:), th(:), sc(:)
        integer         , intent(in ) :: iglobal_gain, ix(:)
        integer         , intent(out) :: iover_l(0:20)
        real (kind = kd), intent(out) :: distortion
        real (kind = kd) ::  dx, ds2(0:20), as2(0:20), bw, a, d
        integer :: i, istart, iend, iscband, ipos
        iover_l = 0
        ds2 = 0.0_kd   
        as2 = 0.0_kd
        distortion = 0.0_kd !
        do iscband = 0, 20
            bw     = real(iscalefactorband_l(iscband, 1), kind = kd)  
            istart =      iscalefactorband_l(iscband, 2) + 1
            iend   =      iscalefactorband_l(iscband, 3) + 1
            do i = istart, iend
                dx = abs(x(i)) - real(abs(ix(i)), kind = kd)**(4.0_kd / 3.0_kd) &
                               * 2.0_kd**(real(iglobal_gain - 210, kind = kd) / 4.0_kd) / sc(i)
                if (abs(dx) > 0.0_kd) ds2(iscband) = ds2(iscband) + log(abs(dx)) 
                as2(iscband) = as2(iscband) + log(th(i))
            end do
            distortion = distortion + (ds2(iscband) - as2(iscband)) / bw ! geometric average (x1...xn)^1/n
        end do
        ipos = maxloc(ds2 - as2, mask = ds2 > as2, dim = 1)
        if (ipos /= 0) iover_l(ipos - 1) = 1
        ! band 21
        bw     = real(iend - istart + 1, kind = kd)  
        istart = iscalefactorband_l(20, 3) + 2
        iend   = 576
        d = 0.0_kd
        a = 0.0_kd
        do i = istart, iend
            dx = abs(x(i)) - real(abs(ix(i)), kind = kd)**(4.0_kd / 3.0_kd) &
                           * 2.0_kd**(real(iglobal_gain - 210, kind = kd) / 4.0_kd) / sc(i)  
            if (abs(dx) > 0.0_kd) d = d + log(abs(dx)) 
            a = a + log(th(i))
        end do
        distortion = distortion + (d - a) / bw 
    end subroutine calc_dist_long
!------------------------------------------------------------------------------------------------
    subroutine calc_dist_short(iglobal_gain, x, ix, th, sc, iover_s, distortion)
        real (kind = kd), intent(in ) :: x(:), th(:), sc(:)
        integer         , intent(in ) :: iglobal_gain, ix(:)
        integer         , intent(out) :: iover_s(0:11, 1:3)
        real (kind = kd), intent(out) :: distortion
        real (kind = kd) :: dx, ds2(0:11, 1:3), as2(0:11, 1:3), bw, a, d
        integer :: i, istart, iend, k, iscband, iwin, ipos(2)
        k = 0
        iover_s = 0
        ds2 = 0.0_kd
        as2 = 0.0_kd
        distortion = 0.0_kd
        do iscband = 0, 11
            bw = real(iscalefactorband_s(iscband, 1), kind = kd)
            istart =  iscalefactorband_s(iscband, 2) + 1
            iend   =  iscalefactorband_s(iscband, 3) + 1
            do iwin = 1, 3
                do i = istart, iend
                    k = k + 1
                    dx = abs(x(k)) - real(abs(ix(k)), kind = kd)**(4.0_kd / 3.0_kd) &
                                   * 2.0_kd**( real(iglobal_gain - 210, kind = kd) / 4.0_kd) / sc(k) 
                    if (abs(dx) > th(k)) ds2(iscband, iwin) = ds2(iscband, iwin) + log(abs(dx)) 
                    as2(iscband, iwin) = as2(iscband, iwin) + log(th(k))
                end do
                distortion = distortion + (ds2(iscband, iwin) - as2(iscband, iwin)) / bw 
            end do
        end do
        ipos = maxloc(ds2 - as2, mask = ds2 > as2) 
        if (any(ipos /= 0)) iover_s(ipos(1) - 1, ipos(2)) = 1
        ! band 12
        bw     = real(iend - istart + 1, kind = kd) 
        istart = iscalefactorband_s(11, 3) + 2
        iend   = 576
        d = 0.0_kd 
        a = 0.0_kd
        do i = istart, iend
            dx = abs(x(i)) - real(abs(ix(i)), kind = kd)**(4.0_kd / 3.0_kd) &
                           * 2.0_kd**(real(iglobal_gain - 210, kind = kd) / 4.0_kd) / sc(i)
            if (abs(dx) > 0.0_kd) d = d + log(abs(dx)) 
            a = a + log(th(i))
        end do
        distortion = distortion + (d - a) / bw 
    end subroutine calc_dist_short
!------------------------------------------------------------------------------------------
    subroutine calc_dist_mixed(iglobal_gain, x, ix, th, sc, iover_l, iover_s, distortion)
        real (kind = kd), intent(in ) :: x(:), th(:), sc(:)
        integer         , intent(in ) :: iglobal_gain, ix(:)
        integer         , intent(out) :: iover_l(0:20), iover_s(0:11, 1:3)
        real (kind = kd), intent(out) :: distortion
        real (kind = kd) :: dx, ds2_l(0:7), as2_l(0:7), ds2_s(0:11, 1:3), as2_s(0:11, 1:3), bw, a, d
        integer :: i, istart, iend, k, iscband, iwin, ipos(2)
        k = 0
        iover_l = 0
        iover_s = 0
        ds2_l = 0.0_kd
        as2_l = 0.0_kd
        ds2_s = 0.0_kd
        as2_s = 0.0_kd
        distortion = 0.0_kd
        do iscband = 0, 7
            bw = real(iscalefactorband_l(iscband, 1), kind = kd)
            istart =  iscalefactorband_l(iscband, 2) + 1
            iend   =  iscalefactorband_l(iscband, 3) + 1
            do i = istart, iend
                k = k + 1
                dx = abs(x(i)) - real(abs(ix(i)), kind = kd)**(4.0_kd / 3.0_kd) &
                               * 2.0_kd**( real(iglobal_gain - 210, kind = kd) / 4.0_kd) / sc(i)
                if (abs(dx) > 0.0_kd) ds2_l(iscband) = ds2_l(iscband) + log(abs(dx))
                as2_l(iscband) = as2_l(iscband) + log(th(i))
            end do
            distortion = distortion + (ds2_l(iscband) - as2_l(iscband)) / bw 
        end do
        do iscband = 3, 11
            do iwin = 1, 3 
                bw = real(iscalefactorband_s(iscband, 1), kind = kd)
                istart =  iscalefactorband_s(iscband, 2) + 1
                iend   =  iscalefactorband_s(iscband, 3) + 1  
                do i = istart, iend
                    k = k + 1
                    dx = abs(x(k)) - real(abs(ix(k)), kind = kd)**(4.0_kd / 3.0_kd) & 
                                   * 2.0_kd**( real(iglobal_gain - 210, kind = kd) / 4.0_kd) / sc(k)
                    if (abs(dx) > 0.0_kd) ds2_s(iscband, iwin) = ds2_s(iscband, iwin) + log(abs(dx))
                    as2_s(iscband, iwin) = as2_s(iscband, iwin) + log(th(k))
                end do
                distortion = distortion + (ds2_s(iscband, iwin) - as2_s(iscband, iwin)) / bw 
            end do
        end do
        if ( maxval(ds2_l - as2_l, mask = ds2_l > as2_l) > maxval(ds2_s - as2_s, mask = ds2_s > as2_s) ) then
            ipos(1) = maxloc(ds2_l - as2_l, mask = ds2_l > as2_l, dim = 1)
            if (ipos(1) /= 0) iover_l(ipos(1) - 1) = 1
        else
            ipos = maxloc(ds2_s - as2_s, mask = ds2_s > as2_s) 
            if (any(ipos /= 0)) iover_s(ipos(1) - 1, ipos(2)) = 1 
        end if
        ! band 12
        bw     = real(iend - istart + 1, kind = kd)  
        istart = iscalefactorband_s(11, 3) + 2
        iend   = 576
        d = 0.0_kd 
        a = 0.0_kd
        do i = istart, iend
            dx = abs(x(i)) - real(abs(ix(i)), kind = kd)**(4.0_kd / 3.0_kd) &
                           * 2.0_kd**( real(iglobal_gain - 210, kind = kd) / 4.0_kd) / sc(i)
            if (abs(dx) > 0.0_kd) d = d + log(abs(dx)) 
            a = a + log(th(i))
        end do
        distortion = distortion + (d - a) / bw 
    end subroutine calc_dist_mixed
!------------------------------------------------------------------------------------------------
end module mod_layer3