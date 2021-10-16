module mod_mdct
    implicit none
    private
    public :: mdct_initialize, sub_mdct
    real (kind = 8):: pi
    integer :: indx9(9)
    complex (kind = 8), save:: omega9(0:8), omega9s(9), sqrt3_2, omega3(0:2), omega3s(3)
    real(kind = 8), parameter::c(0:7) = [-0.6000d0, -0.5350d0, -0.3300d0, -0.1850d0, -0.0950d0, -0.0410d0, -0.0142d0, -0.0037d0]
    real(kind = 8), save :: ca(0:7), cs(0:7), window_n(36), window_s(12), window_start(36), window_stop(36)
    real(kind = 8), save, allocatable :: subbuff(:, :, :)
contains
!---------------------------------------------------------------------------------------------
    subroutine mdct_initialize()
        integer :: i
        pi = 4.0d0 * atan(1.0d0)
        subbuff = 0.0d0
        call fft9_initialize()
        call fft3_initialize()
        ! iso table b.9 coefficients for aliasing reduction:
        do i = 0, 7
            cs(i) =  sqrt( 1.0d0    / ( 1.0d0 + c(i)**2 ) )
            ca(i) = -sqrt( c(i)**2  / ( 1.0d0 + c(i)**2 ) )
        end do
        ! iso 2.4.3.4.10.3 windowing, c.1.5.3.3
        ! normal window
        do i = 1, 36 
           window_n(i) = sin( real(2 * i - 1, kind = 8) / 72.0d0 * pi)
        end do
        ! short window
        do i = 1, 12 
            window_s(i) = sin( real(2 * i - 1, kind = 8) / 24.0d0 * pi)
        end do
        ! start / stop window
        do i = 1, 18 
            window_start(i     ) = sin( real(2 *  i       - 1, kind = 8) / 72.0d0 * pi)
            window_stop (i + 18) = sin( real(2 * (i + 18) - 1, kind = 8) / 72.0d0 * pi)
        end do
        do i = 1, 6
            window_start(i + 24) = sin( real(2 * (i + 6) - 1, kind = 8) / 24.0d0 * pi)
            window_stop (i +  6) = sin( real(2 *  i      - 1, kind = 8) / 24.0d0 * pi)
        end do
        window_start(19:24) = 1.0d0
        window_start(31:36) = 0.0d0
        window_stop ( 1: 6) = 0.0d0
        window_stop (13:18) = 1.0d0
    end subroutine mdct_initialize
!-------------------------------------------------------------------------------------------------
    subroutine cp_subband(subband, subbuff)
        real (kind = 8), intent(in ) :: subband(:, :, :) ! 32 * 36 * nchannel
        real (kind = 8), intent(out) :: subbuff(:, :, :) ! 32 * 54 * nchannel
        subbuff = eoshift(subbuff, 36, 0.0d0, 2)
        subbuff(:     , 19:54  , :) =  subband(:, :, :)           
        subbuff(2:32:2, 20:54:2, :) = -subbuff(2:32:2, 20:54:2, :) ! iso figure a.4 layer iii decoder diagram
                                                 !    ~~~~~~~~~~~~ not written in the text but only in figure 
    end subroutine cp_subband
!--------------------------------------------------------------------------------------------------
    subroutine sub_mdct(subband, r_mdct, mblock_type, q_alias)
        real (kind = 8), intent(in ) :: subband(:, :, :)
        real (kind = 8), intent(out) :: r_mdct (:, :, :, :)
        integer        , intent(in ) :: mblock_type(:)
        logical        , intent(in ) :: q_alias 
        integer :: igranule, ichannel
        logical, save :: qfirst = .true.
        if (qfirst) then
            qfirst = .false.
            allocate( subbuff(32, 54, size(subband, 3)) )
            subbuff = 0.0d0
            call mdct_initialize()
        end if
        call cp_subband(subband, subbuff)
        do ichannel = 1, size(subband, 3)
            do igranule = 1, 2
                select case (mblock_type(igranule))
                case ( 0) ! long block 
                    call sub_mdct_normal(igranule, 1, 32, subbuff(:, :, ichannel), r_mdct(:, :, igranule, ichannel))
                    call anti_alias(r_mdct(1:32, :, igranule, ichannel))
                case (10) ! start block : long -> short
                    call sub_mdct_start (igranule, 1, 32, subbuff(:, :, ichannel), r_mdct(:, :, igranule, ichannel))
                    call anti_alias(r_mdct(1:32, :, igranule, ichannel))
                case (11) ! start block : long -> mixed
                    call sub_mdct_normal(igranule, 1,  2, subbuff(:, :, ichannel), r_mdct(:, :, igranule, ichannel))
                    call sub_mdct_start (igranule, 3, 32, subbuff(:, :, ichannel), r_mdct(:, :, igranule, ichannel))
                    call anti_alias(r_mdct(1:32, :, igranule, ichannel))
                case (30) ! stop block : short -> long
                    call sub_mdct_stop  (igranule, 1, 32, subbuff(:, :, ichannel), r_mdct(:, :, igranule, ichannel))
                    call anti_alias(r_mdct(1:32, :, igranule, ichannel))
                case (31) ! stop block : mixed -> long
                    call sub_mdct_normal(igranule, 1,  2, subbuff(:, :, ichannel), r_mdct(:, :, igranule, ichannel))
                    call sub_mdct_stop  (igranule, 3, 32, subbuff(:, :, ichannel), r_mdct(:, :, igranule, ichannel))
                    call anti_alias(r_mdct(1:32, :, igranule, ichannel))
                case (20) ! short block : short mode
                    call sub_mdct_short(igranule, 1, 32, subbuff(:, :, ichannel), r_mdct(:, :, igranule, ichannel))
                case (21) ! short block : mixed mode
                    call sub_mdct_normal(igranule, 1,  2, subbuff(:, :, ichannel), r_mdct(:, :, igranule, ichannel))
                call sub_mdct_short (igranule, 3, 32, subbuff(:, :, ichannel), r_mdct(:, :, igranule, ichannel))
                    if (q_alias) call anti_alias(r_mdct(1:2, :, igranule, ichannel))
                case default
                    stop 'error : subroutine sub_mdct : block_type error'
                end select 
            end do
        end do
    end subroutine sub_mdct
!--------------------------------------------------------------------------------------------------
    subroutine sub_mdct_normal(igranule, n0, n1, subbuff, r_mdct)
        integer        , intent(in ) :: igranule, n0, n1 
        real (kind = 8), intent(in ) :: subbuff(:, :)
        real (kind = 8), intent(out) :: r_mdct (:, :)
        real (kind = 8)              :: wk(36)
        integer :: iband, k
        k = 18 * (igranule - 1) + 1
        do iband = n0, n1 
            wk = window_n * subbuff(iband, k:k + 35)
            call mdct36(wk, r_mdct(iband, :))
        end do
    end subroutine sub_mdct_normal
!--------------------------------------------------------------------------------------------------
    subroutine sub_mdct_start(igranule, n0, n1, subbuff, r_mdct)
        integer        , intent(in ) :: igranule, n0, n1 
        real (kind = 8), intent(in ) :: subbuff(:, :)
        real (kind = 8), intent(out) :: r_mdct (:, :)
        real (kind = 8)              :: wk(36)
        integer :: iband, k
        k = 18 * (igranule - 1) + 1
        do iband = n0, n1 
            wk = window_start * subbuff(iband, k:k + 35)
            call mdct36(wk, r_mdct(iband, :))
        end do
    end subroutine sub_mdct_start
!--------------------------------------------------------------------------------------------------
    subroutine sub_mdct_stop(igranule, n0, n1, subbuff, r_mdct)
        integer        , intent(in ) :: igranule, n0, n1 
        real (kind = 8), intent(in ) :: subbuff(:, :)
        real (kind = 8), intent(out) :: r_mdct (:, :)
        real (kind = 8)              :: wk(36)
        integer :: iband, k
        k = 18 * (igranule - 1) + 1
        do iband = n0, n1 
            wk = window_stop * subbuff(iband, k:k + 35)
            call mdct36(wk, r_mdct(iband, :))
        end do
    end subroutine sub_mdct_stop
!--------------------------------------------------------------------------------------------------
    subroutine sub_mdct_short(igranule, n0, n1, subbuff, r_mdct)
        integer        , intent(in ) :: igranule, n0, n1
        real (kind = 8), intent(in ) :: subbuff(:, :)
        real (kind = 8), intent(out) :: r_mdct (:, :)
        real (kind = 8)              :: wk(12)
        integer :: iband, k, m, m0, n 
        m0 = 18 * (igranule - 1) + 1 + 6
        do k = 1, 3
            m = m0 + 6 * (k - 1) 
            n =  1 + 6 * (k - 1)
            do iband = n0, n1 
                wk = window_s * subbuff(iband, m:m + 11)
                call mdct12(wk, r_mdct(iband, n:n + 5))
            end do
        end do
    end subroutine sub_mdct_short
!--------------------------------------------------------------------------------------------------
    subroutine anti_alias(x) ! c.1.5.3.3 aliasing-butterfly, 2.4.3.4.10.1 alias reduction, table b.9  
        real (kind = 8), intent(in out) :: x(:, :) ! x(32, 18)
        real (kind = 8) :: tmp1, tmp2
        integer :: iband, k
        do iband = 1, size(x, 1) - 1
            do k = 0, 7
                tmp1 = x(iband    , 18 - k) 
                tmp2 = x(iband + 1,  1 + k) 
                x(iband    , 18 - k) =  tmp1 * cs(k) + tmp2 * ca(k)
                x(iband + 1,  1 + k) =  tmp2 * cs(k) - tmp1 * ca(k) 
            end do
        end do
    end subroutine anti_alias
!--------------------------------------------------------------------------------------------------
    subroutine mdct36(x_in, x_out) ! iso c.1.5.3.3 mdct:   ! mdct(4n) <=> fft(n) proof: http://members.at.infoseek.co.jp/kitaurawa/mdct.pdf   (jp)
        real (kind = 8), intent(in ) :: x_in (:) ! 36  
        real (kind = 8), intent(out) :: x_out(:) ! 18
        real (kind = 8):: x_shift(36), x_re, x_im
        complex (kind = 8):: fft(9)
        integer:: i
        x_shift( 1: 9) = -x_in(28:36)
        x_shift(10:36) =  x_in( 1:27)
        do i = 1, 9
            x_re = x_shift(2 * i -  1) - x_shift(38 - 2 * i) 
            x_im = x_shift(2 * i + 17) - x_shift(20 - 2 * i) 
            fft(i) = cmplx(x_re, x_im, kind = 8) * omega9s(i)
        end do
        call fft9(fft)
        fft = fft * omega9s
        do i = 1, 9
            x_out(2 * i - 1) =  real(      fft(i)      , kind = 8) 
            x_out(2 * i    ) =  real(aimag(fft(10 - i)), kind = 8)  
        end do
    end subroutine mdct36
!---------------------------------------------------------------------------------------------------
    subroutine fft9(fft) ! fft of 3^n (n=2)
        complex (kind = 8), intent(in out):: fft(:)
        complex (kind = 8) :: c1, c2, c3, c4, tmp1, tmp2, tmp3
        integer, parameter :: nn = 9
        integer :: i, j, iphase1, iphase2, m1, m2, m3, k3, kn3
        fft = fft(indx9) / 9.0d0 ! reorder and normalize
        do k3 = 1, 2 ! 3^n (n=2)
            kn3 = 3**(k3 - 1)
            do i = 1, nn, 3 * kn3 
                do j = 1, kn3
                    iphase1 = 3**(2 - k3) * (j - 1) 
                    iphase2 = 2 * iphase1
                    c1 = omega9( mod(iphase1, nn) )
                    c2 = omega9( mod(iphase2, nn) )
                    m1 = i + j - 1
                    m2 = m1 + kn3
                    m3 = m2 + kn3
                    tmp1 =      fft(m1)
                    tmp2 = c1 * fft(m2)
                    tmp3 = c2 * fft(m3)
                    fft(m1) = tmp1 + tmp2 + tmp3
                    c3 = tmp1 - 0.5d0 * ( tmp2 + tmp3 )
                    c4 =     sqrt3_2  * ( tmp2 - tmp3 ) ! sqrt3_2 = i sqrt(3) / 2
                    fft(m2) = c3 + c4
                    fft(m3) = c3 - c4
                end do
            end do
        end do
    end subroutine fft9
!-----------------------------------------------------------------------
    subroutine fft9_initialize()
        integer :: i, j, k, n
        pi = 4.0d0 * atan(1.0d0)
        sqrt3_2 = cmplx(0.0d0, sqrt(0.75d0), kind = 8) ! isqrt(3) / 2
        do i = 1, 9
            n = 0
            k = i - 1
            do j = 2 - 1, 0, -1
                n = n + mod(k, 3) * 3**j
                k = k / 3
            end do
            indx9(i) = n + 1
            !omega9(i - 1) = exp( cmplx(0.0d0, 2.0d0 * pi /  9.0d0 * real(i - 1, kind = 8), kind = 8) )
            !omega9s(i)    = exp( cmplx(0.0d0, 2.0d0 * pi / 36.0d0 * ( 1.0d0 / 8.0d0 + real(i - 1, kind = 8) ), kind = 8) )
            omega9(i - 1) = cmplx( cos( real( 2 * i - 2, kind = 8) /   9.0d0 * pi ),   &
                                   sin( real( 2 * i - 2, kind = 8) /   9.0d0 * pi ), kind = 8 ) 
            omega9s(i)    = cmplx( cos( real( 8 * i - 7, kind = 8) / 144.0d0 * pi ), & 
                                   sin( real( 8 * i - 7, kind = 8) / 144.0d0 * pi ), kind = 8 ) 
        end do
    end subroutine fft9_initialize
!-----------------------------------------------------------------------
    subroutine mdct12(x_in, x_out)
        real (kind = 8), intent(in ) :: x_in (:) ! 12
        real (kind = 8), intent(out) :: x_out(:) !  6
        real (kind = 8):: x_shift(12), x_re, x_im  
        complex (kind = 8):: fft(3)
        integer:: i  
        x_shift( 1: 3) = -x_in(10:12)
        x_shift( 4:12) =  x_in( 1: 9)
        do i = 1, 3
            x_re = x_shift(2 * i - 1) - x_shift(14 - 2 * i) 
            x_im = x_shift(2 * i + 5) - x_shift( 8 - 2 * i) 
            fft(i) = cmplx(x_re, x_im, kind = 8) * omega3s(i)
        end do
        call fft3(fft)
        fft = fft * omega3s
        do i = 1, 3
            x_out(2 * i - 1) =  real(      fft(i)     , kind = 8) 
            x_out(2 * i    ) =  real(aimag(fft(4 - i)), kind = 8)  
        end do
    end subroutine mdct12
!-----------------------------------------------------------------------
    subroutine fft3(fft)
        complex (kind = 8), intent(in out):: fft(:)
        complex (kind = 8) :: c2, c3, c4
        c2 = fft(2) + fft(3)
        c3 = fft(1) - 0.5d0 * c2
        c4 = sqrt3_2 * ( fft(2) - fft(3) ) ! sqrt3_2 = i sqrt(3) / 2
        fft(1) = fft(1) + c2
        fft(2) = c3 + c4
        fft(3) = c3 - c4
        fft = fft / 3.0d0 ! normalization
    end subroutine fft3
!--------------------------------------------------------------------------------------------------
    subroutine fft3_initialize()
        integer :: i
        pi = 4.0d0 * atan(1.0d0)
        sqrt3_2 = cmplx(0.0d0, sqrt(0.75d0), kind = 8) ! isqrt(3) / 2
        do i = 1, 3
            !omega3(i - 1) = exp( cmplx(0.0d0, 2.0d0 * pi /  3.0d0 * real(i - 1, kind = 8), kind = 8) )
            !omega3s(i)    = exp( cmplx(0.0d0, 2.0d0 * pi / 12.0d0 * ( 1.0d0 / 8.0d0 + real(i - 1, kind = 8) ), kind = 8) )
            omega3(i - 1) = cmplx( cos( real( 2 * i - 2, kind = 8) /  3.0d0 * pi ),   &
                                   sin( real( 2 * i - 2, kind = 8) /  3.0d0 * pi ), kind = 8 ) 
            omega3s(i)    = cmplx( cos( real( 8 * i - 7, kind = 8) / 48.0d0 * pi ),   & 
                                   sin( real( 8 * i - 7, kind = 8) / 48.0d0 * pi ), kind = 8 ) 
        end do
    end subroutine fft3_initialize
!-----------------------------------------------------------------------
end module mod_mdct