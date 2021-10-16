module mod_inner_loop
    use mod_mpg
    use mod_huffman
    implicit none
    private
    public :: inner_loop ! subroutine
contains
!-----------------------------------------------------------------------------------------------
    subroutine inner_loop(ibit, iblock_type, wk, i_mdct, side, itot_bits) ! iso c.1.5.4.4
        integer             , intent(in    ) :: ibit, iblock_type
        real (kind = 8)     , intent(in    ) :: wk(:)
        type (side_info_sub), intent(in out) :: side
        integer             , intent(   out) :: i_mdct(:), itot_bits
        integer :: iq0, iqquant, iquantanf, ibigvalues, icount1, & 
                   n0, n1, itab, iregion0, iregion1, ndiv
        integer, parameter :: ndiv0 = 32, ndiv1 = 32
        logical :: qfirst 
        call calc_quantanf(wk, iquantanf) ! not wk 
        ndiv    = ndiv0
        iqquant = 0
        iq0     = 0
        qfirst  = .true. 
        do while (ndiv /= 0) ! obtain starting quantization step
            call quantization(iqquant, iquantanf, wk, i_mdct)
            if ( qfirst .and. iquantanf > -210 .and. maxval(i_mdct) <= 2**13 - 1) then ! 15 + 2**13 - 1
                iquantanf = iquantanf - 63 ! in iso document 2.4.2.7 "big_values", the maximum absolute value is constrained to 8191 = 2^13 - 1.
                cycle                      ! however it is possible to use up to 2^13 - 1 + 15 = 8206.   
            else                        ! iso sample program dist10 uses 8206.
                qfirst = .false.           ! because iso documet section 2 is normative, 2^13 - 1 is chosen in uzura. 
            end if                              
            if (maxval(i_mdct) > 2**13 - 1) then  ! bi-section method
                iq0 = iqquant                           
                iquantanf = iquantanf + ndiv 
            else
             !  ndiv = ndiv / 2 ! bi-section
                ndiv = int(ndiv * 0.618) ! fibonacchi
                iqquant = iq0 + ndiv
            end if 
        end do
        ndiv    = ndiv1
        iqquant = max(0, -iquantanf - 210)
        iq0     = iqquant
        do while (ndiv > 0)
            call quantization(iqquant, iquantanf, wk, i_mdct)
            itot_bits = 0
            side%iglobal_gain = iqquant + iquantanf + 210 
            call divide(i_mdct, ibigvalues, icount1)
            side%ibig_values = ibigvalues 
            side%icount1     = icount1
            call count_one(i_mdct, ibigvalues, icount1, itab, itot_bits)
            side%icount1table_select = itab
            select case (iblock_type) ! iso 2.4.2.7 window_switching_flag
            case (0)
                call sub_divide(ibigvalues, iregion0, iregion1)  ! iso 
                iregion0 = min(15, iregion0)                     ! iso 2.4.2.7 region0_count 
                iregion1 = min( 7, iregion1, 19 - iregion0) 
                side%iregion0_count = iregion0     !  4 bits
                side%iregion1_count = iregion1     !  3 bits
                n0 = iscalefactorband_l(iregion0, 3) + 1 
                n1 = iscalefactorband_l(iregion0 + iregion1 + 1, 3) + 1
                call select_table2(i_mdct,      1, n0, itab, itot_bits)
                side%itable_select(1) = itab
                call select_table2(i_mdct, n0 + 1, n1, itab, itot_bits)
                side%itable_select(2) = itab
                call select_table2(i_mdct, n1 + 1, ibigvalues * 2, itab, itot_bits)
                side%itable_select(3) = itab
            case (10, 11, 30, 31) ! switching block
                iregion0 = 7                
                iregion1 = 36 
                side%iregion0_count = iregion0      !  4 bits
                side%iregion1_count = iregion1      !  3 bits
                n0 = iscalefactorband_l(iregion0, 3) + 1  
                n1 = ibigvalues * 2
                call select_table2(i_mdct,      1, n0, itab, itot_bits)
                side%itable_select(1) = itab
                call select_table2(i_mdct, n0 + 1, n1, itab, itot_bits)
                side%itable_select(2) = itab
            case (20) ! short 
                iregion0 = 8                
                iregion1 = 36 
                side%iregion0_count = iregion0        
                side%iregion1_count = iregion1        
                n0 = iscalefactorband_s( (iregion0 + 1) / 3, 2) * 3   
                n1 = ibigvalues * 2  
                call select_table2(i_mdct,      1, n0, itab, itot_bits)
                side%itable_select(1) = itab
                call select_table2(i_mdct, n0 + 1, n1, itab, itot_bits)
                side%itable_select(2) = itab
            case (21) ! mixed
                iregion0 = 7                
                iregion1 = 36 
                side%iregion0_count = iregion0        
                side%iregion1_count = iregion1        
                n0 = iscalefactorband_l(iregion0, 3) + 1
                n1 = ibigvalues * 2  
                call select_table2(i_mdct,      1, n0, itab, itot_bits)
                side%itable_select(1) = itab
                call select_table2(i_mdct, n0 + 1, n1, itab, itot_bits)
                side%itable_select(2) = itab
            case default
                stop ' error : subroutine inner_loop ' 
            end select
            if (itot_bits > ibit) then ! bi-section method: speeds up c.1.5.4.4.2
                iq0 = iqquant
                iqquant = iqquant + ndiv 
            else
             !  ndiv = ndiv / 2 ! bi-section
               ndiv = int(ndiv * 0.618) ! fibonacchi
               iqquant = iq0 + ndiv 
            end if 
        end do
        where (wk < 0.0d0) i_mdct = -i_mdct
    end subroutine inner_loop
!---------------------------------------------------------------------------------------------
    subroutine calc_quantanf(wk, iquantanf)     ! iso c.1.5.4.2.1
        real (kind = 8), intent(in ) :: wk(:)
        integer, intent(out)         :: iquantanf
        integer :: i
        real (kind = 8) :: sfm, sum1, sum2, tmp
        sum1 =  0.0d0
        sum2 =  0.01d0 
        do i = 1, size(wk) ! 576
            tmp = wk(i)**2.0d0
            if (tmp > 0.0d0) sum1 = sum1 + log(tmp) 
            sum2 = sum2 + tmp
        end do
        sfm = exp( sum1 / 576.0d0) * 576.0d0 / sum2 
        iquantanf = int( 8.0d0 * log(sfm) ) 
    end subroutine calc_quantanf
!------------------------------------------------------------------------------------------------
    subroutine quantization(iqquant, iquantanf, r_mdct, i_mdct) ! iso c.1.5.4.4.1
        integer        , intent(in ) :: iqquant, iquantanf
        real (kind = 8), intent(in ) :: r_mdct(:)
        integer        , intent(out) :: i_mdct(:)
        real (kind = 8) :: denom, tmp(576)
        denom = 2.0d0 ** ( -real(iqquant + iquantanf, kind = 8) / 4.0d0 )
        tmp    = abs(r_mdct) * denom
        i_mdct = nint( sqrt(tmp * sqrt(tmp)) - 0.0946d0 ) ! nint(tmp**(3/4) - 0.0946)
    end subroutine quantization
!----------------------------------------------------------------------------------------------
    subroutine divide(i_mdct, ibigvalues, icount1) ! iso c.1.5.4.4.3, c.1.5.4.4.4
        integer, intent(in ) :: i_mdct(:)
        integer, intent(out) :: ibigvalues, icount1
        integer :: i, ibig, izero
        do i = 576, 110, -2
            izero = i
            if (i_mdct(i) /= 0 .or. i_mdct(i - 1) /= 0) exit 
        end do
        do i = izero, 110, -4
            ibig = i
            if ( abs(i_mdct(i    )) > 1 .or. abs(i_mdct(i - 1)) > 1 .or. &        !  0, +1, - 1
                 abs(i_mdct(i - 2)) > 1 .or. abs(i_mdct(i - 3)) > 1        ) exit  
        end do
        ibigvalues = ibig / 2 
        icount1 = (izero - ibig) / 4 
    end subroutine divide 
!----------------------------------------------------------------------------------------------
    subroutine count_one(i_mdct, ibigvalues, icount1, itab, isum) ! iso 2.4.2.7, c.1.5.4.4.5
        integer, intent(in    ) :: i_mdct(:), ibigvalues, icount1
        integer, intent(   out) :: itab
        integer, intent(in out) :: isum
        integer :: isum0, isum1, i, k1, k2, k3, k4
        isum0 = 0
        isum1 = 0
        do i = ibigvalues * 2 + 1, ibigvalues * 2 + icount1 * 4, 4 
            k1 = abs(i_mdct(i    ))
            k2 = abs(i_mdct(i + 1))
            k3 = abs(i_mdct(i + 2))
            k4 = abs(i_mdct(i + 3))
            isum0 = isum0 + huff_qa(k1, k2, k3, k4)%leng + k1 + k2 + k3 + k4 
            isum1 = isum1 + huff_qb(k1, k2, k3, k4)%leng + k1 + k2 + k3 + k4
        end do
        if (isum0  <= isum1) then
            itab = 0 ! use table a
            isum = isum + isum0
        else
            itab = 1 ! use table b
            isum = isum + isum1 
        end if
    end subroutine count_one
!----------------------------------------------------------------------------------------------
    subroutine sub_divide(ibigvalues, iregion0, iregion1) ! iso c.1.5.4.4.6
        integer, intent(in ) :: ibigvalues
        integer, intent(out) :: iregion0, iregion1
        integer :: n0, n1, i
        n0 = 2 * ibigvalues * r0  
        n1 = 2 * ibigvalues * r1  
        ! division suggested in iso document c.1.5.4.4.6 is 1/3 : 5/12 : 1/4     
        do i = 0, 20
            if ( n0 >= iscalefactorband_l(i, 3) ) iregion0 = min( 15, max( 0, i               ) ) 
            if ( n1 >= iscalefactorband_l(i, 3) ) iregion1 = min(  7, max( 0, i - iregion0 - 1) )
        end do
    end subroutine sub_divide
!----------------------------------------------------------------------------------------------
    subroutine select_table2(ix, n0, n1, itab, isum) ! c.1.5.4.4.5
        integer, intent(in    ) :: ix(:), n0, n1
        integer, intent(   out) :: itab
        integer, intent(in out) :: isum
        integer :: imax, isum0, isum1, isum2, isum3, i, itab1, itab2
        imax = maxval( abs(ix(n0:n1)) )
        itab = -999
        if (imax <= 15) then 
            do i = 13, 0, -1
                if (imax <= huff(i)%nmax) itab = i
            end do
            select case (itab)
            case (0)
               itab = 0
            case (1)
               itab = 1
            case (2) 
               call count_bits(2, ix, n0, n1, isum1)
               call count_bits(3, ix, n0, n1, isum2)
               if (isum1 < isum2) then 
                  itab = 2
               else
                  itab = 3
               end if
            case (5)
                call count_bits(5, ix, n0, n1, isum1)
                call count_bits(6, ix, n0, n1, isum2)
                if (isum1 < isum2) then 
                    itab = 5
                else
                    itab = 6
                end if
            case (7)
                call count_bits(7, ix, n0, n1, isum1)
                call count_bits(8, ix, n0, n1, isum2)
                call count_bits(9, ix, n0, n1, isum3)
                isum0 = min(isum1, isum2, isum3)
                if      (isum1 == isum0) then 
                    itab = 7
                else if (isum2 == isum0) then
                    itab = 8
                else
                    itab = 9
                end if
            case (10)
                call count_bits(10, ix, n0, n1, isum1)
                call count_bits(11, ix, n0, n1, isum2)
                call count_bits(12, ix, n0, n1, isum3)
                isum0 = min(isum1, isum2, isum3)
                if      (isum1 == isum0) then 
                    itab = 10
                else if (isum2 == isum0) then
                    itab = 11
                else
                    itab = 12
                end if
            case (13)
                call count_bits(13, ix, n0, n1, isum1)
                call count_bits(15, ix, n0, n1, isum2)
                if (isum1 < isum2) then 
                    itab = 13
                else
                    itab = 15
                end if
            case default
                stop 'something wrong : select_table2 : <= 15'
            end select    
        else
            do i = 16, 23
                if (imax <= huff(i)%linmax + 15) then 
                    itab1 = i
                    exit
                end if
            end do
            do i = 24, 31
                if (imax <= huff(i)%linmax + 15) then 
                    itab2 = i
                    exit
                end if
            end do
            call count_bits(itab1, ix, n0, n1, isum1)
            call count_bits(itab2, ix, n0, n1, isum2)
            if (isum1 < isum2) then 
               itab = itab1
            else
               itab = itab2
            end if
        end if
        call count_bits(itab, ix, n0, n1, isum0)
        isum = isum + isum0
    end subroutine select_table2
!---------------------------------------------------------------------------------------------
    subroutine count_bits(nhuff, ix, n0, n1, isum) ! c.1.5.4.4.8
        integer   , intent(in ) :: nhuff, ix(:), n0, n1
        integer   , intent(out) :: isum
        integer :: i, k1, k2
        isum = 0
        if (nhuff ==  0) return
        if (nhuff <= 15) then 
            do i = n0, n1, 2
                k1 = abs( ix(i    ) )
                k2 = abs( ix(i + 1) )
                isum = isum + huff(nhuff)%leng(k1, k2)
                if (k1 /= 0) isum = isum + 1 ! one bit for sign
               if (k2 /= 0) isum = isum + 1 ! one bit for sign
            end do
        else
            do i = n0, n1, 2
                k1 = abs( ix(i    ) )
                k2 = abs( ix(i + 1) )
                if (k1 > 14) then 
                    k1 = 15
                    isum = isum + huff(nhuff)%linbits
                end if
                if (k2 > 14) then 
                    k2 = 15
                    isum = isum + huff(nhuff)%linbits
                end if
                isum = isum + huff(nhuff)%leng(k1, k2)
                if (ix(i    ) /= 0) isum = isum + 1 ! one bit for sign
                if (ix(i + 1) /= 0) isum = isum + 1 ! one bit for sign
            end do
        end if
    end subroutine count_bits
!--------------------------------------------------------------------------------------------
end module mod_inner_loop