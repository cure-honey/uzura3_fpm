module mod_encode
    use bit_io
    use mod_mpg
    use mod_layer3
    use mod_huffman
    use mod_crc
    implicit none
    private
    public encode_all
contains
!-------------------------------------------------------------------
    subroutine encode_all(mpg, ix, nchannel, ianc) ! iso 2.4.1
        type (mpeg_parameters), intent(in) :: mpg
        integer               , intent(in) :: ix(:, :, :), nchannel, ianc
        call clear_bit_buff  ()  ! bit strings are written to bit_string (private) in module bit_io
        call encode_header   (mpg) 
        call encode_crc      (mpg, nchannel)
        call encode_side_info(     nchannel) ! side_info, scalefactor
        call encode_part2_3  (ix,  nchannel)
        call encode_ancillary(ianc)
    end subroutine encode_all
!-------------------------------------------------------------------
    subroutine encode_header(mpg)         ! iso 2.4.1.3 
        type (mpeg_parameters), intent(in) :: mpg
        call put_bits_c('11111111111'      )  !sync word
        call put_bits(2, mpg%mtype         )  !mpeg1
        call put_bits(2, mpg%layer         )  !layer 
        call put_bits(1, mpg%icrc          )  !crc check 
        call put_bits(4, mpg%ibit_rate     )  !bitrate 
        call put_bits(2, mpg%isample_rate  )  !sampling frequency 44.1
        call put_bits(1, mpg%ipadding      )  !ipadding
        call put_bits(1, mpg%iprivate      )  !private bit : unused
        call put_bits(2, mpg%mode          )  !stereo
        call put_bits(2, mpg%mode_extension)  !mode
        call put_bits(1, mpg%icopyright    )
        call put_bits(1, mpg%ioriginal     )
        call put_bits(2, mpg%iemphasis     )
    end subroutine encode_header
!-------------------------------------------------------------------
    subroutine encode_crc(mpg, nchannel) ! iso 2.4.3.1, table a.9, table b.5 
        type (mpeg_parameters), intent(in) :: mpg
        integer, intent(in) :: nchannel
        integer :: ichannel, igranule, icrc
        if (mpg%icrc == 0) then ! if crc is on
            icrc = 65535           ! z'0000ffff' crc16 initial value
            call crc16(4, mpg%ibit_rate     , icrc)
            call crc16(2, mpg%isample_rate  , icrc)
            call crc16(1, mpg%ipadding      , icrc)
            call crc16(1, mpg%iprivate      , icrc)
            call crc16(2, mpg%mode          , icrc)
            call crc16(2, mpg%mode_extension, icrc)
            call crc16(1, mpg%icopyright    , icrc)
            call crc16(1, mpg%ioriginal     , icrc)
            call crc16(2, mpg%iemphasis     , icrc)
        ! 
            call crc16(9, side_info%main_data_begin, icrc)
            select case (nchannel)
            case (1) ! mono
                call crc16(5, side_info%iprivate_bits, icrc)
            case (2) ! stereo
                call crc16(3, side_info%iprivate_bits, icrc)
            case default
                stop 'illeagal input nchannel: subroutine encode_crc '
            end select
            do ichannel = 1, nchannel
               call crc16dim(1, side_info%iscfsi(1:4, ichannel), icrc)
            end do
            do igranule = 1, 2
                do ichannel = 1, nchannel
                    call  crc16(12, side_info%sub(igranule, ichannel)%ipart2_3_length       , icrc)
                    call  crc16( 9, side_info%sub(igranule, ichannel)%ibig_values           , icrc)
                    call  crc16( 8, side_info%sub(igranule, ichannel)%iglobal_gain          , icrc)
                    call  crc16( 4, side_info%sub(igranule, ichannel)%iscalefac_compress    , icrc)
                    call  crc16( 1, side_info%sub(igranule, ichannel)%iwindow_switching_flag, icrc)
                    if (side_info%sub(igranule, ichannel)%iwindow_switching_flag == 1) then ! short/mixed-block 
                        call crc16( 2, side_info%sub(igranule, ichannel)%iblock_type           , icrc)
                        call crc16( 1, side_info%sub(igranule, ichannel)%mixied_block_flag     , icrc)
                        call crc16dim( 5, side_info%sub(igranule, ichannel)%itable_select(1:2) , icrc)
                        call crc16dim( 3, side_info%sub(igranule, ichannel)%isubblock_gain(1:3), icrc)
                    else ! long-block  
                        call crc16dim( 5, side_info%sub(igranule, ichannel)%itable_select(1:3) , icrc)
                        call crc16( 4, side_info%sub(igranule, ichannel)%iregion0_count        , icrc)
                        call crc16( 3, side_info%sub(igranule, ichannel)%iregion1_count        , icrc)
                    end if
                    call  crc16( 1, side_info%sub(igranule, ichannel)%ipreflag              , icrc)
                    call  crc16( 1, side_info%sub(igranule, ichannel)%iscalefac_scale       , icrc)
                    call  crc16( 1, side_info%sub(igranule, ichannel)%icount1table_select   , icrc)
                end do
            end do
            call put_bits(16, icrc) ! write result icrc :16bits
        end if
    end subroutine encode_crc
!-------------------------------------------------------------------
    subroutine encode_side_info(nchannel) ! iso 2.4.1.7
        integer, intent(in) :: nchannel
        integer :: ichannel, igranule
        call put_bits(9, side_info%main_data_begin)
        select case (nchannel)
        case (1) ! mono
            call put_bits(5, side_info%iprivate_bits)
        case (2) ! stereo
            call put_bits(3, side_info%iprivate_bits)
        case default
            stop 'illeagal input nchannel: subroutine encode_side_info '
        end select
        do ichannel = 1, nchannel
            call put_bits_dim(1, side_info%iscfsi(1:4, ichannel))
        end do
        do igranule = 1, 2
            do ichannel = 1, nchannel
                call  put_bits(12, side_info%sub(igranule, ichannel)%ipart2_3_length        )
                call  put_bits( 9, side_info%sub(igranule, ichannel)%ibig_values            )
                call  put_bits( 8, side_info%sub(igranule, ichannel)%iglobal_gain           )
                call  put_bits( 4, side_info%sub(igranule, ichannel)%iscalefac_compress     )
                call  put_bits( 1, side_info%sub(igranule, ichannel)%iwindow_switching_flag )
                if (side_info%sub(igranule, ichannel)%iwindow_switching_flag == 1) then 
                    call put_bits( 2, side_info%sub(igranule, ichannel)%iblock_type            )
                    call put_bits( 1, side_info%sub(igranule, ichannel)%mixied_block_flag      )
                    call put_bits_dim( 5, side_info%sub(igranule, ichannel)%itable_select(1:2) )
                    call put_bits_dim( 3, side_info%sub(igranule, ichannel)%isubblock_gain(1:3))
                else  
                    call put_bits_dim( 5, side_info%sub(igranule, ichannel)%itable_select(1:3) )
                    call put_bits( 4, side_info%sub(igranule, ichannel)%iregion0_count         )
                    call put_bits( 3, side_info%sub(igranule, ichannel)%iregion1_count         )
                end if
                call put_bits( 1, side_info%sub(igranule, ichannel)%ipreflag               )
                call put_bits( 1, side_info%sub(igranule, ichannel)%iscalefac_scale        )
                call put_bits( 1, side_info%sub(igranule, ichannel)%icount1table_select    )
            end do
        end do
    end subroutine encode_side_info 
!-------------------------------------------------------------------------------------
    subroutine encode_part2_3(ix, nchannel) ! iso 2.4.1.7
        integer, intent(in) :: ix(:, :, :), nchannel
        integer :: ichannel, igranule, n
        do igranule = 1, 2
            do ichannel = 1, nchannel
            ! scale factors
                if (side_info%sub(igranule, ichannel)%iwindow_switching_flag == 1 .and. &
                    side_info%sub(igranule, ichannel)%iblock_type            == 2         ) then
                    if ( side_info%sub(igranule, ichannel)%mixied_block_flag    == 1) then ! mixed-block
                        n = len_scale_compress(side_info%sub(igranule, ichannel)%iscalefac_compress, 1) 
                        call put_bits_dim(n, scfct(igranule, ichannel)%long(0:7))
                        n = len_scale_compress(side_info%sub(igranule, ichannel)%iscalefac_compress, 1) 
                        call put_bits_dim2(n, scfct(igranule, ichannel)%ishort(3:5, 1:3))
                        n = len_scale_compress(side_info%sub(igranule, ichannel)%iscalefac_compress, 2) 
                        call put_bits_dim2(n, scfct(igranule, ichannel)%ishort(6:11, 1:3))
                    else ! short-block
                        n = len_scale_compress(side_info%sub(igranule, ichannel)%iscalefac_compress, 1) 
                        call put_bits_dim2(n, scfct(igranule, ichannel)%ishort(0:5, 1:3))
                        n = len_scale_compress(side_info%sub(igranule, ichannel)%iscalefac_compress, 2) 
                        call put_bits_dim2(n, scfct(igranule, ichannel)%ishort(6:11, 1:3))
                    end if
                else ! long block
                    if (side_info%iscfsi(1, ichannel) == 0 .or. igranule == 1) then !??????????
                        n = len_scale_compress(side_info%sub(igranule, ichannel)%iscalefac_compress, 1) 
                        call put_bits_dim(n, scfct(igranule, ichannel)%long( 0: 5))
                    end if 
                    if (side_info%iscfsi(2, ichannel) == 0 .or. igranule == 1) then 
                        n = len_scale_compress(side_info%sub(igranule, ichannel)%iscalefac_compress, 1) 
                        call put_bits_dim(n, scfct(igranule, ichannel)%long(6:10))
                    end if
                    if (side_info%iscfsi(3, ichannel) == 0 .or. igranule == 1) then
                        n = len_scale_compress(side_info%sub(igranule, ichannel)%iscalefac_compress, 2) 
                        call put_bits_dim(n, scfct(igranule, ichannel)%long(11:15))
                    end if 
                    if (side_info%iscfsi(4, ichannel) == 0 .or. igranule == 1) then 
                        n = len_scale_compress(side_info%sub(igranule, ichannel)%iscalefac_compress, 2) 
                        call put_bits_dim(n, scfct(igranule, ichannel)%long(16:20))
                    end if
                end if
                ! huffman codes
                call encode_big_region(ix(:, igranule, ichannel), igranule, ichannel)
                call encode_quadruples(ix(:, igranule, ichannel), igranule, ichannel)
             end do
        end do
    end subroutine encode_part2_3
!-------------------------------------------------------------------
    subroutine encode_big_region(ix, igranule, ichannel) ! iso 2.4.1.7
        integer, intent(in) :: ix(:), igranule, ichannel
        integer :: n0, n1, itab, iregion0, iregion1
        select case (side_info%sub(igranule, ichannel)%iblock_type)  
        case (0) ! long-block
            iregion0 = side_info%sub(igranule, ichannel)%iregion0_count
            iregion1 = side_info%sub(igranule, ichannel)%iregion1_count
            n0   = 1
            n1   = iscalefactorband_l(iregion0, 3) + 1 
            itab = side_info%sub(igranule, ichannel)%itable_select(1)
            call encode_big(itab, n0, n1, ix)
            n0   = n1 + 1
            n1   = iscalefactorband_l(iregion0 + iregion1 + 1, 3) + 1
            itab = side_info%sub(igranule, ichannel)%itable_select(2)
            call encode_big(itab, n0, n1, ix)
            n0   = n1 + 1
            n1   = side_info%sub(igranule, ichannel)%ibig_values * 2 
            itab = side_info%sub(igranule, ichannel)%itable_select(3) 
            call encode_big(itab, n0, n1, ix)
        case (1, 3) ! long-block start/stop
            iregion0 = side_info%sub(igranule, ichannel)%iregion0_count
            iregion1 = side_info%sub(igranule, ichannel)%iregion1_count
            n0   = 1
            n1   = iscalefactorband_l(iregion0 + 1, 2) 
            itab = side_info%sub(igranule, ichannel)%itable_select(1)
            call encode_big(itab, n0, n1, ix)
            n0   = n1 + 1
            n1   = side_info%sub(igranule, ichannel)%ibig_values * 2 
            itab = side_info%sub(igranule, ichannel)%itable_select(2)
            call encode_big(itab, n0, n1, ix)
        case (2) ! short-block or mixed-block
            if (side_info%sub(igranule, ichannel)%mixied_block_flag == 0) then ! short block
                iregion0 = side_info%sub(igranule, ichannel)%iregion0_count
                iregion1 = side_info%sub(igranule, ichannel)%iregion1_count
                n0   = 1
                n1   = iscalefactorband_s( (iregion0 + 1) / 3, 2) * 3 
                itab = side_info%sub(igranule, ichannel)%itable_select(1)
                call encode_big(itab, n0, n1, ix)
                n0   = n1 + 1
                n1   = side_info%sub(igranule, ichannel)%ibig_values * 2 
                itab = side_info%sub(igranule, ichannel)%itable_select(2)
                call encode_big(itab, n0, n1, ix)
            else if (side_info%sub(igranule, ichannel)%mixied_block_flag == 1) then ! mixed block
                iregion0 = side_info%sub(igranule, ichannel)%iregion0_count
                iregion1 = side_info%sub(igranule, ichannel)%iregion1_count
                n0   = 1
                n1   = iscalefactorband_l(iregion0, 3) + 1 
                itab = side_info%sub(igranule, ichannel)%itable_select(1)
                call encode_big(itab, n0, n1, ix)
                n0   = n1 + 1
                n1   = side_info%sub(igranule, ichannel)%ibig_values * 2 
                itab = side_info%sub(igranule, ichannel)%itable_select(2)
                call encode_big(itab, n0, n1, ix)
            else
                stop 'error : encode_big_region : mixed_block_flag '
            end if
        case default
            write(*, *) 'error : encode_big_region : iblock_type ', side_info%sub(igranule, ichannel)%iblock_type 
            stop 
        end select
    end subroutine encode_big_region
!-------------------------------------------------------------------
    subroutine encode_quadruples(ix, igranule, ichannel) ! iso 2.4.1.7, 2.4.2.7 huffmancodebits()
        integer, intent(in) :: ix(:), igranule, ichannel
        integer :: i, n0, n1, itab, nbits, k1, k2, k3, k4, is1, is2, is3, is4
        n0 =      1 + side_info%sub(igranule, ichannel)%ibig_values * 2 
        n1 = n0 - 1 + side_info%sub(igranule, ichannel)%icount1 * 4 
        itab = side_info%sub(igranule, ichannel)%icount1table_select
        do i = n0, n1, 4
            k1 = abs( ix(i    ) )
            k2 = abs( ix(i + 1) )
            k3 = abs( ix(i + 2) )
            k4 = abs( ix(i + 3) )
            is1 = iand( 1, ishftc( ix(i    ), 1 ) ) ! get sign bit 
            is2 = iand( 1, ishftc( ix(i + 1), 1 ) ) ! bit 31 of [31....0]
            is3 = iand( 1, ishftc( ix(i + 2), 1 ) ) ! cyclic shift + and b'00...001' 
            is4 = iand( 1, ishftc( ix(i + 3), 1 ) ) ! positive 0 / negative 1
            if (itab == 0) then  ! table a      iso table b.7
                nbits = huff_qa(k1, k2, k3, k4)%leng 
                call put_bits(nbits, huff_qa(k1, k2, k3, k4)%icod)
            else if (itab == 1) then  ! table b      iso table b.7
                nbits = huff_qb(k1, k2, k3, k4)%leng 
                call put_bits(nbits, huff_qb(k1, k2, k3, k4)%icod)
            else
                stop 'error '
            end if 
            if (k1 /= 0) call put_bits(1, is1)
            if (k2 /= 0) call put_bits(1, is2)
            if (k3 /= 0) call put_bits(1, is3)
            if (k4 /= 0) call put_bits(1, is4)
        end do
    end subroutine encode_quadruples
!-------------------------------------------------------------------
    subroutine encode_big(itab, n0, n1, ix)  ! iso 2.4.1.7, 2.4.2.7 huffmancodebits(), c.1.5.3.7
        integer, intent(in) :: itab, n0, n1, ix(:)
        integer :: i, is1, is2, k1, k2, linbitsx, linbitsy
        do i = n0, n1, 2
            if (itab == 0) exit
            k1 = abs( ix(i    ) )
            k2 = abs( ix(i + 1) )
            is1 = iand( 1, ishftc(ix(i    ), 1) ) ! get sign bit ! bit 31 of [31..0]  ! positive = 0  
            is2 = iand( 1, ishftc(ix(i + 1), 1) ) ! cyclic shift + and b'00...001'    ! negative = 1
            if (itab <=15) then
                call put_bits( huff(itab)%leng(k1, k2), huff(itab)%icod(k1, k2) )   
                if (k1 /= 0) call put_bits(1, is1)
                if (k2 /= 0) call put_bits(1, is2)
            else
                if (k1 > 14) then 
                    linbitsx = k1 - 15
                    k1 = 15   
                end if
                if (k2 > 14) then 
                    linbitsy = k2 - 15
                    k2 = 15   
                end if
                call put_bits( huff(itab)%leng(k1, k2), huff(itab)%icod(k1, k2) )   
                if (k1 == 15       ) call put_bits(huff(itab)%linbits, linbitsx)  
                if (ix(i)     /=  0) call put_bits(1, is1)
                if (k2 == 15       ) call put_bits(huff(itab)%linbits, linbitsy)  
                if (ix(i + 1) /=  0) call put_bits(1, is2)
            end if
        end do
    end subroutine encode_big
!-------------------------------------------------------------------
    subroutine encode_ancillary(ianc) ! iso c.1.5.3.6
        integer, intent(in) :: ianc
        integer :: i
        do i = 1, ianc        ! fill remaining bits with 0 
            call put_bits(1, 0)
        end do
    end subroutine encode_ancillary
!-------------------------------------------------------------------
end module mod_encode