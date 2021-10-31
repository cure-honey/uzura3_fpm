module mod_mpg
    use kind_m
    implicit none
    public
    type:: mpeg_parameters
        integer :: mtype
        integer :: layer
        integer :: ibit_rate
        integer :: isample_rate
        integer :: ipadding
        integer :: iprivate
        integer :: icrc
        integer :: mode
        integer :: mode_extension
        integer :: icopyright
        integer :: ioriginal
        integer :: iemphasis
    end type mpeg_parameters
    !mpeg1 / audio
    integer, parameter :: mpeg_frame_size(3)      = (/1152, 1152, 384/)         ! iso 2.4.2.1 frame
    integer, parameter :: mpeg_sample_rates(0:3)  = (/44100, 48000, 32000, 0/)  ! iso 2.4.2.3 sampling_frequency
    integer, parameter :: mpeg_bit_rates(0:14, 3) = &                           ! iso 2.4.2.3 bitrate_index
      reshape( (/ 0, 32, 40, 48,  56,  64,  80,  96, 112, 128, 160, 192, 224, 256, 320,    &
                  0, 32, 48, 56,  64,  80,  96, 112, 128, 160, 192, 224, 256, 320, 384,    &
                  0, 32, 64, 96, 128, 160, 192, 224, 256, 288, 320, 352, 384, 414, 448 /), &
               (/15, 3/) )
    character (len = 8) :: mpeg_mode_names(4)      = (/'stereo  ', 'j-stereo', 'dual-ch ', 'mono    '/)  ! iso 2.4.2.3 mode
    character (len = 3) :: mpeg_layer_names(3)     = (/'iii', 'ii ', 'i  '/)                             ! iso 2.4.2.3 layer
    character (len = 7) :: mpeg_version_names(0:3) = (/'mpeg2.5', '       ', 'mpeg-ii', 'mpeg-i '/)      ! iso 2.4.2.3 id ! mpeg2.5 non-standard 
    character (len = 7) :: mpeg_demp_names(4)      = (/'none   ', '50/15us', '       ', 'citt   '/)      ! iso 2.4.2.3 emphasis
!-------------------------------------------------------------------------------------------
!mpeg1 / layer3:   iso 2.4.1.7, 2.4.2.7 
    type :: side_info_sub
        integer :: ipart2_3_length        ! 12 bits
        integer :: ibig_values            !  9 bits
        integer :: iglobal_gain           !  8 bits
        integer :: iscalefac_compress     !  4 bits
        integer :: iwindow_switching_flag !  1 bit
       ! if short block 
        integer ::  iblock_type           !  2 bits
        integer ::  mixied_block_flag     !  1 bit
        integer ::  itable_select(3)      !  5 bits
        integer ::  isubblock_gain(3)     !  3 bits
       ! if long block
       !integer ::  itable_select(3, 2, 2) 
        integer ::  iregion0_count        !  4 bits
        integer ::  iregion1_count        !  3 bits
       ! 
        integer :: ipreflag               !  1 bit
        integer :: iscalefac_scale        !  1 bit
        integer :: icount1table_select    !  1 bit
    !........
        integer :: icount1  ! local use (not real side_info)
    end type side_info_sub
    !
    type :: side_info_                 !136 bits for mono / 256 bits for stereo
        integer :: main_data_begin        !  9 bits
        integer :: iprivate_bits          !  5 bits for mono /   3 bits for stereo
        integer :: iscfsi(4, 2)           !  1 bit 
        type (side_info_sub) :: sub(2, 2)
    end type side_info_
    !
    integer, parameter :: len_scale_compress(0:15, 2) = & ! iso 2.4.2.7 scalefac_compress[gr][ch]
             reshape( (/ 0, 0, 0, 0, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, &
                         0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 3  /), (/16, 2/) )
    integer, save :: iscalefactorband_l(0:20, 3), iscalefactorband_s(0:11, 3) ! public
!-------------------------------------------------------------------------------------------
    !global variables
    !debug variables
    real (kind = kd), save :: tot_sc_l(0:20) = 0.0_kd, tot_sc_s(0:11, 3) = 0.0_kd
    integer, save :: m1 = 0, m3 = 0, mix = 0, long = 0, nshort = 0, ntable(0:31) = 0, ntab_ab(0:1) = 0
    integer, save :: ms = 0, ns = 0, ns1 = 0, ns2 = 0, nn1 = 0, nn2 = 0, nbits(14)
    integer, save :: n_sc_l = 0, n_sc_s = 0, n_sub_gain = 0, n_scale = 0, n_emph = 0
    !system parameters
    ! switches
    logical, save :: qms_stereo = .true., q_alias = .true., q_mask = .true., q_sm = .true.
    logical, save :: q_vbr = .false., q_rio500 = .false., q_info = .true.
    real (kind = kd), save :: cut_factor = 1.0_kd, distortion_max = 0.0_kd, skip = 1.3_kd
    ! parameters
    integer         , save :: icut = 26                 ! cut at 24 16.5khz; 25 17.2khz; 26 17.9khz; 27 18.6khz
    integer         , save :: mblock_type_param = 20    ! default short block! long = 0, short = 20, mixed = 21
    real (kind = kd), save :: ath_min   = -125.0_kd     ! offset  for ath  90.3db = 2^-15 16bit wav pcm assumed   
    real (kind = kd), save :: ath_max   =   0.0_kd      ! ceiling for ath                    (see init_absolute_threshold inpsycho.f90) 
    real (kind = kd), save :: switch    =   1.1_kd      ! long/short window switching factor (see switch_q in psycho.f90)
    real (kind = kd), save :: xms       =   0.8_kd      ! ns/ms switching factor             (see mid_side in layer3.f90)
    real (kind = kd), save :: xsm       =   1.5_kd      ! short/mixed switching factor       (see mid_side in layer3.f90)
    real (kind = kd), save :: offset    =  40.0_kd ![db]! offset for masking                 (see psycho in psycho.f90)
    real (kind = kd), save :: tempo     =   0.85_kd     ! temporal masking parameter         (see psycho in psycho.f90)
    real (kind = kd), save :: pm_factor =   1.0_kd      ! factor for psychoacoustic moment   (see psycho in psycho.f90) 
    real (kind = kd), save :: factor    =   0.4_kd      ! distribute bits between 2 granules * n channels by total intensity (see av_bits in layer3.f90) 
    real (kind = kd), save :: r0 = 0.33_kd, r1 = 0.75_kd! iso suggests r0 = 0.33_kd, r1 = 0.75_kd (see layer3.f90)
!-------------------------------------------------------------------------------------------
end module mod_mpg