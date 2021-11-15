module mod_polyphase
    use kind_m
    implicit none
    private
    public polyphase_filter36
    integer, parameter :: nsize = 512, nband = 32
    real (kind = kd), save :: coeff(0:nband - 1, 0:nband - 1), prototype(nsize)
    real (kind = kd), save :: pi
contains
!----------------------------------------------------------------------
    subroutine load_window(w)
        real (kind = kd), intent(out) :: w(:) ! iso table c.1 coefficients ci of the analysis window
        w(  1:  4) = (/  0.000000000d0, -0.000000477d0, -0.000000477d0, -0.000000477d0 /) 
        w(  5:  8) = (/ -0.000000477d0, -0.000000477d0, -0.000000477d0, -0.000000954d0 /) 
        w(  9: 12) = (/ -0.000000954d0, -0.000000954d0, -0.000000954d0, -0.000001431d0 /) 
        w( 13: 16) = (/ -0.000001431d0, -0.000001907d0, -0.000001907d0, -0.000002384d0 /) 
        w( 17: 20) = (/ -0.000002384d0, -0.000002861d0, -0.000003338d0, -0.000003338d0 /) 
        w( 21: 24) = (/ -0.000003815d0, -0.000004292d0, -0.000004768d0, -0.000005245d0 /) 
        w( 25: 28) = (/ -0.000006199d0, -0.000006676d0, -0.000007629d0, -0.000008106d0 /) 
        w( 29: 32) = (/ -0.000009060d0, -0.000010014d0, -0.000011444d0, -0.000012398d0 /) 
        w( 33: 36) = (/ -0.000013828d0, -0.000014782d0, -0.000016689d0, -0.000018120d0 /) 
        w( 37: 40) = (/ -0.000019550d0, -0.000021458d0, -0.000023365d0, -0.000025272d0 /) 
        w( 41: 44) = (/ -0.000027657d0, -0.000030041d0, -0.000032425d0, -0.000034809d0 /) 
        w( 45: 48) = (/ -0.000037670d0, -0.000040531d0, -0.000043392d0, -0.000046253d0 /) 
        w( 49: 52) = (/ -0.000049591d0, -0.000052929d0, -0.000055790d0, -0.000059605d0 /) 
        w( 53: 56) = (/ -0.000062943d0, -0.000066280d0, -0.000070095d0, -0.000073433d0 /) 
        w( 57: 60) = (/ -0.000076771d0, -0.000080585d0, -0.000083923d0, -0.000087261d0 /) 
        w( 61: 64) = (/ -0.000090599d0, -0.000093460d0, -0.000096321d0, -0.000099182d0 /) 
        w( 65: 68) = (/  0.000101566d0,  0.000103951d0,  0.000105858d0,  0.000107288d0 /) 
        w( 69: 72) = (/  0.000108242d0,  0.000108719d0,  0.000108719d0,  0.000108242d0 /) 
        w( 73: 76) = (/  0.000106812d0,  0.000105381d0,  0.000102520d0,  0.000099182d0 /) 
        w( 77: 80) = (/  0.000095367d0,  0.000090122d0,  0.000084400d0,  0.000077724d0 /) 
        w( 81: 84) = (/  0.000069618d0,  0.000060558d0,  0.000050545d0,  0.000039577d0 /) 
        w( 85: 88) = (/  0.000027180d0,  0.000013828d0, -0.000000954d0, -0.000017166d0 /) 
        w( 89: 92) = (/ -0.000034332d0, -0.000052929d0, -0.000072956d0, -0.000093937d0 /) 
        w( 93: 96) = (/ -0.000116348d0, -0.000140190d0, -0.000165462d0, -0.000191212d0 /) 
        w( 97:100) = (/ -0.000218868d0, -0.000247478d0, -0.000277042d0, -0.000307560d0 /) 
        w(101:104) = (/ -0.000339031d0, -0.000371456d0, -0.000404358d0, -0.000438213d0 /) 
        w(105:108) = (/ -0.000472546d0, -0.000507355d0, -0.000542164d0, -0.000576973d0 /) 
        w(109:112) = (/ -0.000611782d0, -0.000646591d0, -0.000680923d0, -0.000714302d0 /) 
        w(113:116) = (/ -0.000747204d0, -0.000779152d0, -0.000809669d0, -0.000838757d0 /) 
        w(117:120) = (/ -0.000866413d0, -0.000891685d0, -0.000915051d0, -0.000935555d0 /) 
        w(121:124) = (/ -0.000954151d0, -0.000968933d0, -0.000980854d0, -0.000989437d0 /) 
        w(125:128) = (/ -0.000994205d0, -0.000995159d0, -0.000991821d0, -0.000983715d0 /) 
        w(129:132) = (/  0.000971317d0,  0.000953674d0,  0.000930786d0,  0.000902653d0 /) 
        w(133:136) = (/  0.000868797d0,  0.000829220d0,  0.000783920d0,  0.000731945d0 /) 
        w(137:140) = (/  0.000674248d0,  0.000610352d0,  0.000539303d0,  0.000462532d0 /) 
        w(141:144) = (/  0.000378609d0,  0.000288486d0,  0.000191689d0,  0.000088215d0 /) 
        w(145:148) = (/ -0.000021458d0, -0.000137329d0, -0.000259876d0, -0.000388145d0 /) 
        w(149:152) = (/ -0.000522137d0, -0.000661850d0, -0.000806808d0, -0.000956535d0 /) 
        w(153:156) = (/ -0.001111031d0, -0.001269817d0, -0.001432419d0, -0.001597881d0 /) 
        w(157:160) = (/ -0.001766682d0, -0.001937389d0, -0.002110004d0, -0.002283096d0 /) 
        w(161:164) = (/ -0.002457142d0, -0.002630711d0, -0.002803326d0, -0.002974033d0 /) 
        w(165:168) = (/ -0.003141880d0, -0.003306866d0, -0.003467083d0, -0.003622532d0 /) 
        w(169:172) = (/ -0.003771782d0, -0.003914356d0, -0.004048824d0, -0.004174709d0 /) 
        w(173:176) = (/ -0.004290581d0, -0.004395962d0, -0.004489899d0, -0.004570484d0 /) 
        w(177:180) = (/ -0.004638195d0, -0.004691124d0, -0.004728317d0, -0.004748821d0 /) 
        w(181:184) = (/ -0.004752159d0, -0.004737377d0, -0.004703045d0, -0.004649162d0 /) 
        w(185:188) = (/ -0.004573822d0, -0.004477024d0, -0.004357815d0, -0.004215240d0 /) 
        w(189:192) = (/ -0.004049301d0, -0.003858566d0, -0.003643036d0, -0.003401756d0 /) 
        w(193:196) = (/  0.003134727d0,  0.002841473d0,  0.002521515d0,  0.002174854d0 /) 
        w(197:200) = (/  0.001800537d0,  0.001399517d0,  0.000971317d0,  0.000515938d0 /) 
        w(201:204) = (/  0.000033379d0, -0.000475883d0, -0.001011848d0, -0.001573563d0 /) 
        w(205:208) = (/ -0.002161503d0, -0.002774239d0, -0.003411293d0, -0.004072189d0 /) 
        w(209:212) = (/ -0.004756451d0, -0.005462170d0, -0.006189346d0, -0.006937027d0 /) 
        w(213:216) = (/ -0.007703304d0, -0.008487225d0, -0.009287834d0, -0.010103703d0 /) 
        w(217:220) = (/ -0.010933399d0, -0.011775017d0, -0.012627602d0, -0.013489246d0 /) 
        w(221:224) = (/ -0.014358521d0, -0.015233517d0, -0.016112804d0, -0.016994476d0 /) 
        w(225:228) = (/ -0.017876148d0, -0.018756866d0, -0.019634247d0, -0.020506859d0 /) 
        w(229:232) = (/ -0.021372318d0, -0.022228718d0, -0.023074150d0, -0.023907185d0 /) 
        w(233:236) = (/ -0.024725437d0, -0.025527000d0, -0.026310921d0, -0.027073860d0 /) 
        w(237:240) = (/ -0.027815342d0, -0.028532982d0, -0.029224873d0, -0.029890060d0 /) 
        w(241:244) = (/ -0.030526638d0, -0.031132698d0, -0.031706810d0, -0.032248020d0 /) 
        w(245:248) = (/ -0.032754898d0, -0.033225536d0, -0.033659935d0, -0.034055710d0 /) 
        w(249:252) = (/ -0.034412861d0, -0.034730434d0, -0.035007000d0, -0.035242081d0 /) 
        w(253:256) = (/ -0.035435200d0, -0.035586357d0, -0.035694122d0, -0.035758972d0 /) 
        w(257:260) = (/  0.035780907d0,  0.035758972d0,  0.035694122d0,  0.035586357d0 /) 
        w(261:264) = (/  0.035435200d0,  0.035242081d0,  0.035007000d0,  0.034730434d0 /) 
        w(265:268) = (/  0.034412861d0,  0.034055710d0,  0.033659935d0,  0.033225536d0 /) 
        w(269:272) = (/  0.032754898d0,  0.032248020d0,  0.031706810d0,  0.031132698d0 /) 
        w(273:276) = (/  0.030526638d0,  0.029890060d0,  0.029224873d0,  0.028532982d0 /) 
        w(277:280) = (/  0.027815342d0,  0.027073860d0,  0.026310921d0,  0.025527000d0 /) 
        w(281:284) = (/  0.024725437d0,  0.023907185d0,  0.023074150d0,  0.022228718d0 /) 
        w(285:288) = (/  0.021372318d0,  0.020506859d0,  0.019634247d0,  0.018756866d0 /) 
        w(289:292) = (/  0.017876148d0,  0.016994476d0,  0.016112804d0,  0.015233517d0 /) 
        w(293:296) = (/  0.014358521d0,  0.013489246d0,  0.012627602d0,  0.011775017d0 /) 
        w(297:300) = (/  0.010933399d0,  0.010103703d0,  0.009287834d0,  0.008487225d0 /) 
        w(301:304) = (/  0.007703304d0,  0.006937027d0,  0.006189346d0,  0.005462170d0 /) 
        w(305:308) = (/  0.004756451d0,  0.004072189d0,  0.003411293d0,  0.002774239d0 /) 
        w(309:312) = (/  0.002161503d0,  0.001573563d0,  0.001011848d0,  0.000475883d0 /) 
        w(313:316) = (/ -0.000033379d0, -0.000515938d0, -0.000971317d0, -0.001399517d0 /) 
        w(317:320) = (/ -0.001800537d0, -0.002174854d0, -0.002521515d0, -0.002841473d0 /) 
        w(321:324) = (/  0.003134727d0,  0.003401756d0,  0.003643036d0,  0.003858566d0 /) 
        w(325:328) = (/  0.004049301d0,  0.004215240d0,  0.004357815d0,  0.004477024d0 /) 
        w(329:332) = (/  0.004573822d0,  0.004649162d0,  0.004703045d0,  0.004737377d0 /) 
        w(333:336) = (/  0.004752159d0,  0.004748821d0,  0.004728317d0,  0.004691124d0 /) 
        w(337:340) = (/  0.004638195d0,  0.004570484d0,  0.004489899d0,  0.004395962d0 /) 
        w(341:344) = (/  0.004290581d0,  0.004174709d0,  0.004048824d0,  0.003914356d0 /) 
        w(345:348) = (/  0.003771782d0,  0.003622532d0,  0.003467083d0,  0.003306866d0 /) 
        w(349:352) = (/  0.003141880d0,  0.002974033d0,  0.002803326d0,  0.002630711d0 /) 
        w(353:356) = (/  0.002457142d0,  0.002283096d0,  0.002110004d0,  0.001937389d0 /) 
        w(357:360) = (/  0.001766682d0,  0.001597881d0,  0.001432419d0,  0.001269817d0 /) 
        w(361:364) = (/  0.001111031d0,  0.000956535d0,  0.000806808d0,  0.000661850d0 /) 
        w(365:368) = (/  0.000522137d0,  0.000388145d0,  0.000259876d0,  0.000137329d0 /) 
        w(369:372) = (/  0.000021458d0, -0.000088215d0, -0.000191689d0, -0.000288486d0 /) 
        w(373:376) = (/ -0.000378609d0, -0.000462532d0, -0.000539303d0, -0.000610352d0 /) 
        w(377:380) = (/ -0.000674248d0, -0.000731945d0, -0.000783920d0, -0.000829220d0 /) 
        w(381:384) = (/ -0.000868797d0, -0.000902653d0, -0.000930786d0, -0.000953674d0 /) 
        w(385:388) = (/  0.000971317d0,  0.000983715d0,  0.000991821d0,  0.000995159d0 /) 
        w(389:392) = (/  0.000994205d0,  0.000989437d0,  0.000980854d0,  0.000968933d0 /) 
        w(393:396) = (/  0.000954151d0,  0.000935555d0,  0.000915051d0,  0.000891685d0 /) 
        w(397:400) = (/  0.000866413d0,  0.000838757d0,  0.000809669d0,  0.000779152d0 /) 
        w(401:404) = (/  0.000747204d0,  0.000714302d0,  0.000680923d0,  0.000646591d0 /) 
        w(405:408) = (/  0.000611782d0,  0.000576973d0,  0.000542164d0,  0.000507355d0 /) 
        w(409:412) = (/  0.000472546d0,  0.000438213d0,  0.000404358d0,  0.000371456d0 /) 
        w(413:416) = (/  0.000339031d0,  0.000307560d0,  0.000277042d0,  0.000247478d0 /) 
        w(417:420) = (/  0.000218868d0,  0.000191212d0,  0.000165462d0,  0.000140190d0 /) 
        w(421:424) = (/  0.000116348d0,  0.000093937d0,  0.000072956d0,  0.000052929d0 /) 
        w(425:428) = (/  0.000034332d0,  0.000017166d0,  0.000000954d0, -0.000013828d0 /) 
        w(429:432) = (/ -0.000027180d0, -0.000039577d0, -0.000050545d0, -0.000060558d0 /) 
        w(433:436) = (/ -0.000069618d0, -0.000077724d0, -0.000084400d0, -0.000090122d0 /) 
        w(437:440) = (/ -0.000095367d0, -0.000099182d0, -0.000102520d0, -0.000105381d0 /) 
        w(441:444) = (/ -0.000106812d0, -0.000108242d0, -0.000108719d0, -0.000108719d0 /) 
        w(445:448) = (/ -0.000108242d0, -0.000107288d0, -0.000105858d0, -0.000103951d0 /) 
        w(449:452) = (/  0.000101566d0,  0.000099182d0,  0.000096321d0,  0.000093460d0 /) 
        w(453:456) = (/  0.000090599d0,  0.000087261d0,  0.000083923d0,  0.000080585d0 /) 
        w(457:460) = (/  0.000076771d0,  0.000073433d0,  0.000070095d0,  0.000066280d0 /) 
        w(461:464) = (/  0.000062943d0,  0.000059605d0,  0.000055790d0,  0.000052929d0 /) 
        w(465:468) = (/  0.000049591d0,  0.000046253d0,  0.000043392d0,  0.000040531d0 /) 
        w(469:472) = (/  0.000037670d0,  0.000034809d0,  0.000032425d0,  0.000030041d0 /) 
        w(473:476) = (/  0.000027657d0,  0.000025272d0,  0.000023365d0,  0.000021458d0 /) 
        w(477:480) = (/  0.000019550d0,  0.000018120d0,  0.000016689d0,  0.000014782d0 /) 
        w(481:484) = (/  0.000013828d0,  0.000012398d0,  0.000011444d0,  0.000010014d0 /) 
        w(485:488) = (/  0.000009060d0,  0.000008106d0,  0.000007629d0,  0.000006676d0 /) 
        w(489:492) = (/  0.000006199d0,  0.000005245d0,  0.000004768d0,  0.000004292d0 /) 
        w(493:496) = (/  0.000003815d0,  0.000003338d0,  0.000003338d0,  0.000002861d0 /) 
        w(497:500) = (/  0.000002384d0,  0.000002384d0,  0.000001907d0,  0.000001907d0 /) 
        w(501:504) = (/  0.000001431d0,  0.000001431d0,  0.000000954d0,  0.000000954d0 /) 
        w(505:508) = (/  0.000000954d0,  0.000000954d0,  0.000000477d0,  0.000000477d0 /) 
        w(509:512) = (/  0.000000477d0,  0.000000477d0,  0.000000477d0,  0.000000477d0 /)
    end subroutine load_window
    !----------------------------------------------------------------------
    subroutine initialize_polyphase()
        integer :: i, j
        pi = 4.0_kd * atan(1.0_kd)  
        call load_window(prototype)  
        ! c.1.3 analysis subbband filter
        forall (i = 0:31, j = 0:31) coeff(i, j) = cos(mod((2 * i + 1) * j, 128) * pi / 64.0_kd)  
    end subroutine initialize_polyphase
!----------------------------------------------------------------------
    subroutine polyphase_filter(z, s)
        real (kind = kd), intent(in ) :: z(:)
        real (kind = kd), intent(out) :: s(:)
        real (kind = kd) ::  y(64), f(0:31)
        integer :: i
        forall (i = 1:64) y(i) = sum(z(i:512:64)) 
        f( 0)    = y(17)
        f( 1:16) = y(18:33) + y(16: 1:-1)
        f(17:31) = y(34:48) - y(64:50:-1)
        s = matmul(coeff, f)
    end subroutine polyphase_filter
!----------------------------------------------------------------------
    subroutine polyphase_filter36(pcm, subband)
        real (kind = kd), intent(in ) :: pcm(:, :)
        real (kind = kd), intent(out) :: subband(:, :, :)
        real (kind = kd) :: tmp(512)
        integer        :: i, ichannel, istart, iend
        logical, save :: qfirst = .true.
        if (qfirst) then
            qfirst = .false.
            call initialize_polyphase()
        end if
        do i =  1, 36
            istart = 32 * (i - 1) + 512 
            iend   = 32 * (i - 1) +   1
            do ichannel = 1, size(pcm, 2)
                tmp = prototype * pcm(istart:iend:-1, ichannel)
                call polyphase_filter(tmp, subband(:, i, ichannel))
            end do
        end do
    end subroutine polyphase_filter36
!----------------------------------------------------------------------
end module mod_polyphase