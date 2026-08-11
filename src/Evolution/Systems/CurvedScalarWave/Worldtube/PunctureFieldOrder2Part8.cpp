
// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/Worldtube/PunctureFieldOrder2Impl.hpp"

namespace CurvedScalarWave::Worldtube::detail {

// NOLINTNEXTLINE(google-readability-function-size, readability-function-size)
void puncture_field_2_part_8(const std::array<double, order2_n_doubles>& d,
                             const DataVector& Dx, const DataVector& Dy,
                             DynamicBuffer<DataVector>& temps) {
  DataVector& sc_16 = temps.at(3238);
  DataVector& sc_24 = temps.at(3246);
  sc_24 = (-d[344]) * sc_16;
  DataVector& dv_1762 = temps.at(1565);
  DataVector& dv_2293 = temps.at(1886);
  DataVector& dv_3009 = temps.at(2534);
  DataVector& dv_3344 = temps.at(2828);
  DataVector& dv_3580 = temps.at(3053);
  DataVector& dv_740 = temps.at(689);
  DataVector& sc_28 = temps.at(3250);
  sc_28 = (-d[230]) *
              ((6.0 * d[1850] + 162.0) * dv_740 - dv_1762 * (dv_3344 + 9.0) +
               239.0 * dv_2293 - 54.0 * dv_3009) +
          dv_3580;
  DataVector& dv_14 = temps.at(14);
  DataVector& dv_15 = temps.at(15);
  DataVector& dv_2629 = temps.at(2209);
  DataVector& dv_2793 = temps.at(2346);
  DataVector& dv_2948 = temps.at(2476);
  DataVector& dv_565 = temps.at(518);
  sc_28 +=
      (-d[996]) * ((-d[1291]) * (dv_2629 + dv_565) - dv_1762 * dv_2793 +
                   dv_2948 + d[2923] * (d[2344] * dv_14 + d[2345] * dv_15));
  DataVector& dv_1 = temps.at(1);
  DataVector& dv_121 = temps.at(119);
  DataVector& dv_2235 = temps.at(1836);
  DataVector& dv_637 = temps.at(587);
  sc_28 += d[1549] * ((d[1068] * d[2342] + d[1820] - 9.0) * dv_121 -
                      Dy * ((9.0 - d[2294]) * dv_637 + (-d[2343]) * dv_2235 +
                            (83.0 * d[36]))) +
           d[2341] * dv_1;
  DataVector& dv_115 = temps.at(113);
  sc_28 += d[355] * ((d[1076] + 93.0 * d[7] - 25.0) * dv_115 +
                     dv_1 * (d[9] + 261.0 * dv_1));
  DataVector& sc_13 = temps.at(3235);
  sc_13 = d[1250] * sc_28;
  DataVector& dv_3264 = temps.at(2751);
  DataVector& dv_3270 = temps.at(2757);
  DataVector& dv_3325 = temps.at(2811);
  DataVector& dv_3445 = temps.at(2921);
  DataVector& dv_3473 = temps.at(2948);
  DataVector& dv_3540 = temps.at(3013);
  DataVector& dv_3582 = temps.at(3055);
  DataVector& dv_45 = temps.at(45);
  DataVector& dv_802 = temps.at(750);
  DataVector& sc_10 = temps.at(3232);
  sc_10 =
      (-d[1649]) * (dv_3325 + dv_3473) +
      d[1172] * ((d[1680] * d[2301]) * dv_45 +
                 d[36] * (192.0 * dv_3270 - 2.0 * dv_3582 + dv_802 * d[2924]) +
                 dv_3264 + dv_3540) +
      dv_3445 + sc_13;
  DataVector& dv_128 = temps.at(126);
  DataVector& dv_2435 = temps.at(2025);
  DataVector& dv_3040 = temps.at(2480);
  DataVector& dv_3138 = temps.at(2627);
  DataVector& dv_3287 = temps.at(2774);
  sc_10 += d[1549] *
           ((-d[1]) * (21.0 * dv_3138 + d[2924] * (dv_128 - 52.0 * dv_15)) +
            d[2338] * dv_3040 + 42.0 * dv_2435 + 1965.0 * dv_3287);
  DataVector& dv_2786 = temps.at(2339);
  DataVector& dv_5 = temps.at(5);
  sc_10 += d[1640] * (d[1651] * dv_115 - dv_5 * (d[2336] + dv_2786));
  DataVector& dv_1503 = temps.at(1383);
  DataVector& dv_2048 = temps.at(1693);
  DataVector& dv_3308 = temps.at(2795);
  sc_10 +=
      d[1665] *
      (-Dx * ((-d[1]) * (131.0 * dv_1503 + 9.0) + d[1599] + d[2337] * dv_1) +
       d[1240] * dv_2048 + 54.0 * dv_3308);
  DataVector& dv_1552 = temps.at(1426);
  DataVector& dv_1809 = temps.at(1595);
  DataVector& dv_2785 = temps.at(2338);
  DataVector& dv_2801 = temps.at(2354);
  DataVector& dv_3191 = temps.at(2678);
  sc_10 += dv_2785 * (d[2335] * dv_1552 + dv_1809 - dv_2801 + dv_3191);
  DataVector& dv_2667 = temps.at(2230);
  DataVector& dv_2800 = temps.at(2353);
  sc_10 +=
      -dv_2667 * ((-d[273]) * (d[1578] + d[2339] * dv_1552) +
                  (-d[283]) * (d[1567] + dv_2800) +
                  (3.0 * d[50]) * (Dy * d[2337] + d[1712]) +
                  (12.0 * d[48] * d[2921]) * ((132.0 * d[36]) + Dy * d[2340]) +
                  (-d[1659]));
  sc_16 = (d[122] * d[2921]) * sc_10;
  DataVector& dv_2481 = temps.at(2070);
  DataVector& dv_2797 = temps.at(2350);
  DataVector& dv_2828 = temps.at(2380);
  DataVector& dv_3363 = temps.at(2846);
  DataVector& dv_3371 = temps.at(2853);
  DataVector& sc_2 = temps.at(3224);
  sc_2 = (-d[996]) * (Dx * ((366.0 * d[2927] + 7.0) * dv_2481 +
                            (-d[1]) * dv_2797 + d[2350]) -
                      dv_3363 + dv_3371) +
         d[120] * dv_2828;
  DataVector& dv_1676 = temps.at(1510);
  DataVector& dv_2519 = temps.at(2107);
  DataVector& dv_3079 = temps.at(2572);
  DataVector& dv_3465 = temps.at(2940);
  DataVector& dv_3626 = temps.at(3097);
  DataVector& dv_3632 = temps.at(3102);
  sc_2 += d[230] * (-Dx * (d[2170] + dv_3626) + d[1356] * dv_3270 +
                    128.0 * dv_3308) +
          dv_2519 * dv_3465 +
          dv_3632 * ((-166.0 * d[36]) + d[2338] * dv_1676 + d[2339] * dv_3079);
  DataVector& dv_2804 = temps.at(2357);
  DataVector& dv_794 = temps.at(742);
  sc_2 += -dv_2804 * ((-d[2346]) * Dy + (-d[1895]) - 753.0 * dv_794);
  sc_28 = d[1250] * sc_2;
  DataVector& sc_17 = temps.at(3239);
  sc_17 = (d[1401] * d[2343] * d[50] + d[2364] - d[2373] * d[457] +
           7584.0 * d[992] - d[996] * (d[1306] + 55.0)) *
          dv_14;
  DataVector& dv_538 = temps.at(493);
  sc_17 += dv_5 * ((-d[1107]) * ((d[2283] - 130.0) * Dy + (396.0 * d[36])) +
                   (-d[58]) * ((464.0 * d[36]) + d[2373] * dv_538) +
                   d[1254] * ((d[1929] + 109.0) * dv_1552 + d[1724]) +
                   d[391] * (d[1518] + 459.0 * dv_1) + d[879]);
  sc_2 = d[6] * sc_17;
  DataVector& dv_3315 = temps.at(2802);
  DataVector& dv_3461 = temps.at(2936);
  DataVector& dv_3605 = temps.at(3077);
  DataVector& dv_3607 = temps.at(3078);
  DataVector& dv_3608 = temps.at(3079);
  DataVector& dv_3609 = temps.at(3080);
  DataVector& dv_3610 = temps.at(3081);
  DataVector& dv_3611 = temps.at(3082);
  DataVector& dv_3613 = temps.at(3084);
  DataVector& dv_3624 = temps.at(3095);
  DataVector& dv_3627 = temps.at(3091);
  DataVector& dv_3628 = temps.at(3098);
  DataVector& dv_3629 = temps.at(3099);
  sc_13 = -84.0 * dv_3315 - dv_3461 + dv_3605 - dv_3607 + dv_3608 - dv_3609 -
          800.0 * dv_3610 - 804.0 * dv_3611 + dv_3613 - 220.0 * dv_3624 +
          dv_3627 - 800.0 * dv_3628 - 48.0 * dv_3629;
  DataVector& dv_3297 = temps.at(2784);
  DataVector& dv_3614 = temps.at(3085);
  DataVector& dv_3615 = temps.at(3086);
  DataVector& dv_3616 = temps.at(3087);
  DataVector& dv_3618 = temps.at(3089);
  DataVector& dv_3621 = temps.at(3092);
  DataVector& dv_3623 = temps.at(3094);
  DataVector& dv_708 = temps.at(658);
  sc_13 += (-d[2307]) * dv_14 + (-d[2307]) * dv_15 + (-d[2353]) * dv_708 +
           192.0 * dv_3297 + 420.0 * dv_3614 + 804.0 * dv_3615 +
           1260.0 * dv_3616 + 2412.0 * dv_3618 + 3072.0 * dv_3621 +
           1040.0 * dv_3623;
  DataVector& dv_2220 = temps.at(1822);
  DataVector& dv_2307 = temps.at(1899);
  DataVector& dv_2612 = temps.at(2194);
  DataVector& dv_46 = temps.at(46);
  DataVector& dv_829 = temps.at(773);
  sc_13 += (-d[2357]) * dv_829 + (-d[2358]) * dv_2307 + (-d[2368]) * dv_2612 +
           (-d[2368]) * dv_46 + (-108.0 * d[120]) * dv_2220;
  DataVector& dv_1566 = temps.at(1438);
  DataVector& dv_1578 = temps.at(1450);
  DataVector& dv_2211 = temps.at(1814);
  DataVector& dv_258 = temps.at(251);
  DataVector& dv_3557 = temps.at(3030);
  DataVector& dv_635 = temps.at(585);
  sc_13 += (36.0 * d[1792]) *
               ((-d[1443]) * (d[2371] * dv_15 + d[2372] * dv_14) +
                (-d[2370]) * dv_635 + d[1153] * dv_258 + d[1439] * dv_3557) +
           (192.0 * d[151]) * dv_1566 + (208.0 * d[2087]) * dv_1578 +
           d[1767] * dv_2211 + sc_28;
  DataVector& dv_2306 = temps.at(1898);
  DataVector& dv_2730 = temps.at(2283);
  DataVector& dv_3030 = temps.at(2367);
  DataVector& dv_3551 = temps.at(3024);
  DataVector& dv_3630 = temps.at(3100);
  DataVector& dv_3631 = temps.at(3101);
  DataVector& dv_613 = temps.at(564);
  sc_13 += d[2353] * dv_2629 + d[2359] * dv_2306 + d[2365] * dv_2730 +
           d[463] * dv_3631 + dv_3551 * dv_3631 +
           dv_3630 * ((d[1107] - d[2369]) + dv_3030 + 655.0 * dv_613) + sc_2;
  sc_10 = (d[20] * d[55]) * sc_13;
  DataVector& dv_2606 = temps.at(2188);
  DataVector& dv_3265 = temps.at(2752);
  DataVector& dv_3411 = temps.at(2821);
  DataVector& dv_352 = temps.at(335);
  sc_2 = (-d[1088]) * dv_2606 + (-d[416]) * dv_3411 +
         d[1191] * (Dx * (d[2928] * dv_3265 + d[1966] * dv_1552 + d[255]) +
                    d[1240] * dv_352 - dv_3308);
  DataVector& dv_2775 = temps.at(2328);
  DataVector& dv_3267 = temps.at(2754);
  DataVector& dv_3435 = temps.at(2911);
  DataVector& dv_3565 = temps.at(3038);
  DataVector& dv_3599 = temps.at(3072);
  DataVector& dv_833 = temps.at(776);
  sc_2 += d[2922] *
          (d[20] * dv_2775 + d[31] * (d[1106] * dv_3565 + dv_1762) +
           d[395] * (d[1242] * dv_3565 + dv_3267 * dv_833 + dv_3435 + dv_3599));
  DataVector& dv_1683 = temps.at(1516);
  DataVector& dv_2186 = temps.at(1792);
  DataVector& dv_2363 = temps.at(1953);
  DataVector& dv_3519 = temps.at(2993);
  DataVector& dv_972 = temps.at(878);
  sc_2 += Dx * d[1248] *
              ((-d[2921]) * (d[1832] * dv_972 + d[444]) + (-d[524]) - dv_2186) -
          dv_3519 * (d[2322] * dv_1683 + dv_2363);
  sc_13 = d[1581] * sc_2;
  DataVector& dv_154 = temps.at(152);
  DataVector& dv_2585 = temps.at(2171);
  DataVector& dv_2953 = temps.at(2481);
  DataVector& dv_3449 = temps.at(2924);
  DataVector& dv_3603 = temps.at(3075);
  DataVector& dv_482 = temps.at(450);
  DataVector& sc_29 = temps.at(3251);
  sc_29 = (-d[50]) * ((-d[1193]) * (dv_154 + dv_482) + dv_2585 -
                      dv_3449 * (dv_2953 + 2.0) + dv_3603) +
          (116.0 * d[151]) * dv_46;
  DataVector& dv_2441 = temps.at(2031);
  DataVector& dv_558 = temps.at(511);
  DataVector& dv_724 = temps.at(673);
  sc_29 += d[287] * ((-d[2923]) * ((-d[2160]) * dv_15 + d[2331] * dv_14) +
                     d[1242] * dv_724) +
           d[631] * ((d[1140] * d[2332] + 199.0 * d[7] - 12.0) * dv_14 +
                     Dy * ((-d[2244]) * dv_558 + d[1560] + dv_2441));
  sc_29 += d[894] * dv_1;
  sc_17 = d[1625] * sc_29;
  DataVector& dv_2341 = temps.at(1931);
  DataVector& dv_3262 = temps.at(2749);
  DataVector& dv_3296 = temps.at(2783);
  DataVector& dv_3513 = temps.at(2987);
  DataVector& dv_3604 = temps.at(3076);
  DataVector& dv_69 = temps.at(69);
  DataVector& sc_26 = temps.at(3248);
  sc_26 =
      (-d[1251]) * (3.0 * dv_3262 + dv_3604) +
      d[1638] *
          (d[1787] * dv_3270 -
           dv_3296 * ((-d[1]) * (dv_2341 + 2.0) + d[1420] + d[2239] * dv_69) +
           dv_3513);
  DataVector& dv_3397 = temps.at(2874);
  DataVector& dv_486 = temps.at(454);
  DataVector& dv_61 = temps.at(61);
  sc_26 +=
      d[77] * ((-d[1]) * ((-d[2924]) * (dv_486 + dv_61) + 131.0 * dv_3138) +
               Dx * d[2333] * dv_3397 + 262.0 * dv_2435 + 567.0 * dv_3287);
  DataVector& dv_1502 = temps.at(1382);
  sc_26 += Dx * d[104] * ((-d[1598] + d[1846]) * Dy + d[1636] * dv_1502);
  sc_29 = sc_26 * d[2921];
  DataVector& dv_2234 = temps.at(1835);
  sc_28 = d[1652] * (Dy * ((-69.0 * d[2921]) + dv_2234) + 272.0 * dv_14) +
          sc_17 + sc_29;
  DataVector& dv_240 = temps.at(236);
  sc_28 += -dv_2667 * ((-d[259]) * ((-d[1908]) + d[2329] * dv_1683) +
                       d[104] * ((d[1806] + 2.0) * dv_240 + d[1636]) + d[1633] +
                       d[348] * ((169.0 * d[36]) + d[2330] * dv_558));
  sc_2 = d[362] * sc_28;
  DataVector& dv_0 = temps.at(0);
  DataVector& dv_1115 = temps.at(1012);
  DataVector& dv_3421 = temps.at(2897);
  DataVector& dv_3596 = temps.at(3069);
  DataVector& sc_12 = temps.at(3234);
  DataVector& sc_14 = temps.at(3236);
  DataVector& sc_18 = temps.at(3240);
  DataVector& sc_6 = temps.at(3228);
  DataVector& sc_8 = temps.at(3230);
  sc_8 = (-d[2320]) * dv_3421 +
         (3.0 * d[1931]) * ((-d[1229]) * dv_45 + d[554] * dv_1115 +
                            dv_0 * (d[1622] + dv_3596)) +
         (-d[1327] * d[1675]) * dv_45 + sc_10 + sc_12 + sc_14 + sc_16 + sc_18 +
         sc_24 + sc_6;
  DataVector& dv_3154 = temps.at(2643);
  DataVector& dv_3227 = temps.at(2714);
  DataVector& dv_3597 = temps.at(3070);
  DataVector& dv_3598 = temps.at(3071);
  sc_8 += (2.0 * d[2928] * d[340] * d[6]) *
              ((-d[631]) * dv_3597 + d[230] * dv_1 +
               d[50] * ((-d[1106]) * dv_3598 + dv_1762) - dv_3154 +
               348.0 * dv_3227) +
          sc_13 + sc_2;
  DataVector& dv_2617 = temps.at(2199);
  DataVector& dv_2647 = temps.at(1425);
  sc_8 += (d[1618] * d[2922]) * Dx *
          ((-d[77]) * (d[1086] + d[2322] * dv_637 + dv_2647) +
           d[2321] * dv_2617 + d[348] * dv_3596);
  DataVector& sc_19 = temps.at(3241);
  sc_19 = (-d[422]) * sc_8;
  DataVector& dv_2243 = temps.at(143);
  DataVector& dv_2247 = temps.at(1845);
  DataVector& dv_2249 = temps.at(1847);
  DataVector& dv_2251 = temps.at(1849);
  DataVector& dv_2575 = temps.at(2162);
  DataVector& dv_3176 = temps.at(2663);
  DataVector& dv_3177 = temps.at(2664);
  DataVector& dv_34 = temps.at(34);
  DataVector& dv_561 = temps.at(514);
  sc_24 = d[6] * (dv_2243 + dv_2575) + d[7] * (dv_3176 + dv_34 + dv_561) +
          dv_0 * ((2.0 * d[2928]) * (dv_3177 + 8.0) + (64.0 * d[2923]) * Dy +
                  (-d[1501]) - dv_2247 - dv_2249) +
          dv_2251;
  DataVector& dv_591 = temps.at(543);
  sc_24 += d[2925] * ((-d[1362]) * dv_591 + d[29] * dv_561 + dv_3176 * d[2921]);
  sc_16 = (-d[2928]) * sc_24;
  DataVector& dv_1878 = temps.at(1600);
  DataVector& dv_2184 = temps.at(1785);
  DataVector& dv_2244 = temps.at(529);
  DataVector& dv_241 = temps.at(237);
  DataVector& dv_545 = temps.at(499);
  DataVector& dv_615 = temps.at(566);
  sc_10 = d[1240] * (d[1928] * dv_45 + dv_2244) +
          d[1273] * ((-d[2923]) *
                         (d[139] * dv_1878 + d[180] * dv_615 + d[37] * dv_545) +
                     d[1822] * dv_241 + d[6] * dv_2184) +
          sc_16;
  sc_13 = (-d[19]) * sc_10;
  DataVector& dv_1889 = temps.at(517);
  DataVector& dv_531 = temps.at(486);
  sc_6 =
      d[20] * ((d[1] * d[1031]) * dv_1889 + d[1081] * dv_531 + d[1677] * dv_1);
  DataVector& dv_16 = temps.at(16);
  DataVector& dv_2262 = temps.at(1858);
  DataVector& dv_2868 = temps.at(2418);
  DataVector& dv_2892 = temps.at(2438);
  DataVector& dv_2921 = temps.at(2458);
  DataVector& dv_525 = temps.at(480);
  sc_6 +=
      d[77] * ((12.0 * d[147]) * dv_16 + (80.0 * d[2927] * d[2923]) * dv_525 -
               dv_2262 - 27.0 * dv_2868 - dv_2892 - dv_2921);
  DataVector& dv_2258 = temps.at(1854);
  DataVector& dv_2259 = temps.at(1855);
  DataVector& dv_3179 = temps.at(2666);
  DataVector& dv_669 = temps.at(619);
  DataVector& dv_99 = temps.at(97);
  sc_6 += d[91] * ((1.0 - d[1355]) * dv_99 + d[7] * (dv_2259 - dv_3179) +
                   dv_2258 + dv_669 * d[2927]);
  sc_24 = sc_6 * d[2922];
  DataVector& dv_2253 = temps.at(542);
  DataVector& dv_2257 = temps.at(522);
  sc_16 =
      (-d[1402]) * ((-d[2928]) * dv_2253 + (50.0 * d[280] * d[2927]) * dv_16) +
      (d[2928] * d[2924]) * dv_2257 + sc_24;
  DataVector& dv_1631 = temps.at(1476);
  DataVector& dv_2387 = temps.at(1977);
  DataVector& dv_2500 = temps.at(2088);
  DataVector& dv_3178 = temps.at(2665);
  sc_16 += -Dx * ((-d[795]) * (d[1076] * dv_1 + dv_2500) +
                  d[209] * ((d[1929] + 5.0) * dv_2387 + d[88]) +
                  d[91] * ((d[1930] - 1.0) * dv_1631 + d[1403] + dv_3178));
  DataVector& dv_2240 = temps.at(1841);
  DataVector& dv_2254 = temps.at(1851);
  DataVector& dv_3065 = temps.at(2558);
  sc_16 += -dv_3065 * ((-d[169] + 320.0 * d[2927] + 25.0) * dv_5 +
                       d[1] * (d[9] + dv_2240) + dv_2254);
  sc_10 = (-d[2920]) * sc_16;
  DataVector& dv_1962 = temps.at(1678);
  DataVector& dv_3159 = temps.at(2648);
  DataVector& dv_3160 = temps.at(2649);
  DataVector& dv_3162 = temps.at(2651);
  DataVector& dv_3163 = temps.at(2652);
  DataVector& dv_3164 = temps.at(2653);
  DataVector& dv_3167 = temps.at(2655);
  DataVector& dv_3168 = temps.at(1824);
  DataVector& dv_3170 = temps.at(2657);
  DataVector& dv_3171 = temps.at(2658);
  DataVector& dv_3174 = temps.at(2661);
  sc_2 = -60.0 * dv_1962 + dv_3159 - 240.0 * dv_3160 - 45.0 * dv_3162 -
         dv_3163 - 172.0 * dv_3164 - 400.0 * dv_3167 - dv_3168 -
         160.0 * dv_3170 + 33.0 * dv_3171 - 76.0 * dv_3174;
  DataVector& dv_1200 = temps.at(1095);
  DataVector& dv_2212 = temps.at(1815);
  DataVector& dv_2692 = temps.at(2247);
  DataVector& dv_3155 = temps.at(2644);
  DataVector& dv_3161 = temps.at(2650);
  DataVector& dv_3165 = temps.at(2654);
  DataVector& dv_3166 = temps.at(1414);
  DataVector& dv_3173 = temps.at(2660);
  DataVector& dv_756 = temps.at(704);
  sc_2 += (-d[1033]) * dv_3154 + (-d[1079]) * dv_3155 + (-d[1081]) * dv_3165 +
          (-d[1161]) * dv_1200 + (-d[1161]) * dv_756 + (-d[1388]) * dv_2212 +
          (-d[1395]) * dv_2692 + 480.0 * dv_3161 + 86.0 * dv_3166 +
          200.0 * dv_3173;
  DataVector& dv_1116 = temps.at(1013);
  DataVector& dv_1478 = temps.at(1359);
  DataVector& dv_190 = temps.at(187);
  DataVector& dv_1961 = temps.at(1677);
  DataVector& dv_3158 = temps.at(2647);
  DataVector& dv_691 = temps.at(641);
  DataVector& dv_732 = temps.at(681);
  sc_2 += (-d[1434]) * dv_1116 + (-d[1469]) * dv_732 + (-d[151]) * dv_3158 +
          (-d[1655]) * dv_190 + (-d[1919]) * dv_352 + (-d[1921]) * dv_691 +
          (-d[221]) * dv_1478 + (-d[6]) * dv_1961 + (-d[7]) * dv_3154 +
          (-d[751]) * dv_15 + sc_13;
  DataVector& dv_1547 = temps.at(1421);
  DataVector& dv_3157 = temps.at(2646);
  DataVector& dv_3175 = temps.at(2662);
  DataVector& dv_48 = temps.at(48);
  DataVector& dv_594 = temps.at(545);
  sc_2 += (-d[930]) * dv_1547 + (-172.0 * d[604]) * dv_3175 +
          (80.0 * d[1923]) * dv_1200 + (240.0 * d[1923]) * dv_48 +
          (720.0 * d[2927]) * dv_3157 + (d[1030] * d[1543]) * dv_594 +
          (d[1079] * d[276]) * dv_1200 + (d[1242] * d[1292]) * dv_46 +
          (d[1738] * d[273]) * dv_48 + (-d[1030] * d[1926]) * dv_16 + sc_10;
  DataVector& dv_2200 = temps.at(1804);
  DataVector& dv_2210 = temps.at(1813);
  DataVector& dv_263 = temps.at(256);
  DataVector& dv_3151 = temps.at(2640);
  DataVector& dv_3152 = temps.at(2641);
  sc_2 += (-d[1030] * d[604]) * dv_669 + (-d[1033] * d[1161]) * dv_16 +
          (-d[1242] * d[1537]) * dv_756 + (-d[255] * d[367]) * dv_1200 +
          d[1053] * dv_2200 + d[1079] * dv_3174 + d[1081] * dv_3152 +
          d[1388] * dv_263 + d[1434] * dv_2210 + d[1486] * dv_3151;
  DataVector& dv_123 = temps.at(121);
  DataVector& dv_143 = temps.at(141);
  DataVector& dv_1548 = temps.at(1422);
  DataVector& dv_2082 = temps.at(1712);
  DataVector& dv_2918 = temps.at(2457);
  DataVector& dv_3153 = temps.at(2642);
  DataVector& dv_3169 = temps.at(2656);
  DataVector& dv_508 = temps.at(469);
  sc_2 += d[1543] * dv_1548 + d[1595] * dv_3153 + d[1726] * dv_143 +
          d[1918] * dv_2082 + d[1919] * dv_3169 + d[1919] * dv_508 +
          d[1920] * dv_99 + d[1921] * dv_123 + d[1922] * dv_99 +
          d[1923] * dv_2918;
  DataVector& dv_2217 = temps.at(1819);
  DataVector& dv_2219 = temps.at(1821);
  sc_2 += d[1924] * dv_48 + d[1925] * dv_3157 + d[1927] * dv_115 +
          d[273] * dv_2217 + d[273] * dv_2219 + d[371] * dv_1;
  DataVector& dv_1017 = temps.at(917);
  DataVector& dv_1031 = temps.at(931);
  DataVector& dv_1546 = temps.at(1420);
  DataVector& dv_2196 = temps.at(1800);
  DataVector& dv_2228 = temps.at(1830);
  DataVector& dv_2232 = temps.at(513);
  DataVector& dv_2245 = temps.at(516);
  DataVector& dv_2829 = temps.at(2381);
  DataVector& dv_3156 = temps.at(2645);
  sc_2 +=
      (86.0 * d[48]) * dv_1017 * dv_1503 +
      d[52] * ((-d[2928]) * dv_2232 +
               (-d[1081] * d[2922]) * (dv_1200 + dv_2196) + d[1240] * dv_2228) +
      d[807] * dv_1031 + d[930] * dv_1546 + dv_2245 * dv_3156 -
      dv_2829 * dv_3156;
  sc_8 = (-d[75]) * sc_2;
  DataVector& dv_1513 = temps.at(1392);
  DataVector& dv_2302 = temps.at(1895);
  DataVector& dv_2498 = temps.at(2086);
  DataVector& dv_2700 = temps.at(2253);
  DataVector& dv_3205 = temps.at(2692);
  DataVector& dv_3346 = temps.at(2830);
  DataVector& dv_3514 = temps.at(2988);
  sc_16 = d[1250] * ((-d[2921]) * (Dy * (d[1189] + dv_2700) + dv_2302) +
                     dv_3514 * d[2923]) +
          d[2241] * dv_3205 + dv_1513 * dv_2498 + 20.0 * dv_2606 + dv_3346;
  sc_16 += d[2923] * ((d[1739] + d[2240]) * dv_45 + d[27] * dv_3262);
  sc_13 = (-d[1581]) * sc_16;
  DataVector& dv_3543 = temps.at(3016);
  DataVector& dv_3544 = temps.at(3017);
  DataVector& dv_3545 = temps.at(3018);
  DataVector& dv_450 = temps.at(423);
  sc_12 = d[1322] * (Dx * ((-127.0 * d[255]) + d[2928] * (dv_3544 + 29.0) -
                           dv_2730 - dv_3545) +
                     d[1240] * dv_450 + dv_3543);
  DataVector& dv_1112 = temps.at(1009);
  DataVector& dv_3539 = temps.at(3012);
  sc_12 += d[1409] * (-36.0 * dv_1112 + 19.0 * dv_3262 + dv_3539);
  DataVector& dv_2774 = temps.at(2327);
  DataVector& dv_3541 = temps.at(3014);
  DataVector& dv_3542 = temps.at(3015);
  sc_12 +=
      d[20] * ((-d[37]) * (Dx * (dv_2774 - 147.0) + 64.0 * dv_3270 + dv_3542) +
               (546.0 * d[7]) * Dx * Dy - dv_3540 - dv_3541);
  DataVector& dv_3258 = temps.at(2745);
  sc_12 += dv_3258 * (Dy * d[2248] + d[314] + 513.0 * dv_794);
  sc_6 = (-d[86]) * sc_12;
  DataVector& dv_2867 = temps.at(693);
  DataVector& dv_292 = temps.at(284);
  DataVector& dv_3538 = temps.at(3011);
  sc_18 = (-d[252]) * (930.0 * dv_14 + dv_15) +
          d[20] * ((16.0 * d[2928] * d[2925]) * dv_292 - dv_2867 - dv_3538);
  DataVector& dv_2752 = temps.at(2305);
  DataVector& dv_29 = temps.at(29);
  DataVector& dv_2979 = temps.at(2506);
  sc_18 += d[259] * ((104.0 * d[7] - 195.0 * d[2927] + 119.0) * dv_29 +
                     Dy * ((-d[1355] - 7.0) * Dy + (268.0 * d[36]) + dv_2979)) +
           d[50] * dv_2752;
  sc_12 = (d[535] * d[2921]) * sc_18;
  DataVector& dv_3083 = temps.at(2576);
  DataVector& dv_3274 = temps.at(2761);
  DataVector& dv_3350 = temps.at(2834);
  DataVector& dv_3537 = temps.at(3010);
  sc_24 = (-d[1463]) * (Dy * (d[507] + dv_3083) + dv_3350) +
          (20.0 * d[1792]) * (d[1249] * ((-d[2254]) * dv_14 + dv_3537) +
                              d[252] * dv_3274 + d[781] * dv_740) +
          sc_12 + sc_6;
  DataVector& dv_3452 = temps.at(2927);
  DataVector& dv_3479 = temps.at(2954);
  DataVector& dv_3533 = temps.at(3007);
  DataVector& dv_3534 = temps.at(394);
  DataVector& dv_3535 = temps.at(3008);
  DataVector& dv_3536 = temps.at(3009);
  sc_24 += d[1107] * (dv_3479 + dv_3533 * d[2923]) +
           d[1184] * (d[1] * dv_3536 + d[1106] * dv_3534 +
                      dv_1762 * (3.0 - dv_3535) + dv_3452);
  DataVector& dv_2881 = temps.at(2430);
  DataVector& dv_3312 = temps.at(2799);
  DataVector& dv_3436 = temps.at(2912);
  sc_24 += d[1409] * (d[1242] * dv_3534 + d[152] * dv_3436 + dv_2881 -
                      16.0 * dv_3312) +
           d[271] * dv_46;
  DataVector& dv_1821 = temps.at(1607);
  DataVector& dv_2291 = temps.at(1884);
  DataVector& dv_2942 = temps.at(2473);
  sc_24 += dv_2942 * ((-d[20]) * (203.0 * Dy + d[2245]) +
                      (32.0 * d[2928] * d[2921]) * (d[46] + dv_2291) +
                      (38.0 * d[50]) - 56.0 * dv_1821);
  sc_16 = (-d[299]) * sc_24;
  DataVector& dv_2697 = temps.at(1998);
  DataVector& dv_2899 = temps.at(2444);
  DataVector& dv_2903 = temps.at(715);
  DataVector& dv_3516 = temps.at(2990);
  sc_12 = (d[116] * d[1583]) * dv_2220 +
          d[286] * ((-d[631]) * (-Dy * dv_2697 + dv_29) + d[1527] * dv_3516 +
                    d[50] * dv_3514 + dv_2899 + dv_3154) +
          d[325] * dv_2903;
  DataVector& dv_1576 = temps.at(1448);
  DataVector& dv_2350 = temps.at(1940);
  DataVector& dv_2666 = temps.at(2229);
  DataVector& dv_3219 = temps.at(2706);
  DataVector& dv_3515 = temps.at(2989);
  sc_12 += -dv_2350 *
           ((-d[20]) * (d[314] + dv_3219) + (-d[1100] - d[2242]) * dv_2666 +
            (2.0 * d[50]) * dv_1576 + (4.0 * d[2928] * d[2921]) * dv_3515);
  sc_24 = (-d[340]) * sc_12;
  sc_14 = (-d[225]) * ((47.0 * d[7]) * Dx - 8.0 * dv_3262 - dv_3539) +
          (d[2266] * d[283]) * dv_732;
  DataVector& dv_2743 = temps.at(2296);
  DataVector& dv_2744 = temps.at(2297);
  DataVector& dv_3271 = temps.at(2758);
  DataVector& dv_657 = temps.at(607);
  sc_14 += d[1472] * (Dx * ((5.0 * d[2923] * (d[2156] - 83.0)) * Dy +
                            d[2928] * dv_2743 - dv_2744) +
                      d[1240] * dv_657 + dv_3271);
  DataVector& dv_1583 = temps.at(1455);
  DataVector& dv_3294 = temps.at(2781);
  DataVector& dv_3550 = temps.at(3023);
  sc_14 +=
      d[50] *
      ((-d[37]) * (Dx * (128.0 * dv_1503 - 145.0) + 64.0 * dv_1583 + dv_3550) +
       (256.0 * d[2928] * d[147]) * Dx + (942.0 * d[7]) * Dx * Dy -
       705.0 * dv_3294);
  DataVector& dv_2044 = temps.at(1689);
  DataVector& dv_2722 = temps.at(2275);
  sc_14 += dv_2722 * (d[1751] + d[2262] * dv_2044 + 2430.0 * dv_794);
  sc_18 = (-d[2922]) * sc_14;
  sc_29 = dv_5;
  DataVector& dv_3295 = temps.at(2782);
  DataVector& dv_3387 = temps.at(2869);
  sc_29 *= (-d[20]) * (d[2928] * (83.0 - dv_3295) + d[1748] + dv_3387) +
           (-d[49]) * (d[306] + 381.0 * dv_1) +
           (d[2928] * d[2921]) *
               ((31.0 - d[2249]) * dv_2044 + (656.0 * d[36]) + 240.0 * dv_794) +
           (-d[1715]);
  sc_28 = (d[2271] + 15.0 * d[273] * (d[1595] - d[1880] + 14.0) -
           1529.0 * d[277] + d[50] * (d[1265] - 119.0 * d[2923])) *
              dv_14 +
          sc_29;
  sc_14 = d[535] * sc_28;
  DataVector& dv_2732 = temps.at(2285);
  DataVector& dv_2734 = temps.at(2287);
  DataVector& dv_2880 = temps.at(2429);
  DataVector& dv_3277 = temps.at(2764);
  DataVector& dv_3558 = temps.at(3031);
  sc_6 = d[1184] * ((-d[1193]) * dv_3558 + d[1338] * dv_3557 -
                    dv_1762 * dv_2732 + dv_3277) +
         d[1409] * (d[1242] * dv_3557 + 83.0 * dv_2692 + dv_2734 * dv_2880 -
                    57.0 * dv_3312) +
         sc_18;
  DataVector& dv_2301 = temps.at(1894);
  DataVector& dv_2963 = temps.at(2491);
  DataVector& dv_3559 = temps.at(3032);
  sc_6 += d[1431] * (Dy * ((-d[1939]) - 28.0 * dv_1502) - 19.0 * dv_2301) +
          d[1611] * ((-d[147]) * dv_3558 + dv_2963 + dv_3452 + dv_3559 +
                     26.0 * dv_740);
  DataVector& dv_2802 = temps.at(2355);
  DataVector& dv_2914 = temps.at(2453);
  DataVector& dv_3145 = temps.at(2634);
  sc_6 +=
      d[1802] * ((d[1254] * d[2268] + d[2270]) * dv_14 + d[2269] * dv_3145) +
      d[2202] * dv_258 +
      dv_2942 * ((-d[20]) * (d[1751] + dv_2802) +
                 (4.0 * d[2928] * d[2921]) * (d[2928] + dv_2914) + (-d[2173]) -
                 470.0 * dv_1821) +
      sc_14;
  sc_12 = (-d[342]) * sc_6;
  DataVector& dv_2030 = temps.at(196);
  DataVector& dv_3528 = temps.at(3002);
  DataVector& dv_3529 = temps.at(3003);
  sc_28 =
      d[20] * ((-d[1]) * (dv_2030 * d[2924] + dv_3138) +
               d[167] * (Dx + dv_3529) + d[2243] * dv_732 - 90.0 * dv_3287) +
      d[276] * (-20.0 * dv_3138 + dv_3262 + dv_3528);
  DataVector& dv_1817 = temps.at(1603);
  DataVector& dv_2379 = temps.at(1969);
  DataVector& dv_2384 = temps.at(1974);
  DataVector& dv_3428 = temps.at(2904);
  DataVector& dv_3439 = temps.at(2915);
  DataVector& dv_3440 = temps.at(2916);
  DataVector& dv_3523 = temps.at(2997);
  DataVector& dv_3530 = temps.at(3004);
  DataVector& dv_751 = temps.at(699);
  sc_28 += d[77] * ((-d[2246]) * dv_3428 + d[36] * (-dv_3440 - dv_3530) +
                    dv_3439 + dv_3523) -
           dv_1817 * (d[2248] * dv_751 + 28.0 * dv_2379 + dv_2384);
  sc_18 = d[27] * sc_28;
  DataVector& dv_3531 = temps.at(3005);
  sc_29 = (-d[2250]) * dv_3531 +
          (-d[287]) * ((d[1805] + d[2251]) * dv_14 +
                       Dy * ((-d[2252] - 58.0) * Dy + d[314] + 930.0 * dv_794));
  DataVector& dv_2937 = temps.at(2469);
  DataVector& dv_3532 = temps.at(3006);
  sc_29 +=
      d[50] * ((-d[80]) * (-dv_1762 * (dv_3295 - 17.0) + dv_2937 + dv_3532) +
               d[2240] * dv_635) +
      d[59] * (dv_2301 + dv_240 * (d[152] + dv_1502));
  DataVector& dv_3405 = temps.at(2882);
  DataVector& dv_3443 = temps.at(2919);
  DataVector& dv_683 = temps.at(633);
  DataVector& dv_96 = temps.at(94);
  sc_29 += d[724] *
           ((-d[80]) * ((d[1806] - 7.0) * dv_683 + (d[1100] + 3.0) * dv_96) +
            d[1242] * dv_3516 + d[147] * dv_3405 + 256.0 * dv_2293 + dv_3443);
  sc_28 = sc_29 * d[2922];
  DataVector& dv_2713 = temps.at(2266);
  DataVector& dv_2714 = temps.at(2267);
  sc_14 = d[1378] * (Dy * (d[77] * dv_2713 + dv_2714) + d[325] * dv_29) + sc_18;
  DataVector& dv_2115 = temps.at(1732);
  DataVector& dv_2255 = temps.at(1852);
  DataVector& dv_2490 = temps.at(2078);
  sc_14 += dv_2115 * ((-d[50]) * (d[1464] + dv_2255) + (32.0 * d[151]) * Dy +
                      (2.0 * d[2928] * d[20]) *
                          (d[1683] + d[2244] * dv_972 + 64.0 * dv_794) +
                      (-d[1624]) - 1026.0 * dv_2490) +
           sc_28;
  sc_6 = (-d[344]) * sc_14;
  DataVector& dv_2796 = temps.at(2349);
  DataVector& dv_2993 = temps.at(2520);
  DataVector& dv_3525 = temps.at(2999);
  DataVector& dv_3526 = temps.at(3000);
  sc_29 = (-d[1412]) * (dv_2993 + 23.0 * dv_740 + dv_833 * (13.0 - dv_2796)) +
          (-d[50]) * (Dy * ((-d[2247]) + 58.0 * dv_1502) + dv_3526) +
          (d[49] * d[7]) * dv_3525;
  DataVector& dv_155 = temps.at(153);
  DataVector& dv_2320 = temps.at(1912);
  DataVector& dv_3266 = temps.at(2753);
  DataVector& dv_3527 = temps.at(3001);
  sc_29 += d[1865] * ((-d[434]) * dv_3266 + (13.0 * d[2921]) * dv_635) +
           d[77] * (d[1242] * dv_3525 + d[80] * (70.0 * dv_14 + dv_155) +
                    dv_2320 + dv_3527);
  sc_18 = d[1250] * sc_29;
  DataVector& dv_3488 = temps.at(2963);
  DataVector& dv_3521 = temps.at(2995);
  DataVector& dv_3522 = temps.at(2996);
  sc_28 = (d[1100] + d[1209]) * dv_3521 +
          (-d[1411]) * (14.0 * dv_3138 + dv_3488 + dv_3522);
  DataVector& dv_149 = temps.at(147);
  DataVector& dv_2733 = temps.at(2286);
  DataVector& dv_30 = temps.at(30);
  DataVector& dv_3182 = temps.at(2669);
  DataVector& dv_3288 = temps.at(2775);
  DataVector& dv_3503 = temps.at(2977);
  sc_28 +=
      (-d[20]) * ((-d[1484]) * dv_45 + (d[9] * d[2924]) * (-dv_149 - dv_30) +
                  dv_3182 * (15.0 - dv_2733) + dv_3288 - 445.0 * dv_3503);
  DataVector& dv_3427 = temps.at(2903);
  DataVector& dv_3430 = temps.at(2906);
  DataVector& dv_3524 = temps.at(2998);
  sc_28 +=
      (-d[436]) * ((-d[2244]) * dv_3524 +
                   (d[2928] * d[2923]) * ((126.0 * d[2924]) * dv_14 - dv_3430) -
                   dv_3427 - dv_3523) +
      sc_18;
  DataVector& dv_139 = temps.at(137);
  DataVector& dv_2884 = temps.at(2433);
  sc_28 += d[1378] * ((-d[407]) * (Dy * dv_2884 + dv_139) + Dy * d[448] +
                      d[88] * (32.0 * dv_2868 + dv_833));
  DataVector& dv_2487 = temps.at(2075);
  DataVector& dv_2705 = temps.at(2258);
  DataVector& dv_3113 = temps.at(2605);
  sc_28 += dv_2115 *
           ((-d[20]) * (d[1748] + d[88] * (17.0 - dv_2487) + 223.0 * dv_1) +
            (65.0 * d[276]) + d[436] * ((-d[2246]) * Dy + d[2245] + dv_2705) +
            dv_3113);
  sc_14 = (-d[362]) * sc_28;
  DataVector& dv_3241 = temps.at(2728);
  DataVector& dv_3298 = temps.at(2785);
  DataVector& dv_703 = temps.at(653);
  sc_17 =
      (-d[50]) * ((-d[1360]) * ((-d[2259]) * dv_292 + dv_3241 + 119.0 * dv_740 +
                                dv_833 * (145.0 - dv_3551)) +
                  365.0 * dv_3298) +
      (d[1339] * d[2250]) * dv_703;
  DataVector& dv_2799 = temps.at(2352);
  sc_17 += d[287] * ((d[2054] + 299.0 * d[7] - 58.0) * dv_96 +
                     Dy * ((d[1880] + 18.0) * Dy + d[314] + 396.0 * dv_794)) +
           d[59] * (Dy * ((-103.0 * d[7]) + dv_2799) + 38.0 * dv_2301);
  DataVector& dv_100 = temps.at(98);
  DataVector& dv_2559 = temps.at(2147);
  DataVector& dv_2931 = temps.at(2464);
  sc_17 += d[724] *
           ((-64.0 * d[147]) * dv_292 + d[1242] * (305.0 * dv_14 + dv_2931) +
            d[80] * (d[2260] * dv_96 + d[2261] * dv_100) - 266.0 * dv_2293 +
            dv_833 * (dv_2559 + 55.0));
  sc_29 = sc_17 * d[2922];
  DataVector& dv_3247 = temps.at(2734);
  DataVector& dv_449 = temps.at(422);
  DataVector& dv_454 = temps.at(427);
  sc_18 = (-d[1490]) * dv_3247 +
          (-d[50]) * (d[153] * (-83.0 * Dx + dv_3529 + dv_3550) +
                      d[9] * (83.0 * dv_3138 + d[2924] * (dv_449 + dv_454)) -
                      270.0 * dv_3287 + 705.0 * dv_3503);
  DataVector& dv_2503 = temps.at(2091);
  DataVector& dv_2726 = temps.at(2279);
  DataVector& dv_2727 = temps.at(2280);
  sc_18 += d[1378] * (Dy * ((-d[77]) * (d[1557] + dv_2503) +
                            d[116] * (d[1751] + dv_2727) + d[2130] + dv_2726) +
                      d[1602] * dv_115);
  DataVector& dv_3546 = temps.at(3019);
  sc_18 += d[1594] * (54.0 * dv_3138 + 38.0 * dv_3262 - dv_3546);
  DataVector& dv_2358 = temps.at(1948);
  DataVector& dv_3547 = temps.at(3020);
  DataVector& dv_3548 = temps.at(3021);
  DataVector& dv_3549 = temps.at(3022);
  sc_18 +=
      d[724] *
      (d[2257] * dv_3549 +
       d[36] * (-Dx * (dv_2358 - 55.0) + 134.0 * dv_1583 + 305.0 * dv_3270) -
       dv_3547 + dv_3548);
  DataVector& dv_3017 = temps.at(2542);
  DataVector& dv_3447 = temps.at(2923);
  sc_18 +=
      dv_2667 *
      ((-d[77]) * ((-d[2258]) * dv_972 + (657.0 * d[36]) + dv_3447) +
       (-103.0 * d[276]) +
       d[20] * ((-d[88]) * (dv_3017 - 21.0) + (256.0 * d[255]) + 471.0 * dv_1) +
       d[49] * (d[1083] + 1701.0 * dv_1));
  DataVector& dv_2135 = temps.at(1751);
  DataVector& dv_2402 = temps.at(1992);
  sc_18 +=
      -dv_2722 * (d[1461] + d[2256] * dv_2387 + d[88] * dv_2402 + dv_2135) +
      sc_29;
  sc_28 = (d[122] * d[2921]) * sc_18;
  sc_17 = (-d[286]);
  DataVector& dv_2642 = temps.at(2216);
  DataVector& dv_2735 = temps.at(2288);
  DataVector& dv_2750 = temps.at(2303);
  DataVector& dv_3018 = temps.at(2543);
  sc_17 *=
      (d[1388] + d[1450] * (d[1595] + d[1802] + 6.0) + d[1603] -
       396.0 * d[277]) *
          dv_29 +
      Dy * ((-d[1443]) * (d[88] + 897.0 * dv_1) +
            (-d[51]) * ((-d[2928]) * dv_3018 + d[1588] + dv_2642) +
            (d[401] * d[2923]) +
            d[631] * ((41.0 - d[2252]) * dv_538 + (254.0 * d[36]) + dv_2735) +
            dv_2750);
  DataVector& dv_163 = temps.at(161);
  DataVector& dv_2720 = temps.at(2273);
  DataVector& sc_25 = temps.at(3247);
  sc_25 = (-d[1472]) *
          (Dx * ((-d[2928]) * (dv_3544 - 3.0) + d[1639] + dv_2720 + dv_3545) +
           d[1240] * dv_163 + dv_3271);
  DataVector& dv_3362 = temps.at(2845);
  DataVector& dv_3520 = temps.at(2994);
  sc_25 += (-d[50]) * (d[37] * (Dx * (dv_2774 - 13.0) + dv_3542) -
                       446.0 * dv_3247 + 445.0 * dv_3294 + dv_3541) +
           (-d[59]) * (37.0 * dv_1112 + dv_3362 - dv_3520);
  sc_25 += (4.0 * d[48] * d[2921]) * Dx *
               ((-d[2256]) * dv_2044 + d[314] + 1701.0 * dv_794) +
           (32.0 * d[151] * d[2266] * d[2923]) * Dx * Dy;
  sc_26 = sc_25 * d[2922];
  DataVector& dv_3228 = temps.at(2715);
  DataVector& dv_3320 = temps.at(2807);
  DataVector& dv_3556 = temps.at(3029);
  sc_29 = (-d[1100]) *
              ((-d[1527]) * dv_3320 + (-d[2091]) * (d[2265] * dv_15 + dv_3556) +
               d[2264] * dv_635 + d[283] * dv_15) -
          320.0 * dv_3228 + sc_17;
  DataVector& dv_1888 = temps.at(1548);
  DataVector& dv_2343 = temps.at(1933);
  DataVector& dv_2717 = temps.at(2270);
  DataVector& dv_3555 = temps.at(3028);
  DataVector& dv_59 = temps.at(59);
  sc_29 += d[1107] * ((12.0 * d[147]) * dv_3555 + d[1106] * (dv_1888 + dv_96) -
                      dv_2343 * dv_59 + dv_2717);
  DataVector& dv_3554 = temps.at(3027);
  sc_29 += d[1184] * (d[1407] * dv_3555 + d[1504] * dv_3554 +
                      dv_1762 * (dv_3535 + 29.0) - 48.0 * dv_2293);
  DataVector& dv_2448 = temps.at(2037);
  DataVector& dv_2747 = temps.at(2300);
  sc_29 += d[1409] * ((-d[266]) * dv_2220 + (40.0 * d[147]) * dv_635 +
                      (d[2928] * d[2925]) * dv_3554 - dv_2747 * dv_2880) +
           d[1594] * (Dy * ((-d[1368]) + dv_2448) + 20.0 * dv_2301);
  DataVector& dv_2152 = temps.at(1765);
  DataVector& dv_2385 = temps.at(1975);
  sc_29 += dv_2942 * ((-d[436]) * (d[1540] + dv_2152) + (21.0 * d[20]) * Dy +
                      (-d[1541]) - dv_2385) +
           sc_26;
  sc_18 = (d[19] * d[57]) * sc_29;
  DataVector& sc_23 = temps.at(3245);
  sc_23 = d[631];
  DataVector& dv_2390 = temps.at(1980);
  DataVector& dv_282 = temps.at(274);
  sc_23 *= (-96.0 * d[147]) * dv_282 +
           (2.0 * d[2923]) * (d[2260] * dv_163 + d[2261] * dv_139) +
           (d[2928] * d[2925]) * (134.0 * dv_14 + 305.0 * dv_15) -
           657.0 * dv_2293 - dv_833 * (dv_2390 - 55.0);
  sc_25 = (-d[1443]) *
          ((d[1930] - 381.0 * d[7] + 52.0) * dv_29 -
           Dy * ((155.0 * d[2927] - 104.0) * Dy + d[1751] + 1529.0 * dv_794));
  DataVector& dv_2510 = temps.at(2098);
  DataVector& dv_3318 = temps.at(2805);
  sc_25 += (-d[51]) * ((-d[2923]) * ((-d[2259]) * dv_282 + dv_3538 -
                                     dv_833 * (dv_2796 - 147.0)) +
                       d[2018] * dv_635) +
           (40.0 * d[2163]) * dv_3318 +
           d[57] * (Dy * ((-d[2263]) + dv_2510) + dv_3526);
  sc_25 += sc_23;
  sc_17 = (2.0 * d[2922]) * sc_25;
  DataVector& dv_480 = temps.at(448);
  sc_26 = (-d[50]) * (d[153] * (-65.0 * Dx + dv_3542 + dv_703 * d[2924]) +
                      d[9] * (65.0 * dv_3138 + d[2924] * (dv_115 + dv_480)) -
                      406.0 * dv_3287 + 495.0 * dv_3503);
  DataVector& dv_2537 = temps.at(2125);
  sc_26 += (-d[724]) * ((-d[2258]) * dv_3524 +
                        (-d[36]) * (Dx * (dv_2537 + 55.0) + 305.0 * dv_1583 +
                                    134.0 * dv_3270) +
                        d[1557] * dv_2379 + dv_3547) +
           (-d[1492] * (d[152] + d[1802])) * dv_45 + sc_17;
  sc_26 += (2.0 * d[142] * d[2921]) * (Dy * ((-d[436]) * (d[1540] + dv_2245) +
                                             (32.0 * d[20]) * (Dy + d[31]) +
                                             (-7.0 * d[50]) + 184.0 * dv_1821) +
                                       d[1615] * dv_691);
  DataVector& dv_3237 = temps.at(2724);
  DataVector& dv_3506 = temps.at(2980);
  DataVector& dv_3552 = temps.at(3025);
  DataVector& dv_3553 = temps.at(3026);
  sc_26 += (2.0 * d[57] * d[2923]) * (58.0 * dv_3138 + dv_3522 - dv_3552) +
           (4.0 * d[48] * d[2921]) *
               ((-d[9]) * (dv_3506 + dv_3553) + d[2262] * dv_3237 +
                32.0 * dv_2435 + 235.0 * dv_3287);
  DataVector& dv_2749 = temps.at(2302);
  DataVector& dv_2754 = temps.at(2307);
  DataVector& dv_3464 = temps.at(2939);
  DataVector& dv_3499 = temps.at(2973);
  sc_26 +=
      -dv_2115 *
      ((-d[287]) * (d[9] + 1215.0 * dv_1) +
       (-d[50]) * (d[1748] + d[88] * dv_2754 + dv_3499) + (d[834] * d[2923]) +
       d[724] * ((-d[2257]) * dv_240 + (133.0 * d[36]) + dv_3464) +
       256.0 * dv_2749);
  sc_29 = (d[50] * d[52]) * sc_26;
  DataVector& dv_2992 = temps.at(2519);
  DataVector& dv_3422 = temps.at(2898);
  DataVector& dv_3517 = temps.at(2991);
  DataVector& dv_3518 = temps.at(2992);
  sc_17 = (-d[286]) * (d[31] * (dv_3422 + dv_3518) + d[497] * dv_1 +
                       d[2921] * (dv_1762 - dv_2937 + dv_2992)) +
          (-d[37]) * (dv_2692 + dv_3517 * d[2923]);
  DataVector& dv_2938 = temps.at(2470);
  DataVector& dv_3462 = temps.at(2937);
  sc_17 += (-d[2921]) *
           (d[1247] * dv_3517 + d[2240] * dv_740 + dv_2717 - dv_2938 + dv_3462);
  DataVector& dv_231 = temps.at(228);
  sc_17 += d[1250] *
           (d[20] * (-15.0 * dv_1112 + 20.0 * dv_3262 + dv_3520) +
            dv_231 * ((-d[2243]) * Dy + (-d[31]) * (32.0 * dv_1502 - 17.0) +
                      (25.0 * d[7]) * Dy) +
            dv_3515 * dv_3519);
  DataVector& dv_1989 = temps.at(102);
  DataVector& dv_3328 = temps.at(2814);
  sc_17 += d[1412] * (Dy * (d[284] + dv_2448) - dv_3350) +
           dv_3328 * ((d[2078] + d[444]) + dv_1989);
  sc_26 = d[1466] * sc_17;
  DataVector& dv_3323 = temps.at(1417);
  sc_10 = (2.0 * d[1931]) * dv_3323 + sc_12 + sc_13 + sc_14 + sc_16 + sc_18 +
          sc_24 + sc_26 + sc_28 + sc_29 + sc_6;
  sc_2 = (d[2272] * d[465]) * sc_10;
  DataVector& dv_2476 = temps.at(2065);
  DataVector& dv_3299 = temps.at(2786);
  DataVector& dv_351 = temps.at(334);
  DataVector& dv_402 = temps.at(382);
  sc_14 = d[1549] * ((-5.0 * d[2076]) * dv_740 +
                     (12.0 * d[2928] * d[2925]) * (dv_351 + dv_402) +
                     (18.0 * d[2928]) * Dy * (dv_2476 + 3.0) - dv_3299) +
          240.0 * dv_3297;
  DataVector& dv_1229 = temps.at(1124);
  DataVector& dv_3233 = temps.at(2720);
  sc_14 +=
      d[355] * ((-d[2923]) * (d[2067] * dv_15 + d[2068] * dv_14) + dv_3233) +
      d[57] * (287.0 * dv_1229 + 528.0 * dv_3298);
  DataVector& dv_2404 = temps.at(1994);
  sc_14 += d[996] * ((d[2000] * d[2075] + 360.0 * d[7] - 27.0) * dv_29 +
                     Dy * ((d[1217] - 27.0) * dv_240 + d[1498] + dv_2404));
  sc_28 = (-d[2922]) * sc_14;
  DataVector& dv_3139 = temps.at(2628);
  DataVector& dv_3260 = temps.at(2747);
  DataVector& dv_3293 = temps.at(2780);
  DataVector& dv_645 = temps.at(595);
  sc_6 = d[287] * (d[1240] * (dv_100 + dv_645) + d[2073] * dv_732 + dv_3139 +
                   432.0 * dv_3287) +
         d[50] * ((-73.0 * d[2928]) * dv_3262 + (141.0 * d[2928] * d[7]) * Dx -
                  141.0 * dv_3260 - dv_3293);
  DataVector& dv_3263 = temps.at(2750);
  DataVector& dv_649 = temps.at(599);
  sc_6 += d[724] * ((-d[2074]) * dv_3247 +
                    d[484] * (dv_149 * d[2924] + dv_3296 * (dv_3295 + 3.0) +
                              dv_649 * d[2924]) -
                    144.0 * dv_3263 + 132.0 * dv_3294);
  sc_6 += Dx * d[807] *
          ((d[1209] * d[2000] - 82.0 * d[7]) * Dy + d[1341] * dv_1502);
  sc_14 = (-d[2921]) * sc_6;
  DataVector& dv_1709 = temps.at(1534);
  DataVector& dv_3292 = temps.at(2779);
  sc_18 = (d[142] * d[1926]) * (Dy * ((-d[1626]) - dv_1709) - dv_3292) + sc_14 +
          sc_28;
  DataVector& dv_2159 = temps.at(1771);
  sc_18 += d[638] * dv_231 *
           (d[1480] + d[20] * ((714.0 * d[36]) + d[2070] * dv_2159) +
            d[77] * ((-d[2071]) * dv_1 + (66.0 * d[2928])) +
            d[91] * ((d[1217] - 15.0) * Dy + d[1341]));
  sc_29 = (-d[122]) * sc_18;
  DataVector& dv_2431 = temps.at(2021);
  DataVector& dv_2650 = temps.at(2152);
  DataVector& dv_3098 = temps.at(2591);
  sc_14 =
      d[2928] * (d[1106] * ((7.0 * d[2923]) * dv_635 - dv_2650) +
                 d[2921] * (21.0 * dv_2301 + dv_240 * ((-d[1814]) + dv_2431))) +
      d[2050] * dv_3098;
  DataVector& dv_3068 = temps.at(2561);
  DataVector& dv_3254 = temps.at(2741);
  sc_14 += dv_0 *
           ((-d[354]) * ((17.0 * d[36]) + 66.0 * dv_3068) + d[2928] * dv_3254);
  sc_18 = (-d[1466]) * sc_14;
  DataVector& dv_3268 = temps.at(2755);
  DataVector& dv_3273 = temps.at(2760);
  sc_12 = d[1472] * (Dx * ((-126.0 * d[255]) - dv_3273) + d[1240] * dv_143 +
                     78.0 * dv_3271) +
          d[2051] * dv_3268 +
          dv_2722 * (Dy * d[2065] + d[304] + 207.0 * dv_794);
  DataVector& dv_3269 = temps.at(2756);
  DataVector& dv_850 = temps.at(789);
  sc_12 += dv_850 * (d[2064] + dv_3269);
  sc_6 = (-d[2922]) * sc_12;
  DataVector& dv_3275 = temps.at(2762);
  DataVector& dv_3276 = temps.at(2763);
  DataVector& dv_3278 = temps.at(2765);
  DataVector& dv_3279 = temps.at(2766);
  DataVector& dv_3280 = temps.at(2767);
  DataVector& dv_3281 = temps.at(2768);
  sc_24 = (-d[1256]) * (d[1193] * dv_3276 + dv_1762 * dv_3275 + dv_3278) +
          d[2066] * dv_3274 +
          d[436] * ((-6.0 * d[147]) * dv_3276 + d[1242] * dv_3274 - dv_3279 +
                    dv_3280 + dv_3281);
  sc_24 += d[50] * (Dy * ((-136.0 * d[7]) + 99.0 * dv_1502) + 68.0 * dv_2301);
  sc_12 = d[2928] * sc_24;
  DataVector& dv_1731 = temps.at(1547);
  sc_28 = (-88.0 * d[1792]) * ((d[2928] * d[99] * d[2921]) * dv_635 +
                               d[252] * dv_1731 + d[348] * dv_740) +
          sc_12 + sc_6;
  DataVector& dv_2641 = temps.at(1754);
  DataVector& dv_503 = temps.at(464);
  sc_28 +=
      d[1882] * (d[27] * (Dy * (Dy * d[2061] + d[2062]) + d[2061] * dv_14) +
                 d[9] * ((-d[1360]) * (d[2060] * dv_14 + dv_503) + dv_2641) -
                 99.0 * dv_613);
  DataVector& dv_1637 = temps.at(1482);
  DataVector& dv_2349 = temps.at(1939);
  sc_28 += d[162] * dv_2349 * ((26.0 * d[2921]) - dv_1637);
  sc_14 = (-d[299]) * sc_28;
  sc_12 = (d[1085] * d[1470]) * dv_45 + d[1469] * dv_45;
  DataVector& dv_2430 = temps.at(2020);
  DataVector& dv_297 = temps.at(289);
  DataVector& dv_3256 = temps.at(2743);
  DataVector& dv_3257 = temps.at(2744);
  DataVector& dv_957 = temps.at(864);
  sc_12 +=
      d[6] * ((-d[1472]) * (d[2053] * dv_3257 + dv_833) + (-d[1496]) * dv_15 +
              (3.0 * d[50]) * (Dy * dv_2430 + dv_297) +
              (4.0 * d[48] * d[2921]) * dv_3256 - dv_957);
  DataVector& dv_3255 = temps.at(2742);
  sc_12 += dv_2350 * ((-d[2051]) * dv_2617 +
                      (-d[436]) * ((-d[2052]) * Dy + d[31] + d[7] * dv_3255) +
                      (-d[1569]) + d[20] * dv_3254);
  sc_28 = (-d[2928] * d[120]) * sc_12;
  DataVector& dv_2512 = temps.at(2100);
  sc_16 = (d[1217] + 75.0 * d[7] - 8.0) * dv_3154 +
          (-d[1443]) * ((-d[1242]) * dv_669 + d[255] * dv_2512 +
                        d[2923] * (d[2067] * dv_14 + d[2068] * dv_15)) +
          d[451] * dv_1;
  DataVector& dv_3283 = temps.at(2770);
  sc_16 +=
      d[50] * ((-d[2069] - 10.0) * dv_2937 + (-6.0 * d[2928]) * Dy * dv_3283 +
               (24.0 * d[2928] * d[2925]) * dv_3257 - 156.0 * dv_2293);
  DataVector& dv_3140 = temps.at(2629);
  sc_16 +=
      d[631] * ((d[1345] - d[2056] - 9.0) * dv_29 +
                Dy * ((d[2054] - 9.0) * dv_240 + d[2060] * dv_3140 + d[304]));
  sc_24 = d[94] * sc_16;
  DataVector& dv_3286 = temps.at(2773);
  DataVector& dv_3289 = temps.at(2776);
  sc_13 = (-d[287]) * ((-d[2065]) * dv_732 + d[1240] * dv_3256 - dv_3286 -
                       dv_3288 + dv_3289);
  DataVector& dv_3290 = temps.at(2777);
  DataVector& dv_3291 = temps.at(2778);
  sc_13 += (-d[724]) *
           ((-d[2056]) * dv_45 +
            (-d[484]) * (Dx * dv_3291 + dv_594 * d[2924] + dv_99 * d[2924]) +
            (4.0 * d[2055] * d[7]) * Dx * Dy - dv_3290);
  DataVector& dv_2497 = temps.at(2085);
  DataVector& dv_3261 = temps.at(2748);
  DataVector& dv_3284 = temps.at(2771);
  DataVector& dv_3285 = temps.at(2772);
  sc_13 += d[50] * ((-d[1221]) * dv_3262 + (26.0 * d[2928] * d[7]) * Dx -
                    dv_3261 - dv_3285) +
           dv_3284 * ((30.0 * d[2923]) * (dv_2497 + dv_751) +
                      (-d[2000] * d[532]) * Dx);
  sc_16 = sc_13 * d[2921];
  DataVector& dv_2450 = temps.at(2039);
  sc_6 = (d[142] * d[363]) * (Dy * dv_2450 - dv_29) + sc_24;
  DataVector& dv_1823 = temps.at(1609);
  DataVector& dv_3282 = temps.at(2769);
  sc_6 += dv_3065 * ((d[2000] - 82.0) * dv_1823 + d[1471] +
                     d[50] * ((-d[494]) + Dy * d[2057]) +
                     d[631] * (d[9] + 207.0 * dv_1) + dv_3282) +
          sc_16;
  sc_12 = d[123] * sc_6;
  DataVector& dv_2432 = temps.at(2022);
  sc_25 = (-d[1172]) * ((36.0 - 143.0 * d[2927]) * Dy + (-d[2082]) * dv_794 +
                        (63.0 * d[36])) +
          (-d[1549]) * ((-d[166]) * (dv_2432 + 9.0) + (714.0 * d[255]) +
                        d[2076] * dv_2387);
  DataVector& dv_2470 = temps.at(2059);
  sc_25 += (-d[1669]) * ((517.0 * d[2927] - 26.0) * dv_1 +
                         (-d[1]) * (dv_2470 + 8.0) + d[1400]) +
           (720.0 * d[384] * d[7]) * Dy + d[57] * (d[2064] + 264.0 * dv_3068);
  sc_17 = Dy * sc_25;
  sc_13 = (d[1] * (d[2091] * (d[2090] - 12.0) +
                   d[287] * (d[2092] + d[266] * d[2925]) +
                   d[50] * (78.0 * d[1242] - 175.0 * d[2923]) +
                   d[807] * (d[2089] - 16.0)) +
           d[2000] * d[2088]) *
              dv_29 +
          1152.0 * dv_3315 + sc_17;
  sc_24 = sc_13 * d[2922];
  DataVector& dv_3313 = temps.at(2800);
  DataVector& dv_3314 = temps.at(2801);
  DataVector& dv_608 = temps.at(559);
  sc_16 = (-d[996]) * ((-d[2083]) * dv_732 + d[1240] * (dv_3314 + dv_608) -
                       64.0 * dv_3260 - 765.0 * dv_3287 + dv_3313) +
          (d[1660] * d[2085]) * dv_732;
  DataVector& dv_600 = temps.at(551);
  sc_16 += (-d[2041] * d[258]) * (Dy * ((-d[303]) - dv_600) - dv_482);
  DataVector& dv_129 = temps.at(127);
  DataVector& dv_656 = temps.at(606);
  sc_16 += d[1549] * ((-d[1580] * d[2070]) * dv_45 +
                      d[403] * (dv_129 * d[2924] + dv_3296 * (dv_2432 + 3.0) +
                                dv_656 * d[2924]) -
                      306.0 * dv_3263 + 264.0 * dv_3294);
  DataVector& dv_343 = temps.at(326);
  DataVector& dv_729 = temps.at(678);
  sc_16 += d[355] * ((-d[2086]) * dv_3247 +
                     d[37] * (dv_343 * d[2924] + dv_352 * d[2924] + dv_729) +
                     462.0 * dv_3294);
  sc_16 += d[57] * ((-d[1587]) * dv_3262 + (99.0 * d[2928] * d[7]) * Dx -
                    99.0 * dv_3260 - dv_3293) +
           sc_24;
  DataVector& dv_3259 = temps.at(2746);
  sc_16 +=
      -dv_3259 * ((-d[724]) * (d[2084] * dv_69 + d[42]) +
                  (d[48] * d[2921]) * (Dy * d[2086] + d[1752]) + (-d[451]) +
                  d[50] * ((300.0 * d[36]) + Dy * d[2074]) - dv_3282);
  sc_6 = d[124] * sc_16;
  DataVector& dv_3309 = temps.at(2796);
  DataVector& dv_3310 = temps.at(2797);
  sc_17 =
      (-d[57]) * (Dy * ((-146.0 * d[7]) + 141.0 * dv_1502) + 73.0 * dv_2301) +
      d[1494] * (d[290] * dv_2220 + dv_3310 * d[2923]) + dv_3309;
  DataVector& dv_2449 = temps.at(2038);
  sc_17 += d[276] *
           ((48.0 * d[2928]) * Dy * dv_2449 +
            (600.0 * d[2928] * d[2925]) * dv_635 - dv_3299 - 365.0 * dv_740);
  DataVector& dv_3311 = temps.at(2798);
  sc_17 += d[724] * (d[1242] * dv_3310 + d[1356] * dv_2220 + dv_3311 +
                     300.0 * dv_3312 - 72.0 * dv_740);
  sc_13 = d[2928] * sc_17;
  DataVector& dv_2781 = temps.at(2334);
  sc_25 =
      (-d[1472]) * ((-d[2923]) * ((d[2081] + 345.0) * dv_15 + d[2082] * dv_14) +
                    dv_2781) +
      (-d[287]) *
          ((d[1068] - 15.0) * dv_128 + Dy * ((d[2056] - 25.0) * Dy + d[1752]));
  sc_25 +=
      (-d[50]) * (Dy * ((1416.0 * d[36]) + Dy * d[2080]) + d[2080] * dv_14) +
      (141.0 * d[57]) * Dy + (1200.0 * d[151] * d[2923]) * dv_14;
  sc_17 = d[1053] * sc_25;
  DataVector& dv_112 = temps.at(110);
  DataVector& dv_2474 = temps.at(2063);
  DataVector& dv_3304 = temps.at(2791);
  sc_23 = d[1549] * (Dx * ((-902.0 * d[2927] - 141.0) * dv_2387 +
                           (24.0 * d[2928]) * dv_2474 + (-708.0 * d[255])) +
                     300.0 * dv_3271 + 300.0 * dv_3308) +
          dv_112 * ((287.0 * d[36]) + dv_3269) + 960.0 * dv_3304;
  DataVector& dv_2471 = temps.at(2060);
  DataVector& dv_2524 = temps.at(2112);
  DataVector& dv_3305 = temps.at(2792);
  DataVector& dv_3306 = temps.at(2793);
  DataVector& dv_3307 = temps.at(2794);
  sc_23 += -dv_2524 * ((-d[88]) * dv_2471 + dv_3305) +
           dv_3307 * ((-126.0 * d[36]) + Dy * d[2083] + d[2084] * dv_3306);
  sc_25 = sc_23 * d[2922];
  DataVector& dv_2027 = temps.at(90);
  DataVector& dv_3149 = temps.at(2638);
  DataVector& dv_3300 = temps.at(2787);
  DataVector& dv_3301 = temps.at(2788);
  DataVector& dv_3303 = temps.at(2790);
  sc_24 = (-60.0 * d[142]) * dv_3300 * ((-d[2078]) - dv_2027) +
          d[2056] * ((-d[1549] * d[2079]) * dv_635 + d[1554] * dv_3301 +
                     d[384] * dv_3149 + d[683] * dv_3303 + d[780] * dv_740) +
          sc_13 + sc_17 + sc_25;
  sc_16 = d[127] * sc_24;
  DataVector& dv_2626 = temps.at(2206);
  sc_13 = (-d[1339]) * dv_2626 +
          (-d[57]) * (Dy * ((-52.0 * d[7]) + 89.0 * dv_1502) + 26.0 * dv_2301);
  DataVector& dv_2261 = temps.at(1857);
  DataVector& dv_2765 = temps.at(2318);
  DataVector& dv_3321 = temps.at(2808);
  sc_13 += d[1411] * ((12.0 * d[2928]) * Dy * dv_3275 +
                      (24.0 * d[2928] * d[2925]) * dv_3321 - 312.0 * dv_2293 -
                      91.0 * dv_740) +
           d[630] * (d[1338] * dv_3320 + dv_2261 * (dv_2765 + 4.0));
  DataVector& dv_3322 = temps.at(2809);
  sc_13 += d[724] * ((-d[246]) * dv_2220 + (36.0 * d[147]) * dv_3321 +
                     (d[2928] * d[2925]) * dv_3320 - dv_3322 - 54.0 * dv_740);
  sc_17 = d[2928] * sc_13;
  DataVector& dv_294 = temps.at(286);
  DataVector& dv_3319 = temps.at(2806);
  sc_23 = d[1409] * (Dy * ((-d[2062]) + Dy * d[2094]) + d[2094] * dv_14) +
          d[287] * ((d[2000] - 54.0) * dv_14 +
                    Dy * ((d[2000] + 93.0) * Dy + (-d[1752]))) +
          89.0 * dv_294 + dv_3319;
  DataVector& dv_2335 = temps.at(1926);
  sc_23 += d[724] * (d[1360] * ((d[2059] + 180.0) * dv_15 + dv_2048) + dv_2335);
  sc_13 = d[1053] * sc_23;
  DataVector& sc_9 = temps.at(3231);
  sc_9 = (-d[1653]) *
             (Dx * (d[2095] + dv_3273) + d[1240] * dv_96 - 69.0 * dv_3308) +
         (d[1660] * (d[1366] + d[2085])) * dv_45 +
         dv_3307 * (Dy * d[2073] + d[1498] + d[2071] * dv_794);
  DataVector& dv_2331 = temps.at(1922);
  sc_9 += Dx * d[230] * (Dy * d[1307] + d[1318]) -
          dv_2524 * (-120.0 * dv_2331 + dv_3305);
  sc_23 = sc_9 * d[2922];
  DataVector& dv_2837 = temps.at(2389);
  DataVector& dv_3133 = temps.at(2624);
  DataVector& dv_3150 = temps.at(2639);
  DataVector& dv_3316 = temps.at(2803);
  DataVector& dv_3317 = temps.at(2804);
  sc_25 = d[162] * dv_2837 * (dv_3133 + d[2921]) +
          d[2056] * ((-d[683]) * (d[269] * dv_15 + dv_3317) +
                     (-d[1536] * d[2093]) * dv_635 + d[1554] * dv_3316 +
                     d[223] * dv_740 + d[384] * dv_3150) +
          sc_13 + sc_17 + sc_23;
  sc_24 = d[128] * sc_25;
  sc_17 = (66.0 * d[1792] * d[2058]) * dv_635;
  DataVector& dv_2547 = temps.at(2135);
  DataVector& dv_2843 = temps.at(2395);
  sc_17 += d[2928] *
           (d[20] * dv_2843 + d[314] * (d[2053] * dv_3266 + dv_833) +
            d[354] * ((-d[255]) * dv_2547 + (3.0 * d[2928]) * Dy * dv_3267 +
                      (12.0 * d[2928] * d[2925]) * dv_3266 - 55.0 * dv_740));
  sc_13 = sc_17 * d[2922];
  sc_23 = (-d[162]) * dv_2606 +
          d[20] * ((-d[1099]) * dv_3262 + (89.0 * d[2928] * d[7]) * Dx -
                   89.0 * dv_3260 - dv_3261);
  DataVector& dv_527 = temps.at(482);
  sc_23 +=
      d[77] * (d[1307] * dv_45 + d[2057] * dv_3247 +
               d[403] * (Dx * dv_3265 + dv_123 * d[2924] + dv_527 * d[2924]) +
               dv_3264) +
      dv_3258 * ((-d[2052]) * dv_1 - dv_2363);
  sc_23 += -dv_3259 * ((2.0 * d[2921]) * ((39.0 * d[36]) + d[2055] * dv_240) +
                       (-d[388]) - 102.0 * dv_1229) +
           sc_13;
  sc_25 = d[362] * sc_23;
  DataVector& dv_265 = temps.at(258);
  sc_26 = d[1864] * (dv_1112 + dv_265 - dv_3138) + sc_12 + sc_14 + sc_16 +
          sc_18 + sc_24 + sc_25 + sc_28 + sc_29 + sc_6;
  sc_10 = (d[25] * d[807]) * sc_26;
  DataVector& dv_148 = temps.at(146);
  DataVector& dv_3400 = temps.at(2877);
  DataVector& dv_3575 = temps.at(3048);
  sc_12 =
      (-d[1319]) * (dv_148 + dv_3400) +
      (-d[58]) * ((-d[2110] + 64.0 * d[2927] - 75.0) * dv_14 -
                  Dy * ((125.0 - d[2036]) * Dy + (374.0 * d[36]) + dv_3575));
  sc_12 += (-d[631]) * (d[1338] * (251.0 * dv_14 + 287.0 * dv_15) + dv_59) +
           (6.0 * d[57]) * ((-d[1717]) * dv_538 + dv_3257 * d[2925]);
  sc_12 += (4.0 * d[48] * d[2921]) *
           (d[2193] * dv_14 + dv_1709 * ((-d[2290]) * Dy + d[484]));
  sc_6 = (-d[638]) * sc_12;
  DataVector& dv_2862 = temps.at(2414);
  DataVector& dv_3014 = temps.at(2539);
  DataVector& dv_3369 = temps.at(2851);
  DataVector& dv_3594 = temps.at(3067);
  DataVector& dv_3595 = temps.at(3068);
  DataVector& dv_694 = temps.at(644);
  sc_28 = (-d[273]) * ((-d[1344]) * dv_3595 + (d[9] * d[2925]) * dv_258 +
                       d[1106] * (dv_3014 + dv_694) + dv_2963 + dv_3369) +
          (-d[57]) * (d[2280] * dv_3594 + dv_2862 * dv_558);
  DataVector& dv_3302 = temps.at(2789);
  DataVector& dv_3588 = temps.at(3061);
  sc_28 += d[1153] * (-dv_3302 + dv_3317) +
           d[276] * ((-d[2198]) * dv_3594 + (d[2928] * d[2925]) * dv_3595 +
                     Dy * d[2928] * (dv_3588 + 103.0) - 115.0 * dv_2293);
  DataVector& dv_446 = temps.at(419);
  sc_28 +=
      d[630] * ((-d[1338]) * dv_258 + d[1242] * dv_446 + dv_2293 + dv_3436);
  sc_12 = d[1] * sc_28;
  DataVector& dv_3577 = temps.at(3050);
  sc_14 = (d[1902] + d[2310]) * dv_3577 +
          (-d[1669]) * (Dx * (d[2928] * dv_2754 + d[1406] + d[2316] * dv_1631) +
                        dv_3271 - dv_3543);
  DataVector& dv_2853 = temps.at(2405);
  DataVector& dv_3584 = temps.at(3057);
  DataVector& dv_3585 = temps.at(3058);
  DataVector& dv_671 = temps.at(621);
  sc_14 += d[1549] * ((-d[1240]) * dv_671 +
                      Dx * ((-1139.0 * d[255]) + d[2928] * (dv_3584 + 47.0) -
                            dv_2853 - dv_3585) +
                      110.0 * dv_3308);
  DataVector& dv_2170 = temps.at(1778);
  DataVector& dv_2968 = temps.at(2496);
  DataVector& dv_3569 = temps.at(3042);
  DataVector& dv_3593 = temps.at(3066);
  sc_14 += d[198] * ((-d[36]) * (dv_2170 * (dv_2968 - 20.0) + dv_3593) +
                     (9.0 * d[2928] * d[147]) * Dx - dv_3569) +
           dv_3300 * ((-d[2304]) * dv_240 + d[2302] + 13269.0 * dv_794);
  sc_28 = sc_14 * d[2922];
  sc_16 = sc_12 + sc_6;
  DataVector& dv_1655 = temps.at(1499);
  DataVector& dv_3581 = temps.at(3054);
  sc_16 +=
      d[1035] * (d[1669] * ((-d[2298]) * dv_15 + d[2291] * dv_14) +
                 d[2293] * (d[2318] * dv_14 + d[2319] * dv_15) +
                 d[2317] * dv_446 + d[762] * (dv_1655 + dv_3292) - dv_3581);
  DataVector& dv_2429 = temps.at(2019);
  DataVector& dv_2759 = temps.at(2312);
  DataVector& dv_2831 = temps.at(2383);
  DataVector& dv_793 = temps.at(741);
  sc_16 += dv_2759 * (d[1300] * ((-d[1540]) - dv_2429) + d[151] * dv_793 +
                      d[2311] + 1125.0 * dv_2831) +
           sc_28;
  sc_24 = (-d[1852]) * sc_16;
  DataVector& dv_2778 = temps.at(2331);
  sc_18 = (d[1048] + d[1593] + 3.0) * dv_2778 +
          (-d[1669]) * (d[2928] * (3.0 - dv_2358) + d[1406] + d[2316] * dv_69);
  DataVector& dv_2776 = temps.at(2329);
  DataVector& dv_2830 = temps.at(2382);
  DataVector& dv_3561 = temps.at(3034);
  sc_18 += d[308] * ((4.0 * d[2928]) * (dv_2830 + 25.0) +
                     (120.0 * d[2289] * d[2923]) * Dy + (-2267.0 * d[255]) -
                     dv_2776) +
           d[59] * ((-d[37]) * (dv_2968 - 40.0) + (-d[1888]) - dv_3561);
  sc_18 += d[996] * ((-d[1273] - 103.0) * Dy + d[1701] + 1879.0 * dv_794);
  sc_14 = -Dy * sc_18;
  DataVector& dv_1571 = temps.at(1443);
  sc_6 = (d[1] * (d[101] * (d[1242] + d[1516]) + d[1031] * d[223] -
                  d[151] * (d[1580] + 12.0) + d[273] * (47.0 - 568.0 * d[7]) +
                  d[58] * (-d[1258] + d[1906] + d[2281])) +
          d[2314] * d[2927]) *
             dv_29 +
         d[2313] * dv_1571 + sc_14;
  sc_12 = (-d[2922]) * sc_6;
  DataVector& dv_3459 = temps.at(2934);
  DataVector& dv_3591 = temps.at(3064);
  DataVector& dv_3592 = temps.at(3065);
  DataVector& dv_679 = temps.at(629);
  sc_28 =
      (d[1368] + d[2310]) * dv_3591 +
      (-d[355]) * ((-d[1345]) * dv_3459 +
                   d[36] * ((-d[2924]) * dv_679 + dv_2170 * dv_2754 + dv_3270) +
                   d[88] * dv_2379 + dv_3592);
  DataVector& dv_3572 = temps.at(3045);
  sc_28 += (-d[59]) * (d[167] * (-35.0 * Dx + dv_3593) +
                       d[2286] * (dv_14 + dv_154) + 70.0 * dv_3260 + dv_3572) +
           sc_12;
  DataVector& dv_2456 = temps.at(2045);
  DataVector& dv_2501 = temps.at(2089);
  sc_28 += d[265] *
           ((-d[1711]) * dv_691 + Dy * (d[2311] + d[50] * (d[1426] - dv_2501) +
                                        dv_2456 + 440.0 * dv_2831));
  DataVector& dv_2242 = temps.at(1843);
  DataVector& dv_3573 = temps.at(3046);
  sc_28 += d[308] *
           ((-d[2282]) * dv_2379 +
            (2.0 * d[2928] * d[2923]) * (Dx * (251.0 * dv_1503 + 50.0) +
                                         504.0 * dv_1583 + dv_2242 * d[2924]) +
            (4.0 * d[2284] * d[7]) * Dx * Dy - dv_3573);
  DataVector& dv_1835 = temps.at(1620);
  DataVector& dv_281 = temps.at(273);
  sc_28 += d[604] * ((-d[88]) * (dv_281 * d[2924] + dv_3553) +
                     d[2303] * dv_1835 + 1757.0 * dv_3287 + dv_3313);
  DataVector& dv_2467 = temps.at(2056);
  DataVector& dv_2468 = temps.at(2057);
  DataVector& dv_2846 = temps.at(2398);
  DataVector& dv_2859 = temps.at(2411);
  DataVector& dv_3086 = temps.at(2579);
  sc_28 += dv_3086 *
           ((-d[50]) * ((-d[2287]) * dv_2468 + (481.0 * d[36]) + dv_2467) +
            d[1450] * (d[1518] + 2645.0 * dv_1) + d[151] * dv_2846 +
            d[1611] * ((-d[2177]) * dv_240 + (-d[1399])) + d[225] * dv_2859);
  sc_16 = (-d[1863]) * sc_28;
  DataVector& dv_2443 = temps.at(2033);
  DataVector& dv_2552 = temps.at(2140);
  DataVector& dv_2680 = temps.at(2239);
  DataVector& dv_3096 = temps.at(2589);
  DataVector& dv_3115 = temps.at(2606);
  DataVector& dv_3189 = temps.at(2676);
  DataVector& dv_3562 = temps.at(3035);
  DataVector& dv_3564 = temps.at(3037);
  sc_14 = (-d[142]) * dv_833 + (d[2273] * d[535]) * dv_45 +
          d[1746] * (dv_2443 + dv_2680 * d[2922]) +
          d[648] * (Dy * dv_2552 + dv_3562 * d[2925]) +
          dv_3189 * ((247.0 * d[648]) + dv_3096 - dv_3115) + dv_3564;
  DataVector& dv_2753 = temps.at(2306);
  sc_14 += d[2923] *
           ((d[1140] * d[2922] + d[1240]) * dv_14 +
            (-245.0 * d[1240] + 12.0 * d[2922] * (-d[2162] + d[286] + 20.0)) *
                dv_15 +
            dv_2443 * ((-d[1269]) - dv_2753));
  sc_6 = (-d[308]) * sc_14;
  DataVector& dv_3144 = temps.at(2633);
  sc_18 = (2.0 * d[142]) * dv_3321 - 985.0 * dv_3144 +
          d[2922] * ((d[1048] + d[1813]) * dv_14 +
                     Dy * ((256.0 * d[2927] + 53.0) * dv_240 + d[1980] -
                           2453.0 * dv_794));
  DataVector& dv_2376 = temps.at(1966);
  DataVector& dv_2455 = temps.at(2044);
  sc_18 += -Dy * ((-d[2277]) * dv_751 + 24.0 * dv_2376 + dv_2455);
  sc_14 = (-d[604]) * sc_18;
  DataVector& dv_142 = temps.at(140);
  DataVector& dv_2208 = temps.at(1811);
  DataVector& dv_2369 = temps.at(1959);
  DataVector& dv_3566 = temps.at(3039);
  DataVector& dv_3567 = temps.at(3040);
  sc_12 =
      (-d[384]) * dv_2208 * ((-d[2251] - d[2274]) * Dy + dv_3566 + dv_3567) -
      dv_142 * ((-d[1020]) * (d[1686] + dv_2369) + d[1686] * dv_2376 +
                d[1960] * dv_751) +
      sc_14 + sc_6;
  DataVector& dv_3100 = temps.at(2087);
  DataVector& dv_3221 = temps.at(2708);
  sc_12 += (4.0 * d[151] * d[2921]) * Dy *
           (dv_3096 - dv_3221 + d[2922] * (d[2276] * dv_1631 + dv_3100));
  sc_28 = (-d[1872]) * sc_12;
  DataVector& dv_2639 = temps.at(2215);
  sc_14 = (24.0 * d[3]) * dv_2759 +
          d[1053] * (-dv_2639 + d[2923] * ((-d[2923]) * dv_527 - dv_3422));
  DataVector& dv_2820 = temps.at(2372);
  DataVector& dv_2822 = temps.at(2374);
  DataVector& dv_3560 = temps.at(3033);
  sc_14 += -dv_0 * (d[2928] * dv_3560 + d[27] * (d[31] * dv_2822 + dv_3561)) -
           dv_2443 * dv_2820;
  sc_12 = (-d[1931]) * sc_14;
  DataVector& dv_3574 = temps.at(3047);
  sc_29 =
      (-d[273]) * ((d[2230] * d[2287]) * dv_45 +
                   d[31] * (Dx * (50.0 - dv_3574) + 504.0 * dv_3270 + dv_3530) +
                   23.0 * dv_3263 - dv_3573);
  DataVector& dv_2183 = temps.at(1790);
  DataVector& dv_286 = temps.at(278);
  DataVector& dv_3490 = temps.at(2965);
  sc_29 += d[51] * (d[1435] * (dv_1583 - dv_2183 + dv_3490) +
                    d[2286] * (dv_286 + dv_352) + 90.0 * dv_3260 + dv_3572);
  DataVector& dv_2227 = temps.at(1829);
  DataVector& dv_2246 = temps.at(1844);
  DataVector& dv_2300 = temps.at(1893);
  DataVector& dv_2445 = temps.at(2035);
  DataVector& dv_2772 = temps.at(2325);
  DataVector& dv_2844 = temps.at(2396);
  sc_29 += d[807] * dv_751 * (-dv_2227 + dv_2772) +
           dv_2445 * (d[1828] + d[2285] * dv_1631 +
                      d[306] * (dv_2300 + dv_2844) + 69.0 * dv_2246);
  sc_18 = (-d[2921]) * sc_29;
  DataVector& dv_2805 = temps.at(2358);
  sc_23 =
      (-d[2288]) * dv_99 +
      (-d[59]) * (d[37] * (Dy * ((-d[1189] - 40.0) + dv_2805) + 3.0 * dv_3536) +
                  81.0 * dv_3298);
  sc_23 += (-d[996]) * ((-d[2054] - 365.0 * d[7] + 53.0) * dv_29 +
                        Dy * (d[1362] + d[1966] * dv_240 - 187.0 * dv_794));
  sc_23 +=
      (d[2928] * d[50]) * ((-d[1714]) * dv_3266 +
                           d[1343] * (d[2289] * dv_139 + d[2290] * dv_163) +
                           dv_1762 * (251.0 * dv_1502 + 50.0) -
                           481.0 * dv_2293 + 252.0 * dv_3559) +
      (-d[683] * (d[1242] - d[1360] * d[2276])) * dv_99;
  sc_29 = sc_23 * d[2922];
  DataVector& dv_3486 = temps.at(2961);
  sc_6 =
      d[1493] *
          (Dy * ((-d[2921]) * (d[2282] + dv_2503) + (72.0 * d[319]) + dv_2781) +
           32.0 * dv_3486) +
      sc_18;
  DataVector& dv_2838 = temps.at(2390);
  DataVector& dv_3116 = temps.at(2607);
  sc_6 +=
      dv_3116 *
          ((-d[20]) * ((-d[2284]) * dv_1709 + (2267.0 * d[36]) + dv_2467) +
           d[1153] + d[1333] * dv_2838 + d[259] * (d[1356] + 4729.0 * dv_1)) +
      sc_29;
  sc_14 = (-d[362]) * sc_6;
  DataVector& dv_2611 = temps.at(2193);
  DataVector& dv_2823 = temps.at(2375);
  DataVector& dv_345 = temps.at(328);
  DataVector& dv_3482 = temps.at(2957);
  DataVector& dv_3571 = temps.at(3044);
  sc_23 = d[319] * (d[1193] * (dv_15 + dv_345) + d[2281] * dv_3482 + dv_3452 +
                    dv_833 * (6.0 - dv_3571)) +
          d[51] * (d[2280] * dv_3482 + dv_240 * dv_2823) - dv_2611;
  DataVector& dv_3483 = temps.at(2958);
  sc_23 += d[77] * ((-d[1242]) * dv_527 + dv_3483 * d[2923]);
  sc_18 = d[2928] * sc_23;
  sc_13 = d[243] * ((-d[1678]) * Dy + (3.0 * d[2925]) * dv_3266) +
          d[36] * (-2453.0 * dv_14 - dv_30);
  DataVector& dv_3120 = temps.at(2611);
  sc_13 += d[2921] * ((-d[2017] + d[2230] + 105.0) * dv_115 +
                      Dy * ((-d[1076]) * Dy + (505.0 * d[36]) + dv_3120));
  sc_23 = d[1882] * sc_13;
  DataVector& dv_2114 = temps.at(1731);
  DataVector& dv_3501 = temps.at(2975);
  sc_17 =
      (-d[248]) *
      ((-d[36]) * (dv_121 * d[2924] + dv_2170 * (dv_2114 - 20.0) + dv_3501) +
       (3.0 * d[2928] * d[147]) * Dx - dv_3569);
  DataVector& dv_3218 = temps.at(2705);
  DataVector& dv_3457 = temps.at(2932);
  DataVector& dv_3570 = temps.at(3043);
  sc_17 +=
      (-d[259]) * (Dx * ((12.0 * d[2279] * d[2923]) * Dy + (-471.0 * d[255]) +
                         d[2928] * (dv_3570 + 106.0) - dv_3218) +
                   248.0 * dv_3271 + dv_3457);
  sc_17 += Dx * d[48] * (Dy * d[2277] + d[1980] - 985.0 * dv_794);
  sc_13 = d[86] * sc_17;
  DataVector& dv_2272 = temps.at(1868);
  DataVector& dv_2589 = temps.at(2164);
  DataVector& dv_3451 = temps.at(2926);
  DataVector& dv_3568 = temps.at(3041);
  sc_29 = (d[1035] * d[2921]) * (d[396] * ((-d[2278]) * dv_14 + dv_3537) +
                                 d[416] * dv_740 + d[563] * dv_3568) +
          dv_3451 * ((2.0 * d[2921]) * (d[1495] + dv_2589) +
                     (-d[239] * d[2923]) - dv_2272) +
          sc_13 + sc_18 + sc_23;
  sc_6 = d[1466] * sc_29;
  DataVector& dv_2816 = temps.at(2368);
  DataVector& dv_3477 = temps.at(2952);
  sc_18 =
      d[2928] * (d[253] * (Dy * dv_2816 + dv_115 * d[2925]) + d[255] * dv_3477 +
                 d[2921] * (d[1242] * dv_3477 + d[147] * dv_450 +
                            240.0 * dv_2868 + dv_3527));
  sc_18 += d[1798] * ((-d[31]) * dv_3565 + d[135] * dv_635);
  sc_23 = sc_18 * d[2922];
  DataVector& dv_3470 = temps.at(2945);
  sc_13 = d[1534] * (d[2053] * dv_3470 + dv_833) +
          d[259] * (d[2273] * dv_3549 + d[36] * (-245.0 * dv_3270 + dv_3430) +
                    dv_3427 + dv_3564);
  DataVector& dv_2156 = temps.at(1768);
  DataVector& dv_2825 = temps.at(2377);
  DataVector& dv_2954 = temps.at(2482);
  DataVector& dv_3487 = temps.at(2962);
  DataVector& dv_3563 = temps.at(3036);
  sc_13 += dv_3065 * (dv_2825 + dv_3563 +
                      d[2921] * ((-d[1140]) * Dy + (247.0 * d[36]) + dv_2647)) +
           dv_3487 * (d[167] * dv_2954 + d[1960] * dv_1 +
                      d[257] * (dv_2156 + dv_3083)) +
           sc_23;
  DataVector& dv_1083 = temps.at(981);
  sc_13 += (d[1048] + d[1209]) * dv_1 * dv_1083;
  sc_29 = d[1581] * sc_13;
  DataVector& dv_3496 = temps.at(2970);
  sc_17 = (-d[223]) * ((-d[1691]) * Dy + (3.0 * d[2925]) * dv_635) +
          d[287] * (d[2301] * dv_99 + dv_1229) + dv_3496;
  DataVector& dv_2662 = temps.at(768);
  DataVector& dv_2746 = temps.at(2299);
  sc_17 +=
      d[50] *
          ((-d[2210] + 134.0 * d[2927] - 175.0) * dv_96 -
           Dy * ((-d[2161]) * dv_2746 + (1139.0 * d[36]) + 132.0 * dv_794)) +
      d[631] * (dv_2662 + d[2923] * (1879.0 * dv_14 + 568.0 * dv_15));
  sc_18 = (-d[638]) * sc_17;
  DataVector& dv_2849 = temps.at(2401);
  DataVector& dv_3444 = temps.at(2920);
  sc_9 = (-d[1595]) * dv_3444 +
         d[1549] * (Dx * ((-1122.0 * d[255]) + d[2928] * (dv_3584 + 103.0) -
                          dv_2849 - dv_3585) +
                    d[1698] * dv_3270 + d[1704] * dv_1583);
  DataVector& dv_1881 = temps.at(275);
  DataVector& dv_2848 = temps.at(2400);
  DataVector& dv_3583 = temps.at(3056);
  sc_9 +=
      d[780] * ((-d[31]) * (dv_3582 + dv_3583) - dv_3290 - dv_3569) +
      dv_2804 * (d[2928] * (dv_2390 - dv_2848) + d[167] + d[2295] * dv_1881);
  sc_9 += dv_3300 * (d[2302] + d[2303] * dv_1709 + 7935.0 * dv_794);
  sc_17 = (-d[2922]) * sc_9;
  DataVector& dv_2673 = temps.at(2235);
  DataVector& dv_2841 = temps.at(2393);
  DataVector& dv_3586 = temps.at(3059);
  DataVector& dv_3587 = temps.at(3060);
  DataVector& dv_3589 = temps.at(3062);
  DataVector& sc_15 = temps.at(3237);
  sc_15 = d[1209] * dv_3586 + d[230] * (d[2280] * dv_3587 + dv_240 * dv_2841) +
          d[276] * ((75.0 * d[2923]) * dv_3587 +
                    (10.0 * d[2928] * d[2925]) * dv_3589 - dv_2673 -
                    dv_833 * (dv_3588 + 47.0));
  DataVector& dv_291 = temps.at(283);
  DataVector& dv_3352 = temps.at(2836);
  DataVector& dv_3475 = temps.at(2950);
  DataVector& dv_668 = temps.at(618);
  sc_15 += d[630] * (d[1106] * dv_291 + dv_2692 + dv_3475) +
           d[631] * ((10.0 * d[147]) * dv_3589 + d[1193] * dv_291 + dv_2650 +
                     dv_3352 + d[2923] * (47.0 * dv_15 + dv_668));
  sc_9 = d[1] * sc_15;
  sc_23 = d[1035] * ((-d[2293]) * (d[2299] * dv_14 + d[2300] * dv_15) +
                     (d[2298] * d[683]) * dv_691 +
                     d[762] * (67.0 * dv_15 + dv_671) + dv_3580 + dv_3581) +
          sc_17 + sc_18 + sc_9;
  DataVector& dv_1068 = temps.at(967);
  DataVector& dv_2696 = temps.at(2250);
  DataVector& dv_2798 = temps.at(2351);
  sc_23 += dv_2798 * ((2.0 * d[20]) * (d[1698] + dv_2696) + (-d[2297]) -
                      1757.0 * dv_1068);
  sc_13 = d[1841] * sc_23;
  DataVector& dv_2651 = temps.at(758);
  DataVector& dv_3102 = temps.at(2594);
  DataVector& dv_3576 = temps.at(3049);
  sc_18 = (-d[1472]) * (d[1746] * dv_3576 + dv_1762 * dv_3102 + dv_2651 +
                        d[2923] * (-53.0 * dv_15 + dv_96)) +
          (-d[1339] * d[532]) * dv_503;
  DataVector& dv_2855 = temps.at(2407);
  DataVector& dv_2856 = temps.at(2408);
  DataVector& dv_2891 = temps.at(2437);
  sc_18 += d[276] * ((-d[1407]) * dv_3576 +
                     Dy * ((-d[2928]) * dv_3571 + (630.0 * d[2923]) * Dy +
                           (248.0 * d[2928] * d[7] - d[2127]))) +
           dv_142 * (d[284] * dv_2856 + dv_2855 + dv_2891);
  DataVector& dv_2418 = temps.at(2008);
  DataVector& dv_2857 = temps.at(2409);
  DataVector& dv_2858 = temps.at(2410);
  sc_18 += d[1107] * dv_1 * (d[2928] * (dv_2857 + dv_2858) + d[153] + dv_2418);
  sc_17 = d[2928] * sc_18;
  DataVector& dv_1708 = temps.at(116);
  DataVector& dv_3197 = temps.at(2684);
  sc_15 =
      (-d[1107]) * (Dy * ((-d[1086]) - dv_3197) + dv_121) + d[1389] * dv_1708 +
      d[58] * ((d[1035] + d[1686]) * dv_29 +
               Dy * ((45.0 - d[2294]) * dv_240 + (157.0 * d[36]) + dv_3575));
  DataVector& dv_2977 = temps.at(2504);
  sc_15 +=
      Dy * d[225] * ((-d[155]) - dv_2977) +
      d[724] * ((-d[2923]) * (187.0 * dv_14 + 730.0 * dv_15) + 12.0 * dv_833);
  sc_18 = d[1053] * sc_15;
  DataVector& dv_2784 = temps.at(2337);
  DataVector& dv_3062 = temps.at(2555);
  DataVector& sc_30 = temps.at(3252);
  sc_30 = (-d[1680] - d[2274] - 9.0) * dv_3577 +
          (-d[1669]) * (d[2295] * dv_3062 + dv_2435 - dv_2784 * dv_833);
  DataVector& dv_3466 = temps.at(2941);
  DataVector& dv_3579 = temps.at(3052);
  DataVector& dv_582 = temps.at(534);
  sc_30 += d[308] * (Dx * ((-d[2279]) * dv_582 + (505.0 * d[255]) +
                           d[2928] * (6.0 - dv_3570) + dv_3579) +
                     d[1240] * dv_649 + dv_3466);
}
}  // namespace CurvedScalarWave::Worldtube::detail
