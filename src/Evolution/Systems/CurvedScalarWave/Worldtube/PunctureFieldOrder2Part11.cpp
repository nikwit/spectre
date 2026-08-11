
// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/Worldtube/PunctureFieldOrder2Impl.hpp"

namespace CurvedScalarWave::Worldtube::detail {

// NOLINTNEXTLINE(google-readability-function-size, readability-function-size)
void puncture_field_2_part_11(const std::array<double, order2_n_doubles>& d,
                              const DataVector& Dx, const DataVector& Dy,
                              const DataVector& z,
                              DynamicBuffer<DataVector>& temps,
                              const gsl::not_null<Order2Vars*> result) {
  DataVector& dv_0 = temps.at(0);
  DataVector& dv_3781 = temps.at(997);
  DataVector& dv_3953 = temps.at(3144);
  DataVector& sc_30 = temps.at(3252);
  sc_30 += d[1858] * ((42.0 * d[1182] - d[2699] * d[2700] + d[2701]) * Dy +
                      (-d[2308]) * dv_0 + (-d[2336]) * dv_3953) +
           d[2320] * dv_3781;
  DataVector& dv_2544 = temps.at(2132);
  DataVector& dv_3948 = temps.at(1365);
  DataVector& dv_3949 = temps.at(1177);
  DataVector& dv_752 = temps.at(700);
  sc_30 += d[2668] * ((-d[2293]) * ((d[2299] + 67.0 * d[6]) * Dx + dv_3948) +
                      d[1550] * (dv_2544 - 107.0 * dv_752) + dv_3949);
  DataVector& dv_2170 = temps.at(1778);
  DataVector& dv_3944 = temps.at(3112);
  DataVector& dv_3945 = temps.at(3212);
  DataVector& dv_751 = temps.at(699);
  sc_30 +=
      d[2692] * ((-d[1249]) * ((d[2278] + 34.0 * d[6]) * dv_2170 + dv_3945) +
                 d[91] * (64.0 * dv_751 + dv_752) + dv_3944) +
      d[2698] * ((-214.0 * d[1182] - d[1735] * d[42] * d[50] - d[2697] +
                  54.0 * d[57] * ypdot) *
                     Dy +
                 (d[104] + d[2696]) * dv_3953 + d[2314] * dv_0);
  DataVector& sc_28 = temps.at(3250);
  sc_28 = d[422] * sc_30;
  DataVector& dv_2469 = temps.at(2058);
  DataVector& dv_3907 = temps.at(1196);
  DataVector& sc_14 = temps.at(3236);
  sc_14 = (-d[1581]) * dv_2469 +
          (-d[1834]) * ((-d[1818] - d[1875] + 13.0 * d[20] * ypdot) * Dy +
                        d[1470] * dv_0 + dv_3907) +
          (d[2622] * d[361]) * dv_752;
  DataVector& dv_3866 = temps.at(3127);
  DataVector& dv_3906 = temps.at(1146);
  sc_14 += d[1841] * ((-d[20]) * (136.0 * dv_751 + 99.0 * dv_752) +
                      d[1478] * ((d[2254] + 26.0 * d[6]) * Dx + dv_3906) +
                      d[91] * (-dv_3866 - 46.0 * dv_751));
  DataVector& dv_2227 = temps.at(1829);
  DataVector& dv_3904 = temps.at(925);
  DataVector& dv_3905 = temps.at(663);
  sc_14 += d[2623] * ((-yp) * (dv_3904 + dv_3905) + d[648] * dv_2227);
  DataVector& dv_2419 = temps.at(2009);
  DataVector& dv_2737 = temps.at(2290);
  DataVector& dv_3882 = temps.at(3198);
  DataVector& dv_3908 = temps.at(3157);
  sc_14 += d[341] * ((-d[1491]) * (-dv_3906 + dv_3908) +
                     (-d[50]) * (dv_2737 + 89.0 * dv_752) +
                     d[287] * (dv_2419 - 41.0 * dv_752) + 480.0 * dv_3882);
  DataVector& dv_1674 = temps.at(1508);
  DataVector& dv_3884 = temps.at(1042);
  sc_14 += d[342] * ((d[1254] * (d[2267] + 25.0 * d[6]) - d[2270]) * dv_2170 +
                     (d[2625] * xpdot) * Dy) +
           d[343] * ((36.0 * M * d[20] * d[2626] - d[2627]) * Dy +
                     (d[1488] + 78.0 * d[274]) * dv_1674 + 576.0 * dv_3884);
  DataVector& dv_1564 = temps.at(1437);
  DataVector& dv_263 = temps.at(256);
  DataVector& dv_3899 = temps.at(3145);
  DataVector& dv_5 = temps.at(5);
  sc_14 +=
      d[345] *
          ((d[1497] + 408.0 * d[274]) * dv_1564 +
           (-164.0 * d[252] + 36.0 * d[259] * d[2624] - 141.0 * d[319]) * dv_5 +
           612.0 * dv_3884) +
      d[362] * ((d[1476] + 72.0 * d[36]) * dv_3899 +
                (-d[1299] * d[502] - d[1637] - 89.0 * d[319]) * Dy +
                d[558] * dv_263);
  sc_30 = d[443] * sc_14;
  DataVector& dv_1696 = temps.at(1527);
  DataVector& sc_15 = temps.at(3237);
  sc_15 = (-d[123]) * ((d[101] * d[2186] + d[2634]) * Dy +
                       (-d[2099]) * dv_1564 + d[1430] * dv_263) +
          (-d[124]) *
              ((d[101] * d[2635] + d[133] * d[36] + d[1395] + d[1985]) * Dy +
               (d[1710] + d[2104] + 74.0 * d[277]) * dv_1696 +
               (d[1590] + d[62]) * dv_263);
  DataVector& dv_2525 = temps.at(2113);
  DataVector& dv_265 = temps.at(258);
  DataVector& dv_3891 = temps.at(3165);
  DataVector& dv_3912 = temps.at(3189);
  sc_15 +=
      (-d[127]) * ((-d[807]) * dv_3912 +
                   d[101] * ((d[2180] + d[2596]) * dv_2170 + 221.0 * dv_265) +
                   d[273] * (-dv_2525 - dv_3866) + d[346] * dv_3891);
  DataVector& dv_2018 = temps.at(1523);
  DataVector& dv_3851 = temps.at(1198);
  DataVector& dv_3883 = temps.at(891);
  DataVector& dv_3914 = temps.at(3172);
  sc_15 += (-d[128]) *
           (d[101] * ((-d[2633]) * dv_2170 + (221.0 * xpdot * ypdot) * Dy) +
            d[273] * (dv_3883 + 30.0 * dv_751) + d[313] * (dv_2018 + dv_3851) +
            d[346] * dv_3914);
  DataVector& dv_1819 = temps.at(1605);
  DataVector& dv_2124 = temps.at(1740);
  DataVector& dv_3909 = temps.at(1108);
  DataVector& dv_780 = temps.at(728);
  sc_15 += (-d[299]) * ((-d[1421] - d[1426]) * dv_3909 + d[248] * dv_780 +
                        d[2630] * dv_752) +
           (-4.0 * d[1466]) * dv_265 + (d[120] * d[2099] * xpdot) * Dy +
           d[122] * ((d[1984] + d[2631] + d[36] * d[448]) * Dy +
                     (-d[380]) * dv_2124 +
                     (-d[1256] - d[1440] - 73.0 * d[252]) * dv_1819);
  DataVector& dv_1 = temps.at(1);
  DataVector& dv_3896 = temps.at(3191);
  sc_15 += d[362] * (d[2149] * dv_1696 + d[2628] * dv_1 - dv_3896);
  sc_14 = d[550] * sc_15;
  DataVector& dv_3850 = temps.at(1144);
  DataVector& dv_3854 = temps.at(647);
  DataVector& sc_4 = temps.at(3226);
  sc_4 = -dv_3850 - dv_3854;
  DataVector& dv_1530 = temps.at(1407);
  DataVector& dv_1552 = temps.at(1426);
  DataVector& dv_3131 = temps.at(2622);
  DataVector& dv_3846 = temps.at(1625);
  DataVector& dv_3861 = temps.at(3220);
  DataVector& sc_34 = temps.at(3256);
  DataVector& sc_7 = temps.at(3229);
  sc_4 += (-d[2577]) * (xp * (d[2576] * dv_0 + dv_1552 + dv_3846) +
                        yp * ((3.0 * xpdot) * Dy + (3.0 * d[7] * xpdot) * Dy -
                              dv_1530 - dv_3131)) +
          (-d[2594]) * dv_3861 + sc_34 + sc_7;
  DataVector& dv_3845 = temps.at(1622);
  DataVector& dv_3848 = temps.at(1100);
  DataVector& dv_3849 = temps.at(1494);
  DataVector& dv_3856 = temps.at(1649);
  DataVector& dv_3859 = temps.at(3114);
  DataVector& dv_3860 = temps.at(1252);
  DataVector& dv_3877 = temps.at(245);
  DataVector& sc_12 = temps.at(3234);
  DataVector& sc_17 = temps.at(3239);
  DataVector& sc_33 = temps.at(3255);
  DataVector& sc_36 = temps.at(3258);
  sc_4 += (d[427] * d[475] * d[852]) * dv_3860 + d[124] * dv_3845 +
          d[171] * (dv_3848 * yp + dv_3849 * xp) +
          d[260] * ((-d[20]) * dv_3856 + dv_3859) + d[2607] * dv_3877 + sc_12 +
          sc_17 + sc_33 + sc_36;
  DataVector& dv_3898 = temps.at(1044);
  DataVector& dv_3903 = temps.at(1599);
  sc_4 += d[2621] * (d[120] * dv_3898 +
                     d[299] * (d[2106] * dv_752 + d[2617] * dv_2170) - dv_3903);
  DataVector& dv_3871 = temps.at(3160);
  DataVector& dv_3872 = temps.at(3131);
  DataVector& dv_3873 = temps.at(979);
  DataVector& dv_3875 = temps.at(3129);
  DataVector& dv_3876 = temps.at(3207);
  sc_4 += d[331] * (d[121] * (Dx * d[2600] - dv_3876) +
                    d[55] * ((-d[2598]) * dv_2170 + d[305] * dv_752) + dv_3871 -
                    dv_3872 - dv_3873 - dv_3875);
  DataVector& dv_3855 = temps.at(1283);
  DataVector& sc_13 = temps.at(3235);
  DataVector& sc_6 = temps.at(3228);
  sc_4 += d[63] * dv_3855 + sc_13 + sc_14 + sc_28 + sc_30 + sc_6;
  DataVector& dv_1496 = temps.at(1376);
  DataVector& dv_229 = temps.at(226);
  DataVector& sc_11 = temps.at(3233);
  sc_11 = (2.0 * d[0]) * dv_1496 * dv_229 * sc_4;
  DataVector& dv_289 = temps.at(281);
  DataVector& dv_3182 = temps.at(2669);
  DataVector& dv_4 = temps.at(4);
  DataVector& dv_627 = temps.at(578);
  sc_28 = (-d[306]) * dv_289 + (-d[2231] - d[2702]) * dv_3182 +
          (-d[1244] * d[1374]) * dv_1 + (12.0 * d[147] * d[92] * yp) * Dx +
          xpdot * ((-4.0 * d[2586]) * dv_4 + (d[1766] * d[2585]) * Dx +
                   Dy * d[1120] + d[1118] * dv_627);
  DataVector& dv_119 = temps.at(117);
  DataVector& dv_2803 = temps.at(2356);
  DataVector& dv_726 = temps.at(675);
  sc_28 +=
      ypdot * (d[1121] * dv_119 + d[164] * dv_726 + d[795] * dv_119 - dv_2803);
  sc_30 = (-d[168]) * sc_28;
  sc_28 = (-d[237]);
  DataVector& dv_3862 = temps.at(368);
  DataVector& dv_3863 = temps.at(993);
  DataVector& dv_3867 = temps.at(3153);
  DataVector& dv_3868 = temps.at(3118);
  DataVector& dv_3869 = temps.at(3199);
  DataVector& dv_3870 = temps.at(3117);
  sc_28 *=
      (-d[55]) * dv_3862 +
      d[50] * ((-d[20]) * dv_3867 + d[77] * ((-d[2125]) * dv_2170 + dv_3868) +
               d[91] * dv_3863) +
      d[61] * (Dy * d[2034] + d[2597] * dv_0 + dv_3869) + dv_3870;
  DataVector& dv_146 = temps.at(144);
  DataVector& dv_34 = temps.at(34);
  DataVector& dv_3853 = temps.at(3176);
  sc_6 = d[153] * ((-d[133]) * Dx + Dx * d[112] + dv_34 * xp) + dv_3853 * xp +
         xpdot * ((-d[2702]) * dv_5 + (d[2580] * xp) * Dx +
                  (2.0 * M * d[137]) * Dy - 10.0 * dv_146);
  DataVector& dv_851 = temps.at(790);
  sc_6 += ypdot * ((-d[1367]) * dv_119 + (-d[697] - d[72]) * dv_726 +
                   d[1096] * dv_119 + dv_851);
  sc_13 = (-d[2583]) * sc_6;
  DataVector& dv_190 = temps.at(187);
  DataVector& dv_2947 = temps.at(718);
  DataVector& dv_3955 = temps.at(777);
  DataVector& dv_3956 = temps.at(1059);
  DataVector& dv_868 = temps.at(804);
  DataVector& dv_973 = temps.at(879);
  sc_33 =
      (-d[1170]) * ((-d[149]) * dv_190 + (-d[202]) * dv_5 + (-d[811]) * dv_146 +
                    d[120] * dv_2947 + 16.0 * dv_3955 + 40.0 * dv_3956) +
      (-d[290]) * dv_973 + (-d[510]) * dv_868;
  DataVector& dv_386 = temps.at(369);
  DataVector& dv_875 = temps.at(796);
  DataVector& dv_876 = temps.at(118);
  DataVector& dv_902 = temps.at(834);
  sc_33 += (-d[563]) *
               ((-d[194]) * dv_119 + (-d[195]) * dv_876 + (22.0 * d[122]) * Dy +
                (66.0 * d[55] * yp) * Dx - 46.0 * dv_875 - dv_902) +
           (-d[232] * d[273]) * Dx + (-M * d[1292]) * dv_876 +
           (-d[1] * d[188]) * dv_386;
  DataVector& dv_3954 = temps.at(588);
  DataVector& dv_896 = temps.at(828);
  DataVector& dv_897 = temps.at(829);
  DataVector& dv_970 = temps.at(876);
  sc_33 += (-d[112] * d[2723]) * dv_3954 + (-d[120] * d[278]) * dv_119 +
           (-d[1774] * d[3]) * dv_3954 + (2.0 * M * d[187] * d[50] * xp) * Dy +
           (2.0 * M * d[7] * d[92] * yp) *
               (77.0 * dv_876 + 148.0 * dv_896 + 86.0 * dv_897 + dv_970);
  DataVector& dv_153 = temps.at(151);
  DataVector& dv_853 = temps.at(792);
  sc_33 += (2.0 * M * d[92] * xpdot * ypdot) *
               ((d[2077] * (d[162] + d[195])) * dv_4 + Dy * d[202] +
                Dy * d[204] + 76.0 * dv_153 + 52.0 * dv_853) +
           (4.0 * M * d[187] * d[19] * d[20]) * Dx;
  DataVector& dv_613 = temps.at(564);
  DataVector& dv_854 = temps.at(793);
  sc_33 += (2.0 * M * d[6] * d[92] * xp * yp) *
           (d[1226] * dv_4 + 77.0 * dv_613 + 86.0 * dv_854);
  sc_6 = (d[207] * d[48]) * sc_33;
  sc_14 = (-d[1134]) * dv_3861 +
          (-d[171]) * ((-yp) * dv_3848 + (-xp) * dv_3849) - dv_3850 - dv_3854 +
          sc_28 + sc_30;
  DataVector& dv_126 = temps.at(124);
  DataVector& dv_231 = temps.at(228);
  sc_14 +=
      (-d[2577]) * (dv_3846 * xp + xpdot * (d[154] * dv_126 + d[2576] * dv_4) +
                    ypdot * ((-d[1189] - 2.0) * dv_231 + (3.0 * xp) * Dy)) +
      sc_13;
  sc_14 += (-d[260]) * (d[20] * dv_3856 - dv_3859) +
           (-d[2621]) *
               ((-d[120]) * dv_3898 +
                d[299] * ((-d[2617]) * dv_2170 + d[2724] * dv_752) + dv_3903);
  DataVector& dv_2161 = temps.at(1772);
  DataVector& dv_2163 = temps.at(1774);
  sc_14 += (-d[331]) * (d[121] * ((-d[2600]) * Dx + dv_3876) +
                        d[55] * (d[1126] * dv_752 + d[2598] * dv_2170) -
                        dv_3871 + dv_3872 + dv_3873 + dv_3875) +
           (-xp) * dv_2161 + (-xp) * dv_2163;
  DataVector& dv_1906 = temps.at(753);
  DataVector& dv_1907 = temps.at(710);
  DataVector& dv_1908 = temps.at(714);
  DataVector& dv_2165 = temps.at(1776);
  DataVector& dv_6 = temps.at(6);
  DataVector& dv_977 = temps.at(831);
  sc_14 +=
      (-xp) * dv_2165 +
      (-d[197] * (d[383] + d[48] * (d[102] + d[238])) +
       d[2436] * (-d[1172] + 2.0 * d[1614] * d[19] - d[2704] - d[401]) +
       d[2703] * d[781] - d[76] * xp * (d[383] + d[48] * (d[480] + d[541])) -
       xpdot * (M * (2.0 * d[19] * d[20] * d[548] + 13.0 * d[20] * d[55] -
                     d[300] - 21.0 * d[340] - d[48] * d[59]) +
                d[3] * (d[48] * (d[2705] + d[545] + 57.0 * d[57]) + d[489]))) *
          dv_977 +
      (-d[2709] *
       (-d[425] * (d[144] + d[424]) +
        xp * (d[436] * (d[410] + d[706]) - d[440] * (d[116] + d[437]) +
              ypdot * (d[441] + d[91] * (-d[177] - d[2704] - d[2708]))) -
        xpdot * (d[2706] + d[287] * (-d[2707] + 33.0 * d[55] - d[59]) +
                 d[435] * (d[324] + d[433] + d[59])))) *
          dv_6 +
      (-d[2716] *
       (d[189] * (-d[1443] * (d[338] + d[60] + d[729]) - d[427] * d[446] +
                  32.0 * d[459] * d[591]) +
        xp *
            (-d[125] * d[2665] + d[2662] * d[299] * ypdot +
             d[56] * (d[101] * d[2253] + d[1690] - d[1991] + 41.0 * d[274]) +
             d[60] * (-d[1443] * d[2372] + d[2039] + d[2714] + 62.0 * d[274])) +
        xpdot * (d[2715] * d[461] +
                 d[321] * (d[340] - 134.0 * d[341] - 305.0 * d[342] - d[454]) -
                 d[354] * d[356] * (-d[228] - 17.0 * d[55] - d[729]) +
                 d[427] * d[463] * d[591] + d[460]))) *
          dv_6 +
      (4.0 * d[48] * d[71]) * dv_3877 + (4.0 * d[384] * d[420] * xp) * dv_1908 +
      (8.0 * d[151] * d[390] * xp) * dv_1907 +
      (16.0 * d[375] * d[43] * xp) * dv_1906 +
      (4.0 * d[384] * d[420] * yp *
       (-d[397] * (d[145] * d[994] + d[382] * (d[144] + d[394])) +
        xp * (d[411] * (d[410] + 100.0 * d[603]) -
              d[415] * (d[382] * (d[116] + d[413]) - d[412] * d[994]) +
              ypdot * (d[1104] * (d[2719] + d[2720] + d[57]) +
                       d[2717] * d[373] + d[417])) +
        xpdot *
            (-d[403] * (d[382] * (d[400] + d[458] + d[729]) + d[399] * d[994]) +
             yp * (-d[2717] * d[353] + d[405] +
                   d[95] * (d[2718] - d[408] + d[57]))))) *
          dv_6 +
      (8.0 * d[151] * d[390] * yp *
       (d[1198] * (d[1053] * d[380] + d[2644]) +
        d[126] * (d[1628] + d[252] * d[2639] + d[2711]) + d[1514] * d[2592] +
        d[2088] * d[562] + d[2654] * d[496] + d[2710] +
        d[61] * (-d[1549] * d[2658] + d[1702] + d[2657] * d[683] + d[2713]))) *
          dv_6 +
      sc_6;
  DataVector& dv_1912 = temps.at(1633);
  sc_14 +=
      (16.0 * d[375] * d[43] * yp *
       (d[189] * (12.0 * d[19] * d[352] * ypdot - d[2721] * d[353] -
                  d[439] * (-d[199] + 124.0 * d[55] - d[57])) +
        xp * (d[1127] * d[121] * d[2722] +
              d[387] * (d[373] * d[994] +
                        d[92] * (-d[225] + 13.0 * d[55] - d[841])) -
              d[439] * (d[121] * (d[72] + d[796]) + d[230] - 53.0 * d[55])) +
        xpdot * (d[121] * d[312] * d[365] + d[155] * d[352] * d[63] +
                 d[36] * (-245.0 * d[1466] + d[1618] -
                          d[20] * d[234] * (d[1954] + d[91]) +
                          d[2098] * (-252.0 * d[20] + d[48]) +
                          d[227] * d[57] * (d[2696] + d[48]))))) *
          dv_6 +
      (16.0 * d[130] * d[19] * d[375] * d[92] * yp) * dv_1912 +
      (384.0 * d[20] * d[4] * d[52] * d[538] * d[551]) * dv_6 +
      (16.0 * d[130] * d[375] * d[92] * xp * yp *
       (-d[101] * d[541] * d[92] + d[1130] * d[121] * d[180] +
        d[121] * d[35] * d[479] +
        d[1382] * (-d[1348] * d[485] * yp - d[2206] * d[50] * d[92] +
                   d[491] * ypdot) -
        d[484] * (d[481] * d[994] + d[92] * (d[483] + d[57])))) *
          dv_6 +
      (32.0 * d[1133] * d[19] * d[382] * d[40] * d[475] * yp) * dv_6 +
      (128.0 * d[0] * d[1129] * d[19] * d[382] * d[475] * yp) * dv_6;
  DataVector& dv_10 = temps.at(10);
  DataVector& dv_224 = temps.at(221);
  sc_4 = (6.0 * d[130] * d[16]) * dv_10 * dv_224 * sc_14;
  DataVector& dv_1604 = temps.at(196);
  DataVector& dv_1844 = temps.at(1548);
  DataVector& dv_1845 = temps.at(269);
  DataVector& dv_3766 = temps.at(1206);
  DataVector& dv_3767 = temps.at(1150);
  DataVector& dv_3776 = temps.at(3163);
  DataVector& dv_3815 = temps.at(1637);
  DataVector& dv_3819 = temps.at(105);
  DataVector& dv_3821 = temps.at(259);
  DataVector& dv_3835 = temps.at(130);
  DataVector& dv_4062 = temps.at(916);
  DataVector& sc_21 = temps.at(3243);
  DataVector& sc_5 = temps.at(3227);
  sc_5 = (-d[1097]) * dv_3819 - dv_1604 * dv_3821 - dv_1844 * dv_3835 -
         dv_1845 * dv_3835 - dv_3766 * dv_4062 - dv_3767 * dv_4062 -
         dv_3776 * dv_3815 + sc_21;
  DataVector& dv_3783 = temps.at(3151);
  DataVector& dv_3787 = temps.at(3166);
  DataVector& dv_3812 = temps.at(1518);
  DataVector& dv_3816 = temps.at(1511);
  DataVector& dv_3817 = temps.at(1542);
  DataVector& dv_3822 = temps.at(66);
  DataVector& dv_3828 = temps.at(249);
  DataVector& dv_3829 = temps.at(271);
  DataVector& dv_3830 = temps.at(293);
  DataVector& dv_3831 = temps.at(266);
  DataVector& dv_3837 = temps.at(279);
  DataVector& dv_4002 = temps.at(1011);
  sc_5 += -dv_3783 * dv_3812 - dv_3787 * dv_3816 - dv_3817 * dv_3837 -
          dv_3822 * dv_3828 - dv_3822 * dv_4002 - dv_3829 * dv_3830 -
          dv_3829 * dv_3831;
  DataVector& dv_1867 = temps.at(305);
  DataVector& dv_3839 = temps.at(1608);
  DataVector& dv_3841 = temps.at(111);
  DataVector& dv_3842 = temps.at(250);
  DataVector& dv_3844 = temps.at(1541);
  DataVector& dv_3957 = temps.at(1029);
  DataVector& dv_4003 = temps.at(1131);
  DataVector& dv_4004 = temps.at(54);
  DataVector& dv_4005 = temps.at(1318);
  sc_5 += (6.0 * d[1060] * d[130] * d[16]) * dv_224 * dv_3957 -
          dv_3844 * ((3.0 * M * d[12]) * dv_1867 * dv_3817 - dv_3839 -
                     dv_3841 * dv_3842) -
          dv_4003 * dv_4004 - dv_4003 * dv_4005 + sc_11;
  DataVector& dv_11 = temps.at(11);
  DataVector& dv_19 = temps.at(19);
  DataVector& dv_3813 = temps.at(965);
  DataVector& dv_3820 = temps.at(233);
  DataVector& dv_4012 = temps.at(976);
  sc_5 += (6.0 * d[0] * d[1060]) * dv_11 * dv_229 * dv_4012 +
          (8.0 * d[0]) * dv_1496 * dv_1604 * dv_19 * dv_4012 +
          (4.0 * M * d[16] * d[74]) * dv_229 * dv_3817 * dv_3820 +
          (8.0 * M * d[16] * d[74]) * dv_1867 * dv_229 * dv_3813 +
          (12.0 * M * d[75]) * dv_11 * dv_19 * dv_3817 * dv_3820 + sc_4;
  sc_5 += (24.0 * M * d[75]) * dv_11 * dv_1604 * dv_1867 * dv_3820 +
          (24.0 * M * d[75]) * dv_11 * dv_1867 * dv_19 * dv_3813 +
          (36.0 * d[130] * d[16]) * dv_10 * dv_1604 * dv_229 * dv_3957 +
          (16.0 * M * d[16] * d[74]) * dv_1604 * dv_1867 * dv_19 * dv_3820 +
          (24.0 * M * d[1060] * d[75]) * dv_10 * dv_1867 * dv_19 * dv_3820;
  DataVector& dv_1611 = temps.at(1461);
  DataVector& sc_0 = temps.at(3222);
  sc_0 = ((1.0 / 48.0) * M * d[1024] * d[1025]) * dv_1611 * sc_5;
  DataVector& dv_225 = temps.at(222);
  DataVector& dv_3768 = temps.at(3186);
  DataVector& dv_3769 = temps.at(1065);
  DataVector& dv_3771 = temps.at(1139);
  DataVector& dv_3773 = temps.at(1047);
  DataVector& dv_3778 = temps.at(679);
  DataVector& dv_43 = temps.at(43);
  DataVector& sc_1 = temps.at(3223);
  DataVector& sc_3 = temps.at(3225);
  sc_3 = ((1.0 / 2.0) * M * d[23]) *
             ((-d[1060]) * dv_3773 + dv_3778 + dv_43 * xp) +
         ((7.0 / 48.0) * M * d[1024] * d[1025]) * dv_1604 * dv_225 * dv_3771 -
         dv_1604 * dv_3768 - dv_1604 * dv_3769 - dv_1604 + sc_0 + sc_1;
  DataVector& dv_1499 = temps.at(1379);
  get<0>(get<::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>,
                           Frame::Inertial>>(*result)) = dv_1499 * sc_3;
  DataVector& dv_4077 = temps.at(1204);
  DataVector& dv_4079 = temps.at(1238);
  DataVector& dv_4082 = temps.at(3195);
  DataVector& dv_4083 = temps.at(1614);
  DataVector& dv_4092 = temps.at(3192);
  sc_21 = -dv_4077 - dv_4079 + dv_4082 + dv_4083 + dv_4092;
  DataVector& dv_2181 = temps.at(1788);
  DataVector& dv_3808 = temps.at(3174);
  DataVector& dv_4084 = temps.at(108);
  DataVector& dv_4086 = temps.at(3140);
  DataVector& dv_4087 = temps.at(1356);
  DataVector& dv_4088 = temps.at(2947);
  DataVector& dv_4089 = temps.at(1658);
  DataVector& dv_4090 = temps.at(383);
  sc_21 +=
      d[54] * ((-d[50]) * dv_4087 + d[53] * ((xpdot * yp) * dv_3808 - dv_4088) +
               d[63] * (-dv_2181 + yp * (-dv_4089 + dv_4090 * ypdot)) +
               dv_4084 + dv_4086);
  DataVector& dv_1507 = temps.at(1387);
  sc_11 = dv_1507 * sc_21;
  DataVector& dv_1597 = temps.at(145);
  DataVector& dv_4093 = temps.at(1485);
  DataVector& dv_4094 = temps.at(1513);
  sc_4 = (d[1062] * d[5] + d[68] * d[77] + ypdot) * dv_1597 +
         (-d[17]) * dv_4093 + dv_4094 + sc_11;
  sc_5 = (2.0 * rp) * sc_4;
  DataVector& dv_1606 = temps.at(205);
  DataVector& dv_1607 = temps.at(198);
  DataVector& dv_1616 = temps.at(1466);
  DataVector& dv_3784 = temps.at(3209);
  DataVector& dv_3785 = temps.at(3123);
  DataVector& dv_3788 = temps.at(1465);
  DataVector& dv_4066 = temps.at(884);
  DataVector& dv_4074 = temps.at(1339);
  DataVector& dv_4075 = temps.at(277);
  DataVector& dv_4076 = temps.at(1163);
  DataVector& dv_79 = temps.at(77);
  DataVector& dv_86 = temps.at(84);
  sc_1 = (9.0 * M) * dv_1607 * dv_1611 * dv_1616 - dv_1606 * dv_3784 -
         dv_1607 * dv_3785 -
         dv_3788 * ((-d[17]) * dv_4076 + dv_4066 + dv_4075) - dv_4074 * dv_79 -
         dv_4074 * dv_86 + sc_5;
  sc_0 = (-d[2572]) * sc_1;
  DataVector& dv_122 = temps.at(120);
  DataVector& dv_138 = temps.at(136);
  DataVector& dv_2009 = temps.at(1515);
  DataVector& dv_2019 = temps.at(95);
  DataVector& dv_2024 = temps.at(1550);
  DataVector& dv_210 = temps.at(207);
  DataVector& dv_4158 = temps.at(1338);
  DataVector& dv_543 = temps.at(497);
  DataVector& dv_586 = temps.at(538);
  DataVector& dv_96 = temps.at(94);
  sc_30 = dv_138 * ((-d[2727]) * dv_2019 + (-d[54]) * dv_2024 + d[48] * dv_96) +
          dv_153 * (d[110] * (dv_2009 + dv_586) + d[2729] * (dv_122 + dv_543) +
                    d[54] * (dv_210 + dv_4158));
  DataVector& dv_121 = temps.at(119);
  DataVector& dv_280 = temps.at(272);
  DataVector& dv_294 = temps.at(286);
  DataVector& dv_328 = temps.at(311);
  DataVector& dv_380 = temps.at(363);
  DataVector& dv_51 = temps.at(51);
  DataVector& dv_680 = temps.at(630);
  DataVector& dv_771 = temps.at(719);
  sc_30 += dv_294 * (d[2727] * (dv_280 + dv_51) + d[48] * (dv_380 + dv_771) +
                     d[54] * (dv_121 + dv_328 + dv_680));
  DataVector& dv_154 = temps.at(152);
  DataVector& dv_1644 = temps.at(1489);
  DataVector& dv_1651 = temps.at(1496);
  DataVector& dv_2025 = temps.at(1543);
  DataVector& dv_2035 = temps.at(1250);
  DataVector& dv_340 = temps.at(323);
  DataVector& dv_38 = temps.at(38);
  DataVector& dv_3972 = temps.at(417);
  DataVector& dv_3973 = temps.at(392);
  sc_30 += dv_340 * (d[110] * (dv_38 + dv_3973) + d[2728] * dv_1651 +
                     d[54] * dv_2035) +
           dv_853 * (d[2728] * dv_1644 + d[48] * (dv_154 + dv_3972 + dv_51) +
                     d[54] * dv_2025);
  sc_28 = (-d[119]) * sc_30;
  DataVector& dv_2053 = temps.at(402);
  DataVector& dv_2095 = temps.at(635);
  DataVector& dv_4156 = temps.at(163);
  sc_33 = -dv_875 * (d[114] * dv_2053 + d[72] * dv_4156 + d[968] * dv_2095);
  DataVector& dv_2039 = temps.at(1685);
  DataVector& dv_2072 = temps.at(1708);
  DataVector& dv_2086 = temps.at(1313);
  DataVector& dv_2201 = temps.at(1805);
  DataVector& dv_3015 = temps.at(2540);
  DataVector& dv_4160 = temps.at(398);
  DataVector& dv_4164 = temps.at(1156);
  DataVector& dv_4166 = temps.at(1659);
  DataVector& dv_496 = temps.at(461);
  DataVector& dv_877 = temps.at(809);
  DataVector& dv_880 = temps.at(812);
  sc_33 += -dv_877 * (d[114] * (dv_3015 + dv_4164 + dv_496) +
                      d[118] * (dv_2072 + dv_4166 + dv_496) + d[72] * dv_4160) -
           dv_880 * ((-d[114]) * dv_2039 + (-d[118]) * dv_2086 + dv_2201);
  DataVector& dv_2041 = temps.at(396);
  DataVector& dv_2088 = temps.at(425);
  DataVector& dv_4157 = temps.at(133);
  DataVector& dv_918 = temps.at(848);
  sc_33 += -dv_918 * (d[114] * dv_2041 + d[118] * dv_2088 + d[72] * dv_4157);
  DataVector& dv_130 = temps.at(128);
  DataVector& dv_185 = temps.at(182);
  DataVector& dv_2045 = temps.at(1690);
  DataVector& dv_4155 = temps.at(333);
  DataVector& dv_4165 = temps.at(1650);
  DataVector& dv_426 = temps.at(401);
  DataVector& dv_455 = temps.at(428);
  DataVector& dv_486 = temps.at(454);
  DataVector& dv_919 = temps.at(849);
  sc_33 += (d[120] * rp) * Dx * (d[578] * dv_185 + d[9] * dv_2045) -
           dv_919 * (d[114] * (dv_130 + dv_426 + dv_486) +
                     d[118] * (dv_4165 + dv_455 + dv_680) + d[72] * dv_4155);
  sc_30 = sc_33 * xpdot;
  DataVector& dv_131 = temps.at(129);
  DataVector& dv_1996 = temps.at(1531);
  DataVector& dv_2059 = temps.at(1701);
  DataVector& dv_2068 = temps.at(451);
  DataVector& dv_2073 = temps.at(1709);
  DataVector& dv_3981 = temps.at(194);
  DataVector& dv_4161 = temps.at(920);
  DataVector& dv_4162 = temps.at(407);
  DataVector& dv_4163 = temps.at(3143);
  DataVector& dv_487 = temps.at(455);
  DataVector& dv_594 = temps.at(545);
  DataVector& dv_625 = temps.at(576);
  sc_36 = dv_2059 * dv_4161 - dv_2068 * dv_3956 +
          dv_3981 * (-dv_2073 - dv_4164) + dv_4163 * (-dv_1996 - dv_4162) +
          dv_625 * (dv_131 + dv_487 + dv_594);
  DataVector& dv_2069 = temps.at(1706);
  DataVector& dv_3980 = temps.at(371);
  sc_36 += -dv_2069 * dv_3980;
  sc_17 = d[114] * sc_36;
  DataVector& dv_2075 = temps.at(142);
  DataVector& dv_508 = temps.at(469);
  DataVector& dv_556 = temps.at(509);
  sc_34 = dv_2075 * dv_3955 + dv_3981 * (-dv_4166 - dv_508 - dv_51) +
          dv_4163 * (-dv_4162 - dv_556);
  DataVector& dv_2078 = temps.at(467);
  DataVector& dv_2083 = temps.at(472);
  DataVector& dv_463 = temps.at(408);
  DataVector& dv_503 = temps.at(464);
  sc_34 += -dv_2078 * dv_3956 - dv_2083 * dv_3980 +
           dv_625 * (-dv_4165 + dv_463 + dv_503);
  sc_36 = d[118] * sc_34;
  DataVector& dv_161 = temps.at(159);
  DataVector& dv_168 = temps.at(166);
  DataVector& dv_2386 = temps.at(1976);
  DataVector& dv_601 = temps.at(552);
  sc_12 = (-d[1975]) * dv_6 *
              (d[50] * dv_4155 + d[63] * (dv_168 + dv_601) + dv_161 -
               dv_2386 * dv_4) +
          sc_17 + sc_36;
  sc_33 = sc_12 * ypdot;
  sc_13 = sc_28 + sc_30 + sc_33;
  DataVector& dv_3982 = temps.at(370);
  sc_6 = dv_3982 * sc_13;
  DataVector& dv_2000 = temps.at(330);
  DataVector& dv_2001 = temps.at(101);
  DataVector& dv_2985 = temps.at(2512);
  DataVector& dv_333 = temps.at(316);
  DataVector& dv_381 = temps.at(364);
  sc_33 = -dv_294 * (d[1109] * (dv_2985 + dv_51) +
                     d[118] * (dv_2001 + dv_333 + dv_556) + d[49] * dv_4155) -
          dv_340 * (d[114] * dv_2000 + d[2726] * dv_381 + d[49] * dv_4156);
  DataVector& dv_106 = temps.at(104);
  DataVector& dv_16 = temps.at(16);
  DataVector& dv_1988 = temps.at(337);
  DataVector& dv_2007 = temps.at(1522);
  DataVector& dv_572 = temps.at(525);
  DataVector& dv_676 = temps.at(626);
  sc_33 += (d[19] * d[20]) * Dy *
               ((-d[114]) * (dv_106 + dv_4158) + (-d[49]) * dv_4160 +
                (9.0 * d[12]) * (dv_16 + dv_572 + dv_676)) -
           dv_853 * (d[114] * dv_1988 + d[2725] * dv_2007 + d[49] * dv_4157);
  DataVector& dv_125 = temps.at(123);
  DataVector& dv_2004 = temps.at(212);
  DataVector& dv_2215 = temps.at(1817);
  sc_33 += Dy * d[55] * (d[114] * dv_125 + d[968] * dv_2004 - dv_2215);
  DataVector& dv_1781 = temps.at(1446);
  sc_13 = -dv_1781 * sc_33;
  DataVector& dv_3962 = temps.at(348);
  DataVector& dv_3963 = temps.at(313);
  DataVector& dv_3970 = temps.at(1202);
  DataVector& dv_3971 = temps.at(443);
  DataVector& dv_4065 = temps.at(1017);
  DataVector& dv_4104 = temps.at(264);
  sc_14 = (-yp) * dv_3970 + d[1063] * dv_3963 + dv_3962 * dv_4104 +
          dv_3971 * dv_4065 + sc_13 + sc_6;
  DataVector& dv_1864 = temps.at(200);
  sc_21 = dv_1864 * sc_14;
  DataVector& dv_3958 = temps.at(3175);
  DataVector& dv_3959 = temps.at(1507);
  DataVector& dv_3960 = temps.at(697);
  DataVector& dv_3964 = temps.at(342);
  DataVector& dv_4096 = temps.at(1661);
  DataVector& dv_4097 = temps.at(1535);
  DataVector& dv_4102 = temps.at(56);
  DataVector& dv_4103 = temps.at(3180);
  DataVector& dv_4105 = temps.at(116);
  sc_11 = d[1186] * dv_4097 - dv_3958 * dv_4103 - dv_3959 * dv_4096 -
          dv_3960 * dv_4105 + dv_3964 * dv_4102 + sc_21;
  sc_4 = (-d[108]) * sc_11;
  sc_6 = xpdot;
  DataVector& dv_2168 = temps.at(485);
  DataVector& dv_2385 = temps.at(1975);
  DataVector& dv_3296 = temps.at(2783);
  sc_6 *= d[116] * (d[1] * dv_751 - 5.0 * dv_119) +
          xp * (Dy * d[697] + d[37] * dv_4 - dv_2385) +
          yp * ((-xp) * dv_2168 + d[849] * dv_3296);
  DataVector& dv_4114 = temps.at(1611);
  DataVector& dv_4115 = temps.at(240);
  sc_13 = d[1803] * dv_4114 + sc_6 + ypdot * (d[2582] * dv_4 + dv_4115);
  sc_14 = (-d[2583]) * sc_13;
  DataVector& dv_2350 = temps.at(1940);
  DataVector& dv_3637 = temps.at(3107);
  sc_6 = d[121] * ((-d[2599] + d[313] + d[724] + d[751]) * Dy +
                   d[2829] * dv_3637 + d[310] * dv_2350) +
         d[129] * dv_752;
  DataVector& dv_2949 = temps.at(2477);
  DataVector& dv_3902 = temps.at(1646);
  DataVector& dv_3951 = temps.at(3218);
  sc_6 += d[191] *
          (d[102] * dv_3912 + d[259] * ((d[2646] + d[6] + 4.0) * Dx + dv_2949) +
           d[49] * ((d[535] + 23.0) * dv_751 + dv_3902 + dv_3951));
  DataVector& dv_1115 = temps.at(1012);
  DataVector& dv_4136 = temps.at(885);
  DataVector& dv_728 = temps.at(677);
  sc_6 += d[20] * ((M * d[20] - d[1393] - d[313] + 11.0 * d[48] * yp * ypdot) *
                       dv_728 +
                   (d[259] * d[307]) * dv_0 + 30.0 * dv_3884) +
          d[55] * (d[1126] * dv_0 + d[1529] * dv_1115 + dv_4136);
  DataVector& dv_1700 = temps.at(1530);
  DataVector& dv_1727 = temps.at(1546);
  DataVector& dv_3880 = temps.at(3111);
  DataVector& dv_3913 = temps.at(3168);
  sc_6 += d[61] * (d[20] * (dv_1727 + dv_3880) +
                   d[259] * ((d[2019] + d[2641]) * Dx - dv_3913) +
                   d[49] * (d[2819] * dv_1700 + dv_3866));
  sc_13 = (-d[331]) * sc_6;
  DataVector& dv_3552 = temps.at(3025);
  DataVector& dv_790 = temps.at(738);
  sc_28 = (-d[50]) * (36.0 * dv_265 + dv_3552 + 7.0 * dv_780) +
          d[1443] * ((-d[2670] + d[286] + 55.0) * Dx - 610.0 * dv_265) +
          d[1490] * dv_790;
  DataVector& dv_3937 = temps.at(3155);
  sc_28 += d[631] * ((d[1825] + 65.0) * dv_751 + dv_3937 + 82.0 * dv_752);
  sc_30 = d[1863] * sc_28;
  DataVector& dv_2208 = temps.at(1811);
  DataVector& dv_240 = temps.at(236);
  sc_33 =
      (-d[1800]) * ((d[2666] * d[6] - d[2667]) * Dx +
                    (d[1682] + d[259] + d[2848]) * dv_2208) +
      (-d[2004]) *
          ((260.0 * d[252] - d[259] * (d[463] + 65.0) + 28.0 * d[319]) * dv_0 +
           (d[104] + d[2579] + d[471]) * dv_1115 +
           (-d[1293] - d[1856] - d[2849]) * dv_240) +
      (d[1583] * d[1618]) * dv_0;
  DataVector& dv_3932 = temps.at(1070);
  DataVector& dv_69 = temps.at(69);
  sc_33 += d[1466] * ((-d[402]) * dv_1115 + d[2664] * dv_0 + d[27] * dv_3932) +
           d[1581] * ((-d[2640]) * dv_231 + d[2588] * dv_751 + d[86] * dv_69);
  DataVector& dv_1709 = temps.at(1534);
  DataVector& dv_3858 = temps.at(3119);
  sc_33 += d[1845] * ((d[271] - 17.0 * d[276] - 520.0 * d[277] +
                       d[631] * (d[463] + 15.0)) *
                          dv_0 +
                      (-d[101] * d[2367] - d[2366] + d[272] + 34.0 * d[274]) *
                          dv_1709 +
                      (d[2578] + d[2850]) * dv_3858) +
           sc_30;
  DataVector& dv_3924 = temps.at(3181);
  DataVector& dv_3933 = temps.at(1154);
  DataVector& dv_3934 = temps.at(1086);
  DataVector& dv_4138 = temps.at(1307);
  sc_33 += d[1872] * (d[1443] * (-252.0 * dv_265 + dv_3924) +
                      d[50] * (dv_3913 + dv_4138) +
                      d[631] * (dv_3934 - dv_751 + 18.0 * dv_752) + dv_3933);
  DataVector& dv_2850 = temps.at(2402);
  DataVector& dv_3546 = temps.at(3019);
  DataVector& dv_4139 = temps.at(983);
  sc_33 +=
      d[1905] *
      ((-d[50]) * (76.0 * dv_265 + dv_3546 + 27.0 * dv_780) +
       d[1443] * ((d[2672] - 47.0 * d[6]) * Dx - 268.0 * dv_265) +
       d[631] * ((64.0 * d[6] + 83.0) * dv_751 + dv_3934 + dv_4139) + dv_2850);
  DataVector& dv_1579 = temps.at(1451);
  sc_33 +=
      d[1931] * dv_1579 +
      d[2669] * ((d[1153] - d[1443] * d[2371] + d[2714] + 84.0 * d[274]) * Dy +
                 (d[273] * (d[2171] + 83.0) - 28.0 * d[276] - 610.0 * d[277] +
                  d[283]) *
                     dv_0 +
                 (-d[2022] - d[243] - d[2579]) * dv_3858);
  sc_6 = (-d[467]) * sc_33;
  DataVector& dv_4117 = temps.at(680);
  DataVector& dv_4119 = temps.at(246);
  sc_30 = d[19] * (d[1928] * dv_0 + dv_4119) +
          d[52] * ((-d[2823]) * dv_4117 + d[2822] * dv_751);
  DataVector& dv_4120 = temps.at(458);
  DataVector& dv_4121 = temps.at(830);
  DataVector& dv_4122 = temps.at(956);
  sc_30 +=
      xp * (d[20] * (d[2824] * dv_751 - dv_4122) + d[91] * dv_4121 + dv_4120) +
      yp * ((d[138] + d[158] * d[20] - d[213]) * dv_0 +
            (-d[1244] * d[2825]) * dv_5 + d[2826] * dv_1709);
  sc_33 = d[168] * sc_30;
  sc_30 = d[208];
  DataVector& dv_1237 = temps.at(1132);
  DataVector& dv_1531 = temps.at(1408);
  DataVector& dv_3409 = temps.at(2886);
  DataVector& dv_3885 = temps.at(3108);
  DataVector& dv_4141 = temps.at(1350);
  sc_30 *=
      (-d[59]) * ((44.0 * M * yp - 120.0 * d[252] - d[448] * ypdot) * dv_1237 +
                  d[2836] * dv_1531) +
      d[56] * (d[2837] * dv_5 + dv_3409 + dv_3885) +
      d[60] * ((-d[2838]) * dv_240 + d[2137] * dv_0 + 296.0 * dv_3884) +
      dv_4141;
  DataVector& dv_1630 = temps.at(1475);
  DataVector& dv_3889 = temps.at(3109);
  DataVector& dv_3895 = temps.at(3138);
  DataVector& dv_3897 = temps.at(1079);
  DataVector& dv_4142 = temps.at(1111);
  sc_12 =
      d[120] * (-dv_3889 - dv_4142) +
      d[123] * ((-d[1161] - d[2615] + 324.0 * d[48] * yp * ypdot +
                 16.0 * d[50] * ypdot) *
                    dv_0 +
                d[1993] * dv_1709 + dv_3897) +
      d[124] * ((d[1161] + d[1769] - 43.0 * d[273] + 288.0 * d[277]) * dv_1630 +
                (-d[1987] - d[2013] - d[2839]) * dv_1709 + dv_3895);
  DataVector& dv_3023 = temps.at(2507);
  DataVector& dv_3878 = temps.at(1067);
  DataVector& dv_3886 = temps.at(3120);
  DataVector& dv_3892 = temps.at(1209);
  DataVector& dv_3894 = temps.at(52);
  DataVector& dv_3910 = temps.at(1141);
  sc_12 += d[127] * ((-d[1450]) * (dv_3023 + dv_3910) - dv_3892) +
           d[128] * (dv_3894 - dv_4142) +
           d[1841] * ((-d[1425]) * dv_3910 + d[2614] * dv_751 + dv_3886) +
           d[2613] * dv_3878;
  DataVector& dv_3485 = temps.at(2960);
  DataVector& dv_3887 = temps.at(1234);
  sc_12 += d[345] * ((d[1442] - d[259] + d[2629]) * dv_1630 +
                     (-d[1421] - d[2206]) * dv_3485 + dv_3887);
  sc_28 = d[2616] * sc_12;
  DataVector& dv_3940 = temps.at(1083);
  DataVector& dv_4148 = temps.at(1648);
  sc_36 =
      (d[1466] * d[20]) *
          ((4.0 * d[1681] * yp - 497.0 * d[36]) * dv_0 +
           (d[1676] + d[2516]) * dv_240 + (-d[2096]) * dv_4148) +
      (d[20] * d[362]) *
          (d[259] * ((d[2680] - 251.0 * d[6]) * Dx - 504.0 * dv_265) + dv_3940);
  DataVector& dv_3567 = temps.at(3040);
  DataVector& dv_3939 = temps.at(3214);
  sc_36 +=
      d[1797] * (d[273] * (-490.0 * dv_265 + dv_3924) + dv_3939) +
      d[1804] * ((d[1153] - d[1699] + d[2679] - 497.0 * d[274]) * dv_0 +
                 (-d[2213] * d[273] + d[2691]) * dv_240 + d[1713] * dv_3567);
  DataVector& dv_3151 = temps.at(2640);
  sc_36 +=
      d[1821] *
          ((d[2313] * d[6] + d[2689]) * Dx + (-d[2681] - d[2685]) * dv_2208) +
      d[1931] * dv_3151 +
      d[1951] *
          ((d[2683] * d[6] + d[2686]) * Dx + (d[2682] - d[2684]) * dv_2208);
  DataVector& dv_3941 = temps.at(3202);
  sc_36 += d[2004] * ((-d[1704] + 36.0 * yp * ypdot) * dv_3941 +
                      d[2232] * dv_5 + d[2688] * dv_1564) +
           d[2011] * ((d[2234] - d[2238] * d[273]) * Dy + (-d[1719]) * dv_3567 +
                      d[2690] * dv_1564) +
           d[2673] * dv_3151;
  sc_36 += d[2678] * ((-d[2676] + d[2677]) * Dx + d[1] * dv_265);
  sc_12 = d[377] * sc_36;
  DataVector& dv_833 = temps.at(776);
  sc_17 = (d[116] * d[55]) * ((-179.0 * d[1339] - d[2178] - d[2653]) * dv_0 +
                              (-d[1533] * d[535]) * dv_833 + d[2847] * dv_240) +
          (d[1466] * d[2058]) * dv_1630;
  DataVector& dv_3927 = temps.at(1151);
  DataVector& dv_3930 = temps.at(1219);
  DataVector& dv_3950 = temps.at(3159);
  sc_17 +=
      d[1177] * ((-d[308]) * ((d[1596] + d[2661]) * dv_2170 + dv_3927) +
                 (-d[683]) * ((d[2660] + d[6]) * dv_3296 + 188.0 * dv_265) +
                 d[223] * dv_3950 + dv_3930);
  DataVector& dv_3865 = temps.at(958);
  DataVector& dv_3920 = temps.at(3188);
  DataVector& dv_3926 = temps.at(1077);
  DataVector& dv_3929 = temps.at(3158);
  DataVector& dv_4154 = temps.at(940);
  sc_17 += d[2114] * (d[223] * dv_751 + d[308] * (dv_3926 + dv_4154) +
                      d[683] * (-170.0 * dv_265 + dv_3924) +
                      d[996] * (dv_3865 + dv_3929 + dv_751) + dv_3920);
  DataVector& dv_1549 = temps.at(1423);
  DataVector& dv_3108 = temps.at(2600);
  DataVector& dv_3918 = temps.at(3161);
  DataVector& dv_3921 = temps.at(1015);
  sc_17 +=
      d[2619] * (d[102] * dv_3921 + d[259] * (dv_3918 + dv_4154) - dv_3108) +
      d[2649] * ((d[1172] * (d[1620] + 2.0) - d[151] * d[1866] -
                  122.0 * d[1543] + d[223] + d[466]) *
                     dv_0 +
                 (d[2650] + d[380]) * dv_1549 +
                 (d[1124] * d[683] + d[1536] * d[2643] + d[1665] + d[2712] +
                  d[319] * d[930]) *
                     dv_240);
  DataVector& dv_3919 = temps.at(1743);
  sc_17 += d[2656] * ((-d[2651] * d[638] + d[2655]) * Dx + d[2846] * dv_752) +
           d[300] * ((d[1292] + d[1339] - d[2610] + d[287]) * dv_0 +
                     (-d[2834]) * dv_1115 + d[2637] * dv_34) +
           d[340] * dv_3919;
  sc_36 = d[392] * sc_17;
  sc_34 = (-d[2004]) *
          ((214.0 * d[101] + d[1633] + d[2693] + 312.0 * d[274]) * dv_0 +
           (d[2300] * d[396] - 268.0 * d[252] - 81.0 * d[319]) * dv_5 +
           210.0 * dv_3884);
  DataVector& dv_3946 = temps.at(3187);
  DataVector& dv_3952 = temps.at(1361);
  sc_34 += (-d[2698]) *
           (d[1550] * (107.0 * dv_751 - 20.0 * dv_752) +
            d[2700] * ((d[1735] + 21.0 * d[6]) * dv_3296 + 100.0 * dv_265) -
            dv_3946 + dv_3952);
  DataVector& dv_2712 = temps.at(2265);
  DataVector& dv_3852 = temps.at(1120);
  sc_34 += (d[116] * d[362]) *
               ((-d[136]) * ((d[2694] + 23.0 * d[6]) * Dx + 12.0 * dv_265) +
                d[416] * dv_3852 + d[49] * (dv_2208 + dv_2712)) +
           (d[1931] * d[554]) * dv_1630;
  DataVector& dv_1827 = temps.at(1612);
  DataVector& dv_3942 = temps.at(929);
  sc_34 +=
      d[1581] * (d[396] * (Dx * d[503] + dv_3913) + dv_1827 + dv_3944) +
      d[1845] *
          ((-198.0 * d[1543] + d[1551] + d[1641] + 524.0 * d[604] -
            1392.0 * d[992]) *
               dv_0 +
           (-d[393] + d[91]) * dv_3942 +
           (-220.0 * d[1182] + d[1641] * ypdot - d[2292] * d[2293] - d[2697]) *
               dv_240);
  sc_34 += d[1858] * ((-d[2700]) * ((d[2699] + d[2844]) * Dx + 108.0 * dv_265) +
                      d[1550] * (dv_3904 - 104.0 * dv_752) + dv_3949);
  DataVector& dv_3947 = temps.at(1094);
  DataVector& dv_4068 = temps.at(994);
  sc_34 += d[1872] * (d[1172] * (dv_751 + 64.0 * dv_752) + d[1641] * dv_751 +
                      d[2293] * (-30.0 * dv_265 + dv_3908) + dv_3947) +
           d[2320] * dv_4068;
  DataVector& dv_3943 = temps.at(1169);
  DataVector& dv_4130 = temps.at(1520);
  sc_34 +=
      d[2669] *
          ((d[1661] * d[319] - d[2293] * d[2319] + d[2701]) * Dy +
           (-156.0 * d[1543] - d[1551] + d[1641] + d[2037] - 768.0 * d[992]) *
               dv_1564 +
           (-d[248] - d[48]) * dv_4130) +
      d[2692] * ((-99.0 * d[159] + d[409] + d[91]) * dv_0 + dv_3943);
  sc_17 = d[422] * sc_34;
  sc_7 =
      (-d[1581]) * dv_3904 +
      (-d[1841]) * ((-456.0 * d[159] + 99.0 * d[20] - d[559]) * dv_0 +
                    (d[1249] * d[502] + d[2663] + d[321]) * dv_728 - dv_4130);
  DataVector& dv_2008 = temps.at(171);
  DataVector& dv_4150 = temps.at(1149);
  DataVector& dv_4151 = temps.at(1061);
  sc_7 += (-d[2645]) * (d[1443] * (-dv_2008 - dv_4151) +
                        d[2091] * (-24.0 * dv_265 + dv_3908) +
                        d[50] * (dv_3951 + dv_4150) - 240.0 * dv_3882);
  DataVector& dv_2374 = temps.at(1964);
  sc_7 += (-d[362]) * ((-d[1299]) * ((-d[503]) * Dx + 8.0 * dv_265) +
                       d[20] * (dv_2374 + 52.0 * dv_752) + dv_3108);
  DataVector& dv_2761 = temps.at(2314);
  DataVector& dv_2812 = temps.at(2365);
  DataVector& dv_2851 = temps.at(2403);
  sc_7 += d[1905] * ((-d[20]) * (dv_2851 + 146.0 * dv_752) +
                     d[1299] * ((d[2624] + 17.0 * d[6]) * dv_3296 + dv_3948) +
                     d[91] * ((10.0 * xpdot) * Dy - dv_2812)) +
          d[2622] * dv_2761;
  sc_7 +=
      d[2623] * ((-yp) * dv_2469 + d[2050] * dv_0) +
      d[341] *
          ((-164.0 * d[101] + d[2136] + 456.0 * d[274] - 89.0 * d[50]) * dv_0 +
           (d[2091] * d[2265] - d[2264] + 66.0 * d[277] - d[283]) * dv_1709 +
           552.0 * dv_3884);
  DataVector& dv_4152 = temps.at(918);
  sc_7 += d[342] * ((-d[2269]) * dv_1531 + d[2625] * dv_0 + 600.0 * dv_3884) +
          d[343] * ((36.0 * M * d[20] * (d[2241] + d[2626]) - d[2627]) * Dx +
                    (31.0 * d[101] + d[1487] + 102.0 * d[274]) * dv_4152);
  sc_34 = d[443] * sc_7;
  DataVector& dv_753 = temps.at(701);
  sc_15 = (-d[123]) *
          (d[101] * ((-d[2633]) * Dx + 146.0 * dv_265) + d[1409] * dv_3914 +
           d[273] * (dv_4150 + dv_4151) + d[807] * (dv_1530 + dv_753));
  DataVector& dv_1514 = temps.at(1393);
  DataVector& dv_4153 = temps.at(985);
  sc_15 += (-d[124]) * (d[101] * ((d[2635] + d[2844]) * Dx + 296.0 * dv_265) +
                        d[273] * (dv_2008 - 16.0 * dv_752) +
                        d[313] * (dv_1514 + dv_4152) + dv_4153);
  DataVector& dv_3566 = temps.at(3039);
  sc_15 += (-d[127]) *
               ((d[1709] + d[2312] + d[2845]) * dv_0 +
                (-d[1009] - d[2631] - d[2839]) * dv_240 + d[2219] * dv_3896) +
           (-d[128]) * ((d[1348] + d[1769] - d[1924] + d[2845]) * dv_0 +
                        (d[101] * d[2187] + d[2634]) * dv_240 +
                        (d[243] + 61.0 * d[48]) * dv_3858) +
           (-d[1466]) * dv_3566;
  DataVector& dv_1712 = temps.at(28);
  DataVector& dv_2528 = temps.at(2116);
  sc_15 += (-d[299]) * (d[248] * dv_1115 + d[2630] * dv_0 + d[2843] * dv_4136) +
           (-d[362]) * (d[2843] * dv_751 + d[86] * dv_1712 + dv_2528) +
           (d[120] * d[2099] * xpdot) * Dx + (4.0 * M * d[361] * ypdot) * Dy;
  sc_15 += d[122] * (d[101] * ((d[2181] + d[2695]) * Dx - 114.0 * dv_265) +
                     d[273] * (29.0 * dv_751 - 30.0 * dv_752) +
                     d[313] * dv_751 - dv_4153);
  sc_7 = d[550] * sc_15;
  DataVector& dv_4113 = temps.at(1349);
  DataVector& dv_4116 = temps.at(824);
  DataVector& dv_4127 = temps.at(838);
  sc_21 = (-d[2594] * d[527]) * dv_4127 + (d[476] * d[839] * xp) * dv_3860 +
          d[128] * dv_3845 - dv_4113 - dv_4116 + sc_13 + sc_14 + sc_33 + sc_6;
  DataVector& dv_4111 = temps.at(988);
  DataVector& dv_4112 = temps.at(1512);
  sc_21 += d[171] * (dv_4112 + xpdot * (d[1122] * dv_231 - dv_4111)) + sc_30;
  DataVector& dv_4128 = temps.at(169);
  DataVector& dv_4129 = temps.at(158);
  DataVector& dv_4131 = temps.at(1644);
  DataVector& dv_4133 = temps.at(138);
  DataVector& dv_4134 = temps.at(149);
  DataVector& dv_4135 = temps.at(3148);
  sc_21 +=
      d[237] * ((-d[122]) * dv_4128 + (-d[55]) * dv_4129 + (-d[61]) * dv_4133 +
                d[50] * (d[1566] * dv_0 + dv_4131) - dv_4134 - dv_4135);
  DataVector& dv_4106 = temps.at(1615);
  DataVector& dv_4108 = temps.at(1618);
  DataVector& dv_4109 = temps.at(253);
  DataVector& dv_4110 = temps.at(504);
  DataVector& dv_4123 = temps.at(877);
  DataVector& dv_4124 = temps.at(380);
  DataVector& dv_4125 = temps.at(670);
  DataVector& dv_4126 = temps.at(1490);
  sc_21 += d[2577] *
               ((-d[2821]) * dv_4109 + (-xpdot) * dv_4110 + dv_4106 - dv_4108) +
           d[260] * (d[227] * dv_4126 + d[2827] * dv_3912 + d[52] * dv_4123 +
                     dv_4124 + dv_4125);
  DataVector& dv_1762 = temps.at(1565);
  DataVector& dv_4137 = temps.at(3146);
  DataVector& dv_4140 = temps.at(1573);
  sc_21 += d[2607] * (d[55] * (d[2604] * dv_0 + d[2830] * dv_1762 - dv_4137) +
                      dv_4140) +
           sc_28;
  DataVector& dv_4143 = temps.at(1287);
  DataVector& dv_4144 = temps.at(1643);
  DataVector& dv_4146 = temps.at(1657);
  DataVector& dv_4147 = temps.at(1200);
  DataVector& dv_4149 = temps.at(1147);
  sc_21 +=
      d[2621] *
          (d[120] * (d[2841] * dv_0 + dv_4147) +
           d[122] * (d[259] * ((-d[2840]) * Dx + dv_4144) + dv_4143 - dv_4146) -
           dv_4149) +
      sc_12 + sc_17 + sc_34 + sc_36;
  sc_21 += d[53] * dv_3855 + sc_7;
  sc_11 = (2.0 * d[0]) * dv_1496 * dv_229 * sc_21;
  sc_17 = (-d[2826]) * dv_34 + (d[2825] * d[62]) * dv_1115 + d[1120] * dv_0 +
          d[19] * ((d[1118] * xpdot * yp) * Dx - dv_4119) +
          d[52] * ((-d[2822]) * dv_751 + d[2823] * dv_4117);
  sc_17 += xp * ((-d[91]) * dv_4121 + d[20] * ((-d[2824]) * dv_751 + dv_4122) -
                 dv_4120);
  sc_34 = (-d[168]) * sc_17;
  DataVector& dv_3284 = temps.at(2771);
  DataVector& dv_626 = temps.at(577);
  sc_36 = (d[332] * d[452]) * dv_4 + d[373] * dv_3284 +
          d[638] * ((d[317] * d[61]) * Dx + Dy * d[401] + d[2829] * dv_627 +
                    d[315] * dv_626 + 26.0 * dv_138);
  DataVector& dv_162 = temps.at(160);
  DataVector& dv_2111 = temps.at(1728);
  DataVector& dv_2298 = temps.at(1891);
  sc_36 += d[80] * (d[312] * (d[395] * dv_4 - dv_162 + dv_854) +
                    d[49] * (d[190] * dv_5 + d[423] * dv_4 + d[52] * dv_2298 +
                             22.0 * dv_190)) +
           d[88] * dv_2111;
  DataVector& dv_284 = temps.at(276);
  DataVector& dv_3632 = temps.at(3102);
  DataVector& dv_923 = temps.at(853);
  DataVector& dv_926 = temps.at(835);
  sc_36 +=
      xpdot * ((-d[209] - d[432]) * dv_926 +
               (d[243] * (-d[102] + d[1938] * d[48] + d[2833])) * dv_119 +
               d[1126] * dv_284 + d[307] * dv_3632 + d[311] * dv_231 + dv_923);
  sc_17 = (-d[331]) * sc_36;
  sc_36 = (d[207] * d[48]);
  DataVector& dv_1838 = temps.at(1623);
  sc_36 *= (d[631] * (-d[1578] * d[50] + d[1776] + d[194] * ypdot)) * dv_0 +
           (-d[120] * d[2836]) * dv_1709 +
           d[19] * ((240.0 * d[151] * d[78] - d[1645] + 28.0 * d[2087] +
                     210.0 * d[2355]) *
                        dv_0 +
                    (-d[2838]) * dv_1838 + (296.0 * d[1053]) * dv_294) +
           d[219] * ((105.0 * d[159] + 58.0 * d[48] + d[940]) * dv_0 +
                     Dy * d[2837] + 154.0 * dv_3637) +
           dv_4141;
  sc_7 = (-d[1225]) * dv_4127 +
         (-d[171]) * (-dv_4112 + xpdot * (d[170] * dv_231 + dv_4111)) -
         dv_4113 - dv_4116 + sc_34;
  sc_7 += (-d[237]) *
          (d[122] * dv_4128 + d[50] * ((d[214] * xpdot) * Dx - dv_4131) +
           d[55] * dv_4129 + d[61] * dv_4133 + dv_4134 + dv_4135);
  sc_7 +=
      (-d[2577]) * (d[2821] * dv_4109 - dv_4106 + dv_4108 + dv_4110 * xpdot);
  DataVector& dv_3519 = temps.at(2993);
  DataVector& dv_3636 = temps.at(3106);
  DataVector& dv_4107 = temps.at(262);
  sc_7 += (-d[2583]) *
          (d[1803] * dv_4114 + d[407] * dv_3636 +
           xpdot * ((-d[2702]) * dv_231 +
                    (-2.0 * d[133] - 2.0 * d[1612] - 2.0 * d[72]) * dv_119 +
                    d[137] * dv_3519 + 5.0 * dv_876) +
           ypdot * (d[836] * dv_4107 + dv_4115));
  sc_7 += (-d[260]) * ((-d[227]) * dv_4126 + (-d[2827]) * dv_3912 +
                       (-d[52]) * dv_4123 - dv_4124 - dv_4125);
  sc_7 += (-d[2621]) *
          (d[120] * ((-d[2841]) * dv_0 - dv_4147) +
           d[122] * (d[259] * (Dx * d[2840] - dv_4144) - dv_4143 + dv_4146) +
           dv_4149);
  DataVector& dv_1909 = temps.at(702);
  DataVector& dv_1910 = temps.at(1631);
  DataVector& dv_976 = temps.at(882);
  sc_7 +=
      (-d[2709]) * dv_1909 + (-d[2716]) * dv_1910 +
      (-d[180] * d[19] * (d[383] + d[48] * (-d[112] - d[919])) -
       d[2] * (d[259] * (-d[1550] + 8.0 * d[19] * d[329] - d[429] - d[925]) +
               ypdot * (d[48] * (d[2705] + d[546] + 57.0 * d[55]) + d[489])) +
       d[36] * (d[1833] + 21.0 * d[299] + d[56] * (d[348] + d[48]) -
                d[60] * (d[104] + d[388])) +
       d[540] * d[603] * d[92] -
       d[6] * d[63] * (d[383] + d[48] * (d[241] + d[897]))) *
          dv_977 +
      (-d[1669] * d[549]) * dv_976 +
      (-d[443] *
       (-d[425] * (d[121] + d[423]) +
        xp * (d[436] * (d[371] + d[410]) - d[440] * (-d[19] + 26.0 * d[20]) +
              ypdot * (d[441] + d[91] * (d[2707] - d[442] + d[56]))) -
        xpdot * (d[2706] + d[287] * (d[172] + d[2708] - d[429]) +
                 d[435] * (d[19] * d[430] + d[223] + d[56] + d[866])))) *
          dv_289 +
      (-d[467] * (d[357] * (32.0 * M * d[382] * yp * ypdot - d[447] -
                            d[49] * (d[20] * d[437] + d[229] + d[521])) +
                  xp * (d[180] * (d[453] + d[49] * (-d[339] - d[60] - d[685])) +
                        d[37] * (d[198] * d[48] + d[2852] + 17.0 * d[340] +
                                 d[348] * d[55] + d[547] * d[60]) +
                        d[449] * (-d[112] + d[448])) +
                  xpdot * (d[1606] * d[382] * d[50] +
                           d[321] * (-d[2852] - 305.0 * d[341] - d[455] -
                                     d[456] * d[55]) +
                           d[411] * (d[457] + 32.0 * d[60] + d[685]) + d[460] +
                           d[461] * d[678]))) *
          dv_289 +
      sc_17 + sc_36;
  sc_7 +=
      (4.0 * d[48] * d[71]) *
          (d[55] * ((-d[2603]) * dv_0 + (2.0 * M * d[2830]) * Dy - dv_4137) +
           dv_4140) +
      (4.0 * d[384] * d[420] * yp) * dv_1908 +
      (8.0 * d[151] * d[390] * yp) * dv_1907 +
      (16.0 * d[375] * d[43] * yp) * dv_1906 +
      (4.0 * d[384] * d[420] * xp *
       (-d[397] * (d[145] * d[996] + d[382] * (d[121] + d[393])) +
        xp * (d[411] * (d[410] + 100.0 * d[604]) +
              d[415] * (d[382] * (d[19] - d[414]) + d[412] * d[996]) +
              ypdot * (d[370] * d[418] + d[417] +
                       d[95] * (d[2718] - d[419] + d[55]))) +
        xpdot *
            (d[104] * d[312] * (d[2720] + d[2787] + d[55]) - d[115] * d[406] -
             d[1161] * d[319] * d[399] +
             d[382] * d[403] * (-d[401] - d[458] - d[685]) + d[405] * yp))) *
          dv_6 +
      (8.0 * d[151] * d[390] * xp *
       (d[1199] * (-d[1533] * d[638] + d[2847]) + d[126] * (d[2636] + d[2711]) +
        d[191] * (d[1549] * d[2648] + d[1665] + d[2713] + d[2832] * d[683]) +
        d[2589] * (-d[1334] - 188.0 * d[1339] - d[2653]) + d[2710] +
        d[2846] * d[496] +
        d[50] * d[648] *
            (d[271] - 39.0 * d[276] - 170.0 * d[277] +
             d[724] * (d[1544] + 5.0)))) *
          dv_6;
  sc_7 +=
      (16.0 * d[375] * d[43] * xp *
       (d[206] * (d[1128] * d[261] +
                  d[255] * (d[373] * d[996] + d[92] * (d[200] - d[374])) -
                  d[356] * (d[226] + d[363] - 53.0 * d[57] - d[835])) +
        d[357] * (-d[151] * d[243] * d[353] + 12.0 * d[352] * yp * ypdot -
                  d[356] * (d[229] - d[322] + d[841])) +
        xpdot * (d[1395] * d[352] +
                 d[36] * (d[359] * d[996] + d[92] * (-d[2852] - 763.0 * d[341] -
                                                     253.0 * d[342] - d[360])) +
                 d[366] * d[51]))) *
          dv_6 +
      (16.0 * d[130] * d[20] * d[375] * d[92] * xp) * dv_1912 +
      (384.0 * d[19] * d[4] * d[50] * d[538] * d[551]) * dv_6 +
      (16.0 * d[130] * d[375] * d[92] * xp * yp *
       (d[1130] * d[2851] * xp + d[116] * d[189] * d[479] -
        d[1293] * xp * (d[481] * d[49] + d[810] * d[92]) - d[2703] * d[542] +
        xpdot * (-d[486] * d[62] - d[488] * (d[1471] + d[55]) +
                 2.0 * d[491] * yp * ypdot))) *
          dv_6 +
      (32.0 * d[1133] * d[20] * d[382] * d[40] * d[475] * xp) * dv_6 +
      (128.0 * d[0] * d[1129] * d[20] * d[382] * d[475] * xp) * dv_6;
  sc_21 = (6.0 * d[130] * d[16]) * dv_10 * dv_224 * sc_7;
  DataVector& dv_4101 = temps.at(1140);
  DataVector& dv_4202 = temps.at(3133);
  sc_5 = (-d[1097]) * dv_4097 - dv_1607 * dv_3821 - dv_1844 * dv_4101 -
         dv_1845 * dv_4101 - dv_3766 * dv_4202 - dv_3767 * dv_4202 -
         dv_3812 * dv_4074 + sc_4;
  DataVector& dv_4098 = temps.at(1642);
  DataVector& dv_4099 = temps.at(282);
  sc_5 += -dv_3815 * dv_4065 - dv_3816 * dv_4076 - dv_3828 * dv_4098 -
          dv_3830 * dv_4099 - dv_3831 * dv_4099 - dv_3837 * dv_4096;
  DataVector& dv_4167 = temps.at(356);
  sc_5 += (6.0 * d[1063] * d[130] * d[16]) * dv_224 * dv_3957 -
          dv_3844 * ((3.0 * M * d[12]) * dv_1867 * dv_4096 - dv_3842 * dv_4105 -
                     dv_4103) -
          dv_4002 * dv_4098 - dv_4004 * dv_4167 - dv_4005 * dv_4167;
  DataVector& dv_4095 = temps.at(1168);
  sc_5 += (6.0 * d[0] * d[1063]) * dv_11 * dv_229 * dv_4012 +
          (8.0 * d[0]) * dv_1496 * dv_1607 * dv_19 * dv_4012 +
          (4.0 * M * d[16] * d[74]) * dv_229 * dv_3820 * dv_4096 +
          (8.0 * M * d[16] * d[74]) * dv_1867 * dv_229 * dv_4095 +
          (12.0 * M * d[75]) * dv_11 * dv_19 * dv_3820 * dv_4096 + sc_11 +
          sc_21;
  sc_5 += (24.0 * M * d[75]) * dv_11 * dv_1607 * dv_1867 * dv_3820 +
          (24.0 * M * d[75]) * dv_11 * dv_1867 * dv_19 * dv_4095 +
          (36.0 * d[130] * d[16]) * dv_10 * dv_1607 * dv_229 * dv_3957 +
          (16.0 * M * d[16] * d[74]) * dv_1607 * dv_1867 * dv_19 * dv_3820 +
          (24.0 * M * d[1063] * d[75]) * dv_10 * dv_1867 * dv_19 * dv_3820;
  sc_1 = ((1.0 / 48.0) * M * d[1024] * d[1025]) * dv_1611 * sc_5;
  sc_3 = ((1.0 / 2.0) * M * d[23]) *
             ((-d[1063]) * dv_3773 + dv_4066 + dv_43 * yp) +
         ((7.0 / 48.0) * M * d[1024] * d[1025]) * dv_1607 * dv_225 * dv_3771 -
         dv_1607 * dv_3768 - dv_1607 * dv_3769 - dv_1607 + sc_0 + sc_1;
  get<1>(get<::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>,
                           Frame::Inertial>>(*result)) = dv_1499 * sc_3;
  sc_11 = (2.0 * d[17]) * dv_10;
  DataVector& dv_1527 = temps.at(1404);
  DataVector& dv_1626 = temps.at(1471);
  DataVector& dv_1984 = temps.at(157);
  sc_11 *= d[356] * dv_1626 + d[65] * dv_6 + d[66] * dv_6 +
           d[93] * ((d[395] * xp) * dv_790 + d[19] * (dv_1527 + dv_1984) +
                    d[20] * ((5.0 * ypdot) * Dy - dv_0));
  DataVector& dv_33 = temps.at(33);
  DataVector& dv_42 = temps.at(42);
  DataVector& dv_4205 = temps.at(185);
  sc_21 = (-d[17] * d[2870]) * (d[4] * dv_33 + dv_42) - dv_4205 + sc_11;
  sc_5 = d[69] * sc_21;
  DataVector& dv_1622 = temps.at(1468);
  DataVector& dv_4203 = temps.at(1335);
  sc_0 = (d[1474] * d[92]) * dv_1622 * dv_4203 + d[2869] * dv_79 +
         d[2869] * dv_86 + d[42] * dv_1611 * dv_1616 - dv_3785 + sc_5;
  sc_1 = (-d[2572]) * sc_0;
  DataVector& dv_1717 = temps.at(1538);
  DataVector& dv_3497 = temps.at(2971);
  DataVector& dv_708 = temps.at(658);
  sc_12 =
      d[19] * (d[29] * (-62.0 * dv_1717 + ypdot * (dv_131 - dv_3497 + dv_708)) +
               d[9] * dv_2024);
  DataVector& dv_2006 = temps.at(1521);
  DataVector& dv_2031 = temps.at(1662);
  DataVector& dv_2076 = temps.at(1702);
  DataVector& dv_2917 = temps.at(2456);
  DataVector& dv_703 = temps.at(653);
  DataVector& dv_801 = temps.at(749);
  sc_12 += d[20] * (d[29] * (dv_2076 * ypdot + dv_2917) + d[9] * dv_2031) +
           d[206] * ((-d[1578] - 93.0 * d[3]) * dv_801 +
                     (3.0 * xpdot * yp) * (dv_131 - dv_2006 + dv_703));
  DataVector& dv_3590 = temps.at(3063);
  sc_12 += d[2570] * (dv_2086 * xpdot + dv_3590);
  sc_34 = (-d[12]) * sc_12;
  DataVector& dv_114 = temps.at(112);
  DataVector& dv_1835 = temps.at(1620);
  DataVector& dv_2140 = temps.at(1260);
  DataVector& dv_4159 = temps.at(1652);
  DataVector& dv_679 = temps.at(629);
  sc_28 = (-d[52]) * (dv_1835 + dv_2039 * xpdot) +
          d[19] * (dv_2181 +
                   yp * (Dy * dv_2140 + ypdot * (-dv_114 - dv_4159 - dv_679)));
  DataVector& dv_1818 = temps.at(1604);
  DataVector& dv_2062 = temps.at(154);
  DataVector& dv_2965 = temps.at(2493);
  DataVector& dv_3968 = temps.at(397);
  DataVector& dv_691 = temps.at(641);
  sc_28 += d[20] * ((-yp) * (dv_1818 + dv_2062 * ypdot) + dv_2965) +
           d[206] * ((d[166] + d[1823]) * dv_801 +
                     d[86] * (-dv_114 - dv_3968 - dv_691));
  sc_12 = d[114] * sc_28;
  DataVector& dv_4078 = temps.at(3142);
  DataVector& dv_542 = temps.at(496);
  sc_17 = (-24.0 * d[0]) * (d[19] * dv_2019 + d[20] * dv_1651 - dv_542) +
          12.0 * dv_4078 + sc_12 + sc_34;
  sc_36 = d[92] * dv_10 * sc_17;
  DataVector& dv_1796 = temps.at(177);
  DataVector& dv_4215 = temps.at(3126);
  DataVector& dv_81 = temps.at(79);
  sc_7 = (d[2868] * d[323]) * dv_1796 + d[10] * dv_4215 * dv_81 + sc_36;
  DataVector& dv_157 = temps.at(155);
  DataVector& dv_174 = temps.at(172);
  DataVector& dv_2012 = temps.at(106);
  DataVector& dv_3484 = temps.at(2959);
  DataVector& dv_375 = temps.at(358);
  DataVector& dv_4204 = temps.at(179);
  sc_7 +=
      d[92] * dv_375 *
      (d[114] * (d[19] * dv_125 + d[20] * dv_157 - dv_3484 * dv_4) +
       d[968] * (d[19] * dv_2004 + d[20] * dv_2012 - dv_4204) - 4.0 * dv_174);
  sc_4 = -dv_1864 * sc_7;
  DataVector& dv_1741 = temps.at(134);
  DataVector& dv_1791 = temps.at(1581);
  DataVector& dv_1797 = temps.at(86);
  DataVector& dv_1863 = temps.at(211);
  DataVector& dv_4207 = temps.at(1486);
  DataVector& dv_4209 = temps.at(923);
  DataVector& dv_4216 = temps.at(50);
  sc_11 = (-d[1096] * d[12]) * dv_4209 + (-4.0 * d[107]) * dv_1791 * dv_1863 +
          (M * d[0] * d[16]) * dv_1741 * dv_1797 +
          (2.0 * M * d[92]) * dv_1791 * dv_1797 * dv_4207 +
          (4.0 * M * rp) * dv_1741 * dv_1791 * dv_4216 + sc_4;
  sc_21 = d[108] * sc_11;
  DataVector& dv_1497 = temps.at(1377);
  DataVector& dv_1748 = temps.at(1553);
  DataVector& dv_1753 = temps.at(1557);
  DataVector& dv_1764 = temps.at(1566);
  DataVector& dv_1765 = temps.at(1567);
  DataVector& dv_1789 = temps.at(1579);
  DataVector& dv_1915 = temps.at(867);
  DataVector& dv_1918 = temps.at(1638);
  DataVector& dv_1922 = temps.at(1641);
  DataVector& dv_3770 = temps.at(1569);
  DataVector& dv_4212 = temps.at(3132);
  sc_5 =
      (-d[2869]) * dv_1764 + (-36.0 * d[1003]) * dv_4209 +
      (d[2912] *
       (d[2428] * (d[19] * (d[2896] - d[80] * (d[135] * d[7] + d[2888])) +
                   d[20] * (d[2897] + ypdot * (d[1473] - d[2898] * yp)) -
                   d[2894] * d[609] - d[2895]) +
        d[2882] * (d[2881] + ypdot * (-d[2334] * yp + d[314])) + d[2911] +
        d[43] * (-d[19] * (-d[2820] * d[284] + d[6] * (d[1361] + d[88])) -
                 d[2479] * d[508] +
                 xp * xpdot * (d[2640] * d[280] - ypdot * (d[1876] + d[2871])) +
                 yp * (d[286] * (d[1358] + d[2872] + d[77]) -
                       d[7] * (d[1322] + d[2873]))) +
        d[465] * (-d[286] * d[98] + d[2874] * xp * xpdot - d[2875] * ypdot) +
        d[845] * (-d[19] * (d[2878] + d[2879] * d[7]) - d[2] * d[2880] +
                  d[2877] * d[52] * xpdot -
                  yp * (d[2879] * d[35] +
                        ypdot * (d[162] + d[213] - d[430] * d[7]))) +
        d[873] * (-d[2899] * d[2900] + d[2903]) +
        d[899] * (-d[2883] * d[2884] + d[2885] * d[50] * xp * xpdot +
                  d[2886] * d[52] * xpdot * yp - d[2892] * d[60] +
                  d[55] * (-d[2887] * d[6] + 2.0 * d[2889] * ypdot) -
                  d[57] * (d[2890] + ypdot * (d[1401] - d[2891] * yp))))) *
          dv_1497 +
      (d[558] * d[83]) * dv_1789 +
      (6.0 * d[16] * d[75] *
       (d[2428] * (d[19] * (d[154] * d[88] - d[261] * d[2904] + d[2896]) +
                   d[20] * (d[1616] + d[2897] - d[2898] * d[3] + d[88]) -
                   d[2894] * d[609] - d[2895]) +
        d[2882] * (d[1460] + d[2675] + d[2881] - d[2913] + d[9]) + d[2911] +
        d[43] *
            (4.0 * M * d[412] * d[7] - 5.0 * d[1116] - d[1117] * d[2913] -
             d[391] + d[6] * (d[150] * d[1721] + d[2736] * d[9]) +
             xp * xpdot * ypdot * (-d[2871] + 10.0 * ypdot * (d[19] + d[48]))) +
        d[465] * (d[1281] * d[286] + d[2] * d[2874] - d[2875] * ypdot) +
        d[845] * (-d[19] * (-d[148] + d[1955] + d[2878] + d[9]) -
                  d[2] * d[2880] + d[2877] * d[609] +
                  yp * (-d[1295] + d[20] * ypdot * (d[1580] + d[6]) -
                        d[77] * (d[1868] + d[2247] + 2.0))) +
        d[873] * (-d[2899] * d[2900] + d[2903]) +
        d[899] *
            (-d[2883] * d[2884] + d[2885] * d[716] + d[2886] * d[601] -
             d[2892] * d[60] + d[55] * (-d[2887] * d[6] + d[2889] * d[80]) -
             d[57] * (d[2890] + ypdot * (d[1401] - d[2891] * yp))))) *
          dv_1915 +
      dv_1748 * dv_4212 + dv_1753 * dv_1765 + 18.0 * dv_1918 + 6.0 * dv_1922 -
      dv_3770 + sc_21;
  DataVector& dv_1746 = temps.at(1551);
  DataVector& dv_1755 = temps.at(1559);
  DataVector& dv_1768 = temps.at(1570);
  DataVector& dv_1777 = temps.at(231);
  DataVector& dv_1785 = temps.at(1555);
  DataVector& dv_1786 = temps.at(229);
  DataVector& dv_1974 = temps.at(1673);
  DataVector& dv_1976 = temps.at(498);
  DataVector& dv_250 = temps.at(243);
  DataVector& dv_307 = temps.at(299);
  DataVector& dv_4210 = temps.at(1656);
  DataVector& dv_4211 = temps.at(1068);
  DataVector& dv_4217 = temps.at(387);
  sc_5 += -dv_1746 * dv_1786 - dv_1755 * dv_1785 - dv_1768 * dv_4211 -
          dv_1777 * dv_4211 + dv_1974 * dv_4217 - dv_1976 * dv_4217 +
          dv_250 * dv_4210 + dv_307 * dv_4210;
  DataVector& dv_1778 = temps.at(1568);
  DataVector& dv_1799 = temps.at(1585);
  DataVector& dv_1874 = temps.at(49);
  DataVector& dv_306 = temps.at(298);
  DataVector& dv_4206 = temps.at(1500);
  DataVector& dv_4214 = temps.at(1586);
  sc_5 += (-d[1] * d[207]) * dv_1741 *
              ((-d[1110]) * dv_1874 + (-d[1088] * d[93]) * dv_1867 * dv_4206 -
               4.0 * dv_1791 * dv_4216) +
          (-d[108] * d[323]) * dv_1799 * dv_4207 + d[2912] * dv_306 * dv_4214 +
          d[75] * dv_3766 * dv_4214 + dv_1778 * dv_19 * dv_4212;
  DataVector& dv_1747 = temps.at(1552);
  sc_5 += (d[2868] * d[92]) * dv_10 * dv_1746 * dv_1747;
  sc_0 = (-d[1026] * d[70]) * dv_1611 * sc_5;
  sc_3 = (-d[183] * d[23]) * dv_4203 +
         ((7.0 / 48.0) * M * d[1024] * d[1025]) * dv_225 * dv_3771 - dv_3768 -
         dv_3769 + sc_0 + sc_1 - 1.0;
  get<2>(get<::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>,
                           Frame::Inertial>>(*result)) = dv_1499 * sc_3 * z;
}
}  // namespace CurvedScalarWave::Worldtube::detail
