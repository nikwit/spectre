
// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/Worldtube/PunctureFieldOrder2Impl.hpp"

namespace CurvedScalarWave::Worldtube::detail {

// NOLINTNEXTLINE(google-readability-function-size, readability-function-size)
void puncture_field_2_part_10(const std::array<double, order2_n_doubles>& d,
                              const DataVector& Dx, const DataVector& Dy,
                              DynamicBuffer<DataVector>& temps,
                              const gsl::not_null<Order2Vars*> result) {
  DataVector& sc_17 = temps.at(3239);
  DataVector& sc_6 = temps.at(3228);
  sc_17 = d[28] * sc_6;
  DataVector& dv_2872 = temps.at(2421);
  DataVector& sc_31 = temps.at(3253);
  sc_31 = -dv_2872;
  DataVector& dv_2551 = temps.at(2139);
  DataVector& dv_2873 = temps.at(2422);
  DataVector& dv_2877 = temps.at(2426);
  DataVector& dv_3068 = temps.at(2561);
  DataVector& dv_3238 = temps.at(2725);
  DataVector& dv_3239 = temps.at(2726);
  DataVector& dv_794 = temps.at(742);
  sc_31 +=
      Dx * ((-d[2928]) * ((-d[2928]) * (dv_2551 - 16.0) + dv_2877 + dv_3239) +
            dv_2873 +
            d[2921] * ((-d[36]) * (dv_3238 - 127.0) + d[2030] +
                       884.0 * dv_3068 + 33.0 * dv_794));
  DataVector& sc_12 = temps.at(3234);
  sc_12 = (-d[1250]) * sc_31;
  DataVector& dv_1229 = temps.at(1124);
  DataVector& dv_2137 = temps.at(1753);
  DataVector& dv_2883 = temps.at(2432);
  DataVector& dv_2886 = temps.at(703);
  DataVector& dv_545 = temps.at(499);
  DataVector& dv_711 = temps.at(661);
  DataVector& dv_789 = temps.at(737);
  DataVector& sc_18 = temps.at(3240);
  sc_18 = d[1] * (d[1048] * dv_711 + d[7] * (d[1049] * dv_545 + dv_789) -
                  dv_1229 * dv_2883 + dv_2137) +
          dv_2886 + sc_12;
  DataVector& dv_2220 = temps.at(1822);
  DataVector& dv_2571 = temps.at(652);
  DataVector& dv_2870 = temps.at(2420);
  DataVector& dv_2882 = temps.at(2431);
  DataVector& dv_3241 = temps.at(2728);
  DataVector& dv_696 = temps.at(646);
  DataVector& dv_811 = temps.at(757);
  sc_18 += d[27] * ((-d[2031]) * (dv_696 + dv_811) + d[1723] * dv_2220 +
                    dv_2882 - dv_3241) +
           d[6] * (d[1048] * dv_2571 + dv_2870);
  sc_6 = d[55] * sc_18;
  DataVector& dv_1903 = temps.at(702);
  DataVector& dv_2865 = temps.at(2417);
  DataVector& dv_3134 = temps.at(1861);
  DataVector& dv_3237 = temps.at(2724);
  DataVector& dv_678 = temps.at(628);
  DataVector& sc_13 = temps.at(3235);
  DataVector& sc_16 = temps.at(3238);
  DataVector& sc_24 = temps.at(3246);
  DataVector& sc_29 = temps.at(3251);
  sc_13 =
      (-d[122]) * (d[1048] * ((-d[2922]) * dv_678 + 36.0 * dv_1903 + dv_3237) +
                   dv_2865) +
      d[300] * dv_3134 + sc_16 + sc_17 + sc_24 + sc_29 + sc_6;
  DataVector& sc_28 = temps.at(3250);
  sc_28 = d[301] * sc_13;
  DataVector& dv_0 = temps.at(0);
  DataVector& dv_1478 = temps.at(1359);
  DataVector& dv_1963 = temps.at(1679);
  DataVector& dv_2212 = temps.at(1815);
  DataVector& dv_2586 = temps.at(2172);
  DataVector& dv_263 = temps.at(256);
  DataVector& dv_3060 = temps.at(2553);
  DataVector& dv_46 = temps.at(46);
  DataVector& dv_48 = temps.at(48);
  DataVector& dv_5 = temps.at(5);
  DataVector& dv_700 = temps.at(650);
  sc_6 = (-d[1030]) * dv_1963 + (-d[1525]) * dv_5 + (-d[1576]) * dv_2212 +
         (-d[1724]) * dv_46 + (-d[1724]) * dv_48 + (-d[1786]) * dv_263 +
         (-d[1788]) * dv_1478 + (-d[1789]) * dv_0 + (-d[1790]) * dv_2586 +
         (-d[1791]) * dv_700 + dv_3060;
  DataVector& dv_106 = temps.at(104);
  DataVector& dv_1200 = temps.at(1095);
  DataVector& dv_2141 = temps.at(1258);
  DataVector& dv_2142 = temps.at(740);
  DataVector& dv_297 = temps.at(289);
  DataVector& dv_3061 = temps.at(2554);
  DataVector& dv_615 = temps.at(566);
  DataVector& dv_756 = temps.at(704);
  DataVector& dv_839 = temps.at(782);
  sc_6 += (-d[1792]) * dv_3061 + (d[2928] * d[1307]) * dv_1200 +
          (d[1033] * d[1053]) * dv_106 + (d[1088] * d[6]) * dv_46 +
          (d[255] * d[6]) * dv_615 + (-d[2928] * d[1033]) * dv_297 +
          (-d[1033] * d[1468]) * dv_756 + (-d[1795] * d[3]) * dv_839 +
          d[2928] * dv_2141 + d[2928] * dv_2142;
  DataVector& dv_2143 = temps.at(1756);
  DataVector& dv_2145 = temps.at(1758);
  DataVector& dv_2146 = temps.at(1759);
  DataVector& dv_2147 = temps.at(1760);
  DataVector& dv_2148 = temps.at(1761);
  DataVector& dv_3058 = temps.at(2551);
  DataVector& dv_3059 = temps.at(2552);
  DataVector& dv_567 = temps.at(520);
  DataVector& dv_608 = temps.at(559);
  sc_6 += d[2928] * dv_2143 + d[2928] * dv_2145 + d[2928] * dv_2146 +
          d[1033] * dv_3058 + d[1053] * dv_567 + d[1053] * dv_608 +
          d[1242] * dv_2148 + d[1307] * dv_3059 + d[1406] * dv_48 +
          d[1474] * dv_2147;
  DataVector& dv_1545 = temps.at(1419);
  DataVector& dv_2191 = temps.at(1797);
  DataVector& dv_3038 = temps.at(446);
  DataVector& dv_503 = temps.at(464);
  DataVector& dv_602 = temps.at(553);
  DataVector& dv_810 = temps.at(756);
  sc_6 += d[1786] * dv_2212 + d[1787] * dv_756 + d[1790] * dv_810 +
          d[1791] * dv_1545 + d[1793] * dv_3038 + d[1794] * dv_503 +
          d[1794] * dv_615 + d[255] * dv_567 + d[255] * dv_602 +
          d[266] * dv_2191;
  DataVector& dv_2154 = temps.at(546);
  DataVector& dv_2160 = temps.at(549);
  DataVector& dv_2443 = temps.at(2033);
  DataVector& dv_2950 = temps.at(2478);
  DataVector& dv_3062 = temps.at(2555);
  DataVector& dv_548 = temps.at(502);
  DataVector& dv_555 = temps.at(508);
  DataVector& dv_732 = temps.at(681);
  DataVector& dv_833 = temps.at(776);
  sc_6 += (-d[1447]) * dv_0 * dv_833 - dv_2443 * dv_2950 +
          d[2920] * ((-d[2928]) * dv_2160 +
                     (-d[1217]) * (d[1248] * dv_732 + dv_3062 - dv_548 +
                                   dv_555 * d[2922]) +
                     (d[2928] * d[2924]) * dv_2154);
  DataVector& dv_788 = temps.at(736);
  sc_6 += (-d[1788]) * dv_5 * dv_788;
  sc_13 = d[83] * sc_6;
  DataVector& dv_1842 = temps.at(1626);
  DataVector& dv_2136 = temps.at(1752);
  DataVector& dv_3074 = temps.at(2567);
  DataVector& dv_69 = temps.at(69);
  sc_24 = (-d[416] * d[2927]) * dv_1842 +
          d[49] * ((d[1805] + d[298] + 36.0) * dv_0 +
                   (-d[1778] - d[284]) * dv_69 + 48.0 * dv_2136 + dv_3074);
  DataVector& dv_1115 = temps.at(1012);
  DataVector& dv_2739 = temps.at(2292);
  DataVector& dv_3075 = temps.at(2568);
  DataVector& dv_3080 = temps.at(2573);
  DataVector& dv_3081 = temps.at(2574);
  DataVector& dv_3082 = temps.at(2575);
  DataVector& dv_3084 = temps.at(2577);
  sc_24 +=
      d[77] * ((d[1] * d[147]) + d[1070] * dv_1115 +
               d[648] * (dv_2739 + dv_3081) - dv_3075 - dv_3080 +
               d[2923] * ((-d[2928]) * dv_3084 + (-d[1668] * d[6]) + dv_3082));
  sc_16 = (-d[1807]) * sc_24;
  DataVector& dv_1817 = temps.at(1603);
  DataVector& dv_2170 = temps.at(1778);
  DataVector& dv_3067 = temps.at(2560);
  DataVector& dv_3069 = temps.at(2562);
  DataVector& dv_3087 = temps.at(2580);
  DataVector& dv_3088 = temps.at(2581);
  DataVector& dv_3092 = temps.at(2585);
  sc_29 =
      (-d[1819] - d[567]) * dv_3067 +
      (-d[2922]) * ((-d[1808] - d[1820] - 36.0) * dv_1817 +
                    d[77] * (d[2928] * dv_3092 + d[1667] - dv_3088) + dv_3069) +
      (d[1248] * (d[1519] - d[1792])) * dv_2443 +
      d[321] * ((d[1368] + d[1739] + 15.0) * dv_2170 + dv_3087);
  DataVector& dv_2267 = temps.at(1863);
  sc_29 += d[91] * dv_2267;
  sc_24 = (-d[1821]) * sc_29;
  DataVector& dv_3073 = temps.at(2566);
  sc_18 = d[416] * dv_3073 + d[49] * ((-d[1036] - d[1801]) * dv_0 +
                                      (-d[169] + d[1802] - 18.0) * dv_69 +
                                      (4.0 * d[142]) * Dx - dv_3074);
  DataVector& dv_2344 = temps.at(1934);
  DataVector& dv_3076 = temps.at(2569);
  DataVector& dv_3077 = temps.at(2570);
  DataVector& dv_3078 = temps.at(2571);
  sc_18 += d[77] *
           ((-d[7]) * dv_3077 + d[1579] + d[648] * dv_2344 + dv_3075 + dv_3076 +
            d[2923] * ((-d[1803]) + d[2928] * dv_3078 + d[1029] * dv_0));
  sc_29 = d[1804] * sc_18;
  DataVector& dv_1541 = temps.at(1415);
  sc_12 = (-d[49]) * ((d[1808] + d[1809] + 48.0) * dv_0 +
                      (-d[1070] + d[1620] + 24.0) * dv_69 + 101.0 * dv_1541 +
                      14.0 * dv_2136);
  DataVector& dv_1828 = temps.at(1613);
  DataVector& dv_3085 = temps.at(2578);
  DataVector& dv_558 = temps.at(511);
  sc_12 += (-d[77]) * ((-d[1306]) * Dy +
                       (-d[2923]) * ((-d[2928]) * dv_3085 + d[1811] - dv_3082) +
                       d[1810] + d[648] * (dv_3081 + dv_558 * d[2924]) +
                       dv_3076 + dv_3080) +
           (27.0 * d[20] * d[2927]) * dv_1828;
  sc_18 = d[1812] * sc_12;
  DataVector& dv_3086 = temps.at(2579);
  DataVector& dv_787 = temps.at(735);
  sc_31 = (-d[252]) * ((d[1188] + d[152] + 18.0) * dv_2170 + dv_3087) +
          (-d[1796] - 109.0 * d[36]) * dv_3086 +
          (-d[142] * d[49]) * (d[30] + dv_787) +
          (3.0 * d[2927] * (d[1682] - d[1813] * d[259] + d[563])) * Dx;
  DataVector& dv_1821 = temps.at(1607);
  DataVector& dv_3091 = temps.at(2584);
  sc_31 += d[2922] * ((-d[1805] - d[1814] - 30.0) * dv_1821 +
                      d[259] * (-dv_3088 - dv_3091) + dv_3069);
  sc_12 = d[1816] * sc_31;
  DataVector& dv_1 = temps.at(1);
  DataVector& dv_1762 = temps.at(1565);
  DataVector& dv_2331 = temps.at(1922);
  DataVector& dv_3063 = temps.at(2556);
  DataVector& dv_3064 = temps.at(2557);
  sc_17 =
      (d[1] * d[1797]) * (d[1330] * dv_1 + d[142] * dv_1762 + d[36] * dv_3063 +
                          d[2922] * ((d[1035] + d[298] - 6.0) * dv_833 +
                                     d[2921] * (d[1054] + dv_2331 + dv_3064))) +
      (d[1618] * d[49]) * dv_1541 + sc_16 + sc_24;
  DataVector& dv_1801 = temps.at(1587);
  sc_17 += d[1466] * ((d[1037] * d[49] + d[1796] * d[554]) * dv_0 +
                      (-d[6] * d[630]) + d[49] * dv_1801);
  DataVector& dv_3065 = temps.at(2558);
  DataVector& dv_3066 = temps.at(2559);
  DataVector& dv_3072 = temps.at(2565);
  sc_17 +=
      d[1800] *
          ((-d[1798] + 5.0 * d[36]) * dv_3065 + (d[1449] * d[221]) +
           d[1799] * dv_3067 + d[321] * ((d[1037] - 3.0) * dv_2170 + dv_3066) +
           d[2922] * (d[72] * dv_794 + d[77] * (d[1509] + dv_3064 + dv_3072) +
                      dv_3069)) +
      sc_29;
  sc_17 += sc_12 + sc_18;
  DataVector& dv_3094 = temps.at(2587);
  sc_6 = (-d[2926]) * dv_3094 * sc_17;
  DataVector& dv_2787 = temps.at(2340);
  DataVector& dv_3120 = temps.at(2611);
  sc_29 = (-d[151]) * dv_3120 +
          (-d[1882]) * ((-254.0 * d[20]) + d[166] * ((-d[1418]) - dv_2787) +
                        d[29] * ((d[1079] + 159.0) * Dy + d[1881]));
  DataVector& dv_1502 = temps.at(1382);
  DataVector& dv_2445 = temps.at(2035);
  DataVector& dv_3121 = temps.at(2612);
  sc_29 +=
      (-d[50]) * ((-d[1879]) + d[1704] * dv_1502 + 95.0 * dv_2331 + dv_3121) +
      (360.0 * d[142]) * dv_2445;
  DataVector& dv_2421 = temps.at(2011);
  DataVector& dv_2767 = temps.at(2320);
  DataVector& dv_3123 = temps.at(2614);
  DataVector& dv_751 = temps.at(699);
  DataVector& dv_850 = temps.at(789);
  sc_29 += d[1250] * ((d[1809] - 90.0) * dv_2421 + (120.0 * d[2927]) * dv_850 +
                      (d[1348] * d[1869]) * dv_751 +
                      d[273] * ((-d[1886]) * dv_751 + dv_2767 + dv_3123));
  DataVector& dv_2246 = temps.at(1844);
  DataVector& dv_2431 = temps.at(2021);
  DataVector& dv_2766 = temps.at(2319);
  DataVector& dv_3101 = temps.at(2593);
  DataVector& dv_3122 = temps.at(2613);
  sc_29 += d[273] * ((d[1884] - 517.0) * dv_794 +
                     d[1236] * (dv_2431 + dv_2766 + 6.0) + d[1883] + dv_3122) +
           d[287] * ((-d[1880] - 27.0) * dv_69 +
                     (3.0 * d[2928]) * (dv_1502 + dv_3101) + (-d[1671]) -
                     18.0 * dv_2246);
  sc_18 = (-d[1841]) * sc_29;
  DataVector& dv_126 = temps.at(124);
  DataVector& dv_1503 = temps.at(1383);
  DataVector& dv_240 = temps.at(236);
  DataVector& dv_2704 = temps.at(2257);
  DataVector& dv_2782 = temps.at(2335);
  sc_24 =
      (d[1868] * d[514]) * Dx +
      d[1250] * ((-d[1869]) * dv_2704 +
                 (-d[259]) * ((-d[1870] - 27.0) * dv_240 + d[1362] + dv_2782) +
                 (-d[1569]) +
                 d[20] * ((d[1871] + 179.0) * dv_1 + (-d[1221]) * dv_1503 +
                          (42.0 * d[255]))) +
      d[1268] * dv_126;
  DataVector& dv_2606 = temps.at(2188);
  sc_24 += d[1474] * dv_2606;
  sc_29 = (d[2928] * d[1872]) * sc_24;
  DataVector& dv_2666 = temps.at(2229);
  DataVector& dv_3124 = temps.at(2615);
  DataVector& dv_3125 = temps.at(2616);
  sc_16 = (d[1707] + d[1802] * d[532]) * dv_3125 +
          (-d[1896]) * ((d[1802] - 2.0) * dv_2666 + d[1563] +
                        d[20] * ((d[1273] - 109.0) * Dy + d[1895]) +
                        d[77] * (d[9] + 101.0 * dv_1)) +
          dv_3124;
  DataVector& dv_2459 = temps.at(2048);
  DataVector& dv_2762 = temps.at(2315);
  DataVector& dv_2772 = temps.at(2325);
  DataVector& dv_2808 = temps.at(2361);
  DataVector& dv_3083 = temps.at(2576);
  DataVector& dv_3100 = temps.at(2087);
  sc_16 += (-d[273]) * ((-d[1873]) * dv_794 + (-72.0 * d[1862]) +
                        d[403] * (dv_2459 + dv_3083 + 18.0) + dv_3122) +
           d[157] * (dv_2762 + dv_2772 + dv_2808 + dv_3100);
  DataVector& dv_1803 = temps.at(1589);
  DataVector& dv_3105 = temps.at(2597);
  DataVector& dv_3128 = temps.at(2619);
  sc_16 +=
      d[50] *
          ((-106.0 * d[255]) + d[1221] * dv_1502 + 53.0 * dv_2331 + dv_3121) +
      d[94] * ((d[1081] - 261.0 * d[7] + 108.0) * dv_3105 +
               (-d[1897] * d[72]) * dv_751 +
               d[20] * ((-d[1221]) * dv_1803 + d[1878] * dv_751 + dv_3128));
  sc_24 = d[1845] * sc_16;
  DataVector& dv_1816 = temps.at(1602);
  DataVector& dv_3106 = temps.at(2598);
  DataVector& dv_793 = temps.at(741);
  sc_31 =
      (d[1870] + d[1887] + 189.0) * dv_1816 +
      (-d[273]) * ((-540.0 * d[2927] - 1037.0) * dv_794 +
                   d[1236] * (dv_2431 + dv_3106) + d[1888] + dv_793 * d[2927]) +
      (-d[1209] * d[1802] - d[507]) * dv_3125 - dv_3124;
  DataVector& dv_2385 = temps.at(1975);
  DataVector& dv_3126 = temps.at(2617);
  DataVector& dv_538 = temps.at(493);
  sc_31 += d[1053] *
           ((d[1273] - 7.0) * dv_2385 + (-190.0 * d[50]) +
            d[185] * ((d[1273] + 29.0) * dv_538 + d[1889]) - 1734.0 * dv_3126);
  DataVector& dv_3127 = temps.at(2618);
  sc_31 += d[1250] *
           ((d[1081] - 867.0 * d[7] + 378.0) * dv_2445 + (-d[1884]) * dv_850 +
            (-d[273]) * ((-d[1891]) * dv_751 + dv_3123 + dv_3128) +
            d[1890] * dv_3127);
  sc_31 += d[50] * ((-254.0 * d[255]) + d[1704] * dv_1503 + 95.0 * dv_2772 -
                    240.0 * dv_3073);
  sc_16 = d[1852] * sc_31;
  DataVector& dv_3129 = temps.at(2620);
  DataVector& sc_32 = temps.at(3254);
  sc_32 =
      (-1110.0 * d[159] + d[162] - d[1884] * d[280] + 517.0 * d[20]) * dv_3086 +
      (d[1489] * d[631] - d[1710] + d[1894]) * dv_3129;
  DataVector& dv_2371 = temps.at(1961);
  DataVector& dv_2437 = temps.at(2027);
  DataVector& dv_2495 = temps.at(2083);
  DataVector& dv_3117 = temps.at(2608);
  DataVector& dv_3130 = temps.at(2621);
  sc_32 += (-d[2928]) *
           ((d[1242] + d[1447] - 36.0 * d[2923]) * dv_3117 +
            (-d[1682]) * ((53.0 * d[2923]) * Dx - dv_2437 - dv_2495) +
            (-d[50]) * (53.0 * dv_1803 + dv_3130 * d[2924]) + 84.0 * dv_2371);
  DataVector& dv_1769 = temps.at(1571);
  sc_32 += (-d[1893]) * dv_1769;
  DataVector& dv_2477 = temps.at(2066);
  DataVector& dv_3114 = temps.at(2049);
  DataVector& dv_3119 = temps.at(2610);
  DataVector& dv_624 = temps.at(575);
  sc_32 += d[1250] *
           ((-d[1409]) * dv_3119 + (d[1161] * d[1890]) * dv_1 +
            d[101] * ((d[1100] + 27.0) * dv_624 + d[1362] - 1083.0 * dv_794) +
            d[273] * ((-d[42]) * (dv_2477 + dv_3114) + (270.0 * d[255]) +
                      d[1886] * dv_1));
  sc_31 = d[1863] * sc_32;
  DataVector& sc_34 = temps.at(3256);
  sc_34 = (2.0 * d[2922]);
  DataVector& dv_2997 = temps.at(2522);
  DataVector& dv_3109 = temps.at(2601);
  DataVector& dv_3132 = temps.at(2623);
  DataVector& dv_3133 = temps.at(2624);
  sc_34 *= (-d[1897]) * dv_3132 + (-d[51]) * ((127.0 * d[36]) + Dy * d[1884]) +
           (d[2928] * d[20]) * ((-d[42]) * (dv_2477 + dv_3109) +
                                (252.0 * d[255]) + d[1891] * dv_1) +
           (3.0 * d[48] * d[2921]) * ((-d[7]) * dv_2997 + d[37] + dv_3133);
  DataVector& dv_728 = temps.at(677);
  DataVector& sc_33 = temps.at(3255);
  sc_33 = (-d[1893]) * (d[262] + dv_728) +
          (-d[1441] - 2166.0 * d[159] + d[1884] * d[292] + 1037.0 * d[20]) *
              dv_3086 +
          (-d[1452] - d[1899] - d[1901] + d[1903]) * dv_3129 + sc_34;
  DataVector& dv_2376 = temps.at(1966);
  DataVector& dv_2397 = temps.at(1987);
  DataVector& dv_2784 = temps.at(2337);
  DataVector& dv_2826 = temps.at(2378);
  DataVector& dv_3131 = temps.at(2622);
  sc_33 += d[2928] * ((-d[1299]) * (dv_2376 + dv_3131 - 45.0 * dv_751) +
                      (-d[1682]) * (dv_2397 + dv_2437 - 87.0 * dv_751) +
                      d[1904] * dv_2784 + dv_2826);
  sc_32 = d[1905] * sc_33;
  DataVector& dv_231 = temps.at(228);
  DataVector& dv_2469 = temps.at(2058);
  DataVector& dv_2792 = temps.at(2345);
  DataVector& dv_3089 = temps.at(2582);
  DataVector& dv_3116 = temps.at(2607);
  sc_34 = (-d[1570]) * ((-d[259]) * ((-d[46]) * (dv_2792 + dv_3089 + 18.0) +
                                     d[1576] + d[1878] * dv_1) +
                        d[1287] * (d[9] + dv_2469) + d[243] * dv_3119) +
          (d[1873] * d[2921] - 534.0 * d[36]) * dv_3116 +
          (-d[1273] * (d[1874] * d[259] + d[1877])) * dv_231 +
          (72.0 * d[142] * d[20] * d[48]);
  DataVector& dv_1112 = temps.at(1009);
  DataVector& dv_1805 = temps.at(1591);
  DataVector& dv_2396 = temps.at(1986);
  DataVector& dv_3118 = temps.at(2609);
  sc_34 +=
      d[2928] * ((d[1230] + d[1650]) * dv_3117 + d[138] * dv_1112 +
                 d[1628] * ((-d[2924]) * dv_3118 + dv_2396 + 109.0 * dv_751) +
                 d[50] * (95.0 * dv_1803 + 53.0 * dv_1805));
  sc_33 = d[362] * sc_34;
  DataVector& dv_1553 = temps.at(1427);
  DataVector& dv_1564 = temps.at(1437);
  DataVector& dv_1570 = temps.at(1442);
  DataVector& dv_2172 = temps.at(1780);
  DataVector& dv_2341 = temps.at(1931);
  sc_12 =
      (d[361] * d[514] * d[6]) * dv_2172 +
      d[1466] *
          ((-d[2928] * (d[1033] * d[1221] + d[1420] + d[1867]) +
            d[1311] * d[1865]) *
               dv_1564 +
           (d[1] * d[35] * d[495]) +
           d[259] * ((-d[504]) * (dv_1553 + dv_1570) +
                     d[2921] * ((-42.0 * d[7]) + 53.0 * dv_1502 + dv_2341))) +
      sc_16 + sc_18 + sc_24 + sc_29 + sc_31;
  DataVector& dv_2818 = temps.at(2370);
  sc_12 += d[1864] * ((-d[261] * d[2922]) + dv_2818) + sc_32 + sc_33;
  DataVector& dv_6 = temps.at(6);
  sc_17 = (d[1227] * d[376]) * dv_6 * sc_12;
  DataVector& dv_2497 = temps.at(2085);
  DataVector& dv_2554 = temps.at(2142);
  DataVector& dv_3104 = temps.at(2596);
  DataVector& dv_3108 = temps.at(2600);
  sc_31 =
      (-d[1250]) *
          ((d[42] * (d[1068] + d[1785] + 16.0)) * dv_231 + d[1849] * dv_3108 +
           d[348] * ((d[1851] * d[2923]) * Dx - dv_2497 - dv_3104)) +
      dv_2554;
  DataVector& dv_1683 = temps.at(1516);
  DataVector& dv_2359 = temps.at(1949);
  sc_31 += (-d[1638]) * ((d[1035] + 9.0) * dv_1683 +
                         (-d[2928]) * (dv_2359 + dv_3106) + d[167]) +
           (d[1271] + d[1846]) * dv_2666;
  DataVector& dv_2343 = temps.at(1933);
  DataVector& dv_3107 = temps.at(2599);
  DataVector& dv_637 = temps.at(587);
  sc_31 +=
      d[1248] * ((-d[1847]) * dv_3107 + (-d[20]) * (d[1848] * dv_637 + d[314]) +
                 (d[2928] * d[2921]) * (d[88] + dv_2762) + (9.0 * d[50])) +
      d[1251] * ((-d[1428]) - dv_2343);
  DataVector& dv_1504 = temps.at(1384);
  DataVector& dv_1881 = temps.at(275);
  DataVector& dv_2485 = temps.at(2073);
  sc_31 += d[1299] * ((-d[1]) * dv_1504 + (-d[2927] - 4.0) * dv_1881 + dv_2485);
  sc_32 = (-d[1852]) * sc_31;
  DataVector& dv_3099 = temps.at(2592);
  sc_16 = (d[1240] * d[557]) * dv_5 + d[166] * dv_2606 + d[557] * dv_3099;
  DataVector& dv_2508 = temps.at(2096);
  DataVector& dv_2631 = temps.at(2210);
  sc_16 += d[2922] * ((-d[1831]) * dv_2631 +
                      d[259] * (d[1365] + d[1832] * dv_558 + d[7] * dv_2508) +
                      d[348] * ((-d[1829]) * dv_1 - dv_3100));
  sc_31 = d[1834] * sc_16;
  DataVector& dv_2798 = temps.at(2351);
  DataVector& dv_3102 = temps.at(2594);
  sc_24 =
      (-d[1628]) * ((d[1068] - 45.0) * dv_1 + d[1] * (dv_3102 + 6.0) + d[153]) -
      144.0 * dv_2798;
  DataVector& dv_1530 = temps.at(1407);
  DataVector& dv_2383 = temps.at(1973);
  DataVector& dv_2403 = temps.at(1993);
  DataVector& dv_3103 = temps.at(2595);
  sc_24 +=
      d[1250] * ((-d[1831]) * dv_3103 +
                 (-d[348]) * ((-d[1840]) * dv_1530 + dv_2383 + dv_3104) +
                 (3.0 * d[2928] * d[2921] * (d[1068] + d[1620] + 36.0)) * Dx) +
      d[1251] * dv_2403 + d[138] * dv_794;
  DataVector& dv_2135 = temps.at(1751);
  DataVector& dv_2854 = temps.at(2406);
  sc_24 += d[35] * ((-d[1838]) + d[1] * ((-d[278]) + dv_1683) +
                    d[29] * ((d[1278] + 99.0) * Dy + d[1701])) +
           d[436] * ((d[1035] + 3.0) * dv_1881 + d[2928] * (dv_1503 + dv_2854) +
                     d[1054] + dv_2135);
  sc_16 = d[1841] * sc_24;
  DataVector& dv_2821 = temps.at(2373);
  DataVector& dv_2829 = temps.at(2381);
  DataVector& dv_3053 = temps.at(2546);
  sc_29 = (d[1267] + d[1842]) * dv_2666 +
          (-d[436]) * ((-d[2928]) * (dv_1502 + dv_2821) +
                       (-d[2927] - 9.0) * dv_1881 + d[1694] + dv_2829) +
          d[1251] * dv_3053;
  DataVector& dv_2375 = temps.at(1965);
  sc_29 += d[1625] * ((d[1844] + 12.0) * dv_3105 + (-d[1843]) * dv_3103 +
                      d[20] * ((-d[1780]) * dv_1530 + dv_2375 + dv_2497));
  DataVector& dv_2349 = temps.at(1939);
  DataVector& dv_2358 = temps.at(1948);
  sc_29 += d[1628] * ((-d[1]) * (-dv_1502 + dv_2358 + 6.0) + d[1460] +
                      d[1835] * dv_1) +
           d[558] * dv_2349;
  sc_29 += d[6] * ((-d[1831]) * dv_2666 + (-d[348]) * (Dy * d[1836] + d[314]) +
                   (2.0 * d[2928] * d[2921]) * (d[1418] + dv_2469));
  sc_24 = d[1845] * sc_29;
  DataVector& dv_2372 = temps.at(1962);
  DataVector& dv_2819 = temps.at(2371);
  DataVector& dv_3110 = temps.at(2602);
  DataVector& dv_3112 = temps.at(2604);
  DataVector& dv_780 = temps.at(728);
  sc_18 =
      (d[1306] * d[468] - 147.0 * d[159] + 110.0 * d[48] + d[519]) * dv_780 +
      (-d[1009]) * dv_1803 + (-d[1410]) * dv_1683 +
      (-d[1640]) * (d[262] + dv_624) + (d[1306] * d[1857]) * Dx +
      d[138] * dv_2819 + d[1661] * dv_1112 + 189.0 * dv_2372 + 108.0 * dv_3110 +
      dv_3112;
  DataVector& dv_2114 = temps.at(1731);
  DataVector& dv_3079 = temps.at(2572);
  sc_18 += d[431] * dv_3066 +
           d[2922] * ((-d[1843]) * dv_2704 +
                      (-d[259]) * ((-d[1854]) * dv_558 + d[1853] + dv_3079) +
                      (3.0 * d[20]) * ((-d[2928]) * (dv_2114 + dv_3109) +
                                       d[1460] + d[1851] * dv_1) +
                      (-d[1629]));
  sc_29 = d[1858] * sc_18;
  DataVector& dv_2422 = temps.at(2012);
  DataVector& dv_2429 = temps.at(2019);
  sc_34 = (-d[1140] * d[1285] - d[1293] + d[1859] + 135.0 * d[20]) * dv_780 +
          (-d[1410]) * dv_2429 + (-d[1413]) * dv_2422 +
          (-120.0 * d[1862]) * dv_231 + (-d[1469]) + 180.0 * dv_2371 +
          297.0 * dv_2372 + 180.0 * dv_3110 + dv_3112;
  DataVector& dv_2468 = temps.at(2057);
  DataVector& dv_3113 = temps.at(2605);
  sc_34 +=
      d[1250] *
      ((-d[259]) * ((-d[2927] - 3.0) * dv_2468 + d[1853] + 147.0 * dv_794) +
       (-54.0 * d[276]) + d[1849] * dv_3113 +
       d[348] * ((-d[2928]) * (dv_2114 + dv_3114) + d[1694] + d[1840] * dv_69));
  DataVector& dv_3115 = temps.at(2606);
  sc_34 += d[162] * dv_3066 + d[1861] * dv_3115 + d[559] * dv_2819;
  sc_18 = d[1863] * sc_34;
  DataVector& dv_2390 = temps.at(1980);
  DataVector& sc_35 = temps.at(3257);
  sc_35 = (-d[1570]) * ((-d[36]) * (d[1418] + dv_1683) + d[1837] +
                        d[29] * (d[2928] * (dv_2390 + dv_3101 + 6.0) + d[153] +
                                 d[1780] * dv_69)) +
          (d[29] * (d[1236] + d[1835] * d[2921])) * dv_780 +
          (d[9] * (d[1106] * (d[1832] + d[7]) + d[1242])) * dv_231 +
          (24.0 * d[1493]);
  sc_35 += d[1251] * dv_2784 +
           d[1628] * ((-d[1836]) * dv_751 + dv_2397 + dv_2497) +
           d[91] * dv_1112;
  sc_34 = d[362] * sc_35;
  DataVector& dv_3098 = temps.at(2591);
  sc_33 =
      (d[361] * d[557]) * dv_3098 +
      d[1830] * ((d[255] + d[2921] * (d[1242] - d[1829] * d[2923])) * dv_1564 +
                 (-d[1605] * d[554]) + d[1229] * dv_231) +
      sc_16 + sc_18 + sc_24 + sc_29 + sc_31 + sc_32 + sc_34;
  sc_12 = (-d[1077] * d[476]) * dv_6 * sc_33;
  sc_29 = d[1286] * dv_240 + d[1289] * dv_3063 + d[1289] * dv_3066;
  DataVector& dv_2547 = temps.at(2135);
  DataVector& dv_2769 = temps.at(2322);
  DataVector& dv_3184 = temps.at(2671);
  DataVector& dv_3185 = temps.at(2672);
  sc_29 += d[2922] *
           ((d[1068] + d[1189] - 6.0) * dv_3184 +
            (-d[259]) * (d[1934] * dv_69 + dv_3185) +
            d[20] * (d[314] + dv_2547 * d[2927] - dv_2769) + d[50] * dv_3053);
  sc_18 = (-d[1797]) * sc_29;
  DataVector& dv_2837 = temps.at(2389);
  DataVector& dv_3045 = temps.at(2487);
  sc_24 = d[35] * ((-56.0 * d[319] + d[500]) + d[2921] * (d[1418] + dv_3045)) +
          120.0 * dv_2837;
  DataVector& dv_2707 = temps.at(2260);
  DataVector& dv_3186 = temps.at(2673);
  sc_24 += d[2922] * ((-d[110]) * dv_1112 +
                      (-d[1291] + 52.0 * d[2923]) * dv_3105 + d[495] * dv_3186 +
                      d[51] * (dv_2707 + dv_624 * d[2924]) + 40.0 * dv_2372);
  DataVector& dv_1632 = temps.at(1477);
  DataVector& dv_1802 = temps.at(1588);
  DataVector& dv_2484 = temps.at(2072);
  DataVector& dv_2686 = temps.at(2242);
  DataVector& dv_2702 = temps.at(2255);
  sc_24 +=
      d[2921] *
      ((-d[484]) * (dv_1632 + dv_2772) +
       (-d[2921]) * (d[1460] + d[1933] * dv_1 + dv_1802 + dv_2484 + dv_2686) +
       (2.0 * d[20] * d[2923]) * (d[7] + dv_2702 + dv_3101));
  sc_29 = d[1466] * sc_24;
  DataVector& dv_2339 = temps.at(1929);
  DataVector& dv_2379 = temps.at(1969);
  DataVector& dv_3187 = temps.at(2674);
  DataVector& dv_3188 = temps.at(2675);
  sc_16 = (d[1936] * d[77] - d[244] * d[2923] + 296.0 * d[319]) * dv_780 +
          (-d[1190]) * dv_1803 + (-d[1378]) * (d[1935] - 50.0 * dv_5) +
          (-d[1448]) * dv_2339 + (-d[72]) * dv_2379 +
          (-95.0 * d[319] - d[502] * d[511] - d[505]) * dv_3067 +
          (d[115] * d[2925]) * dv_751 + dv_3187 + 30.0 * dv_3188;
  sc_16 +=
      d[1096] * dv_751 + d[116] * dv_2376 + d[1559] * dv_231 + d[367] * dv_2375;
  DataVector& dv_2728 = temps.at(2281);
  DataVector& dv_2834 = temps.at(2386);
  DataVector& dv_2968 = temps.at(2496);
  sc_16 += d[2922] * ((-d[20]) * (d[434] + dv_2834 * d[2927] - 78.0 * dv_794) +
                      (-d[244]) * dv_794 +
                      (2.0 * d[50]) * ((-d[1608]) + dv_2728 + dv_2968) +
                      (3.0 * d[2928] * d[2921]) *
                          ((d[1937] + 8.0) * dv_69 + d[153] - dv_3072));
  sc_24 = d[1800] * sc_16;
  DataVector& dv_1808 = temps.at(1594);
  DataVector& dv_2298 = temps.at(1891);
  DataVector& dv_2842 = temps.at(2394);
  DataVector& dv_3189 = temps.at(2676);
  DataVector& dv_3193 = temps.at(2680);
  sc_31 =
      (-d[1286]) * dv_2298 +
      (-d[259]) * ((-d[1936]) * dv_3189 + (d[2928] * d[1447]) +
                   d[484] * dv_3078 + dv_3193) +
      (-d[276]) * (d[1368] + dv_1502 + dv_2842) +
      (-d[2922]) * ((d[1029] * d[514] - d[110] * (d[1620] + 12.0) +
                     d[1249] * (d[1242] - d[1462]) + d[1297] - d[448] * d[7]) *
                        Dx +
                    d[1285] * dv_1808);
  DataVector& dv_1554 = temps.at(1428);
  DataVector& dv_2748 = temps.at(2301);
  DataVector& dv_3191 = temps.at(2678);
  sc_31 += (6.0 * d[48] * d[2923] * (d[1368] + d[1854])) * Dy +
           d[20] * (d[1509] + dv_1554 + dv_2748 - 53.0 * dv_3073 + dv_3191);
  DataVector& dv_2387 = temps.at(1977);
  sc_31 += d[6] * (d[1433] + d[248] * (d[1540] + dv_2387) + d[490] * dv_1 +
                   d[77] * ((19.0 - d[1806]) * Dy + (-d[1943])));
  sc_16 = d[1804] * sc_31;
  DataVector& dv_2193 = temps.at(732);
  DataVector& dv_2487 = temps.at(2075);
  DataVector& dv_2670 = temps.at(2232);
  sc_32 =
      (10.0 * d[138] + 10.0 * d[781]) * dv_2136 +
      (-d[138] * (d[1209] + d[2927])) * dv_1 +
      d[1411] * ((-d[1276]) + dv_2487 + dv_2728) +
      d[20] * (d[1815] + dv_2193 + 46.0 * dv_2246 + dv_2670 - 95.0 * dv_3073);
  DataVector& dv_3196 = temps.at(2683);
  DataVector& dv_3197 = temps.at(2684);
  sc_32 += d[259] * ((-d[1952]) * dv_3189 - dv_3193 - dv_3196) +
           d[6] * ((-114.0 * d[276]) + d[116] * (d[1464] + 149.0 * dv_1) +
                   d[396] * (d[1953] + dv_3197) + d[48] * dv_2469);
  sc_32 += d[2922] * ((d[1946] + 310.0 * d[196] + d[244] * (d[1189] + 10.0) +
                       d[259] * (d[1258] + 128.0 * d[2923]) +
                       d[2927] * (d[138] + 234.0 * d[159] - d[1954])) *
                          Dx +
                      d[1935] * dv_1808);
  sc_31 = d[1807] * sc_32;
  DataVector& dv_2729 = temps.at(2282);
  DataVector& dv_3011 = temps.at(2536);
  DataVector& dv_3190 = temps.at(2677);
  sc_35 = (d[162] + d[423]) * dv_3011 +
          (-d[259]) * ((-d[1941]) * dv_3189 + 108.0 * dv_3068 + dv_3190) +
          d[1411] * (d[1939] + dv_2702 + dv_2729) + d[1938] * dv_2704;
  DataVector& dv_2591 = temps.at(2175);
  sc_35 += d[20] * ((-d[1940]) * dv_1 + d[1734] + d[622] * dv_1502 +
                    170.0 * dv_2246 + dv_2591);
  DataVector& dv_1014 = temps.at(914);
  sc_35 +=
      d[6] * ((-52.0 * d[276]) + d[367] * (d[1567] + 31.0 * dv_1) +
              d[77] * ((d[1806] + 26.0) * Dy + d[1942]) + 231.0 * dv_1014) +
      d[2922] * ((d[1396] + d[1792] * (234.0 * d[36] - 95.0 * d[2921]) +
                  298.0 * d[196] + d[259] * (d[1258] + d[1685]) +
                  21.0 * d[48] * (d[1580] + 6.0)) *
                     Dx +
                 (d[244] + d[248]) * dv_1808);
  sc_32 = d[1812] * sc_35;
  DataVector& dv_2544 = temps.at(2132);
  DataVector& dv_3194 = temps.at(2681);
  DataVector& dv_3195 = temps.at(2682);
  DataVector& dv_3198 = temps.at(2685);
  DataVector& dv_972 = temps.at(878);
  DataVector& sc_36 = temps.at(3258);
  sc_36 = (d[1817] * d[511] + d[505] - d[506]) * dv_3067 +
          (-d[1952] * d[77] + 75.0 * d[252] + 78.0 * d[319]) * dv_780 +
          d[1297] * dv_2544 +
          d[1359] * ((d[1448] + d[50]) + d[20] * dv_972 + dv_3198) +
          d[1415] * dv_751 + d[162] * dv_3194 + d[1955] * dv_231 + dv_3187 +
          100.0 * dv_3188 + dv_3195;
  DataVector& dv_2758 = temps.at(2311);
  sc_36 += d[367] * dv_2376 + d[48] * dv_2758 + d[631] * dv_1803;
  DataVector& dv_2969 = temps.at(2497);
  DataVector& dv_3199 = temps.at(2686);
  sc_36 += d[2922] * ((-d[20]) * (d[434] + dv_3130 * d[2927] - 296.0 * dv_794) +
                      (-d[259]) * ((-d[1950] - 64.0) * dv_69 + dv_3199) +
                      (-d[51]) * (d[1522] + dv_1502 + dv_2969) +
                      (15.0 * d[48] * (d[1585] + 6.0)) * Dy);
  sc_35 = d[1821] * sc_36;
  DataVector& dv_2726 = temps.at(2279);
  DataVector& dv_613 = temps.at(564);
  DataVector& sc_37 = temps.at(3259);
  sc_37 = (d[1813] * d[510] - 127.0 * d[3]) * dv_3186 +
          (d[1941] * d[77] + d[1947] * d[2923] + 370.0 * d[319]) * dv_780 +
          (52.0 * d[255]) * dv_231 + d[1096] * dv_3194 +
          d[1359] * ((d[1657] - d[1948]) + dv_2726 + 66.0 * dv_613) +
          d[1373] * dv_2379 + d[1472] * dv_1803 + d[1944] * dv_751 + dv_3195;
  DataVector& dv_2521 = temps.at(2109);
  sc_37 +=
      d[1945] * dv_2339 + d[1946] * dv_751 + d[20] * dv_2521 + d[882] * dv_2379;
  DataVector& dv_2536 = temps.at(2124);
  DataVector& dv_3090 = temps.at(2583);
  sc_37 += d[2922] *
           ((d[1068] + d[1949] + 24.0) * dv_3184 +
            (-d[20]) * (Dy * d[1940] + d[444] - 370.0 * dv_794) +
            d[259] * ((d[1950] + 62.0) * dv_69 + d[1499] + d[46] * dv_3090) +
            d[51] * ((-d[1321]) + dv_2536 + dv_2729));
  sc_36 = d[1951] * sc_37;
  DataVector& dv_2392 = temps.at(1982);
  DataVector& dv_2700 = temps.at(2253);
  DataVector& dv_3180 = temps.at(2667);
  DataVector& dv_3183 = temps.at(2670);
  sc_34 = (-d[1289]) * dv_3180 +
          (-d[1581]) * ((-d[1247] + d[1375] + d[1933] * d[2923]) * dv_231 +
                        (-d[1932]) * dv_780 + d[1432] + d[319] * dv_2392 +
                        d[86] * ((-d[2921]) * (dv_2700 + dv_3101) + d[314] +
                                 d[7] * (d[1506] + dv_538)) +
                        dv_3183);
  DataVector& dv_3181 = temps.at(2668);
  sc_34 += (-d[1931]) * dv_3181 + sc_16 + sc_18 + sc_24 + sc_29 + sc_31 +
           sc_32 + sc_35 + sc_36;
  DataVector& dv_2107 = temps.at(1725);
  sc_33 = (-d[1223] * d[92]) * dv_2107 * sc_34;
  DataVector& dv_1579 = temps.at(1451);
  sc_36 =
      (-d[122]) * ((-d[29]) * dv_1502 + (d[395] * d[6]) + d[1778] * dv_1579) +
      (-d[123]) * ((-d[1781] - d[1782]) * Dy + (-d[29]) * dv_3053 +
                   d[1779] * dv_1579 + d[1780] * dv_1115);
  DataVector& dv_2189 = temps.at(1795);
  DataVector& dv_3054 = temps.at(2547);
  DataVector& dv_3055 = temps.at(2548);
  DataVector& dv_3056 = temps.at(2549);
  sc_36 += (d[19] * d[50]) * ((15.0 - d[2927]) * dv_780 + Dx * d[1782] +
                              d[81] * dv_3055 + 33.0 * dv_1112 + dv_3054) +
           (d[20] * d[52]) * ((-d[1783] - d[1785]) * Dy + d[1784] * dv_2189 +
                              d[29] * dv_3056 + d[6] * dv_3055);
  DataVector& dv_265 = temps.at(258);
  sc_36 += (d[55] * d[2921]) * ((d[1068] + 39.0) * dv_780 + (-d[1783]) * Dx +
                                (-d[1092]) * (Dy * d[1779] + d[302]) +
                                9.0 * dv_1112 + dv_3054) +
           (-d[120] * d[1778]) * dv_265;
  DataVector& dv_289 = temps.at(281);
  sc_34 = d[1202] * dv_289 * sc_36;
  sc_16 = d[2922];
  DataVector& dv_2227 = temps.at(1829);
  sc_16 *= (-d[1443]) * (d[1970] * dv_1 + dv_3185) +
           (-d[50]) * ((-d[1968]) + d[1315] * dv_1503 + d[1963] * dv_2227) +
           (d[151] * (d[1070] + d[1813])) * dv_558 +
           d[273] * (Dy * d[1965] + d[1731] - 527.0 * dv_794);
  DataVector& dv_3205 = temps.at(2692);
  sc_31 = (-d[1328]) * dv_3099 + (d[1213] * d[280]) * dv_833 +
          d[259] * ((-d[1328]) * dv_1805 + d[1806] * dv_3205) + sc_16;
  sc_32 = (-d[1797]) * sc_31;
  DataVector& dv_2235 = temps.at(1836);
  DataVector& dv_3071 = temps.at(2564);
  sc_29 = (-d[151]) * dv_2235 +
          (-d[273]) * ((d[1806] + 32.0) * Dy + d[1731] - 676.0 * dv_794) +
          d[101] * ((d[1969] + 80.0) * dv_1 + (-d[46]) * dv_3071 + d[1406]) +
          d[1692];
  DataVector& dv_3215 = temps.at(2702);
  sc_29 += d[51] * ((d[1068] - 13.0) * dv_1881 + (-164.0 * d[255]) +
                    d[2928] * (57.0 * dv_1502 + dv_3215 + 12.0));
  sc_24 = (-d[2922]) * sc_29;
  DataVector& dv_2420 = temps.at(2010);
  DataVector& dv_3111 = temps.at(2603);
  DataVector& dv_3209 = temps.at(2696);
  DataVector& dv_3210 = temps.at(2697);
  DataVector& dv_3211 = temps.at(2698);
  DataVector& dv_3212 = temps.at(2699);
  DataVector& dv_3213 = temps.at(2700);
  DataVector& dv_3214 = temps.at(2701);
  sc_16 = (d[1339] + d[1346] - d[1798] * d[1988] - 496.0 * d[274]) * dv_3063 +
          (d[1984] + d[1987] + 63.0 * d[274]) * dv_3214 + (-d[1331]) * dv_751 +
          (-d[1980]) * dv_3111 + (-d[1982]) * dv_1803 + 90.0 * dv_2420 +
          dv_3209 + dv_3210 + 48.0 * dv_3211 - dv_3212 - 60.0 * dv_3213;
  sc_16 += (-d[273]) * dv_2758 + d[1348] * dv_2379 + d[1348] * dv_3194 +
           d[1640] * (d[1983] - 40.0 * dv_5) + d[1981] * dv_2339 + sc_24;
  sc_31 = (-d[1989]) * sc_16;
  DataVector& dv_2982 = temps.at(2509);
  sc_18 = (-d[101]) * ((-891.0 * d[2927] - 200.0) * dv_1 + dv_3199) +
          (-d[273]) * ((d[1998] + 16.0) * dv_538 + d[1731] - 1488.0 * dv_794) +
          (-d[51]) * ((-d[2928]) * (dv_2982 + dv_3215 + 4.0) +
                      (-d[1971]) * dv_1881 + d[1664]) +
          (-d[1996]);
  sc_18 += (3.0 * d[151] * (d[1140] + d[1997] + 24.0)) * Dy;
  sc_29 = sc_18 * d[2922];
  DataVector& dv_3216 = temps.at(2703);
  DataVector& dv_3217 = temps.at(2704);
  sc_24 = (d[1990] + d[1993]) * dv_3214 +
          (-d[101] * d[1994] + d[1105] * d[1973] + 87.0 * d[1339] +
           676.0 * d[274]) *
              dv_780 +
          (-d[1239]) * dv_850 + (-d[223]) * dv_1803 + (-d[223]) * dv_1805 +
          d[1161] * dv_3194 - 32.0 * dv_3211 + dv_3212 + dv_3216 +
          240.0 * dv_3217;
  DataVector& dv_2386 = temps.at(1976);
  DataVector& dv_2453 = temps.at(2042);
  DataVector& dv_2840 = temps.at(2392);
  sc_24 += d[1327] * ((d[101] + d[50]) + dv_2386 + dv_3107) +
           d[1341] * dv_3111 + d[151] * dv_2453 + d[151] * dv_2840 +
           d[1982] * dv_1805 + sc_29;
  sc_16 = (d[125] * d[52]) * sc_24;
  DataVector& dv_1552 = temps.at(1426);
  DataVector& dv_1667 = temps.at(1458);
  DataVector& dv_2360 = temps.at(1950);
  DataVector& dv_2844 = temps.at(2396);
  DataVector& dv_3207 = temps.at(2694);
  sc_37 = (-d[1412]) * (d[2928] * (dv_2844 + dv_3207 + 4.0) + d[1403] +
                        d[1973] * dv_1552) +
          d[1105] * (dv_1503 + dv_2360) + d[1287] * (-dv_1667 + dv_2772);
  DataVector& dv_2247 = temps.at(1845);
  DataVector& dv_2864 = temps.at(2416);
  sc_37 += d[259] * ((d[1972] + 16.0) * dv_1 + d[1723] + dv_2247 + dv_2864);
  sc_18 = (-d[27]) * sc_37;
  DataVector& dv_1514 = temps.at(1393);
  DataVector& dv_2347 = temps.at(1937);
  DataVector& dv_2442 = temps.at(2032);
  DataVector& dv_3206 = temps.at(2693);
  sc_29 = d[1250] *
              ((-d[151]) * dv_3206 + (-d[1291] + d[1970] * d[2923]) * dv_2445 +
               (4.0 * d[50]) * (d[1971] * dv_1514 + dv_2347 + 35.0 * dv_2375) +
               (10.0 * d[2928] * d[20] * (d[1781] - 4.0)) * Dx) +
          600.0 * dv_2442 + sc_18;
  sc_29 +=
      d[35] * ((-d[243]) * ((-d[1829]) * dv_538 + (170.0 * d[2928] * d[2923])) +
               (d[1389] + 72.0 * d[50]) + 573.0 * dv_3126);
  sc_24 = d[1466] * sc_29;
  DataVector& dv_2545 = temps.at(2133);
  DataVector& dv_3204 = temps.at(2691);
  sc_18 = (d[1505] - d[2923] * (d[155] + d[1965])) * dv_3105 +
          (d[1964] * d[302] + 527.0 * d[36]) * dv_3204 + (-d[1333]) * dv_2784 +
          (-120.0 * d[1493]) +
          d[319] * (d[1966] * dv_2545 + 37.0 * dv_2376 - dv_2396) +
          d[559] * dv_1112;
  sc_18 +=
      d[86] *
      ((-d[2928]) * (d[1731] + dv_3075 - 111.0 * dv_794) + (48.0 * d[319]) +
       d[2921] * ((-62.0 * d[255]) +
                  d[2928] * (120.0 * dv_1502 + 37.0 * dv_1503 + 16.0) +
                  d[1967] * dv_2227));
  sc_29 = d[1581] * sc_18;
  DataVector& dv_2746 = temps.at(2299);
  DataVector& dv_2860 = temps.at(2412);
  DataVector& dv_3208 = temps.at(2695);
  sc_37 = (d[1595] + d[1806] + 12.0) * dv_2860 +
          (-d[1975]) * ((-d[1306]) * dv_794 + d[1579] + d[36] * dv_3078 +
                        dv_2746 * d[2927]) +
          (-d[280]) * dv_3208 + (-d[780]) * dv_3053 +
          (20.0 * d[273]) * (d[1489] * dv_69 + dv_2772);
  DataVector& dv_582 = temps.at(534);
  sc_37 += d[276] *
           ((-d[1974]) + d[2928] * (37.0 * dv_1502 + 120.0 * dv_1503 + 16.0) +
            d[1964] * dv_582);
  DataVector& dv_2491 = temps.at(2079);
  sc_37 +=
      d[286] * ((d[2928] * d[367]) * (d[1567] + dv_2491) +
                d[1448] * ((-d[1976]) * Dy + (-d[1365])) + d[151] * dv_2491 +
                d[50] * ((-d[1498]) + d[1966] * dv_558));
  DataVector& dv_94 = temps.at(92);
  sc_37 +=
      d[2922] *
      ((-d[1070] * d[1979] + d[1331] * d[1874] - d[1443] * (d[1291] - d[1977]) +
        d[273] * (573.0 * d[7] - 32.0) + d[50] * (d[1316] - d[1332])) *
           Dx +
       (-d[1240] * d[280]) * dv_94);
  sc_18 = d[1804] * sc_37;
  DataVector& sc_40 = temps.at(3262);
  sc_40 = (-d[273]) * ((d[2000] + 8.0) * dv_637 + d[1731] - 1845.0 * dv_794) +
          (d[151] * (d[1306] + d[2014] + 18.0)) * dv_538 +
          d[1448] * ((d[1972] + 20.0) * dv_1552 + dv_3091) + d[1692];
  DataVector& dv_2733 = temps.at(2286);
  DataVector& dv_3222 = temps.at(2709);
  sc_40 += d[58] * ((-200.0 * d[255]) + d[2928] * (dv_2733 + dv_3222 + 8.0) +
                    d[1999] * dv_582);
  DataVector& sc_39 = temps.at(3261);
  sc_39 = sc_40 * d[2922];
  DataVector& dv_2464 = temps.at(2053);
  DataVector& dv_3221 = temps.at(2708);
  DataVector& sc_38 = temps.at(3260);
  sc_38 = (d[101] * (d[1593] - 70.0) + d[1389] + d[1395] - d[2013]) * dv_3214 +
          (d[101] * d[2009] + d[1292] * d[1999] + 183.0 * d[1339] +
           1845.0 * d[274]) *
              dv_780 +
          (-d[273]) * dv_2464 + (114.0 * d[1347]) * dv_850 +
          (243.0 * d[50]) * dv_3221 - 126.0 * dv_2420 - dv_3209 - dv_3210 +
          50.0 * dv_3213 + dv_3216 + 342.0 * dv_3217;
  sc_38 += d[1327] * ((d[1657] + d[2012]) + Dy * d[350] + dv_2385) +
           d[1331] * dv_3194 + d[1340] * dv_2379 + d[151] * dv_2464 +
           d[866] * dv_1803 + d[866] * dv_1805 + sc_39;
  sc_37 = d[1816] * sc_38;
  DataVector& dv_3219 = temps.at(2706);
  sc_39 = (-d[885]) * dv_3056 + (-d[1029] - d[284] + 6.0) * dv_3132 +
          d[101] * ((-d[1994]) * dv_794 - dv_3196 - dv_3219);
  DataVector& dv_3017 = temps.at(2542);
  sc_39 += d[1393] * (d[2928] * (81.0 * dv_1502 + dv_3017 + 8.0) + d[1734] +
                      d[1999] * dv_1881) +
           d[1983] * dv_3208;
  DataVector& dv_2589 = temps.at(2164);
  sc_39 +=
      d[6] * ((-d[1105]) * ((-d[1999]) * dv_538 + (100.0 * d[2928] * d[2923])) +
              d[151] * dv_2589 + d[1975] * (d[1953] + dv_3075) + d[235] +
              d[631] * (d[1634] + 661.0 * dv_1));
  DataVector& dv_3218 = temps.at(2705);
  DataVector& dv_3220 = temps.at(2707);
  sc_39 += d[631] * ((-d[2000] - 4.0) * dv_1683 + d[1815] + d[257] * dv_2343 +
                     dv_3218) +
           d[2922] * ((d[101] * (d[1258] + 200.0 * d[2923]) +
                       d[1450] * (469.0 * d[7] - 32.0) +
                       d[1700] * (d[2002] - 28.0 * d[2923]) +
                       d[1806] * (d[1454] + d[2003] + 99.0 * d[277]) +
                       d[2001] * (d[1214] + 24.0)) *
                          Dx +
                      d[1983] * dv_3220);
  sc_38 = d[2004] * sc_39;
  DataVector& dv_2759 = temps.at(2312);
  sc_40 = (-d[101]) * ((-d[2009]) * dv_794 + Dy * d[2008] + dv_3190) +
          (-d[223]) * (d[1736] + dv_1502 + dv_2114) +
          (18.0 * d[2005]) * dv_2759 +
          (d[1331] * (d[1035] + d[1189] + 6.0)) * dv_1;
  DataVector& dv_2402 = temps.at(1992);
  sc_40 += d[1411] * ((-d[2006]) + d[2928] * (57.0 * dv_1503 + dv_3207 + 12.0) +
                      d[1829] * dv_1683) +
           d[273] * ((-d[2007] - 32.0) * dv_1552 + d[1815] + d[622] * dv_2402 +
                     456.0 * dv_2246);
  DataVector& dv_2427 = temps.at(2017);
  sc_40 += d[6] * (d[101] * ((d[1960] + 50.0) * Dy + d[2010]) +
                   d[273] * ((-d[290]) + 1407.0 * dv_1) +
                   d[51] * ((-164.0 * d[36]) + d[1973] * dv_637) + d[780] +
                   159.0 * dv_2427);
  sc_40 += d[2922] *
           ((d[1029] * (d[1388] + d[1769] - 66.0 * d[273] + 297.0 * d[277]) +
             d[1409] * (35.0 * d[1242] - d[1517]) +
             d[1448] * (d[1724] * d[2925] + d[247]) +
             d[2001] * (53.0 * d[7] + 36.0) + d[631] * (661.0 * d[7] - 36.0)) *
                Dx +
            d[2005] * dv_3220);
  sc_39 = d[2011] * sc_40;
  sc_35 = (d[2928] * d[1328]) * dv_3180 +
          d[1931] * ((37.0 * d[255] + d[2921] * (d[1316] + d[1343] * d[1963])) *
                         dv_0 +
                     (d[1311] * d[2924]) * dv_231 + (-d[1311] * d[1605])) +
          sc_16 + sc_18 + sc_24 + sc_29 + sc_31 + sc_32 + sc_37 + sc_38 + sc_39;
  DataVector& dv_3093 = temps.at(2586);
  sc_36 = d[527] * dv_3093 * sc_35;
  DataVector& sc_19 = temps.at(3241);
  DataVector& sc_20 = temps.at(3242);
  DataVector& sc_21 = temps.at(3243);
  DataVector& sc_22 = temps.at(3244);
  DataVector& sc_27 = temps.at(3249);
  sc_21 = sc_19 + sc_20 + sc_22 + sc_27;
  DataVector& dv_2116 = temps.at(1733);
  DataVector& dv_2117 = temps.at(1734);
  DataVector& dv_2119 = temps.at(1736);
  DataVector& dv_2120 = temps.at(1737);
  DataVector& dv_2121 = temps.at(291);
  DataVector& dv_2139 = temps.at(1755);
  DataVector& dv_3057 = temps.at(2550);
  DataVector& dv_547 = temps.at(501);
  DataVector& sc_2 = temps.at(3224);
  DataVector& sc_8 = temps.at(3230);
  sc_21 +=
      (-d[74]) *
          (d[1210] * ((-d[2922]) * dv_2121 + dv_2120) + dv_2139 +
           d[2920] * ((-d[81]) * dv_2117 + dv_2116 + 6.0 * dv_3057 +
                      d[2924] * ((-d[1269]) * dv_547 + dv_2119 + dv_555))) +
      sc_2 + sc_8;
  DataVector& sc_10 = temps.at(3232);
  DataVector& sc_14 = temps.at(3236);
  DataVector& sc_15 = temps.at(3237);
  DataVector& sc_23 = temps.at(3245);
  DataVector& sc_25 = temps.at(3247);
  DataVector& sc_26 = temps.at(3248);
  DataVector& sc_30 = temps.at(3252);
  DataVector& sc_9 = temps.at(3231);
  sc_21 += sc_10 + sc_12 + sc_13 + sc_14 + sc_15 + sc_17 + sc_23 + sc_25 +
           sc_26 + sc_28 + sc_30 + sc_33 + sc_34 + sc_36 + sc_6 + sc_9;
  DataVector& dv_1498 = temps.at(1378);
  DataVector& sc_11 = temps.at(3233);
  sc_11 = -dv_1498 * sc_21;
  DataVector& dv_14 = temps.at(14);
  DataVector& dv_15 = temps.at(15);
  DataVector& dv_1537 = temps.at(1412);
  DataVector& dv_16 = temps.at(16);
  DataVector& dv_1868 = temps.at(255);
  DataVector& dv_1869 = temps.at(267);
  DataVector& dv_1870 = temps.at(280);
  DataVector& dv_1871 = temps.at(285);
  DataVector& dv_1872 = temps.at(248);
  DataVector& dv_1873 = temps.at(13);
  sc_33 = d[1042] * dv_1873 + d[1067] * dv_14 + d[1067] * dv_15 +
          d[1067] * dv_16 + d[107] * dv_0 + dv_1537 - dv_1868 - dv_1869 -
          dv_1870 - dv_1871 - dv_1872;
  sc_33 += d[107] * dv_1;
  DataVector& dv_1874 = temps.at(68);
  sc_34 = 2.0 * dv_1874 * sc_33;
  DataVector& dv_1536 = temps.at(1386);
  DataVector& dv_1589 = temps.at(49);
  DataVector& dv_1595 = temps.at(135);
  DataVector& dv_1601 = temps.at(186);
  DataVector& dv_1623 = temps.at(1469);
  DataVector& dv_1624 = temps.at(1391);
  DataVector& dv_1665 = temps.at(1506);
  DataVector& dv_1704 = temps.at(1533);
  DataVector& dv_1739 = temps.at(3);
  DataVector& dv_1791 = temps.at(1581);
  DataVector& dv_1876 = temps.at(287);
  DataVector& dv_314 = temps.at(300);
  sc_36 = dv_1791 * ((-d[16]) * dv_1665 + (-3.0 * d[1112]) * dv_1623 +
                     (d[1042] * d[105]) * dv_1601 + d[1113] * dv_1589 +
                     d[1114] * dv_1624 + 4.0 * dv_1536 * dv_1595 +
                     dv_314 * (dv_1739 + d[2922] * (dv_1704 + dv_1876)));
  DataVector& dv_10 = temps.at(10);
  DataVector& dv_1612 = temps.at(1462);
  DataVector& dv_1613 = temps.at(1463);
  DataVector& dv_1617 = temps.at(1467);
  DataVector& dv_1618 = temps.at(1396);
  DataVector& dv_1619 = temps.at(32);
  DataVector& dv_1621 = temps.at(1395);
  DataVector& dv_1866 = temps.at(200);
  DataVector& dv_1867 = temps.at(305);
  DataVector& dv_375 = temps.at(358);
  sc_36 += (-d[1109]) * dv_1867 *
               ((-d[1068]) * dv_1866 + (-d[16]) * dv_1617 +
                (-d[21]) * dv_1618 * dv_375 + d[1064] * dv_1612 +
                dv_10 * dv_1621 + dv_1613 * d[2927] + dv_1619) +
           sc_34;
  DataVector& dv_1877 = temps.at(290);
  sc_21 = -dv_1877 * sc_36;
  DataVector& dv_1742 = temps.at(1457);
  DataVector& dv_1778 = temps.at(1568);
  DataVector& dv_1779 = temps.at(54);
  DataVector& dv_1784 = temps.at(229);
  DataVector& dv_1843 = temps.at(1627);
  DataVector& dv_1844 = temps.at(254);
  DataVector& dv_1845 = temps.at(1394);
  DataVector& dv_1865 = temps.at(201);
  DataVector& dv_1916 = temps.at(1636);
  DataVector& dv_1921 = temps.at(1640);
  DataVector& dv_1972 = temps.at(1655);
  DataVector& dv_1976 = temps.at(1319);
  DataVector& sc_5 = temps.at(3227);
  DataVector& sc_7 = temps.at(3229);
  sc_5 = (4.0 * d[1041]) * dv_1921 + (d[1108] * d[549]) * dv_1865 +
         (12.0 * d[1041] * d[465]) * dv_1916 +
         (90.0 * d[208] * d[2927]) * dv_1742 +
         (54.0 * d[16] * d[22] * d[2927]) * dv_1916 + d[1097] * dv_1784 +
         dv_1778 * dv_1779 + dv_1843 * dv_1844 + dv_1843 * dv_1845 +
         dv_1972 * dv_1976 + sc_7;
  DataVector& dv_1525 = temps.at(1402);
  DataVector& dv_1534 = temps.at(1411);
  DataVector& dv_1749 = temps.at(1554);
  DataVector& dv_3765 = temps.at(3221);
  DataVector& dv_3766 = temps.at(3168);
  DataVector& dv_3767 = temps.at(3187);
  DataVector& sc_4 = temps.at(3226);
  sc_5 += -dv_1525 * dv_1749 - dv_1534 * dv_1749 + dv_3765 * dv_3766 +
          dv_3765 * dv_3767 + sc_11 + sc_4;
  DataVector& dv_1539 = temps.at(71);
  DataVector& dv_1575 = temps.at(1447);
  DataVector& dv_1577 = temps.at(1449);
  DataVector& dv_1584 = temps.at(1456);
  DataVector& dv_1761 = temps.at(1564);
  DataVector& dv_1763 = temps.at(53);
  DataVector& dv_1764 = temps.at(1566);
  DataVector& dv_58 = temps.at(58);
  DataVector& dv_75 = temps.at(74);
  sc_5 += -dv_1764 * ((-d[12]) * dv_1539 + (d[10] * d[2927]) * dv_75 +
                      (d[12] * d[9]) * dv_1577 +
                      d[21] * (d[6] * dv_1763 + dv_1575 +
                               d[2922] * (-Dx * dv_1761 + d[1085] * dv_58)) +
                      dv_1584);
  DataVector& dv_1609 = temps.at(1459);
  DataVector& dv_1768 = temps.at(1570);
  DataVector& dv_1776 = temps.at(1578);
  DataVector& dv_1777 = temps.at(1446);
  DataVector& dv_1787 = temps.at(1430);
  DataVector& dv_1918 = temps.at(1638);
  DataVector& dv_1922 = temps.at(1641);
  DataVector& dv_1974 = temps.at(1243);
  sc_5 += -36.0 * dv_1609 * dv_1918 - 12.0 * dv_1609 * dv_1922 -
          dv_1768 * dv_1776 - dv_1768 * dv_1787 - dv_1776 * dv_1777 -
          dv_1777 * dv_1787 - dv_1972 * dv_1974 + sc_21;
  DataVector& dv_1497 = temps.at(1377);
  DataVector& dv_1744 = temps.at(134);
  DataVector& dv_1745 = temps.at(174);
  DataVector& dv_1754 = temps.at(1558);
  DataVector& dv_1755 = temps.at(1559);
  DataVector& dv_1757 = temps.at(295);
  DataVector& dv_1758 = temps.at(1561);
  DataVector& dv_1759 = temps.at(1562);
  DataVector& dv_1789 = temps.at(1579);
  DataVector& dv_1919 = temps.at(1223);
  DataVector& dv_1970 = temps.at(1662);
  DataVector& dv_249 = temps.at(242);
  DataVector& dv_305 = temps.at(297);
  sc_5 += (-d[1041]) * dv_1745 * dv_1755 + (-d[1081]) * dv_1744 * dv_1754 +
          (-d[1098]) * dv_1609 * dv_1789 + (-d[1102]) * dv_1497 * dv_1970 +
          (8.0 * d[107]) * dv_1759 * dv_1919 +
          (216.0 * d[1082]) * dv_1757 * dv_249 +
          (d[1041] * d[83]) * dv_1758 * dv_305;
  DataVector& dv_1752 = temps.at(1556);
  DataVector& dv_1756 = temps.at(1560);
  DataVector& dv_1760 = temps.at(1563);
  DataVector& dv_1782 = temps.at(223);
  DataVector& dv_1785 = temps.at(65);
  DataVector& dv_1799 = temps.at(1585);
  DataVector& dv_1920 = temps.at(1639);
  DataVector& dv_306 = temps.at(298);
  sc_5 += (d[1093] * d[1140]) * dv_1919 * dv_249 +
          (d[1099] * d[549]) * dv_1752 * dv_1799 +
          (168.0 * d[1074] * d[330]) * dv_1757 * dv_306 +
          d[109] * dv_1782 * dv_1799 + d[1141] * dv_1756 * dv_1920 +
          d[306] * dv_1779 * dv_1785 + d[588] * dv_1759 * dv_1760;
  DataVector& dv_1524 = temps.at(1401);
  DataVector& dv_1526 = temps.at(1403);
  DataVector& dv_1743 = temps.at(75);
  DataVector& dv_1746 = temps.at(1551);
  DataVector& dv_1751 = temps.at(241);
  DataVector& dv_1753 = temps.at(1557);
  DataVector& dv_1765 = temps.at(1567);
  DataVector& dv_1790 = temps.at(1580);
  DataVector& dv_1975 = temps.at(1667);
  DataVector& dv_229 = temps.at(226);
  DataVector& dv_82 = temps.at(80);
  sc_5 += (-d[1079]) * dv_1743 * dv_1744 * dv_82 +
          (-228.0 * d[260]) * dv_1751 * dv_1752 * dv_1753 -
          6.0 * dv_1524 * dv_1975 * dv_305 - dv_1526 * dv_1745 * dv_1746 -
          dv_1609 * dv_1760 * dv_1790 + dv_1756 * dv_1765 * dv_229;
  DataVector& dv_1786 = temps.at(923);
  DataVector& dv_19 = temps.at(19);
  DataVector& dv_1914 = temps.at(1635);
  DataVector& dv_224 = temps.at(221);
  sc_5 += (-d[558] * d[75]) * dv_1754 * dv_1756 * dv_19 +
          d[1083] * dv_1609 * dv_1743 * dv_1786 +
          d[1139] * dv_1524 * dv_1914 * dv_224;
  DataVector& dv_1611 = temps.at(1461);
  DataVector& sc_1 = temps.at(3223);
  sc_1 = (-d[1027]) * dv_1611 * sc_5;
  DataVector& dv_1510 = temps.at(1389);
  DataVector& dv_1519 = temps.at(1397);
  DataVector& dv_1608 = temps.at(216);
  DataVector& dv_177 = temps.at(175);
  DataVector& dv_43 = temps.at(43);
  DataVector& dv_80 = temps.at(78);
  DataVector& dv_83 = temps.at(81);
  DataVector& dv_84 = temps.at(82);
  DataVector& sc_0 = temps.at(3222);
  DataVector& sc_3 = temps.at(3225);
  sc_3 = (d[1044] * d[24]) *
             ((-d[1046]) * dv_1525 + (-d[1046]) * dv_1534 +
              (d[39] * d[688]) * dv_1526 + d[1035] * dv_83 + dv_1510 * dv_43 +
              dv_1519 * dv_177 + dv_80 * d[2927]) +
         (-d[1044] * d[1078]) * dv_84 + dv_1608 + sc_0;
  DataVector& dv_1508 = temps.at(8);
  DataVector& dv_1740 = temps.at(288);
  DataVector& dv_225 = temps.at(222);
  DataVector& dv_3768 = temps.at(951);
  DataVector& dv_3769 = temps.at(1418);
  DataVector& dv_3771 = temps.at(1218);
  DataVector& dv_3772 = temps.at(991);
  sc_3 += (-7.0 / 12.0 * d[1047] * d[1074]) * dv_1740 +
          (-7.0 / 48.0 * d[1026]) * dv_1609 * dv_225 * dv_3772 +
          ((1.0 / 8.0) * d[1041] * 1.0 / (d[16] * d[16] * d[16] * d[16]) * 1.0 /
           d[577]) *
              dv_1611 * dv_3772 +
          ((13.0 / 16.0) * d[1025] * d[1074] / pow(d[2926], 31.0)) * dv_1611 *
              dv_3771 +
          d[1028] * dv_1508 + dv_1609 * dv_3768 + dv_1609 * dv_3769 + sc_1;
  DataVector& dv_1499 = temps.at(1379);
  get(get<::Tags::dt<CurvedScalarWave::Tags::Psi>>(*result)) = dv_1499 * sc_3;
  DataVector& dv_3795 = temps.at(3160);
  DataVector& dv_3796 = temps.at(1128);
  DataVector& dv_3799 = temps.at(995);
  DataVector& dv_3800 = temps.at(1115);
  DataVector& dv_3811 = temps.at(1173);
  sc_4 = -dv_3795 - dv_3796 + dv_3799 + dv_3800 + dv_3811;
  DataVector& dv_2181 = temps.at(1788);
  DataVector& dv_3802 = temps.at(1365);
  DataVector& dv_3803 = temps.at(3205);
  DataVector& dv_3804 = temps.at(336);
  DataVector& dv_3805 = temps.at(3154);
  DataVector& dv_3806 = temps.at(3219);
  DataVector& dv_3809 = temps.at(3191);
  DataVector& dv_764 = temps.at(712);
  sc_4 += d[54] * ((-d[53]) * dv_3809 +
                   d[52] * (-dv_2181 + d[2921] * (dv_3804 * d[2923] - dv_764)) +
                   d[63] * ((d[2922] * d[2921]) * dv_3806 - dv_3805) + dv_3802 +
                   dv_3803);
  DataVector& dv_1507 = temps.at(1387);
  sc_11 = dv_1507 * sc_4;
  DataVector& dv_1597 = temps.at(145);
  DataVector& dv_311 = temps.at(93);
  DataVector& dv_3791 = temps.at(1504);
  DataVector& dv_3793 = temps.at(3108);
  DataVector& dv_91 = temps.at(89);
  sc_21 = (d[1059] * d[5] + d[2442] * d[68] + d[2922]) * dv_1597 +
          (-d[17]) * dv_3793 - dv_311 - dv_3791 + dv_91 + sc_11;
  sc_5 = (2.0 * d[2926]) * sc_21;
  DataVector& dv_1603 = temps.at(180);
  DataVector& dv_1604 = temps.at(198);
  DataVector& dv_1616 = temps.at(1466);
  DataVector& dv_3778 = temps.at(1209);
  DataVector& dv_3783 = temps.at(897);
  DataVector& dv_3784 = temps.at(3137);
  DataVector& dv_3785 = temps.at(970);
  DataVector& dv_3786 = temps.at(1086);
  DataVector& dv_3787 = temps.at(3166);
  DataVector& dv_3788 = temps.at(83);
  DataVector& dv_79 = temps.at(77);
  DataVector& dv_86 = temps.at(84);
  sc_0 = (9.0 * d[2928]) * dv_1604 * dv_1611 * dv_1616 - dv_1603 * dv_3784 -
         dv_1604 * dv_3785 - dv_3783 * dv_79 - dv_3783 * dv_86 -
         dv_3788 * ((-d[17]) * dv_3787 + dv_3778 + dv_3786) + sc_5;
  sc_1 = (-d[2572]) * sc_0;
  DataVector& dv_134 = temps.at(132);
  DataVector& dv_2020 = temps.at(1521);
  DataVector& dv_2021 = temps.at(1543);
  DataVector& dv_2023 = temps.at(312);
  DataVector& dv_2029 = temps.at(1523);
  DataVector& dv_284 = temps.at(276);
  DataVector& dv_378 = temps.at(361);
  DataVector& dv_572 = temps.at(525);
  DataVector& dv_585 = temps.at(537);
  DataVector& dv_828 = temps.at(772);
  sc_12 = -dv_134 * (d[110] * (dv_572 + dv_585) + d[2729] * dv_2021 +
                     d[54] * dv_2029) -
          dv_284 *
              (d[2727] * dv_2020 + d[48] * (dv_378 + dv_828) + d[54] * dv_2023);
  DataVector& dv_17 = temps.at(17);
  DataVector& dv_2019 = temps.at(95);
  DataVector& dv_386 = temps.at(369);
  DataVector& dv_3972 = temps.at(179);
  DataVector& dv_404 = temps.at(384);
  sc_12 += -dv_386 * (d[110] * (dv_17 + dv_3972) + d[2728] * dv_2019 +
                      d[54] * (-51.0 * dv_14 + dv_404));
  DataVector& dv_121 = temps.at(119);
  DataVector& dv_3790 = temps.at(1465);
  DataVector& dv_397 = temps.at(378);
  DataVector& dv_3973 = temps.at(332);
  DataVector& dv_398 = temps.at(379);
  DataVector& dv_51 = temps.at(51);
  DataVector& dv_645 = temps.at(595);
  DataVector& dv_975 = temps.at(881);
  sc_12 += -dv_975 * (d[2728] * dv_3790 + d[48] * (dv_121 + dv_3973 + dv_51) +
                      d[54] * (dv_397 + dv_398 + dv_645));
  DataVector& dv_163 = temps.at(161);
  DataVector& dv_1651 = temps.at(1496);
  DataVector& dv_2031 = temps.at(1663);
  sc_12 +=
      Dx * d[57] * ((-d[48]) * dv_163 + d[2727] * dv_1651 + d[54] * dv_2031);
  sc_33 = d[119] * sc_12;
  DataVector& dv_2037 = temps.at(1341);
  DataVector& dv_2051 = temps.at(1695);
  DataVector& dv_2085 = temps.at(1673);
  DataVector& dv_2093 = temps.at(1716);
  DataVector& dv_3955 = temps.at(777);
  DataVector& dv_3975 = temps.at(347);
  DataVector& dv_3977 = temps.at(263);
  DataVector& dv_3980 = temps.at(426);
  sc_17 = -dv_3955 * (d[114] * dv_2037 + d[118] * dv_2085 + dv_3975) -
          dv_3980 * (d[114] * dv_2051 + d[118] * dv_2093 + dv_3977);
  DataVector& dv_2242 = temps.at(1843);
  DataVector& dv_3976 = temps.at(360);
  DataVector& dv_3979 = temps.at(341);
  DataVector& dv_3981 = temps.at(421);
  DataVector& dv_425 = temps.at(400);
  DataVector& dv_426 = temps.at(401);
  DataVector& dv_458 = temps.at(431);
  sc_17 += -dv_3981 * (d[114] * (dv_3979 + dv_425 + dv_426) +
                       d[118] * (dv_2242 + dv_458) + dv_3976);
  DataVector& dv_2317 = temps.at(1909);
  DataVector& dv_3978 = temps.at(406);
  DataVector& dv_418 = temps.at(393);
  DataVector& dv_420 = temps.at(395);
  DataVector& dv_449 = temps.at(422);
  DataVector& dv_451 = temps.at(424);
  DataVector& dv_625 = temps.at(576);
  sc_17 += -dv_625 * (d[114] * (-dv_3979 + dv_418 + dv_420) +
                      d[118] * (-dv_2317 + dv_449 + dv_451) + dv_3978);
  DataVector& dv_114 = temps.at(112);
  DataVector& dv_2046 = temps.at(1691);
  DataVector& dv_2090 = temps.at(408);
  DataVector& dv_3974 = temps.at(195);
  DataVector& dv_61 = temps.at(61);
  DataVector& dv_826 = temps.at(771);
  sc_17 +=
      (d[57] * d[2920]) * Dx * (d[114] * dv_2046 + d[118] * dv_2090 + dv_3974) +
      (d[120] * d[2926]) * Dy *
          (d[578] * (dv_114 + dv_826) + d[9] * (dv_545 + dv_61));
  sc_12 = sc_17 * d[2922];
  DataVector& dv_2070 = temps.at(411);
  DataVector& dv_2081 = temps.at(1711);
  DataVector& dv_2084 = temps.at(1713);
  DataVector& dv_488 = temps.at(456);
  DataVector& dv_510 = temps.at(471);
  DataVector& dv_875 = temps.at(647);
  DataVector& dv_877 = temps.at(809);
  sc_6 = -dv_875 * (d[114] * dv_2070 + d[118] * dv_2081 + dv_3977) -
         dv_877 * (d[114] * (-165.0 * dv_14 + dv_488) +
                   d[968] * (dv_2084 + dv_510) + dv_3978);
  DataVector& dv_2060 = temps.at(148);
  DataVector& dv_2079 = temps.at(115);
  DataVector& dv_918 = temps.at(848);
  sc_6 += -dv_918 * (d[114] * dv_2060 + d[118] * dv_2079 + dv_3975);
  DataVector& dv_184 = temps.at(181);
  DataVector& dv_2058 = temps.at(1700);
  DataVector& dv_481 = temps.at(449);
  DataVector& dv_504 = temps.at(465);
  DataVector& dv_919 = temps.at(849);
  sc_6 += (d[122] * d[2926]) * Dy * (d[578] * dv_184 + d[9] * dv_2058) -
          dv_919 * (d[114] * (75.0 * dv_14 + dv_420 + dv_481) +
                    d[118] * (93.0 * dv_14 + dv_451 + dv_504) + dv_3976);
  DataVector& dv_2062 = temps.at(154);
  DataVector& dv_2076 = temps.at(1702);
  sc_6 += Dx * d[120] * (d[114] * dv_2062 + d[118] * dv_2076 + dv_3974);
  sc_17 = sc_6 * d[2923];
  sc_34 = sc_12 + sc_17 + sc_33;
  DataVector& dv_3982 = temps.at(372);
  sc_36 = dv_3982 * sc_34;
  DataVector& dv_1986 = temps.at(302);
  DataVector& dv_2002 = temps.at(37);
  DataVector& dv_346 = temps.at(329);
  DataVector& dv_351 = temps.at(334);
  DataVector& dv_379 = temps.at(362);
  DataVector& dv_3965 = temps.at(417);
  DataVector& dv_3967 = temps.at(340);
  sc_17 = -dv_284 * (d[1109] * dv_1986 + d[118] * dv_2002 + d[49] * dv_3965) -
          dv_386 * (d[114] * (-45.0 * dv_14 + dv_346 + dv_351) +
                    d[2726] * dv_379 + d[49] * dv_3967);
  DataVector& dv_1993 = temps.at(88);
  DataVector& dv_2010 = temps.at(157);
  DataVector& dv_338 = temps.at(321);
  DataVector& dv_3381 = temps.at(2863);
  DataVector& dv_344 = temps.at(327);
  DataVector& dv_3497 = temps.at(2971);
  DataVector& dv_3966 = temps.at(3164);
  DataVector& dv_3969 = temps.at(392);
  sc_17 +=
      (d[19] * d[20]) * Dx *
          ((-d[114]) * dv_1993 + (-d[49]) * dv_3969 + (9.0 * d[12]) * dv_2010) -
      dv_975 * (d[114] * (dv_3381 + dv_344 + dv_346) +
                d[2725] * (dv_338 + dv_3497) + d[49] * dv_3966);
  DataVector& dv_157 = temps.at(155);
  DataVector& dv_2012 = temps.at(1544);
  DataVector& dv_2905 = temps.at(2447);
  sc_17 += Dx * d[57] * (d[114] * dv_157 + d[968] * dv_2012 + dv_2905);
  DataVector& dv_1781 = temps.at(1435);
  sc_34 = -dv_1781 * sc_17;
  DataVector& dv_3776 = temps.at(1597);
  DataVector& dv_3840 = temps.at(1608);
  DataVector& dv_3962 = temps.at(315);
  DataVector& dv_3963 = temps.at(310);
  DataVector& dv_3970 = temps.at(1670);
  DataVector& dv_3971 = temps.at(373);
  sc_7 = (-d[2920]) * dv_3970 + d[1060] * dv_3963 + dv_3776 * dv_3971 +
         dv_3840 * dv_3962 + sc_34 + sc_36;
  DataVector& dv_1864 = temps.at(459);
  sc_4 = dv_1864 * sc_7;
  DataVector& dv_3817 = temps.at(105);
  DataVector& dv_3819 = temps.at(130);
  DataVector& dv_3838 = temps.at(1622);
  DataVector& dv_3839 = temps.at(3151);
  DataVector& dv_3841 = temps.at(111);
  DataVector& dv_3958 = temps.at(3119);
  DataVector& dv_3959 = temps.at(1507);
  DataVector& dv_3960 = temps.at(253);
  DataVector& dv_3964 = temps.at(194);
  sc_11 = d[1186] * dv_3819 - dv_3817 * dv_3959 + dv_3838 * dv_3964 -
          dv_3839 * dv_3958 - dv_3841 * dv_3960 + sc_4;
  sc_21 = (-d[108]) * sc_11;
  DataVector& dv_3852 = temps.at(3211);
  DataVector& dv_3853 = temps.at(1153);
  DataVector& dv_752 = temps.at(700);
  sc_34 = (-10.0 * d[52]) * dv_0 +
          d[19] * ((-d[139]) * dv_3852 + d[37] * (dv_1514 + dv_752)) +
          d[2920] * (d[2580] * dv_0 + d[2582] * dv_1 + dv_3853);
  DataVector& dv_34 = temps.at(34);
  sc_34 += d[2921] * ((-d[1723]) * dv_231 + d[1096] * dv_752 +
                      d[2923] * ((-d[938]) * Dx + d[648] * dv_34));
  sc_7 = (-d[2583]) * sc_34;
  sc_17 = d[2668];
  DataVector& dv_3934 = temps.at(1177);
  DataVector& dv_3935 = temps.at(886);
  DataVector& dv_3936 = temps.at(3104);
  sc_17 *= (-d[50]) * (8.0 * dv_1112 + dv_3936 + 19.0 * dv_780) +
           d[1443] * ((-d[1627] + d[2253]) * Dx + dv_3935) +
           d[273] * ((d[1825] + 41.0) * dv_1530 + dv_3934 + 65.0 * dv_752) +
           d[391] * dv_751;
  DataVector& dv_2134 = temps.at(1750);
  DataVector& dv_2528 = temps.at(2116);
  DataVector& dv_3468 = temps.at(2943);
  DataVector& dv_3909 = temps.at(1283);
  DataVector& dv_3931 = temps.at(52);
  DataVector& dv_3932 = temps.at(3179);
  sc_36 = d[1466] * (d[1584] * dv_2528 + d[2662] * dv_3909 + d[2664] * dv_752) +
          d[1581] * (dv_2134 + dv_3931 + dv_3932) + d[1583] * dv_3468;
  DataVector& dv_1696 = temps.at(1527);
  sc_36 +=
      d[1797] * ((-d[2665]) * dv_240 + d[1583] * dv_1564 + d[325] * dv_1115) +
      d[1800] * ((-d[2666]) * dv_1115 +
                 (-d[1682] + d[2332] * d[77] - 126.0 * d[252]) * dv_1696 +
                 Dy * d[2667]);
  DataVector& dv_3924 = temps.at(3113);
  DataVector& dv_3933 = temps.at(3112);
  sc_36 +=
      d[1845] * (d[287] * (dv_3924 + dv_3935) +
                 d[50] * ((2.0 * d[6]) * Dx - 20.0 * dv_1112 - 17.0 * dv_265) +
                 d[631] * (-dv_1530 + dv_3934 + 15.0 * dv_752) + dv_3933);
  sc_36 +=
      d[1863] * ((d[1443] * (55.0 - d[2670]) + d[1490] * d[2923] -
                  d[2014] * d[50] + 130.0 * d[274]) *
                     Dy +
                 (d[1592] + d[1763] - 134.0 * d[277]) * dv_1696 +
                 (d[1600] + d[902]) * dv_263) +
      d[1905] *
          ((-d[1251] * d[7] + d[1443] * d[2672] + d[272] + 166.0 * d[274]) *
               Dy +
           (d[2342] * d[631] + d[2671] - 305.0 * d[277]) * dv_1696 +
           (128.0 * d[2928] * d[2921] * d[2923] - d[416] - 94.0 * d[48]) *
               dv_263) +
      d[1931] * dv_265 + sc_17;
  DataVector& dv_2451 = temps.at(2040);
  DataVector& dv_3937 = temps.at(1108);
  DataVector& dv_791 = temps.at(739);
  sc_36 += d[2669] *
           ((-d[1443]) * ((d[2372] + d[6]) * Dx + 305.0 * dv_265) +
            (-d[50]) * (19.0 * dv_1112 + dv_2451 + dv_3936) +
            d[273] * ((d[2241] + 31.0) * dv_1530 + dv_3937 + 83.0 * dv_752) +
            dv_791);
  sc_34 = (-d[467]) * sc_36;
  DataVector& dv_2118 = temps.at(1735);
  DataVector& dv_2479 = temps.at(2068);
  sc_17 =
      d[19] * ((d[1] * (d[1371] + 4.0) - d[3] * (d[1680] + 29.0)) * dv_2170 +
               d[1928] * dv_752) +
      d[52] * ((-d[2585]) * dv_2118 + (12.0 * d[6] * d[2923]) * Dy - dv_2479);
  sc_17 += d[2920] * (Dy * d[1915] + d[2586] * dv_1696 + d[62] * dv_1541);
  DataVector& dv_3024 = temps.at(2503);
  sc_17 += d[2921] * ((d[138] - d[388]) * dv_752 + (-d[62]) * dv_2379 +
                      d[2587] * (dv_126 * d[2922] + 28.0 * dv_2443) +
                      d[80] * ((-d[2581] - d[91]) * Dx + (-d[648]) * dv_3024));
  sc_36 = d[168] * sc_17;
  DataVector& dv_3637 = temps.at(3107);
  DataVector& dv_3878 = temps.at(1151);
  DataVector& dv_3879 = temps.at(3109);
  sc_12 = (-d[122]) * ((d[186] - d[2510] + d[416]) * dv_1564 + Dy * d[2139] -
                       172.0 * dv_3637) +
          (-d[299]) * (d[2138] * dv_752 + d[557] * dv_3879) - 60.0 * dv_3878;
  DataVector& dv_1727 = temps.at(1546);
  DataVector& dv_3880 = temps.at(3203);
  DataVector& dv_3881 = temps.at(1149);
  DataVector& dv_3884 = temps.at(1080);
  sc_12 +=
      (d[2928] * d[59]) * (d[1775] * dv_265 + d[319] * (dv_1727 + dv_3881) +
                           d[436] * (dv_1530 + dv_3880)) +
      d[1199] * ((d[2928] * d[20] * d[2145] - d[2146]) * Dy +
                 (104.0 * d[101] + d[2611] + d[58]) * dv_0 + 163.0 * dv_3884);
  DataVector& dv_3885 = temps.at(958);
  sc_12 += d[337] * ((d[278] * (d[1295] - d[1322] + d[2609])) * dv_0 +
                     Dy * d[2144] + dv_3885);
  DataVector& dv_3036 = temps.at(2524);
  DataVector& dv_3882 = temps.at(3111);
  DataVector& dv_3883 = temps.at(1049);
  sc_12 +=
      d[56] * (d[1443] * (-dv_3883 - 66.0 * dv_751) +
               d[273] * ((d[1261] + d[2608] - 9.0) * dv_2170 + 105.0 * dv_265) +
               dv_3036 + 120.0 * dv_3882);
  sc_12 += d[60] * ((d[1393] + d[1475] * d[3] - d[1490] +
                     d[631] * (d[1673] + d[2608] + 1.0)) *
                        dv_2170 +
                    d[2137] * dv_752);
  sc_17 = d[208] * sc_12;
  DataVector& dv_3869 = temps.at(683);
  DataVector& dv_3870 = temps.at(3189);
  sc_33 = (2.0 * d[52]) * ((-d[2032] - d[2033] * d[259] + d[321]) * Dy +
                           (-d[2597]) * dv_0 - dv_3869) -
          dv_3870;
  DataVector& dv_3862 = temps.at(368);
  DataVector& dv_3864 = temps.at(1175);
  DataVector& dv_3867 = temps.at(1135);
  DataVector& dv_3868 = temps.at(1646);
  sc_33 += d[50] * (d[20] * dv_3867 + d[77] * (d[2125] * dv_2170 - dv_3868) +
                    d[91] * dv_3864) +
           d[55] * dv_3862;
  sc_12 = d[237] * sc_33;
  DataVector& dv_3455 = temps.at(2930);
  DataVector& dv_3563 = temps.at(3036);
  DataVector& dv_3888 = temps.at(1599);
  DataVector& dv_3889 = temps.at(3118);
  DataVector& dv_3892 = temps.at(3127);
  sc_6 = (-d[361]) * dv_3563 + (d[1462] * d[1466]) * dv_0 +
         d[122] * ((-d[1450]) * (dv_3455 + dv_752) - dv_3892) +
         d[123] * (-dv_3888 - dv_3889);
  DataVector& dv_3894 = temps.at(3171);
  DataVector& dv_3895 = temps.at(1231);
  DataVector& dv_3897 = temps.at(3177);
  sc_6 +=
      d[124] * (-dv_3888 + dv_3894) +
      d[127] * ((d[1416] + d[1433] + d[2003]) * dv_2118 +
                (-24.0 * d[1339] - d[1986] * d[287] - d[2615] * d[2923] +
                 d[463] * d[50]) *
                    Dy +
                dv_3895) +
      d[128] *
          ((d[115] * d[7] + 72.0 * d[1339] + d[1990] + d[1992] * d[287]) * Dy +
           (-d[1979]) * dv_1696 + dv_3897);
  DataVector& dv_1770 = temps.at(1572);
  DataVector& dv_3886 = temps.at(1141);
  DataVector& dv_3887 = temps.at(1252);
  sc_6 += d[299] * ((d[1313] * (d[243] + d[490])) * dv_0 + d[2614] * dv_1770 +
                    dv_3887) +
          d[362] * ((-d[1424] - d[309]) * dv_1530 + d[2613] * dv_752 + dv_3886);
  sc_33 = d[2616] * sc_6;
  DataVector& dv_2672 = temps.at(2234);
  DataVector& dv_3915 = temps.at(1198);
  DataVector& dv_3938 = temps.at(1015);
  DataVector& dv_3940 = temps.at(3120);
  sc_13 =
      (d[1466] * d[2921]) *
          (d[259] * ((d[2200] - 124.0 * d[6]) * dv_2170 + dv_3938) + dv_3940) +
      (d[1675] * d[747]) * ((-d[532]) * Dy + dv_1115 + dv_2189) +
      (d[1931] * d[2923]) * (d[2674] * dv_2672 + dv_3915);
  DataVector& dv_3925 = temps.at(1072);
  DataVector& dv_3939 = temps.at(3162);
  DataVector& dv_3941 = temps.at(1188);
  sc_13 += d[1800] * ((-d[2681] + d[2682]) * dv_1696 +
                      (-d[2282] + 72.0 * d[2921] * d[2923]) * dv_3941 +
                      (-d[1485] + d[259] * d[2680] + 180.0 * d[319]) * dv_5) +
           d[1804] * (d[273] * (dv_3925 + dv_3938) + dv_3939);
  sc_13 += d[1821] * ((d[1688] - 126.0 * d[274]) * dv_1696 + Dy * d[2689] +
                      d[2313] * dv_1115) +
           d[1951] * ((-d[2684] - d[2685]) * dv_1696 + Dy * d[2686] +
                      d[2683] * dv_1115);
  DataVector& dv_3245 = temps.at(2732);
  DataVector& dv_3910 = temps.at(1130);
  sc_13 += d[2004] * ((-d[1716] * d[6] - d[2233] * d[273] + d[2234]) * Dx +
                      d[2688] * dv_3910) +
           d[2011] * ((d[1711] * d[286] - d[2237] * d[273] + d[2691]) * Dx +
                      d[2690] * dv_3910) +
           d[2673] * dv_3245;
  sc_13 += d[2678] * ((2.0 * d[1681] * d[2921] - 245.0 * d[36]) * dv_1564 +
                      (-d[2676]) * Dy + Dy * d[2677]);
  sc_6 = d[377] * sc_13;
  sc_28 = (d[1514] * d[340]) * dv_3915 + (d[1830] * d[2058]) * dv_752 +
          d[1177] * ((-d[1549] * d[2661] - d[1981] * d[2660] + d[2317] +
                      d[252] * d[811] + d[885] * d[2923]) *
                         Dy +
                     (56.0 * d[2928] * d[2921] * d[2923] - d[110] - d[370]) *
                         dv_3637 +
                     d[2088] * dv_0);
  DataVector& dv_1700 = temps.at(1530);
  DataVector& dv_3851 = temps.at(1093);
  DataVector& dv_3923 = temps.at(3176);
  DataVector& dv_3927 = temps.at(1027);
  DataVector& dv_3930 = temps.at(3174);
  sc_28 += d[219] * ((-d[308]) * (d[2658] * dv_2298 + dv_3927) +
                     d[223] * (dv_1700 + dv_3851) +
                     d[683] * (d[2657] * dv_2170 + dv_3923) + dv_3930);
  DataVector& dv_3874 = temps.at(1179);
  DataVector& dv_3916 = temps.at(1104);
  DataVector& dv_3919 = temps.at(1077);
  sc_28 +=
      d[2619] * ((-d[278]) * dv_263 + d[2058] * dv_3916 + d[2637] * dv_240) +
      d[2645] * (Dy * d[2644] + d[380] * dv_3874 + dv_3919);
  DataVector& dv_3801 = temps.at(3134);
  DataVector& dv_3917 = temps.at(1070);
  DataVector& dv_3920 = temps.at(663);
  DataVector& dv_3921 = temps.at(3105);
  DataVector& dv_3926 = temps.at(3114);
  sc_28 += d[2649] * (d[1172] * (d[1620] * dv_752 + dv_3921) +
                      d[1549] * (dv_3917 + dv_3926) + d[223] * dv_3801 +
                      d[683] * (dv_3923 + dv_3925) + dv_3920);
  DataVector& dv_3918 = temps.at(1050);
  sc_28 += d[2656] * ((-d[2651]) * dv_3098 + Dy * d[2655] + d[2654] * dv_0) +
           d[300] * (d[1105] * dv_3852 + d[287] * (d[2639] * dv_751 + dv_752) +
                     d[631] * (dv_3917 + dv_3918) + dv_3882);
  sc_13 = d[392] * sc_28;
  sc_30 = (-d[1989]) * ((d[1654] + d[2693] + 300.0 * d[274]) * dv_0 +
                        (d[136] * d[2694] - d[1837] - 262.0 * d[252]) * dv_5 +
                        207.0 * dv_3884);
  DataVector& dv_3948 = temps.at(3156);
  DataVector& dv_3950 = temps.at(3170);
  DataVector& dv_3951 = temps.at(3213);
  DataVector& dv_3952 = temps.at(1254);
  sc_30 += (-d[2669]) *
           ((-d[1641]) * dv_3950 + d[1550] * (-dv_3951 + 104.0 * dv_751) +
            d[2293] * ((d[1825] + d[2318]) * Dx + dv_3948) + dv_3952);
  DataVector& dv_2350 = temps.at(1940);
  DataVector& dv_3832 = temps.at(240);
  DataVector& dv_3907 = temps.at(3144);
  DataVector& dv_3943 = temps.at(3198);
  sc_30 += (d[1618] * d[2920]) * ((d[1637] + d[1819] + d[564]) * Dy +
                                  d[557] * dv_3832 + dv_3907) +
           (3.0 * d[1931] * d[554]) * dv_752 +
           d[1581] * ((-18.0 * d[1706] - 18.0 * d[30]) * dv_2350 + dv_3943);
  DataVector& dv_3908 = temps.at(3172);
  DataVector& dv_3945 = temps.at(1111);
  DataVector& dv_3946 = temps.at(1079);
  DataVector& dv_3947 = temps.at(3145);
  sc_30 += d[1845] * (d[2293] * (2.0 * dv_3908 - dv_3945) +
                      d[996] * (dv_1700 + 131.0 * dv_752) + dv_3946 + dv_3947);
}
}  // namespace CurvedScalarWave::Worldtube::detail
