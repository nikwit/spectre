
// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/Worldtube/PunctureFieldOrder2Impl.hpp"

namespace CurvedScalarWave::Worldtube::detail {

// NOLINTNEXTLINE(google-readability-function-size, readability-function-size)
void puncture_field_2_part_7(const std::array<double, order2_n_doubles>& d,
                             const DataVector& Dx, const DataVector& Dy,
                             DynamicBuffer<DataVector>& temps) {
  DataVector& dv_1803 = temps.at(1589);
  DataVector& dv_1805 = temps.at(1591);
  DataVector& sc_28 = temps.at(3250);
  DataVector& sc_29 = temps.at(3251);
  sc_28 += d[51] * (d[1681] * dv_1805 + 45.0 * dv_1803) + sc_29;
  DataVector& sc_27 = temps.at(3249);
  sc_27 = d[362] * sc_28;
  DataVector& dv_1 = temps.at(1);
  DataVector& dv_1108 = temps.at(1005);
  DataVector& dv_1237 = temps.at(1132);
  DataVector& dv_2760 = temps.at(2313);
  DataVector& dv_2815 = temps.at(2080);
  DataVector& dv_5 = temps.at(5);
  DataVector& sc_20 = temps.at(3242);
  DataVector& sc_21 = temps.at(3243);
  sc_20 = (d[1618] * (d[1207] + d[298])) * dv_1237 +
          (d[340] * d[6]) * ((-d[162]) * dv_5 + d[50] * dv_2815 +
                             d[807] * dv_1 - 245.0 * dv_1108) +
          dv_2760 + sc_21;
  DataVector& dv_2813 = temps.at(2366);
  DataVector& dv_2816 = temps.at(2368);
  DataVector& dv_2817 = temps.at(2369);
  DataVector& dv_2818 = temps.at(2370);
  DataVector& dv_2820 = temps.at(2372);
  DataVector& sc_15 = temps.at(3237);
  DataVector& sc_17 = temps.at(3239);
  DataVector& sc_23 = temps.at(3245);
  DataVector& sc_25 = temps.at(3247);
  DataVector& sc_26 = temps.at(3248);
  DataVector& sc_9 = temps.at(3231);
  sc_20 += d[1581] * ((-d[1213] * d[3]) + d[1244] * dv_2818 + dv_2820 +
                      d[2922] * (d[261] * dv_2816 + dv_2817)) +
           d[1675] * dv_2813 + sc_15 + sc_17 + sc_23 + sc_25 + sc_26 + sc_27 +
           sc_9;
  DataVector& dv_905 = temps.at(837);
  DataVector& sc_24 = temps.at(3246);
  sc_24 = -dv_905 * sc_20;
  DataVector& dv_1502 = temps.at(1382);
  DataVector& dv_2500 = temps.at(2088);
  DataVector& dv_2501 = temps.at(2089);
  DataVector& dv_2502 = temps.at(2090);
  sc_9 = d[1523] +
         d[20] * ((-201.0 * d[255]) +
                  d[2928] * (80.0 * dv_1502 + dv_2502 + 8.0) - dv_2501) +
         d[49] * (44.0 * dv_1 + dv_2500);
  DataVector& dv_2390 = temps.at(1980);
  DataVector& dv_605 = temps.at(556);
  DataVector& dv_794 = temps.at(742);
  sc_9 += d[77] * (d[1524] * (1.0 - dv_2390) - dv_605 + 193.0 * dv_794);
  sc_23 = (-d[2922]) * sc_9;
  DataVector& dv_1700 = temps.at(1530);
  DataVector& dv_2115 = temps.at(1732);
  DataVector& dv_2344 = temps.at(1934);
  DataVector& dv_2454 = temps.at(2043);
  DataVector& dv_2495 = temps.at(2083);
  DataVector& dv_2496 = temps.at(2084);
  DataVector& dv_2497 = temps.at(2085);
  DataVector& dv_2498 = temps.at(2086);
  DataVector& dv_751 = temps.at(699);
  DataVector& dv_763 = temps.at(711);
  sc_25 = (d[1521] + d[1522] * d[48] +
           d[259] * (28.0 * d[1242] - 157.0 * d[2923])) *
              dv_2115 +
          (-d[1449] * d[88]) * (d[1520] + dv_763) + d[1333] * dv_2344 +
          d[319] * (dv_2495 - dv_2496 + 36.0 * dv_751) +
          d[436] * (dv_1700 + dv_2454 + dv_2497 + dv_2498) + sc_23;
  DataVector& dv_2494 = temps.at(2082);
  sc_25 += d[49] * dv_2494;
  sc_26 = d[122] * sc_25;
  DataVector& dv_2375 = temps.at(1965);
  DataVector& dv_2376 = temps.at(1966);
  DataVector& dv_2423 = temps.at(2013);
  DataVector& dv_2446 = temps.at(1398);
  sc_23 = (-d[62]) * dv_2446 +
          (-d[2922]) * (Dx * d[1510] +
                        d[2921] * (39.0 * dv_2375 + 80.0 * dv_2376 + dv_2423));
  DataVector& dv_1229 = temps.at(1124);
  DataVector& dv_2227 = temps.at(1829);
  DataVector& dv_2488 = temps.at(2076);
  sc_23 +=
      d[3] * (d[2928] * (-39.0 * dv_1502 - dv_2488 - 8.0) + d[1509] + dv_2227) +
      d[6] * ((-d[62]) - 313.0 * dv_1229 + d[2921] * (84.0 * Dy + d[1437]));
  DataVector& dv_2486 = temps.at(2074);
  DataVector& dv_69 = temps.at(69);
  sc_23 += d[9] * (dv_2486 + dv_69);
  sc_25 = d[299] * sc_23;
  DataVector& dv_2529 = temps.at(2117);
  DataVector& dv_2534 = temps.at(2122);
  DataVector& dv_2535 = temps.at(2123);
  DataVector& dv_2536 = temps.at(2124);
  DataVector& dv_2537 = temps.at(2125);
  DataVector& dv_538 = temps.at(493);
  sc_17 = (-d[1443]) * ((-d[2928]) * (dv_2536 + dv_2537 + 21.0) +
                        (93.0 * d[255]) + dv_2534 + dv_2535) +
          d[1538] + d[391] * (d[31] + dv_2529 + dv_538);
  DataVector& dv_2234 = temps.at(1835);
  DataVector& dv_2510 = temps.at(2098);
  DataVector& dv_2530 = temps.at(2118);
  DataVector& dv_2531 = temps.at(2119);
  DataVector& dv_2533 = temps.at(2121);
  sc_17 +=
      d[50] * ((-d[1539]) + d[2928] * (72.0 * dv_1502 + dv_2531) - dv_2530) +
      d[724] * ((-d[31]) * (dv_2510 + dv_2533) + (d[147] * d[512]) - dv_2234 +
                241.0 * dv_794);
  sc_15 = (-d[86]) * sc_17;
  DataVector& dv_2392 = temps.at(1982);
  DataVector& dv_2463 = temps.at(2052);
  DataVector& dv_2523 = temps.at(2111);
  DataVector& dv_2524 = temps.at(2112);
  DataVector& dv_2528 = temps.at(2116);
  sc_9 = (d[1247] + d[147] + d[504]) * dv_2524 +
         (-174.0 * d[1339] + d[1443] * (49.0 * d[7] + 26.0) + d[1537] +
          d[273] * (14.0 * d[1242] - 241.0 * d[2923])) *
             dv_2528 +
         (d[1533] * d[1534]) * dv_2463 +
         d[1135] * ((120.0 * d[2923]) * Dx - 72.0 * dv_2375 - dv_2496) +
         d[1532] * dv_2392 + dv_2523 + sc_15;
  DataVector& dv_2455 = temps.at(2044);
  DataVector& dv_2525 = temps.at(2113);
  DataVector& dv_2526 = temps.at(2114);
  DataVector& dv_2527 = temps.at(2115);
  sc_9 += d[1536] * (13.0 * dv_2376 - dv_2455 + dv_2525 + dv_2526 + dv_2527) +
          d[762] * (dv_2375 - 47.0 * dv_2376 + 39.0 * dv_751);
  sc_23 = d[52] * sc_9;
  DataVector& dv_2548 = temps.at(2136);
  DataVector& dv_2550 = temps.at(2138);
  DataVector& dv_728 = temps.at(677);
  sc_21 = (-d[1549]) * ((d[147] * d[1548]) + d[1524] * dv_2550 - dv_728 +
                        157.0 * dv_794) +
          (-d[355]) * (d[37] + dv_2548 + 174.0 * dv_794) + (-d[120] * d[1313]) +
          d[1547] * dv_1;
  DataVector& dv_1503 = temps.at(1383);
  DataVector& dv_2246 = temps.at(1844);
  DataVector& dv_2549 = temps.at(2137);
  DataVector& dv_2551 = temps.at(2139);
  DataVector& dv_2552 = temps.at(2140);
  sc_21 +=
      d[1550] * ((89.0 * d[255]) + d[2928] * (-dv_2551 + dv_2552) +
                 140.0 * dv_2246 + dv_2534) +
      d[57] * (d[2928] * (-39.0 * dv_1503 - dv_2549 - 8.0) + d[1455] + dv_2227);
  sc_17 = sc_21 * d[2922];
  DataVector& dv_2541 = temps.at(2129);
  DataVector& dv_2542 = temps.at(2130);
  DataVector& dv_2543 = temps.at(2131);
  DataVector& dv_2544 = temps.at(2132);
  sc_28 = (-d[1527]) * (85.0 * dv_2376 + dv_2495 - dv_2544) +
          d[276] * ((-d[294]) * dv_1803 + (-d[2924]) * dv_2541 +
                    (84.0 * d[2923]) * Dx) +
          d[391] * (-dv_2542 - dv_2543);
  DataVector& dv_2379 = temps.at(1969);
  DataVector& dv_2439 = temps.at(2029);
  DataVector& dv_2545 = temps.at(2133);
  DataVector& dv_2546 = temps.at(2134);
  sc_28 += d[724] * (-60.0 * dv_2379 + dv_2439 + dv_2527 + dv_2545 + dv_2546) +
           d[780] * dv_1803;
  sc_21 = sc_28 * d[2921];
  DataVector& dv_2538 = temps.at(2126);
  DataVector& dv_2540 = temps.at(2128);
  sc_15 = (-193.0 * d[1543] + d[1546] - 184.0 * d[992]) * dv_2115 +
          d[1542] * (d[259] * dv_2538 + dv_2540) + sc_17 + sc_21;
  sc_9 = d[53] * sc_15;
  DataVector& dv_2517 = temps.at(2105);
  DataVector& dv_2520 = temps.at(2108);
  DataVector& dv_2521 = temps.at(2109);
  DataVector& dv_2522 = temps.at(2110);
  sc_28 = (-d[1443]) * (85.0 * dv_2375 + dv_2521 - dv_2522 - 102.0 * dv_751) +
          (d[2928] * d[20]) * ((d[1531] - 595.0 * d[7] + 40.0) * Dx + dv_2520) -
          170.0 * dv_2517;
  DataVector& dv_2452 = temps.at(2041);
  DataVector& dv_2518 = temps.at(2106);
  sc_28 += d[50] * ((156.0 * d[2923]) * Dx - dv_2452 - dv_2518);
  sc_17 = sc_28 * d[2922];
  DataVector& dv_2157 = temps.at(1769);
  DataVector& dv_2349 = temps.at(1939);
  DataVector& dv_2506 = temps.at(2094);
  DataVector& dv_2507 = temps.at(2095);
  sc_21 = (-d[339]) * dv_2446 + (d[1520] * d[290]) * dv_2349 + d[1515] * dv_1 +
          d[1527] * (d[2928] * (dv_2157 - dv_2506) + dv_2507);
  DataVector& dv_2503 = temps.at(2091);
  DataVector& dv_2505 = temps.at(2093);
  sc_21 +=
      d[276] * ((-d[2928]) * (72.0 * dv_1503 + dv_2505) + d[1525] + dv_2503);
  DataVector& dv_2512 = temps.at(2100);
  DataVector& dv_2513 = temps.at(2101);
  DataVector& dv_2514 = temps.at(2102);
  DataVector& dv_2516 = temps.at(2104);
  DataVector& dv_787 = temps.at(735);
  sc_21 += d[6] * ((-d[273]) * (d[1529] + 1133.0 * dv_1 + dv_2513) +
                   d[1443] * ((89.0 * d[36]) + dv_2514 + dv_787) +
                   d[50] * (d[1528] + dv_2512) + dv_2516);
  DataVector& dv_2114 = temps.at(1731);
  DataVector& dv_2431 = temps.at(2021);
  DataVector& dv_2509 = temps.at(2097);
  DataVector& dv_2511 = temps.at(2099);
  DataVector& dv_582 = temps.at(534);
  sc_21 += d[724] * (d[2928] * (dv_2114 + dv_2431) + d[255] * dv_2511 -
                     dv_2509 + dv_582) +
           sc_17;
  sc_15 = d[55] * sc_21;
  DataVector& dv_231 = temps.at(228);
  sc_29 = (-d[1331] * (d[1555] - 7.0)) * dv_231 + d[1547] * dv_751 +
          d[308] * ((d[1556] - 1133.0 * d[7]) * Dx + dv_2520) +
          d[57] * (-dv_2518 + dv_2521 + 144.0 * dv_751);
  sc_29 += d[996] * ((-d[1557]) * dv_1803 + dv_2376 + dv_2522 + 51.0 * dv_751);
  sc_28 = sc_29 * d[2922];
  DataVector& dv_2247 = temps.at(1845);
  DataVector& dv_2504 = temps.at(2092);
  DataVector& dv_2553 = temps.at(2141);
  DataVector& dv_2554 = temps.at(2142);
  DataVector& dv_2555 = temps.at(2143);
  DataVector& dv_2556 = temps.at(2144);
  DataVector& dv_2557 = temps.at(2145);
  sc_17 = (-d[1532]) * dv_2446 + (-d[1533]) * dv_2554 +
          (-d[1553]) * (dv_2247 - dv_2556 + dv_2557) +
          d[1135] * ((-d[2928]) * (80.0 * dv_1503 + dv_2504 + 8.0) + d[1552] +
                     dv_2555) +
          dv_2553;
  DataVector& dv_1667 = temps.at(1458);
  DataVector& dv_2359 = temps.at(1949);
  DataVector& dv_2477 = temps.at(2066);
  DataVector& dv_2558 = temps.at(2146);
  DataVector& dv_2559 = temps.at(2147);
  DataVector& dv_2560 = temps.at(2148);
  DataVector& dv_2561 = temps.at(2149);
  sc_17 += d[1536] * (d[2928] * (dv_2359 + dv_2477) + d[255] * dv_2561 +
                      dv_1667 - 50.0 * dv_2246) +
           d[1554] * ((-d[2928]) * (dv_2559 + dv_2560) + d[153] + dv_2558);
  DataVector& dv_1685 = temps.at(1517);
  DataVector& dv_2562 = temps.at(2150);
  sc_17 += d[35] * ((-d[273]) * (d[1088] + 595.0 * dv_1 + dv_2513) +
                    (-d[313]) * (d[88] + dv_2562) +
                    (3.0 * d[50]) * ((67.0 * d[36]) + dv_1685) +
                    (2.0 * d[48] * d[2921]) *
                        ((93.0 * d[36]) + 78.0 * Dy + 196.0 * dv_794) +
                    (-d[780]));
  sc_17 += sc_28;
  sc_21 = d[63] * sc_17;
  DataVector& dv_1632 = temps.at(1477);
  DataVector& dv_2482 = temps.at(2002);
  sc_27 = (-3.0 * d[362]) *
              (d[1507] * dv_1805 + d[2922] * (d[1376] + dv_1632 + dv_2482)) +
          sc_25 + sc_26;
  DataVector& dv_2350 = temps.at(1940);
  DataVector& dv_2489 = temps.at(2077);
  DataVector& dv_2490 = temps.at(2078);
  DataVector& dv_2491 = temps.at(2079);
  DataVector& dv_2493 = temps.at(2081);
  sc_27 +=
      d[50] * ((d[1443] * (d[1242] + d[1325]) + d[1515] +
                d[273] * (8.0 - 313.0 * d[7]) +
                d[58] * (-d[1434] * d[2925] - d[1516])) *
                   dv_2350 +
               (d[1] * d[380]) * dv_2349 + (d[1514] * d[273]) * dv_1502 +
               d[1053] * ((d[1517] * d[57]) + d[50] * ((-d[1518]) - dv_2491) +
                          d[631] * dv_2493 + dv_2489 - 170.0 * dv_2490)) +
      sc_23;
  sc_27 += sc_15 + sc_21 + sc_9;
  DataVector& dv_925 = temps.at(846);
  sc_20 = -dv_925 * sc_27;
  DataVector& dv_2389 = temps.at(1979);
  DataVector& dv_2407 = temps.at(1997);
  DataVector& dv_2752 = temps.at(2305);
  DataVector& dv_2794 = temps.at(2347);
  DataVector& dv_2797 = temps.at(2350);
  DataVector& dv_600 = temps.at(551);
  sc_25 = (-d[1172]) *
              ((198.0 * d[255]) + d[2928] * (-dv_2407 - dv_2797) + dv_2752) +
          (-d[1656]) * (d[37] + dv_600 - 189.0 * dv_794) +
          (-d[50] * d[558]) * (dv_2389 + dv_2794) + (108.0 * d[1663]);
  DataVector& dv_2741 = temps.at(2294);
  DataVector& dv_2795 = temps.at(2348);
  DataVector& dv_2796 = temps.at(2349);
  sc_25 += d[230] * ((2.0 * d[2928]) * (50.0 * dv_1503 + dv_2796 + 9.0) +
                     (-d[1664]) - dv_2795) +
           d[384] * dv_2741;
  sc_23 = sc_25 * d[2922];
  DataVector& dv_1112 = temps.at(1009);
  DataVector& dv_1839 = temps.at(1624);
  DataVector& dv_2268 = temps.at(1864);
  DataVector& dv_2373 = temps.at(1963);
  sc_9 = (d[1658] * d[50] + d[1662] + 1512.0 * d[992]) * dv_2268 +
         (-d[1645]) * (dv_1805 + dv_2373) + (192.0 * d[384]) * dv_1112 +
         d[1477] * ((d[101] + d[1657]) + Dy * d[398] + dv_1839);
  DataVector& dv_2780 = temps.at(2333);
  sc_9 += d[1536] *
              (d[1529] * dv_1803 - 10.0 * dv_2376 + 315.0 * dv_2379 + dv_2780) +
          d[1624] * (140.0 * dv_2375 + 100.0 * dv_2376 - 189.0 * dv_751);
  DataVector& dv_2723 = temps.at(2276);
  sc_9 += d[1656] * (6.0 * dv_2379 + dv_2497 + dv_2546 - 25.0 * dv_751) +
          d[704] * (96.0 * dv_2376 + dv_2723 + 65.0 * dv_751) + sc_23;
  sc_15 = (-d[124]) * sc_9;
  DataVector& dv_1504 = temps.at(1384);
  DataVector& dv_2731 = temps.at(2284);
  DataVector& dv_2775 = temps.at(2328);
  DataVector& dv_2798 = temps.at(2351);
  DataVector& dv_2799 = temps.at(2352);
  DataVector& dv_2800 = temps.at(2353);
  sc_23 =
      (-d[1665]) *
          ((-d[2928]) * (dv_2731 + dv_2799 + 9.0) + (67.0 * d[255]) + dv_2775) +
      (30.0 * d[1651]) * dv_2798 +
      d[1549] * ((-d[1668]) * dv_1504 + d[1667] + 1005.0 * dv_2246 - dv_2800) -
      dv_2553;
  DataVector& dv_2369 = temps.at(1959);
  sc_23 += d[1649] * dv_2446 +
           d[1666] * (d[2928] * (192.0 * dv_1502 + dv_2369 + 25.0) + d[1403] -
                      22.0 * dv_1);
  DataVector& dv_1555 = temps.at(1429);
  DataVector& dv_2801 = temps.at(2354);
  sc_23 += d[1669] *
           ((8.0 * d[2928] * d[2924]) * Dx - 75.0 * dv_1 - dv_1555 - dv_2801);
  DataVector& dv_2802 = temps.at(2355);
  DataVector& dv_785 = temps.at(733);
  sc_23 += d[35] * ((-d[1611]) * ((-d[1564]) - dv_785) +
                    (-d[58]) * ((236.0 * d[36]) + dv_2802) + d[1659] +
                    d[391] * (d[1464] + 453.0 * dv_1) +
                    d[724] * (d[1481] + 885.0 * dv_1));
  DataVector& dv_2377 = temps.at(1967);
  DataVector& dv_2709 = temps.at(2262);
  DataVector& dv_2803 = temps.at(2356);
  DataVector& dv_2804 = temps.at(2357);
  sc_23 += d[2922] *
           ((561.0 * d[7] - 125.0) * dv_2804 +
            (d[2928] * (118.0 * d[7] - 9.0)) * dv_2803 +
            d[1172] * (192.0 * dv_2375 + dv_2377 + dv_2709) + d[1646] * dv_751 +
            d[230] * (108.0 * dv_2375 + 256.0 * dv_2376 - 297.0 * dv_751));
  sc_9 = (-d[127]) * sc_23;
  DataVector& dv_2811 = temps.at(2364);
  DataVector& dv_2812 = temps.at(2365);
  sc_17 = (128.0 * d[384]) * dv_751 +
          (d[1161] * (151.0 * d[7] - 25.0)) * dv_231 +
          d[1672] * ((29.0 * d[2928] * d[2925]) * Dx - dv_2497 - dv_2812) +
          d[1674] * dv_2811;
  sc_17 += d[230] * (d[1474] * dv_1803 + 124.0 * dv_2376 - 153.0 * dv_751);
  sc_26 = sc_17 * d[2922];
  DataVector& dv_2250 = temps.at(1848);
  DataVector& dv_2805 = temps.at(2358);
  sc_25 = (-d[1549]) *
              ((110.0 * d[255]) + d[1] * (67.0 * dv_1502 - 55.0 * dv_1503) -
               714.0 * dv_2246 + dv_2556) +
          (-d[1665]) * ((-d[2928]) * (68.0 * dv_1503 + dv_2805 + 3.0) +
                        (34.0 * d[255]) + dv_2250);
  DataVector& dv_1570 = temps.at(1442);
  DataVector& dv_2591 = temps.at(2175);
  DataVector& dv_2808 = temps.at(2361);
  sc_25 += (-d[1669]) * (d[1460] + 125.0 * dv_1 + dv_1570 - dv_2591 + dv_2808) +
           (-72.0 * d[1670]) * dv_2798 + (256.0 * d[384]) * dv_794 +
           d[1645] * dv_2446;
  DataVector& dv_2806 = temps.at(2359);
  sc_25 +=
      d[1666] * (d[2928] * (87.0 * dv_1502 + dv_2806) + d[1671] + 165.0 * dv_1);
  DataVector& dv_2044 = temps.at(1689);
  DataVector& dv_2508 = temps.at(2096);
  DataVector& dv_2810 = temps.at(2363);
  sc_25 += d[6] * ((-d[1672]) * (d[378] + dv_2044) +
                   (-d[230]) * ((98.0 * d[36]) + dv_2508) +
                   d[1549] * (d[1617] + 1101.0 * dv_1) +
                   d[1669] * (d[1] + 561.0 * dv_1) + dv_2810) +
           sc_26;
  sc_23 = (-d[128]) * sc_25;
  DataVector& dv_2616 = temps.at(2198);
  DataVector& dv_2788 = temps.at(2341);
  DataVector& dv_2789 = temps.at(2342);
  DataVector& dv_2791 = temps.at(2344);
  DataVector& dv_2793 = temps.at(2346);
  sc_28 = (-d[1107]) * ((-d[2928]) * (dv_2791 + dv_2793) + d[1655] + dv_2616) +
          (-d[1254]) * (d[1524] + dv_2508 - dv_2789) + (162.0 * d[1135]) +
          d[272] * dv_2788;
  DataVector& dv_2745 = temps.at(2298);
  DataVector& dv_2790 = temps.at(2343);
  sc_28 += d[58] * ((2.0 * d[2928]) * (54.0 * dv_1503 + dv_2745 + 9.0) +
                    (-236.0 * d[255]) - dv_2790);
  sc_17 = (-d[86]) * sc_28;
  DataVector& dv_2008 = temps.at(167);
  DataVector& dv_2719 = temps.at(2272);
  DataVector& dv_2784 = temps.at(2337);
  DataVector& dv_2785 = temps.at(2338);
  sc_26 =
      (d[1607] + d[1650]) * dv_2785 +
      (d[162] * d[319]) * ((-d[1240]) * dv_2719 + dv_2008 + dv_2546) +
      (-d[1248] * (712.0 * d[1339] + d[1654] + 568.0 * d[274])) * dv_231 +
      d[1624] * ((135.0 * d[2923]) * Dx - 256.0 * dv_2375 - 108.0 * dv_2376) +
      d[1649] * dv_2784 + dv_2523 + sc_17;
  DataVector& dv_2437 = temps.at(2027);
  DataVector& dv_2786 = temps.at(2339);
  sc_26 += d[1652] * (d[1651] + dv_2786) +
           d[1653] * (-192.0 * dv_2379 - dv_2437 + dv_2454 + 27.0 * dv_751);
  sc_25 = d[122] * sc_26;
  DataVector& dv_1683 = temps.at(1516);
  DataVector& dv_2331 = temps.at(1922);
  DataVector& dv_2457 = temps.at(2046);
  DataVector& dv_2461 = temps.at(2050);
  DataVector& dv_2470 = temps.at(2059);
  DataVector& dv_2783 = temps.at(2336);
  sc_29 =
      (-d[1172]) * (d[1648] - 87.0 * dv_2331 + dv_2783) +
      (-d[223]) * ((17.0 * d[255]) + d[2928] * (-dv_2461 - dv_2470) + dv_1683) +
      (d[2928] * d[1333]) * (d[1284] + dv_2457 + 94.0 * dv_794) +
      (d[1645] * d[2923]);
  DataVector& dv_2782 = temps.at(2335);
  sc_29 += d[1646] * dv_1 + d[1647] * (d[36] + dv_2782 - dv_600);
  sc_28 = (-d[2922]) * sc_29;
  DataVector& dv_1634 = temps.at(1479);
  DataVector& dv_2779 = temps.at(2332);
  DataVector& sc_30 = temps.at(3252);
  sc_30 = (-d[1644]) * (-dv_1634 - dv_2779) +
          d[1393] * (-30.0 * dv_2376 - dv_2495 - dv_2780);
  DataVector& dv_2399 = temps.at(1989);
  DataVector& dv_2438 = temps.at(2028);
  DataVector& dv_2777 = temps.at(2330);
  DataVector& dv_2778 = temps.at(2331);
  DataVector& dv_2781 = temps.at(2334);
  sc_30 += d[1472] * (dv_2375 - dv_2399 + dv_2438 + dv_2781 * d[2924]) +
           dv_2777 * d[2925] + dv_2778 * d[2924];
  sc_29 = sc_30 * d[2921];
  DataVector& dv_613 = temps.at(564);
  sc_17 = (-84.0 * d[1543] + d[1643]) * dv_2115 +
          d[1640] * ((8.0 * d[48]) * Dy + (-d[50]) - 62.0 * dv_613) + sc_28 +
          sc_29;
  sc_26 = d[123] * sc_17;
  DataVector& dv_1720 = temps.at(1540);
  DataVector& dv_2383 = temps.at(1973);
  DataVector& dv_2393 = temps.at(1983);
  DataVector& dv_2765 = temps.at(2318);
  DataVector& dv_2766 = temps.at(2319);
  sc_29 =
      (-d[1191]) * (d[2928] * (dv_2765 + dv_2766 + 3.0) + d[255] - dv_1683) +
      (-d[1625]) *
          (d[266] * dv_1112 + d[2921] * (-dv_1720 + dv_2383 + dv_2393)) +
      (-d[416]) * dv_2446;
  DataVector& dv_2434 = temps.at(2024);
  DataVector& dv_2763 = temps.at(2316);
  DataVector& dv_2764 = temps.at(2317);
  sc_29 +=
      (2.0 * d[2928]) * dv_2434 +
      (3.0 * d[6]) * (d[524] + dv_2763 + d[2921] * ((34.0 * d[36]) + dv_2764));
  sc_17 = d[1466] * sc_29;
  DataVector& dv_2772 = temps.at(2325);
  DataVector& dv_2776 = temps.at(2329);
  DataVector& sc_31 = temps.at(3253);
  sc_31 = (-d[1632]) * dv_2446 +
          (-d[77]) * (d[1] * (55.0 * dv_1502 - 67.0 * dv_1503) + d[1639] -
                      dv_2775 + dv_2776) +
          d[1637] * ((67.0 * d[2923]) * Dy - 87.0 * dv_2772);
  DataVector& dv_2773 = temps.at(2326);
  DataVector& dv_2774 = temps.at(2327);
  sc_31 += d[1638] *
           ((-d[2928]) * (50.0 * dv_1502 + dv_2774 + 9.0) + d[1606] + dv_2773);
  sc_30 = sc_31 * d[2921];
  DataVector& dv_2442 = temps.at(2032);
  DataVector& dv_2443 = temps.at(2033);
  DataVector& dv_2771 = temps.at(2324);
  sc_28 = (-d[2922]) *
              ((87.0 * d[1242] - 169.0 * d[2923]) * dv_2771 +
               (d[102] * (367.0 * d[7] - 21.0)) * dv_2443 +
               d[58] * (100.0 * dv_2375 + 140.0 * dv_2376 - 171.0 * dv_751) +
               696.0 * dv_2517) -
          1428.0 * dv_2442;
  DataVector& dv_1637 = temps.at(1482);
  DataVector& dv_2770 = temps.at(2323);
  sc_28 += (d[6] * d[2921]) * (d[104] * (d[1636] + dv_1637) + d[1633] +
                               d[1635] * ((-d[1634]) - dv_2556) +
                               d[348] * ((170.0 * d[36]) + dv_2770)) +
           sc_30;
  sc_29 = d[299] * sc_28;
  DataVector& dv_2371 = temps.at(1961);
  sc_30 = 256.0 * dv_2371;
  DataVector& dv_2701 = temps.at(2254);
  DataVector& dv_2768 = temps.at(2321);
  DataVector& dv_2769 = temps.at(2322);
  DataVector& dv_558 = temps.at(511);
  sc_30 +=
      (-d[2922]) * ((-d[348]) * ((-d[1]) * (68.0 * dv_1502 + dv_2701 + 3.0) +
                                 (98.0 * d[255]) + 153.0 * dv_1) +
                    (24.0 * d[2928] * d[2921]) * (d[1631] + dv_2769 - dv_558) +
                    (-d[1630]) - dv_2768);
  DataVector& dv_2345 = temps.at(1935);
  DataVector& dv_2767 = temps.at(2320);
  sc_30 += (d[1627] * (-d[1318] - d[1582])) * dv_231 +
           (d[265] * d[302]) * ((-d[1626]) - dv_728) +
           d[1251] * (dv_1803 + dv_2345) +
           d[1628] * ((27.0 * d[2923]) * Dx - 124.0 * dv_2375 - dv_2767);
  DataVector& dv_2710 = temps.at(2263);
  sc_30 += d[436] * (dv_1720 - 93.0 * dv_2379 + dv_2497 + dv_2710);
  sc_28 = d[362] * sc_30;
  DataVector& dv_2433 = temps.at(2023);
  DataVector& dv_2761 = temps.at(2314);
  sc_21 = (d[1234] + d[1619] + d[1621] * d[396]) * dv_2761 +
          (3.0 * d[1581]) * (d[554] * dv_1805 + d[2922] * (d[1623] + dv_2433)) +
          6.0 * dv_2760 + sc_15 + sc_23 + sc_9;
  DataVector& dv_2762 = temps.at(2315);
  sc_21 += (d[120] * d[638]) *
               ((-d[50]) * (d[1] + dv_2762) + (8.0 * d[151]) * Dy +
                (128.0 * d[2928] * d[20]) * Dy + (-d[1624]) - 348.0 * dv_2490) +
           (d[1618] * d[557]) * dv_1570 + sc_17 + sc_25 + sc_26 + sc_28 + sc_29;
  DataVector& dv_935 = temps.at(855);
  sc_27 = -dv_935 * sc_21;
  DataVector& dv_2112 = temps.at(1729);
  DataVector& dv_6 = temps.at(6);
  DataVector& dv_954 = temps.at(861);
  sc_17 = (d[1170] * d[481]) * dv_6 + (d[333] * (d[1384] + d[2])) * dv_6 +
          d[1065] * dv_954 + d[482] * dv_2112;
  DataVector& dv_0 = temps.at(0);
  DataVector& dv_1995 = temps.at(103);
  sc_17 += d[92] *
           ((-21.0 * d[496]) + d[50] * (d[1385] + dv_1700) +
            d[61] * (d[1386] + 42.0 * dv_0 + 11.0 * dv_1) + d[758] * dv_1995);
  sc_29 = (-d[484]) * sc_17;
  DataVector& dv_1509 = temps.at(1388);
  DataVector& dv_1666 = temps.at(212);
  DataVector& dv_300 = temps.at(292);
  DataVector& dv_871 = temps.at(807);
  DataVector& dv_958 = temps.at(865);
  sc_26 = (d[1343] * (d[1065] * d[382] + d[244] * (d[1391] + d[1394]))) * dv_6 +
          (d[491] * d[80]) * dv_1509 + (d[491] * d[2925]) * dv_871 +
          (-d[1161] * (-d[1387] - d[97])) * dv_300 +
          (-d[1389] * d[485]) * dv_6 + (-d[166] * d[4]) * dv_958 +
          d[487] * dv_1666;
  DataVector& dv_1514 = temps.at(1393);
  DataVector& dv_2208 = temps.at(1811);
  DataVector& dv_2209 = temps.at(1812);
  sc_26 +=
      d[488] *
      ((-d[52]) * dv_2208 + (22.0 * d[53]) * ((-d[1385]) - dv_1514) + d[1390] +
       d[50] * ((21.0 * d[2921] * d[2923]) - 22.0 * dv_0 - dv_2209));
  sc_17 = sc_26 * d[2922];
  DataVector& dv_2206 = temps.at(1802);
  DataVector& dv_2207 = temps.at(1810);
  DataVector& dv_289 = temps.at(281);
  DataVector& dv_884 = temps.at(816);
  DataVector& dv_931 = temps.at(851);
  DataVector& dv_955 = temps.at(862);
  sc_28 = (-d[1291]) * dv_955 + (-140.0 * d[4]) * dv_931 +
          (d[1130] * d[1375]) * dv_289 + (d[1374] * d[477]) * dv_1666 +
          (d[1379] * d[354]) * dv_2207 +
          (d[152] * (-d[1381] * d[244] + 8.0 * d[4] * d[92])) * dv_884 +
          (d[286] * d[85]) * dv_2207 +
          (d[535] * (8.0 * d[1220] + d[1383] * d[244])) * dv_884 +
          (-d[1006] * d[92]) * dv_2206 + (-d[389] * d[85]) * dv_2206 + sc_29;
  DataVector& dv_1911 = temps.at(1632);
  DataVector& dv_2106 = temps.at(1724);
  DataVector& dv_2110 = temps.at(1727);
  sc_28 += (d[1101] * d[1130] * d[1376]) * dv_6 + d[1131] * dv_2106 +
           d[1377] * dv_2106 + d[1377] * dv_2110 + d[1378] * dv_2207 +
           dv_1911 * d[2924] + sc_17;
  DataVector& dv_961 = temps.at(868);
  sc_21 = -dv_961 * sc_28;
  DataVector& dv_2096 = temps.at(841);
  DataVector& dv_2103 = temps.at(820);
  DataVector& dv_2105 = temps.at(1723);
  DataVector& dv_951 = temps.at(819);
  DataVector& dv_965 = temps.at(871);
  DataVector& dv_968 = temps.at(874);
  DataVector& sc_11 = temps.at(3233);
  DataVector& sc_2 = temps.at(3224);
  DataVector& sc_8 = temps.at(3230);
  sc_11 =
      (d[129] * (d[102] * d[2925] + d[1191] * d[1242] - 60.0 * d[1241] +
                 d[1269] * (-d[1318] - d[303]) + d[147] * d[46] +
                 d[155] * d[2921] + d[37] - d[431] * d[2925]) +
       d[191] * (20.0 * d[2928] * d[20] * d[2923] +
                 6.0 * d[2928] * d[280] * d[2924] * d[2922] * d[2921] -
                 d[1161] * d[1303] - d[1297] * d[1341] - 27.0 * d[1339] -
                 d[1340] * d[147] - 150.0 * d[1342] + 36.0 * d[50] * d[7] +
                 6.0 * d[57] * d[2925] -
                 d[6] * (d[1335] + 39.0 * d[1339] + 228.0 * d[274])) +
       d[299] * (d[1312] * d[2924] + d[1317] * d[2922]) +
       d[50] * (d[1327] * d[280] + d[1329] * d[1330] +
                d[2922] * (d[1317] * d[50] + d[1324] * d[49] * d[2921] -
                           d[1331] + d[631] * (4.0 - 141.0 * d[7]))) +
       d[55] * (d[1196] * d[1320] - 300.0 * d[289] +
                d[2922] * (-d[1322] * (d[1321] - 5.0) + 2.0 * d[1324] * d[48] -
                           d[1326] * d[348])) -
       d[61] * (6.0 * d[1241] * d[292] - d[1303] * d[1348] + 9.0 * d[1339] +
                81.0 * d[1342] - d[1343] * d[273] - d[1344] * d[151] -
                d[1345] * d[50] + d[1346] * d[1347] +
                d[1350] * (d[1349] + 189.0 * d[36]) + d[180] * d[559] -
                d[780] * d[2925] + d[866] * d[2925]) -
       d[63] * (d[1056] * d[1336] +
                d[1200] * (d[1107] * (d[1193] + d[1338]) + d[1215] * d[313] +
                           d[1326] * d[50] + d[631] * (d[1337] - 8.0)) +
                d[142] * d[292] * d[510])) *
          dv_965 +
      (d[1205] * d[55] + d[50] * (d[1206] + d[1207] * d[2922] + d[142]) +
       d[52] * (-8.0 * d[1208] + d[2923] * (d[1034] + d[1209])) +
       d[53] * (d[1211] - d[2923] * (d[1212] + d[155] + 13.0 * d[6] + 9.0)) -
       d[63] * (d[1213] + d[256] * d[2924] + d[2922] * (d[1212] + d[1215]))) *
          dv_968 +
      (-d[85]) * dv_2103 + (-d[85]) * dv_2105 + (-d[86]) * dv_2103 +
      (-d[86]) * dv_2105 + (-d[960]) * dv_2096 + (-d[961]) * dv_2096 +
      (-d[1205] * d[299] -
       d[122] * (-20.0 * d[1208] + 2.0 * d[1281] * d[2925] + d[147] -
                 d[535] * d[2923]) -
       d[50] * (d[1056] * d[1289] + d[1286] +
                d[2922] * (d[2928] * d[2921] * (d[1291] - 22.0 * d[2923]) -
                           d[1290] - d[244] + d[50] * d[2925])) +
       d[52] * (d[1031] * d[1292] + 6.0 * d[1208] * d[468] - d[1209] * d[1287] -
                d[1293] * (d[1242] - d[80]) + d[348] * (d[1193] + d[1294]) +
                d[35] * (103.0 * d[3] + d[88])) -
       d[53] * (-d[1031] * d[1300] - d[116] * d[1242] - d[1245] * d[48] -
                d[1269] * (d[430] * d[2923] + d[436] + d[900] * d[2923]) +
                2.0 * d[1285] * d[2924] * d[2922] * d[2921] - d[1301] * d[48] -
                50.0 * d[1302] - d[1303] * d[162] - d[180] * d[88]) -
       d[55] * (d[1283] * d[2924] -
                d[2922] * (-d[110] * d[2925] + d[1284] + 33.0 * d[180] +
                           50.0 * d[35] + d[471] * d[2925])) +
       d[63] * (d[1056] * d[1296] + 9.0 * d[142] * d[468] +
                d[2922] * (d[110] * d[1215] + 19.0 * d[1297] + d[1298] * d[7] +
                           d[1299] * (d[1106] + d[1193])))) *
          dv_951 +
      sc_2 + sc_8;
  DataVector& dv_1879 = temps.at(23);
  DataVector& dv_1886 = temps.at(213);
  DataVector& dv_1892 = temps.at(1390);
  DataVector& dv_1894 = temps.at(1628);
  DataVector& dv_1896 = temps.at(1629);
  DataVector& dv_2098 = temps.at(1718);
  DataVector& dv_2102 = temps.at(1721);
  DataVector& dv_963 = temps.at(869);
  DataVector& dv_978 = temps.at(810);
  DataVector& dv_982 = temps.at(839);
  sc_11 +=
      (d[121] * (d[1056] * d[498] +
                 d[1200] * (d[116] * (d[1264] + d[1265]) + d[1263] +
                            d[259] * (d[1266] + d[1267] - 27.0))) -
       21.0 * d[122] * d[2925] -
       d[20] * (d[1250] * (d[1249] * (d[1269] + d[1271] + 9.0) +
                           d[1259] * d[20] + d[505]) +
                d[1268] * d[29]) -
       d[206] *
           (36.0 * d[2928] * d[20] * d[2924] * d[2922] - d[162] * d[7] +
            4.0 * d[20] * d[2923] * (-d[1272] + 53.0 * d[2923]) -
            d[286] * (d[138] + 165.0 * d[159] - 74.0 * d[20]) -
            d[396] * (d[1193] + d[1246] + d[1257]) + 53.0 * d[50] * d[2925]) -
       d[55] * (d[1132] * d[2924] + d[1250] * d[1259]) -
       d[61] *
           (-36.0 * d[1241] + d[1242] * d[1262] + d[1260] + d[1261] * d[2921] +
            d[147] * d[42] + d[286] * (-51.0 * d[36] - d[493]) +
            d[497] * d[2925] + d[72] * d[2925])) *
          dv_963 +
      (d[1137] * d[1202]) * dv_2098 +
      (d[12] *
       (d[19] * (d[1056] * d[561] +
                 d[1250] * (-d[1249] * (d[1248] + d[298] - 9.0) +
                            d[348] * (-d[1246] - d[1247]) + d[749])) +
        d[20] * (d[1232] * d[557] +
                 d[2922] * (d[1233] + d[1234] + d[396] * (d[1235] + d[6]))) +
        d[206] * (d[1030] * d[1254] + d[1251] * d[2925] + d[1253] -
                  d[1256] * (d[1193] + d[1255]) +
                  d[286] * (-d[266] * d[3] + d[560]) +
                  d[77] * (10.0 * d[1242] + d[1245] - d[1257])) +
        d[226] * (d[1229] + d[1231] * d[2922]) +
        d[52] * (d[1236] + d[1237] + d[1238] + d[1239] * d[2921] -
                 24.0 * d[1241] + d[1242] * d[1243] +
                 d[1244] * (d[135] + d[36]) + d[416] * d[2925]))) *
          dv_982 +
      (d[1217] * d[83]) * dv_1886 + (d[1218] * d[352]) * dv_2102 +
      (d[1222] * d[131]) * dv_1894 + (d[1225] * d[336]) * dv_978 +
      (d[1273] * d[75]) * dv_1879 + (d[1304] * d[1305]) * dv_1892 +
      (d[1351] * d[71]) * dv_1896;
  DataVector& dv_12 = temps.at(12);
  DataVector& dv_1737 = temps.at(1549);
  DataVector& dv_1887 = temps.at(244);
  DataVector& dv_1899 = temps.at(694);
  DataVector& dv_1902 = temps.at(710);
  DataVector& dv_1905 = temps.at(753);
  DataVector& dv_2100 = temps.at(831);
  DataVector& dv_2104 = temps.at(1722);
  sc_11 += (-d[1219] * d[336]) * dv_2104 + (-d[1219] * d[476]) * dv_2102 +
           (-d[1355] * d[208]) * dv_1899 +
           (-256.0 * d[1224] * d[206]) * dv_2104 +
           (32.0 * d[1203] * d[1216]) * dv_2100 +
           (d[2928] * d[1306] * d[330]) * dv_1887 +
           (d[1134] * d[336] * d[528]) * dv_1737 +
           (d[1352] * d[236] * d[2927]) * dv_1902 +
           (d[1354] * d[151] * d[549]) * dv_1905 +
           (-d[1039] * d[552] * d[553]) * dv_12;
  DataVector& dv_1897 = temps.at(1630);
  DataVector& dv_950 = temps.at(822);
  DataVector& dv_966 = temps.at(872);
  DataVector& dv_980 = temps.at(806);
  DataVector& sc_10 = temps.at(3232);
  DataVector& sc_12 = temps.at(3234);
  DataVector& sc_13 = temps.at(3235);
  DataVector& sc_14 = temps.at(3236);
  DataVector& sc_18 = temps.at(3240);
  DataVector& sc_6 = temps.at(3228);
  sc_11 += (-d[108] * d[1353] * d[48]) * dv_1897 +
           (-d[1133] * d[1203] * d[335]) * dv_950 +
           (-d[1201] * d[4] * d[60]) * dv_1737 +
           (-1536.0 * d[14] * d[539] * d[60]) * dv_12 +
           (96.0 * d[206] * d[476] * d[569]) * dv_978 +
           (d[1134] * d[1226] * d[1227] * d[375]) * dv_2098 +
           (-192.0 * d[1129] * d[12] * d[1203] * d[427]) * dv_980 +
           (384.0 * d[1204] * d[206] * d[382] * d[4]) * dv_966 + sc_10 + sc_12 +
           sc_13 + sc_14 + sc_18 + sc_6;
  DataVector& dv_3020 = temps.at(2545);
  sc_11 += d[351] * dv_3020;
  DataVector& dv_1880 = temps.at(294);
  DataVector& dv_1885 = temps.at(479);
  DataVector& dv_2113 = temps.at(1730);
  DataVector& dv_2116 = temps.at(1733);
  DataVector& dv_2117 = temps.at(1734);
  DataVector& dv_2119 = temps.at(1736);
  DataVector& dv_2120 = temps.at(1737);
  DataVector& dv_2121 = temps.at(291);
  DataVector& dv_2139 = temps.at(1755);
  sc_11 +=
      d[74] * (d[1210] * ((-d[2922]) * dv_2121 + dv_2120) + dv_2139 +
               d[2920] * ((-d[81]) * dv_2117 - 6.0 * dv_2113 + dv_2116 +
                          d[2924] * (d[1269] * dv_1880 + dv_1885 + dv_2119)));
  DataVector& dv_2097 = temps.at(873);
  DataVector& dv_2101 = temps.at(1720);
  DataVector& dv_2108 = temps.at(870);
  DataVector& dv_2161 = temps.at(1772);
  DataVector& dv_2163 = temps.at(1774);
  DataVector& dv_2165 = temps.at(1776);
  DataVector& sc_16 = temps.at(3238);
  DataVector& sc_19 = temps.at(3241);
  DataVector& sc_22 = temps.at(3244);
  sc_11 += d[85] * dv_2097 + d[85] * dv_2101 + d[85] * dv_2108 +
           d[86] * dv_2097 + d[86] * dv_2101 + d[86] * dv_2108 +
           dv_1509 * dv_2161 + dv_1509 * dv_2163 + dv_1509 * dv_2165 + sc_16 +
           sc_19 + sc_22;
  DataVector& dv_1626 = temps.at(1471);
  DataVector& dv_1906 = temps.at(695);
  DataVector& dv_1907 = temps.at(707);
  DataVector& dv_1908 = temps.at(730);
  DataVector& dv_1913 = temps.at(1634);
  DataVector& dv_2109 = temps.at(1726);
  DataVector& dv_2164 = temps.at(1775);
  sc_11 += (-d[1050]) * dv_1908 * dv_2164 + (-d[377]) * dv_1509 * dv_1906 +
           (-d[392]) * dv_1509 * dv_1907 + (-d[422]) * dv_1509 * dv_1908 +
           (-d[526]) * dv_1626 * dv_1913 - dv_2109 * dv_2110 + sc_20 + sc_21 +
           sc_24 + sc_27;
  DataVector& dv_1510 = temps.at(1389);
  DataVector& dv_1909 = temps.at(698);
  DataVector& dv_1910 = temps.at(1631);
  DataVector& dv_2099 = temps.at(1719);
  DataVector& dv_2162 = temps.at(1773);
  DataVector& dv_976 = temps.at(882);
  sc_11 += (-d[85]) * dv_2109 * dv_6 + (128.0 * d[206]) * dv_1510 * dv_2099 +
           (d[1278] * d[1279]) * dv_1909 * dv_6 +
           (d[1280] * d[390]) * dv_2162 * dv_976 +
           (d[1309] * d[420]) * dv_1910 * dv_2164 +
           (-d[1307] * d[25]) * dv_1907 * dv_2162 +
           (-d[407] * d[92]) * dv_1913 * dv_2112;
  DataVector& dv_1912 = temps.at(1633);
  DataVector& dv_2107 = temps.at(1725);
  DataVector& dv_2111 = temps.at(1728);
  sc_11 += (-128.0 * d[1129] * d[476]) * dv_2106 * dv_2107 +
           (-96.0 * d[1203] * d[527]) * dv_1912 * dv_2111 +
           (-d[130] * d[1308] * d[375]) * dv_1906 * dv_6;
  DataVector& dv_3052 = temps.at(2485);
  DataVector& sc_4 = temps.at(3226);
  sc_4 = dv_3052 * sc_11;
  DataVector& dv_2168 = temps.at(487);
  DataVector& dv_2171 = temps.at(1779);
  DataVector& dv_2185 = temps.at(1791);
  DataVector& dv_2192 = temps.at(1744);
  DataVector& dv_2674 = temps.at(2236);
  DataVector& dv_344 = temps.at(327);
  DataVector& dv_45 = temps.at(45);
  DataVector& dv_529 = temps.at(484);
  sc_22 = d[1210] * (d[1822] * dv_45 + dv_2185) +
          d[6] * (d[139] * (-dv_2171 - dv_558) - dv_2168 + dv_2674 + dv_344 +
                  dv_529) +
          dv_2192;
  DataVector& dv_1564 = temps.at(1437);
  DataVector& dv_2483 = temps.at(2071);
  sc_22 += dv_1564 * ((15.0 * d[2923]) * Dy + (-d[1734] - d[1823]) - dv_2483);
  sc_24 = (-d[19]) * sc_22;
  sc_22 = d[1222];
  DataVector& dv_1478 = temps.at(1359);
  DataVector& dv_2151 = temps.at(1764);
  DataVector& dv_2155 = temps.at(1767);
  DataVector& dv_3095 = temps.at(2588);
  DataVector& dv_3097 = temps.at(2590);
  DataVector& dv_596 = temps.at(547);
  DataVector& dv_603 = temps.at(554);
  DataVector& dv_610 = temps.at(561);
  DataVector& dv_611 = temps.at(562);
  sc_22 *= (-d[2920]) * (dv_1 * dv_2155 - dv_3097 + dv_596 + dv_611 * d[2922]) +
           (-d[2923]) * (-dv_3095 + d[2921] * (dv_603 + dv_610)) +
           d[170] * dv_1478 + d[6] * dv_2151;
  DataVector& dv_1190 = temps.at(1085);
  DataVector& dv_1461 = temps.at(1343);
  DataVector& dv_1965 = temps.at(1681);
  DataVector& dv_2149 = temps.at(1762);
  DataVector& dv_2205 = temps.at(1809);
  DataVector& dv_48 = temps.at(48);
  sc_20 = (-d[138]) * dv_48 + (-d[1494]) * dv_0 + (-d[159]) * dv_1965 +
          (-d[1761]) * dv_0 + (-d[1824]) * dv_1190 + (-d[1825]) * dv_1108 +
          (-d[1826]) * dv_1461 + (-d[196]) * dv_2149 + (-d[367]) * dv_48 +
          dv_2205 + sc_24;
  DataVector& dv_1200 = temps.at(1095);
  DataVector& dv_1545 = temps.at(1419);
  DataVector& dv_2166 = temps.at(481);
  DataVector& dv_2167 = temps.at(485);
  DataVector& dv_2178 = temps.at(1786);
  DataVector& dv_700 = temps.at(650);
  DataVector& dv_993 = temps.at(893);
  sc_20 += (-d[88]) * dv_1545 +
           (-d[2920]) *
               ((-d[2924]) * dv_2166 +
                dv_2115 * (d[558] * ((-d[1540]) - dv_1) + dv_2167) + dv_2178) +
           (d[1030] * d[276]) * dv_993 + (d[1030] * d[559]) * dv_700 +
           (d[1033] * d[431]) * dv_1200;
  DataVector& dv_131 = temps.at(129);
  DataVector& dv_2179 = temps.at(314);
  DataVector& dv_2190 = temps.at(1796);
  DataVector& dv_450 = temps.at(423);
  DataVector& dv_503 = temps.at(464);
  DataVector& dv_756 = temps.at(704);
  sc_20 += (d[559] * d[7]) * dv_1200 + d[1415] * dv_756 + d[1827] * dv_131 +
           d[1827] * dv_503 + d[1828] * dv_1478 + d[367] * dv_2190 +
           d[431] * dv_2190 + d[464] * dv_0 + d[6] * dv_2179 + d[602] * dv_450 +
           sc_22;
  DataVector& dv_1017 = temps.at(917);
  DataVector& dv_1063 = temps.at(962);
  DataVector& dv_139 = temps.at(137);
  DataVector& dv_2193 = temps.at(732);
  sc_20 += d[76] * dv_139 - 30.0 * dv_1 * dv_1017 + dv_1017 * dv_2591 +
           dv_1063 * dv_2193;
  sc_27 = (-d[131]) * sc_20;
  DataVector& dv_14 = temps.at(14);
  DataVector& dv_1655 = temps.at(1499);
  DataVector& dv_2690 = temps.at(2245);
  DataVector& dv_2694 = temps.at(754);
  sc_16 =
      (-d[319]) * (d[1106] * ((11.0 - d[1076]) * dv_14 - dv_1655) + dv_2690) +
      dv_2694;
  DataVector& dv_2684 = temps.at(784);
  DataVector& dv_2685 = temps.at(2241);
  DataVector& dv_2687 = temps.at(2243);
  DataVector& dv_2689 = temps.at(2217);
  DataVector& dv_3202 = temps.at(2689);
  sc_16 += d[1250] *
           (Dx * ((-d[20]) * (d[1962] * dv_2250 + dv_2685) +
                  (-d[259]) * ((213.0 * d[2927] - 2.0) * dv_558 + dv_2689) +
                  d[1364] + d[49] * (dv_2687 - dv_3202)) +
            dv_2684);
  DataVector& dv_1558 = temps.at(1432);
  DataVector& dv_2137 = temps.at(1753);
  DataVector& dv_545 = temps.at(499);
  DataVector& dv_719 = temps.at(668);
  DataVector& dv_831 = temps.at(775);
  sc_16 += d[49] * (d[1806] * dv_719 + d[7] * (d[1192] * dv_545 + dv_831) +
                    dv_1558 + dv_2137);
  DataVector& dv_2285 = temps.at(1879);
  DataVector& dv_2677 = temps.at(60);
  DataVector& dv_2678 = temps.at(766);
  DataVector& dv_2681 = temps.at(2226);
  DataVector& dv_2683 = temps.at(2237);
  sc_16 +=
      d[6] * ((-d[348]) * dv_2678 + (-d[51]) * dv_2677 + (-d[77]) * dv_2683 +
              (9.0 * d[2927]) * dv_2285 + d[48] * dv_2681);
  DataVector& dv_123 = temps.at(121);
  DataVector& dv_2220 = temps.at(1822);
  DataVector& dv_2277 = temps.at(1871);
  DataVector& dv_2693 = temps.at(2248);
  DataVector& dv_2717 = temps.at(2270);
  DataVector& dv_3012 = temps.at(2537);
  DataVector& dv_3061 = temps.at(2554);
  DataVector& dv_725 = temps.at(674);
  DataVector& dv_843 = temps.at(786);
  sc_16 += d[77] * ((-d[1806] * d[2923]) * (dv_3012 + dv_725) +
                    d[1460] * dv_2220 + d[147] * dv_843 + dv_123 * d[2923] +
                    dv_2277 + dv_2693 - dv_2717 + dv_3061);
  sc_24 = d[19] * sc_16;
  DataVector& dv_2638 = temps.at(2214);
  DataVector& dv_2640 = temps.at(1233);
  DataVector& dv_3200 = temps.at(2687);
  sc_19 = (-d[6]) * (-Dy * (Dy * d[1957] + d[303]) + dv_2638) +
          dv_0 * ((-d[2928]) * (41.0 * dv_1503 + 8.0) + (d[1502] - d[1958]) -
                  dv_3200) +
          dv_2640;
  DataVector& dv_815 = temps.at(761);
  DataVector& dv_835 = temps.at(778);
  sc_19 += d[2924] * ((26.0 * d[2928] * d[2922]) * dv_835 - dv_815);
  sc_16 = d[55] * sc_19;
  sc_18 = (-d[86]);
  DataVector& dv_2248 = temps.at(1846);
  DataVector& dv_2653 = temps.at(781);
  DataVector& dv_2655 = temps.at(2222);
  DataVector& dv_2656 = temps.at(2223);
  sc_18 *= Dx * (d[1411] + d[20] * (dv_2655 + dv_3200) +
                 d[259] * ((369.0 * d[2927] - 8.0) * Dy + (-d[37]) * dv_2656 +
                           (d[147] * d[306]) + 184.0 * dv_794) +
                 d[91] * (d[2928] + dv_2248 + dv_3202)) -
           dv_2653;
  DataVector& dv_2661 = temps.at(785);
  DataVector& dv_2665 = temps.at(1745);
  DataVector& dv_669 = temps.at(619);
  sc_6 = (-d[276]) * ((-d[2923]) * (d[1957] * dv_14 - dv_669) + dv_2661) +
         dv_2665 + sc_18;
  DataVector& dv_2175 = temps.at(1783);
  DataVector& dv_2652 = temps.at(613);
  DataVector& dv_3148 = temps.at(2637);
  DataVector& dv_455 = temps.at(428);
  DataVector& dv_723 = temps.at(672);
  DataVector& dv_727 = temps.at(676);
  DataVector& dv_772 = temps.at(720);
  DataVector& dv_830 = temps.at(774);
  sc_6 +=
      d[101] * ((-d[1306]) * dv_772 + d[1189] * (d[1140] * dv_727 + dv_830) +
                d[1274] * dv_455 + dv_2175) +
      d[273] * ((d[1930] * d[2923]) * dv_723 + dv_2652 - 549.0 * dv_3148);
  DataVector& dv_2325 = temps.at(1917);
  DataVector& dv_2647 = temps.at(1425);
  DataVector& dv_2649 = temps.at(1431);
  DataVector& dv_2658 = temps.at(2225);
  DataVector& dv_2660 = temps.at(2213);
  sc_6 +=
      d[35] * ((-d[77]) * dv_2658 + (d[2928] * d[1192]) * dv_2325 + dv_2660) +
      d[57] * ((-d[2925]) * dv_2649 + dv_2647);
  sc_19 = sc_6 * d[2921];
  DataVector& dv_2646 = temps.at(2220);
  DataVector& dv_3201 = temps.at(2688);
  DataVector& dv_716 = temps.at(665);
  DataVector& dv_732 = temps.at(681);
  DataVector& dv_738 = temps.at(687);
  sc_22 = (-d[52]) * (d[1806] * ((61.0 * d[142]) * dv_716 + d[1959] * dv_732 +
                                 dv_3201 + dv_738 * d[2922]) +
                      dv_2646);
  DataVector& dv_1631 = temps.at(1476);
  DataVector& dv_2634 = temps.at(1799);
  DataVector& dv_2636 = temps.at(2212);
  DataVector& dv_2676 = temps.at(2238);
  DataVector& dv_3203 = temps.at(2690);
  DataVector& dv_672 = temps.at(622);
  DataVector& dv_742 = temps.at(690);
  sc_22 += (-d[2920]) * ((d[135] * d[2927]) *
                             (d[1961] * dv_1631 * dv_231 + d[248] * dv_3203 +
                              d[289] * dv_672 + dv_742 * d[2922]) +
                         dv_2676) +
           d[122] * ((-d[2924]) * dv_2636 + (2.0 * d[2922]) * dv_2634) + sc_16 +
           sc_24;
  sc_22 += sc_19;
  sc_20 = (-d[331]) * sc_22;
  DataVector& dv_15 = temps.at(15);
  DataVector& dv_3432 = temps.at(2908);
  DataVector& dv_3435 = temps.at(2911);
  DataVector& dv_3436 = temps.at(2912);
  DataVector& dv_676 = temps.at(626);
  sc_18 = d[1443] * ((-d[147]) * dv_676 + d[1242] * dv_3432 +
                     d[80] * (d[2164] * dv_14 + d[2165] * dv_15) + dv_3435 +
                     dv_3436) +
          d[2102] * dv_3432;
  DataVector& dv_1709 = temps.at(1534);
  DataVector& dv_2261 = temps.at(1857);
  DataVector& dv_2293 = temps.at(1886);
  DataVector& dv_3275 = temps.at(2762);
  DataVector& dv_3433 = temps.at(2909);
  DataVector& dv_3434 = temps.at(2910);
  sc_18 += d[273] * ((-d[1556] + 591.0 * d[7] + 156.0 * d[2927]) * dv_14 +
                     3.0 * Dy * (d[1473] + d[2163] * dv_1709 + 197.0 * dv_794) -
                     dv_3434) +
           d[50] * (dv_2261 * dv_3275 - 442.0 * dv_2293 + dv_3433);
  sc_18 += d[57] * dv_2209;
  sc_6 = (-d[2922]) * sc_18;
  DataVector& dv_2429 = temps.at(2019);
  DataVector& dv_2487 = temps.at(2075);
  DataVector& dv_3358 = temps.at(2841);
  DataVector& dv_3425 = temps.at(2901);
  DataVector& dv_729 = temps.at(678);
  sc_24 = (d[1140] + d[1209]) * dv_3358 +
          (-d[276]) * (dv_3425 - dv_729 * ((-d[2928]) * (dv_2487 + 3.0) +
                                           (-d[1963]) * dv_2429 + (-d[1509]))) +
          sc_6;
  DataVector& dv_3426 = temps.at(2902);
  DataVector& dv_3427 = temps.at(2903);
  DataVector& dv_3428 = temps.at(2904);
  DataVector& dv_3429 = temps.at(2905);
  DataVector& dv_3430 = temps.at(2906);
  sc_24 +=
      d[1443] * (d[2157] * dv_3428 + d[36] * ((-d[2924]) * dv_3429 + dv_3430) +
                 dv_3426 + dv_3427);
  DataVector& dv_352 = temps.at(335);
  DataVector& dv_454 = temps.at(427);
  DataVector& dv_833 = temps.at(776);
  DataVector& dv_972 = temps.at(878);
  sc_24 += d[1534] * ((-d[354]) * (Dy * (d[1519] + dv_972) + dv_352) +
                      Dy * d[471] + d[1] * (dv_454 * d[2923] + dv_833));
  DataVector& dv_1552 = temps.at(1426);
  DataVector& dv_26 = temps.at(26);
  DataVector& dv_3083 = temps.at(2576);
  DataVector& dv_3345 = temps.at(2829);
  DataVector& dv_3362 = temps.at(2845);
  DataVector& dv_3424 = temps.at(2900);
  DataVector& dv_345 = temps.at(328);
  sc_24 += d[724] * ((-d[1447]) * dv_45 + Dx * d[2159] * dv_1552 +
                     d[1] * dv_1112 * (dv_3083 - 4.0) +
                     d[1240] * (dv_26 + dv_345) + dv_3345) +
           d[780] * (dv_3362 + dv_3424);
  DataVector& dv_3431 = temps.at(2907);
  sc_24 += dv_2115 * ((-d[631]) * ((-d[1548]) * dv_1503 +
                                   (d[1099] + 56.0 * d[255]) + 423.0 * dv_1) +
                      d[1443] * ((-d[2160]) * dv_558 + d[1636] + dv_2514) +
                      d[1700] * ((-d[2154]) - 4.0 * dv_3431) + dv_2516);
  sc_16 = (-d[122]) * sc_24;
  DataVector& dv_1583 = temps.at(1455);
  DataVector& dv_3057 = temps.at(2550);
  DataVector& dv_3339 = temps.at(2514);
  DataVector& dv_3439 = temps.at(2915);
  DataVector& dv_3440 = temps.at(2916);
  sc_13 =
      (-d[1443]) * ((-d[2160]) * dv_3339 + d[36] * (-85.0 * dv_1583 - dv_3440) +
                    dv_3426 + dv_3439) +
      (-d[780]) * dv_3057;
  DataVector& dv_1881 = temps.at(275);
  DataVector& dv_2671 = temps.at(2233);
  DataVector& dv_3412 = temps.at(2888);
  sc_13 +=
      d[276] * (dv_3412 + dv_729 * (d[2928] + d[1967] * dv_1881 + dv_2671));
  DataVector& dv_102 = temps.at(100);
  DataVector& dv_3138 = temps.at(2627);
  DataVector& dv_3287 = temps.at(2774);
  DataVector& dv_3441 = temps.at(2917);
  sc_13 +=
      d[724] * ((-d[1]) * dv_3138 + d[1240] * dv_102 + d[167] * (Dx - dv_3441) +
                d[2155] * dv_732 + 70.0 * dv_3287);
  DataVector& dv_2354 = temps.at(1944);
  sc_13 += dv_2354 * (d[2166] * dv_1514 + 20.0 * dv_2379 - dv_2521);
  sc_18 = (-d[2921]) * sc_13;
  DataVector& dv_96 = temps.at(94);
  sc_12 = (-d[1532]) * dv_1 +
          (-d[1553]) * ((-d[1068] + d[99]) * dv_96 +
                        Dy * ((-d[1848]) * dv_558 + d[1473] + 433.0 * dv_794));
  DataVector& dv_2307 = temps.at(1899);
  DataVector& dv_3019 = temps.at(2544);
  DataVector& dv_3416 = temps.at(2892);
  DataVector& dv_3443 = temps.at(2919);
  sc_12 +=
      d[1550] * ((-d[80]) * (d[2164] * dv_15 + d[2165] * dv_14) + Dy * d[2170] +
                 d[1242] * dv_3416 + d[147] * dv_3019 + dv_3443) +
      d[2167] * dv_2307;
  DataVector& dv_3089 = temps.at(2582);
  DataVector& dv_3442 = temps.at(2918);
  sc_12 += d[308] * ((-d[2168] + 252.0 * d[2927] + 16.0) * dv_14 -
                     Dy * ((-d[1473]) * (dv_3089 - 3.0) + (-d[2152]) * dv_1709 +
                           d[2169] + dv_3442));
  DataVector& dv_3414 = temps.at(2890);
  DataVector& dv_3415 = temps.at(2891);
  sc_12 += d[57] * ((-d[1962]) * dv_3414 + dv_2261 * (dv_1502 - 1.0) +
                    118.0 * dv_2293 + dv_3415);
  sc_13 = sc_12 * d[2922];
  DataVector& dv_29 = temps.at(29);
  sc_6 = d[1542] * (Dy * (d[259] * dv_2538 + dv_2540) + d[380] * dv_29) + sc_18;
  DataVector& dv_240 = temps.at(236);
  DataVector& dv_2736 = temps.at(2289);
  DataVector& dv_3437 = temps.at(2913);
  DataVector& dv_3438 = temps.at(2914);
  sc_6 += dv_2115 * ((-d[1549]) * (d[1464] + dv_2736) +
                     d[57] * ((-d[1962]) * dv_1685 + d[1310]) +
                     d[604] * (Dy * d[1565] + d[1683] + d[2157] * dv_240) +
                     dv_3437 - 353.0 * dv_3438) +
          sc_13;
  sc_24 = (-d[337]) * sc_6;
  DataVector& dv_3270 = temps.at(2757);
  DataVector& dv_3457 = temps.at(2932);
  DataVector& dv_3458 = temps.at(2933);
  sc_12 = d[101] * (Dx * ((2.0 * d[2185] * d[2923]) * Dy + (-d[1500]) +
                          d[2928] * (dv_3458 + 42.0) - 168.0 * dv_2246) +
                    d[88] * dv_3270 + dv_3457);
  DataVector& dv_3400 = temps.at(2877);
  DataVector& dv_3455 = temps.at(2930);
  DataVector& dv_3456 = temps.at(2931);
  sc_12 += d[50] * (Dx * ((-d[1539]) + d[2928] * dv_2531 + d[2183] * dv_2250) +
                    d[1240] * dv_3400 + dv_3456) +
           d[57] * dv_3455;
  DataVector& dv_2378 = temps.at(1968);
  DataVector& dv_3247 = temps.at(2734);
  DataVector& dv_3446 = temps.at(2922);
  DataVector& dv_3459 = temps.at(2934);
  DataVector& dv_3460 = temps.at(2935);
  sc_12 += d[631] * ((-d[37]) * (Dx * dv_2533 + dv_3446) + 450.0 * dv_3247 +
                     15.0 * dv_3459 + dv_3460) +
           dv_2378 * (d[1473] + d[2166] * dv_538 + dv_3442);
  sc_18 = (-d[1570]) * sc_12;
  DataVector& dv_163 = temps.at(161);
  DataVector& dv_3454 = temps.at(2929);
  sc_10 =
      (-d[1443]) * ((d[1565] - d[2008] + 65.0) * dv_29 +
                    Dy * ((187.0 * d[36]) + 4.0 * dv_3454 + 84.0 * dv_794)) +
      Dy * d[235] + d[1319] * (433.0 * dv_14 + dv_163);
  DataVector& dv_3009 = temps.at(2534);
  DataVector& dv_527 = temps.at(482);
  DataVector& dv_740 = temps.at(689);
  sc_10 +=
      d[273] * (Dy * d[216] + Dy * d[2184] - 56.0 * dv_3009 + 1169.0 * dv_740) +
      d[51] * (Dy * ((-d[1528]) + d[2183] * dv_1685) + d[2183] * dv_527);
  sc_12 = (-d[35]) * sc_10;
  DataVector& dv_3277 = temps.at(2764);
  DataVector& dv_3301 = temps.at(2788);
  DataVector& dv_3312 = temps.at(2799);
  DataVector& dv_3350 = temps.at(2834);
  DataVector& dv_3450 = temps.at(2925);
  sc_13 = d[1532] * (Dy * (d[1428] + dv_2359) + dv_3350) +
          d[1536] * (d[1099] * dv_2220 + d[1242] * dv_3301 + dv_2511 * dv_3277 -
                     36.0 * dv_3312 + 16.0 * dv_740) +
          dv_3450 + sc_12 + sc_18;
  DataVector& dv_2301 = temps.at(1894);
  DataVector& dv_3311 = temps.at(2798);
  DataVector& dv_3337 = temps.at(2823);
  DataVector& dv_3419 = temps.at(2895);
  DataVector& dv_3452 = temps.at(2927);
  DataVector& dv_3453 = temps.at(2928);
  sc_13 +=
      d[1554] *
          (d[2053] * dv_3301 + dv_3419 + dv_3452 + dv_833 * (6.0 - dv_3453)) +
      d[1594] * ((-d[1474]) * dv_2301 - dv_2505 * dv_833 + dv_3311 + dv_3337);
  DataVector& dv_2888 = temps.at(1619);
  DataVector& dv_3303 = temps.at(2790);
  DataVector& dv_3353 = temps.at(2837);
  DataVector& dv_635 = temps.at(585);
  sc_13 += d[2182] * ((-d[2179]) * dv_635 +
                      d[101] * ((-d[2180]) * dv_14 + d[2181] * dv_15) +
                      d[1439] * dv_3353 + d[313] * dv_740) +
           d[355] * (dv_2888 + dv_3303 * d[2923]);
  DataVector& dv_2654 = temps.at(2221);
  DataVector& dv_3130 = temps.at(2621);
  DataVector& dv_3451 = temps.at(2926);
  sc_13 +=
      dv_3451 * ((-d[431]) * Dy + (-d[20]) * (d[378] + dv_3130) +
                 (4.0 * d[2928] * d[2921]) * (d[2928] + dv_2654) + (-d[2178]));
  sc_6 = (-d[55]) * sc_13;
  DataVector& dv_2664 = temps.at(2228);
  DataVector& dv_2898 = temps.at(2443);
  DataVector& dv_3413 = temps.at(2889);
  DataVector& dv_691 = temps.at(641);
  DataVector& dv_773 = temps.at(721);
  sc_12 = d[1053] * ((-d[631]) * (-Dy * dv_2493 + dv_691) + d[1527] * dv_3416 +
                     d[50] * dv_3413 + 39.0 * dv_2898 + dv_773) +
          d[380] * dv_2664;
  DataVector& dv_2831 = temps.at(2383);
  DataVector& dv_3410 = temps.at(2887);
  DataVector& dv_3417 = temps.at(2893);
  sc_12 += dv_2350 * ((d[1140] + d[1817]) * dv_2354 + (-d[1443]) * dv_3417 +
                      (-d[58]) * dv_3410 +
                      d[724] * (Dy * d[2147] + d[37] - 88.0 * dv_794)) +
           dv_2831 * ((-d[1513]) * dv_1502 +
                      (24.0 * d[2927] * d[2921] * d[2923]) * Dy);
  sc_13 = (-d[57]) * sc_12;
  DataVector& dv_2667 = temps.at(2230);
  sc_10 = dv_2667;
  DataVector& dv_2751 = temps.at(2304);
  DataVector& dv_3133 = temps.at(2624);
  DataVector& dv_3447 = temps.at(2923);
  DataVector& dv_3448 = temps.at(2120);
  sc_10 *= (-d[101]) * ((-d[2172]) * dv_240 + (191.0 * d[36]) + dv_3447) +
           d[151] * (d[306] + 725.0 * dv_1) +
           d[50] * ((-221.0 * d[36]) + d[1963] * dv_3133) +
           d[724] * ((-d[2928]) * dv_3448 + (28.0 * d[255]) + dv_2751) + d[780];
  DataVector& dv_2306 = temps.at(1898);
  DataVector& dv_3240 = temps.at(2727);
  DataVector& dv_3370 = temps.at(2852);
  sc_8 = (-d[2167]) * dv_2306 +
         (-d[308]) * ((d[1531] + d[1956] - d[2176] + 64.0) * dv_14 -
                      Dy * ((-d[1473]) * (dv_3089 - 13.0) +
                            (-d[2177]) * dv_3240 + Dy * d[2176] + d[2169]) +
                      dv_3434) +
         d[120] * dv_3370;
  DataVector& dv_2334 = temps.at(1925);
  DataVector& dv_3449 = temps.at(2924);
  sc_8 +=
      d[1550] * ((111.0 * d[2927] - 16.0) * dv_2334 - 191.0 * dv_2293 +
                 94.0 * dv_3009 - 112.0 * dv_3312 + dv_3449 * (dv_1502 + 7.0));
  sc_8 += d[355] * ((d[1076] + d[2174] - 42.0) * dv_14 +
                    Dy * (Dy * d[2174] + d[2175] * dv_558 + d[314])) +
          d[57] * (dv_2261 * (dv_2536 + 3.0) - 486.0 * dv_2293 + dv_3433);
  sc_14 = sc_8 * d[2922];
  DataVector& dv_2183 = temps.at(1790);
  DataVector& dv_2229 = temps.at(1831);
  DataVector& dv_3262 = temps.at(2749);
  DataVector& dv_3285 = temps.at(2772);
  DataVector& dv_3445 = temps.at(2921);
  DataVector& dv_351 = temps.at(334);
  sc_18 =
      (-d[1532]) * (dv_3262 + dv_3424) +
      (-d[1536]) * ((-d[1967]) * dv_1 * dv_2183 + d[1240] * (dv_351 + dv_96) +
                    d[167] * (-dv_2229 + dv_3446) + dv_3285 - 95.0 * dv_3287) +
      dv_3445;
  DataVector& dv_2230 = temps.at(1832);
  DataVector& dv_3294 = temps.at(2781);
  DataVector& dv_3296 = temps.at(2783);
  sc_18 += (-d[1550]) *
           ((-d[2172]) * dv_3428 +
            (-d[31]) * (47.0 * dv_1583 + dv_2230 * dv_3296 + 47.0 * dv_3270) +
            d[46] * dv_2379 + 204.0 * dv_3294);
  DataVector& dv_2298 = temps.at(1891);
  DataVector& dv_3272 = temps.at(2759);
  sc_18 += d[1135] * (-dv_2298 * ((-d[1]) * dv_3272 + (-d[1966]) * dv_2429 +
                                  (19.0 * d[2928] * d[7])) +
                      dv_3425);
  DataVector& dv_2519 = temps.at(2107);
  DataVector& dv_2603 = temps.at(2185);
  sc_18 += d[1542] * (-Dy * ((-d[248]) * dv_2794 +
                             (d[2928] * d[2921]) * (d[46] + dv_2519) +
                             (-d[2173]) - dv_2603) +
                      d[1533] * dv_691) +
           sc_10 + sc_14;
  DataVector& dv_2343 = temps.at(1933);
  sc_18 += d[313] * dv_231 *
           ((-d[88]) * dv_2343 + d[1616] + d[1999] * dv_2429 + dv_2801);
  sc_12 = d[205] * sc_18;
  DataVector& dv_1083 = temps.at(981);
  DataVector& dv_3423 = temps.at(2899);
  DataVector& dv_646 = temps.at(596);
  DataVector& dv_653 = temps.at(603);
  sc_8 = (30.0 * d[50]) * dv_751 +
         d[116] * (Dx * ((-59.0 * d[255]) + dv_3423) + d[1240] * dv_646 +
                   d[1240] * dv_653) +
         dv_1083 * dv_3417;
  DataVector& dv_2827 = temps.at(2379);
  sc_8 +=
      dv_2827 * ((-d[37]) * (dv_3083 - 3.0) + Dy * d[2155] + 115.0 * dv_794);
  sc_10 = d[1250] * sc_8;
  DataVector& dv_2468 = temps.at(2057);
  DataVector& dv_3422 = temps.at(2898);
  DataVector& dv_589 = temps.at(541);
  sc_2 = (-d[321]) * (dv_3422 + dv_589 * d[2923]) +
         d[259] * ((-d[1242]) * dv_676 + dv_2261 + 353.0 * dv_740) +
         d[50] * dv_2468;
  sc_2 += d[62] * (Dy * (Dy * d[2153] + d[2154]) + d[2153] * dv_14);
  sc_8 = d[6] * sc_2;
  DataVector& dv_1731 = temps.at(1547);
  DataVector& dv_2692 = temps.at(2247);
  DataVector& dv_3351 = temps.at(2835);
  sc_14 = (-d[1333]) * (dv_2301 + dv_3351) +
          (-d[321]) * (d[80] * dv_1731 + dv_2692);
  DataVector& dv_2651 = temps.at(758);
  DataVector& dv_3421 = temps.at(2897);
  sc_14 +=
      (-d[436]) * (d[1242] * dv_1731 + d[1375] * dv_635 + dv_2651 + dv_3421 +
                   d[2923] * ((d[1998] + 4.0) * dv_15 + d[2152] * dv_14)) +
      sc_10;
  DataVector& dv_2746 = temps.at(2299);
  DataVector& dv_3418 = temps.at(2894);
  DataVector& dv_3420 = temps.at(2896);
  sc_14 += d[253] * (d[2151] * dv_3420 + dv_3277 + dv_3418 - dv_3419) +
           d[88] * dv_2349 * ((d[1524] + d[2078]) + dv_2746) + sc_8;
  sc_18 = d[299] * sc_14;
  DataVector& dv_2606 = temps.at(2188);
  DataVector& dv_3411 = temps.at(2821);
  sc_8 =
      (-d[290]) * dv_2606 + (-d[62]) * dv_3411 +
      (d[2921] * d[2923]) *
          (dv_2298 * ((-d[1]) * (dv_1503 - 1.0) + d[1962] * dv_1881 + d[255]) +
           dv_3412);
  DataVector& dv_1039 = temps.at(939);
  DataVector& dv_2891 = temps.at(2437);
  DataVector& dv_2966 = temps.at(2494);
  sc_8 += d[2922] *
          ((-d[1076]) * (d[2149] * dv_14 - dv_1039) + (-d[36]) * dv_3413 +
           (-d[2921]) *
               (Dy * d[1968] - dv_2261 * (dv_2891 + 1.0) + dv_3414 + dv_3415) +
           dv_2966);
  sc_8 +=
      (-d[9]) * Dx * (d[2147] * dv_1 + dv_2486) +
      (2.0 * d[6]) * Dx *
          (d[62] + 176.0 * dv_1229 + d[2921] * ((-d[2148]) + d[1967] * dv_605));
  sc_14 = d[362] * sc_8;
  DataVector& dv_3467 = temps.at(2942);
  sc_28 = (-d[120]) * dv_2438 +
          (-d[1549]) * ((-d[37]) * (Dx * dv_3448 + dv_3441) +
                        (-3.0 * d[2159]) * dv_45 + 423.0 * dv_3247 + dv_3460) +
          d[59] * (Dx * ((81.0 * d[2928] * d[7]) - dv_3423) + dv_3467);
  DataVector& dv_3465 = temps.at(2940);
  sc_28 += d[604] * (Dx * ((-d[2185]) * dv_69 + (187.0 * d[255]) +
                           d[2928] * (6.0 - dv_3458) + 280.0 * dv_2246) +
                     dv_3467) +
           dv_2503 * dv_3465;
  DataVector& dv_1676 = temps.at(1510);
  DataVector& dv_3249 = temps.at(2736);
  sc_28 += -dv_3249 * (d[1473] + d[1999] * dv_1676 + 725.0 * dv_794);
  sc_2 = (-d[1250]) * sc_28;
  DataVector& dv_3464 = temps.at(2939);
  sc_26 = (-d[1409]) * (d[2148] + d[2188] * dv_538) +
          (-d[1443]) *
              ((13.0 - 61.0 * d[2927]) * dv_558 + (177.0 * d[36]) + dv_3464) +
          d[273] * (d[2184] + d[88] * (4.0 - dv_3089) + 591.0 * dv_1) + d[780];
  sc_26 += d[807] * (d[88] + 139.0 * dv_1);
  sc_29 = -dv_5 * sc_26;
  sc_17 = (-591.0 * d[1543] + d[1646] + d[2188] * d[780] +
           d[363] * (d[1029] + d[1545]) - 556.0 * d[992]) *
              dv_14 +
          sc_29;
  sc_28 = (-d[6]) * sc_17;
  DataVector& dv_2900 = temps.at(2445);
  DataVector& dv_3461 = temps.at(2936);
  sc_10 =
      (-d[1431]) * ((-d[278]) * dv_2301 + dv_2900 - dv_3418 + 21.0 * dv_740) +
      (-d[1532]) * (Dy * (d[507] + dv_1502) + dv_2301) - dv_3461 + sc_2;
  DataVector& dv_2963 = temps.at(2491);
  DataVector& dv_3316 = temps.at(2803);
  DataVector& dv_3349 = temps.at(2833);
  sc_10 += (-d[1536]) * (d[1242] * dv_3316 + dv_2561 * dv_3277 + dv_2963 -
                         60.0 * dv_3312 + dv_3349) +
           (-d[1554]) * ((-d[88]) * dv_2301 + d[1755] * dv_3316 + dv_2651 -
                         dv_833 * (dv_3453 + 42.0)) +
           sc_28;
  DataVector& dv_1762 = temps.at(1565);
  DataVector& dv_2402 = temps.at(1992);
  DataVector& dv_2880 = temps.at(2429);
  DataVector& dv_3341 = temps.at(2825);
  DataVector& dv_3366 = temps.at(2849);
  DataVector& dv_3462 = temps.at(2937);
  DataVector& dv_3463 = temps.at(2938);
  DataVector& dv_609 = temps.at(560);
  sc_10 += (8.0 * d[151] * d[2921]) * (d[1106] * (dv_14 + dv_609) -
                                       dv_1762 * dv_2402 + dv_2880 + dv_3462) +
           (12.0 * d[2927] * d[2921]) *
               (d[101] * (d[2186] * dv_14 + d[2187] * dv_15) +
                d[274] * dv_3366 + dv_3341 + dv_3463);
  sc_10 +=
      (2.0 * d[2928] * d[142] * d[2921]) * Dx *
      ((-d[77]) * (d[1540] + dv_2209) + Dy * d[244] + d[1409] + 48.0 * dv_613);
  sc_8 = d[60] * sc_10;
  DataVector& dv_1115 = temps.at(1012);
  sc_19 = d[1830] * ((-d[1507]) * dv_1115 + d[1507] * dv_2220 +
                     dv_0 * (d[1376] + dv_3410)) +
          sc_12 + sc_13 + sc_14 + sc_16 + sc_18 + sc_24 + sc_6 + sc_8;
  sc_22 = (-d[392]) * sc_19;
  sc_12 = (-d[2922]);
  DataVector& dv_3271 = temps.at(2758);
  DataVector& dv_3521 = temps.at(2995);
  DataVector& dv_3601 = temps.at(3073);
  DataVector& dv_3602 = temps.at(3074);
  DataVector& dv_637 = temps.at(587);
  sc_12 *= d[2321] * dv_3521 +
           d[348] * (Dx * ((-d[1879]) - dv_3602) + 136.0 * dv_3271 + dv_3457) +
           d[50] * dv_3601 +
           dv_2827 * ((-d[2328]) * dv_637 + d[2327] + 171.0 * dv_794);
  DataVector& dv_3600 = temps.at(466);
  sc_18 = (-d[1628]) * (d[1247] * dv_3600 + d[2323] * dv_3420 + dv_2880 +
                        dv_3449 * (11.0 * dv_1502 + 2.0));
  DataVector& dv_101 = temps.at(99);
  DataVector& dv_59 = temps.at(59);
  sc_18 +=
      (-d[35]) *
          ((-d[29]) * (Dy * (d[2133] + d[2326] * dv_538) + d[2326] * dv_96) +
           Dy * d[1838] +
           d[1] * (d[1106] * (d[2325] * dv_14 + dv_101) + dv_59)) +
      sc_12;
  DataVector& dv_2302 = temps.at(1895);
  DataVector& dv_3330 = temps.at(2816);
  DataVector& dv_3568 = temps.at(3041);
  sc_18 += d[1251] * (-dv_2302 - dv_240 * dv_3330) + d[2066] * dv_3568;
  sc_18 += d[436] * ((-d[1344]) * dv_3600 +
                     d[1230] * ((-d[2324]) * dv_14 + d[2175] * dv_15) +
                     d[1242] * dv_3568 + dv_3421 + dv_3452) +
           dv_2554 * ((34.0 * d[2921]) - dv_637);
  sc_14 = (-d[1466]) * sc_18;
  DataVector& dv_3308 = temps.at(2795);
  DataVector& dv_3626 = temps.at(3097);
  sc_6 = d[57] * dv_3601 + d[58] * (Dx * ((-232.0 * d[255]) - dv_3626) +
                                    d[1617] * dv_3270 + 70.0 * dv_3308);
  DataVector& dv_2406 = temps.at(1996);
  DataVector& dv_2444 = temps.at(2034);
  DataVector& dv_2850 = temps.at(2402);
  DataVector& dv_3625 = temps.at(3096);
  sc_6 +=
      -dv_2444 * ((-d[2347]) * dv_1685 + (-d[2348]) * dv_794 + d[2123]) +
      dv_2771 * ((-d[1655]) + d[2928] * (dv_2406 + dv_2791 + 25.0) - dv_3625) +
      dv_2788 * dv_2850;
  sc_13 = (-d[1570]) * sc_6;
  DataVector& dv_2868 = temps.at(2418);
  DataVector& dv_294 = temps.at(286);
  sc_24 = (-4968.0 * d[151]) * dv_2868 +
          d[1700] * (Dy * ((116.0 * d[36]) + Dy * d[2361]) + d[2361] * dv_14) -
          162.0 * dv_294;
  sc_24 +=
      d[287] * ((432.0 * d[2927] - 385.0) * dv_14 +
                dv_240 * ((d[1192] + 7.0) * Dy + (174.0 * d[36]))) +
      d[724] * ((-d[1106]) * ((d[1280] + 133.0) * dv_15 + d[2349] * dv_96) +
                107.0 * dv_833);
  sc_6 = d[35] * sc_24;
  DataVector& dv_3315 = temps.at(2802);
  DataVector& dv_3605 = temps.at(3077);
  DataVector& dv_3607 = temps.at(3078);
  DataVector& dv_3608 = temps.at(3079);
  DataVector& dv_3609 = temps.at(3080);
  DataVector& dv_3613 = temps.at(3084);
  DataVector& dv_3614 = temps.at(3085);
  DataVector& dv_3615 = temps.at(3086);
  DataVector& dv_3616 = temps.at(3087);
  DataVector& dv_3617 = temps.at(3088);
  DataVector& dv_3618 = temps.at(3089);
  DataVector& dv_3622 = temps.at(3093);
  sc_12 = -428.0 * dv_3315 + dv_3450 - dv_3605 + dv_3607 - dv_3608 + dv_3609 -
          dv_3613 - 768.0 * dv_3614 - 384.0 * dv_3615 - 1536.0 * dv_3616 -
          536.0 * dv_3617 - 768.0 * dv_3618 + dv_3622;
  DataVector& dv_3533 = temps.at(3007);
  DataVector& dv_3534 = temps.at(394);
  DataVector& dv_3586 = temps.at(3059);
  DataVector& dv_3610 = temps.at(3081);
  DataVector& dv_3611 = temps.at(3082);
  DataVector& dv_3619 = temps.at(3090);
  DataVector& dv_3623 = temps.at(3094);
  DataVector& dv_3624 = temps.at(3095);
  sc_12 +=
      (-d[2360]) * (d[101] * dv_3533 + d[274] * dv_3534 - dv_3463 + dv_3586) +
      400.0 * dv_3610 + 384.0 * dv_3611 + 268.0 * dv_3619 + 120.0 * dv_3623 +
      804.0 * dv_3624 + sc_13;
  DataVector& dv_2612 = temps.at(2194);
  DataVector& dv_3001 = temps.at(2526);
  DataVector& dv_46 = temps.at(46);
  DataVector& dv_629 = temps.at(580);
  DataVector& dv_703 = temps.at(653);
  sc_12 += d[2341] * dv_2220 + d[2352] * dv_14 + d[2352] * dv_15 +
           d[2353] * dv_703 + d[2354] * dv_2612 + d[2354] * dv_46 +
           d[2357] * dv_629 + d[2358] * dv_3001 + d[2359] * dv_3001 + sc_6;
  DataVector& dv_2759 = temps.at(2312);
  sc_12 += d[102] * dv_2759 * ((d[221] + d[456]) - 483.0 * dv_5);
  sc_18 = (-d[299]) * sc_12;
  sc_13 = (-d[6]);
  sc_13 *= (92.0 * d[1543] + d[1642] + d[230] * d[2323] - d[466] +
            d[996] * (d[1070] - 67.0)) *
               dv_96 +
           Dy * ((-d[1672]) * (d[1889] + dv_3454) +
                 d[1536] * ((d[1280] + 199.0) * dv_1552 + (131.0 * d[2928])) +
                 d[1647] * (d[1] + dv_2562) +
                 d[230] * ((-190.0 * d[36]) + d[2323] * dv_538) + dv_2810);
  DataVector& dv_99 = temps.at(97);
  sc_16 = (-d[230]) *
          (Dx * (d[1968] + dv_3602) + d[1240] * dv_99 - 124.0 * dv_3308);
  DataVector& dv_1806 = temps.at(1592);
  DataVector& dv_2172 = temps.at(1780);
  sc_16 += (-d[371]) * (Dx * d[1648] + Dx * dv_3625 +
                        dv_2172 * (-29.0 * dv_1803 - dv_1806)) +
           (224.0 * d[384]) * dv_732 + d[1649] * dv_751;
  sc_16 += dv_2785 * (d[2335] * dv_538 + d[31] + 540.0 * dv_794) +
           dv_2811 * (d[1572] + d[2329] * dv_794 + d[2333] * dv_1685);
  sc_24 = (-d[2922]) * sc_16;
  DataVector& dv_3606 = temps.at(190);
  DataVector& dv_3612 = temps.at(3083);
  DataVector& dv_3627 = temps.at(3091);
  DataVector& dv_3629 = temps.at(3099);
  sc_6 = 27.0 * dv_3606 + 408.0 * dv_3611 - 198.0 * dv_3612 - 408.0 * dv_3615 -
         1632.0 * dv_3618 - 220.0 * dv_3619 + dv_3622 + 56.0 * dv_3623 -
         1540.0 * dv_3624 - dv_3627 - 96.0 * dv_3629;
  DataVector& dv_2975 = temps.at(2502);
  DataVector& dv_3628 = temps.at(3098);
  sc_6 += (-d[235]) * dv_1229 + (-d[2353]) * dv_2975 + (-d[2356]) * dv_2306 +
          (-d[2364]) * dv_2612 + (-288.0 * d[2358]) * dv_15 +
          (96.0 * d[2365]) * dv_2246 + 524.0 * dv_3315 + 440.0 * dv_3617 +
          1200.0 * dv_3628 + sc_13 + sc_24;
  DataVector& dv_1578 = temps.at(1450);
  DataVector& dv_2006 = temps.at(1514);
  DataVector& dv_2212 = temps.at(1815);
  DataVector& dv_3192 = temps.at(2679);
  sc_6 += (128.0 * d[2288]) * dv_5 + (-d[120] * d[7]) * dv_3192 +
          (-d[2134] * d[36]) * dv_450 + d[1551] * dv_2212 + d[1645] * dv_2220 +
          d[220] * dv_1578 + d[2351] * dv_2006;
  DataVector& dv_3317 = temps.at(2804);
  DataVector& dv_3531 = temps.at(3005);
  DataVector& dv_3554 = temps.at(3027);
  sc_6 += d[2360] * (d[101] * (d[2367] * dv_15 + dv_3317) + d[2366] * dv_635 +
                     d[274] * dv_3554 - dv_3531) +
          d[2362] * dv_14 + d[2362] * dv_15 + d[2363] * dv_2612 +
          d[2363] * dv_46;
  DataVector& dv_3630 = temps.at(3100);
  sc_6 += dv_3630 * ((-d[20]) * dv_2770 + (24.0 * d[48]) * Dy + (-d[51]));
  sc_12 = (-d[341]) * sc_6;
  DataVector& dv_3079 = temps.at(2572);
  sc_16 = (-d[120]) * dv_2775 + (-96.0 * d[384]) * dv_740 +
          d[2293] * ((-d[1140] * d[1591] - 133.0 * d[7] + 18.0) * dv_14 +
                     Dy * ((3.0 - d[2255]) * dv_558 + (-d[2349]) * dv_3079 +
                           (29.0 * d[36])));
  DataVector& dv_121 = temps.at(119);
  DataVector& dv_2352 = temps.at(1942);
  DataVector& dv_3603 = temps.at(3075);
  DataVector& dv_418 = temps.at(393);
  sc_16 += d[230] * ((-d[1193]) * (dv_121 + dv_418) + 169.0 * dv_2293 -
                     dv_3449 * (dv_2352 + 2.0) + dv_3603);
  sc_16 +=
      d[355] * ((-d[1192] - 459.0 * d[7] + 100.0) * dv_14 +
                dv_1709 * ((25.0 - d[1192]) * Dy + d[314] - 237.0 * dv_794));
  DataVector& dv_2861 = temps.at(2413);
  DataVector& dv_724 = temps.at(673);
  sc_16 += d[996] *
           ((-d[1291]) * (dv_2975 + dv_724) + Dy * d[2350] + dv_1762 * dv_2861 +
            d[2923] * (d[2344] * dv_15 + d[2345] * dv_14));
  sc_13 = d[1250] * sc_16;
  DataVector& dv_2809 = temps.at(2362);
  sc_10 = (-d[1672]) * ((-d[2301]) * Dy + d[378]) +
          (-d[457]) * ((33.0 * d[36]) + d[2239] * dv_240) + d[1645] +
          d[1669] * (d[1] + 753.0 * dv_1) + 256.0 * dv_2809;
  sc_10 += d[308] * ((268.0 * d[2928]) + d[2348] * dv_1552);
  sc_16 = -dv_2115 * sc_10;
  DataVector& dv_3004 = temps.at(2529);
  DataVector& dv_3304 = temps.at(2791);
  DataVector& dv_3604 = temps.at(3076);
  sc_24 = (-d[1640]) * ((-d[1670]) * dv_703 + Dy * ((d[1107] - 63.0 * d[50]) +
                                                    dv_3004 + 256.0 * dv_613)) +
          (-d[1645]) * (-dv_3262 - dv_3604) - 512.0 * dv_3304;
  DataVector& dv_2170 = temps.at(1778);
  DataVector& dv_3540 = temps.at(3013);
  DataVector& dv_3578 = temps.at(3051);
  DataVector& dv_774 = temps.at(722);
  sc_24 +=
      d[1172] *
          ((-d[36]) * (192.0 * dv_1583 + dv_2170 * dv_2806 + dv_774 * d[2924]) +
           (3.0 * d[2340] * d[7]) * Dx * Dy - dv_3540 - dv_3578) +
      sc_13;
  DataVector& dv_2920 = temps.at(752);
  DataVector& dv_3260 = temps.at(2747);
  DataVector& dv_3347 = temps.at(2831);
  DataVector& dv_644 = temps.at(594);
  DataVector& dv_683 = temps.at(633);
  sc_24 += d[1549] * ((-d[2347]) * Dx * dv_2250 +
                      (2.0 * d[2928] * d[2924]) * (dv_2920 + dv_683) +
                      (214.0 * d[2928] * d[7]) * Dx - 214.0 * dv_3260 -
                      1449.0 * dv_3287) +
           d[1647] * ((-d[2346]) * dv_732 + d[1240] * dv_644 + dv_3347);
  DataVector& dv_3509 = temps.at(2983);
  DataVector& dv_3574 = temps.at(3047);
  sc_24 +=
      d[1665] *
          ((-d[1787]) * dv_1583 +
           3.0 * Dx *
               ((-d[1]) * (dv_3574 + 2.0) + (23.0 * d[255]) + d[2330] * dv_69) -
           dv_3509) +
      sc_16;
  sc_6 = (-d[343]) * sc_24;
  DataVector& dv_3241 = temps.at(2728);
  DataVector& dv_3283 = temps.at(2770);
  DataVector& dv_3348 = temps.at(2832);
  DataVector& dv_3598 = temps.at(3071);
  DataVector& dv_3599 = temps.at(3072);
  DataVector& dv_565 = temps.at(518);
  sc_10 =
      (-d[1672]) * ((-d[1242]) * dv_565 + d[255] * dv_3348 +
                    d[2923] * ((-d[2160]) * dv_14 + d[2331] * dv_15)) +
      (-d[223]) * (d[1242] * dv_3598 + dv_3241 + dv_3283 * dv_833 - dv_3599);
  DataVector& dv_2991 = temps.at(2518);
  sc_10 += d[1549] * ((-d[1076] - d[2334]) * dv_96 +
                      Dy * (d[2324] * dv_1676 + d[2325] * dv_3079 + d[2327])) +
           d[1551] * dv_2991 + d[1645] * dv_1;
  sc_10 += d[391] * dv_5 * ((d[1306] - 25.0) * dv_240 + d[37] + 621.0 * dv_794);
  sc_13 = (-d[2922]) * sc_10;
  DataVector& dv_2481 = temps.at(2070);
  DataVector& dv_3291 = temps.at(2778);
  DataVector& dv_343 = temps.at(326);
  DataVector& dv_3529 = temps.at(3003);
  sc_28 = d[1710] * ((-d[1240]) * dv_343 +
                     Dx * ((-d[2928]) * dv_3291 + d[167] + d[1832] * dv_2481) +
                     dv_3271) +
          d[384] * dv_3529;
  DataVector& dv_1542 = temps.at(1416);
  DataVector& dv_3183 = temps.at(2670);
  DataVector& dv_3597 = temps.at(3070);
  sc_28 += d[631] * (Dx * d[2328] * dv_1683 + d[2101] * dv_3597 +
                     d[9] * dv_3138 + dv_3183 - 54.0 * dv_3287) +
           dv_1542 * dv_2777;
  DataVector& dv_2018 = temps.at(139);
  sc_28 += d[104] * dv_5 * ((-d[2923]) * (-dv_2018 - dv_2779) + Dx * d[1842]);
  sc_10 = sc_28 * d[2921];
  DataVector& dv_2944 = temps.at(2474);
  DataVector& dv_3314 = temps.at(2801);
  sc_16 = (6.0 * d[2928] * d[142] * d[2921]) *
              (d[116] * (dv_14 + dv_3314) + dv_2944) +
          sc_10 + sc_13;
  sc_16 += -dv_2115 *
           ((-d[91] * (131.0 - d[1806])) * dv_613 + d[308] * (d[9] + dv_2795) +
            d[457] * (d[36] + dv_3431) - dv_3437 + 1392.0 * dv_3438);
}
}  // namespace CurvedScalarWave::Worldtube::detail
