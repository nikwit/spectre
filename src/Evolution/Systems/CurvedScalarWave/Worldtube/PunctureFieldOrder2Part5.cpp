
// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/Worldtube/PunctureFieldOrder2Impl.hpp"

namespace CurvedScalarWave::Worldtube::detail {

// NOLINTNEXTLINE(google-readability-function-size, readability-function-size)
void puncture_field_2_part_5(const std::array<double, order2_n_doubles>& d,
                             const DataVector& Dx, const DataVector& Dy,
                             DynamicBuffer<DataVector>& temps,
                             const gsl::not_null<Order2Vars*> result) {
  DataVector& dv_2611 = temps.at(2193);
  DataVector& dv_4026 = temps.at(1344);
  DataVector& dv_4030 = temps.at(1650);
  DataVector& dv_4038 = temps.at(1645);
  DataVector& dv_4045 = temps.at(996);
  DataVector& dv_4046 = temps.at(3204);
  DataVector& dv_4060 = temps.at(975);
  DataVector& dv_4061 = temps.at(1060);
  DataVector& dv_4062 = temps.at(1084);
  dv_4062 += (-d[53] * d[599]) * dv_2611 + (-d[53] * d[770]) * dv_4026 +
             (-d[627] * d[2920]) * dv_2611 + (-d[767] * d[770]) * dv_4030 +
             (d[1181] * d[186] * d[319]) * dv_4046 +
             (d[119] * d[36] * d[628]) * dv_4061 +
             (d[1233] * d[2] * d[626]) * dv_4046 +
             (d[13] * d[2759] * d[2921]) * dv_4060 +
             (d[19] * d[2369] * d[6]) * dv_4038 +
             (d[2424] * d[319] * d[770]) * dv_4045;
  DataVector& dv_115 = temps.at(113);
  DataVector& dv_1398 = temps.at(1287);
  DataVector& dv_1426 = temps.at(1310);
  DataVector& dv_1968 = temps.at(1683);
  DataVector& dv_3670 = temps.at(3140);
  DataVector& dv_4016 = temps.at(1658);
  DataVector& dv_4023 = temps.at(1656);
  DataVector& dv_4035 = temps.at(1356);
  DataVector& dv_4041 = temps.at(1307);
  DataVector& dv_534 = temps.at(489);
  dv_4062 += (d[2715] * d[2815] * d[2816]) * dv_1426 +
             (d[2719] * d[2794] * d[319]) * dv_1968 +
             (d[273] * d[2731] * d[643]) * dv_4016 +
             (d[2797] * d[2798] * d[619]) * dv_1398 +
             (d[2806] * d[290] * d[420]) * dv_4041 +
             (d[52] * d[615] * d[83]) * dv_534 +
             (d[595] * d[616] * d[2920]) * dv_3670 +
             (-d[1018] * d[682] * d[760]) * dv_4023 +
             (-d[168] * d[2761] * d[2778]) * dv_115 +
             (-d[1921] * d[2807] * d[768]) * dv_4035;
  DataVector& dv_1084 = temps.at(982);
  DataVector& dv_1404 = temps.at(1293);
  DataVector& dv_4031 = temps.at(1271);
  DataVector& dv_4036 = temps.at(1277);
  DataVector& dv_4056 = temps.at(904);
  DataVector& dv_4059 = temps.at(1174);
  DataVector& dv_700 = temps.at(650);
  dv_4062 += (-d[2043] * d[2781] * d[643]) * dv_1084 +
             (-d[2741] * d[621] * d[648]) * dv_700 +
             (-d[2805] * d[598] * d[609]) * dv_4036 +
             (-d[370] * d[676] * d[725]) * dv_4023 +
             (-d[54] * d[598] * d[771]) * dv_4061 +
             (d[2798] * d[682] * d[764] * d[992]) * dv_1404 +
             (-d[231] * d[2761] * d[2926] * d[2922]) * dv_4023 +
             d[1155] * dv_4056 + d[1155] * dv_4059 + d[2446] * dv_4031;
  DataVector& dv_1085 = temps.at(983);
  DataVector& dv_1091 = temps.at(989);
  DataVector& dv_3797 = temps.at(3202);
  DataVector& dv_4051 = temps.at(1661);
  DataVector& dv_4052 = temps.at(3216);
  DataVector& dv_4053 = temps.at(1105);
  DataVector& dv_4054 = temps.at(928);
  DataVector& dv_4055 = temps.at(1180);
  dv_4062 += d[2464] * dv_3797 + d[2589] * dv_1091 + d[2760] * dv_4052 +
             d[2764] * dv_4053 + d[2764] * dv_4055 + d[2768] * dv_4053 +
             d[2768] * dv_4055 + d[277] * dv_1085 + d[2775] * dv_4051 +
             d[2776] * dv_4054;
  DataVector& dv_1048 = temps.at(947);
  DataVector& dv_1090 = temps.at(988);
  DataVector& dv_1166 = temps.at(1061);
  DataVector& dv_3645 = temps.at(3115);
  DataVector& dv_3663 = temps.at(3133);
  DataVector& dv_3984 = temps.at(398);
  DataVector& dv_4018 = temps.at(1643);
  DataVector& dv_4057 = temps.at(1010);
  DataVector& dv_999 = temps.at(899);
  dv_4062 += d[2778] * dv_1048 + d[2779] * dv_4018 + d[2780] * dv_1090 +
             d[2780] * dv_3663 + d[2786] * dv_1166 + d[2786] * dv_3984 +
             d[2807] * dv_4057 + d[52] * dv_999 + d[53] * dv_3645 +
             d[53] * dv_999;
  DataVector& dv_1018 = temps.at(918);
  DataVector& dv_1106 = temps.at(1003);
  DataVector& dv_3731 = temps.at(3195);
  DataVector& dv_3983 = temps.at(964);
  DataVector& dv_3989 = temps.at(443);
  DataVector& dv_3996 = temps.at(306);
  DataVector& dv_4024 = temps.at(1339);
  DataVector& dv_4025 = temps.at(1653);
  DataVector& dv_4028 = temps.at(1647);
  dv_4062 += d[645] * dv_3731 + d[649] * dv_3989 + d[653] * dv_1106 +
             d[653] * dv_3996 + d[677] * dv_4024 + d[727] * dv_4025 +
             d[747] * dv_3983 + d[753] * dv_4028 + d[760] * dv_1018 +
             d[769] * dv_4030;
  DataVector& dv_1020 = temps.at(920);
  DataVector& dv_1040 = temps.at(940);
  DataVector& dv_1057 = temps.at(956);
  DataVector& dv_1120 = temps.at(1017);
  DataVector& dv_1122 = temps.at(1019);
  DataVector& dv_3472 = temps.at(2947);
  DataVector& dv_3673 = temps.at(3143);
  DataVector& dv_3987 = temps.at(371);
  DataVector& dv_4050 = temps.at(1053);
  DataVector& dv_5 = temps.at(5);
  dv_4062 += (-d[580]) * dv_3472 * dv_5 + d[847] * dv_4050 + dv_1020 * d[2920] +
             dv_1040 * d[2922] + dv_1057 * d[2920] + dv_1120 * d[2922] +
             dv_1122 * d[2920] + dv_3673 * d[2920] + dv_3987 * d[2920];
  DataVector& dv_1631 = temps.at(1476);
  DataVector& dv_1771 = temps.at(1573);
  dv_4062 += d[580] * dv_1631 * dv_1771;
  DataVector& dv_4 = temps.at(4);
  DataVector& dv_4063 = temps.at(956);
  DataVector& dv_613 = temps.at(564);
  dv_4063 = (-d[121]) * Dy + d[29] * dv_4 + dv_613;
  DataVector& dv_4064 = temps.at(940);
  DataVector& dv_854 = temps.at(793);
  dv_4064 = d[27] * dv_4 + dv_613 - dv_854;
  DataVector& dv_4065 = temps.at(909);
  dv_4065 = d[21] * dv_4064 + d[4] * dv_4063;
  DataVector& dv_3777 = temps.at(3124);
  DataVector& dv_4066 = temps.at(3124);
  DataVector& dv_871 = temps.at(807);
  dv_4066 = -dv_3777 * dv_4065 + dv_4063 * dv_871;
  DataVector& dv_0 = temps.at(0);
  DataVector& dv_4067 = temps.at(1650);
  dv_4067 = d[1061] * dv_0;
  DataVector& dv_4068 = temps.at(904);
  dv_4068 = d[1] * dv_0;
  DataVector& dv_4069 = temps.at(975);
  dv_4069 = d[78] * dv_4068;
  DataVector& dv_1115 = temps.at(1012);
  DataVector& dv_4070 = temps.at(1683);
  dv_4070 = d[631] * dv_1115;
  DataVector& dv_1531 = temps.at(1408);
  DataVector& dv_263 = temps.at(256);
  DataVector& dv_4071 = temps.at(1271);
  dv_4071 = d[1216] * (dv_1531 + dv_263 - dv_4068);
  DataVector& dv_4072 = temps.at(953);
  dv_4072 = (d[1] * d[2819] + d[282] * d[2926]) * dv_854;
  DataVector& dv_3912 = temps.at(1063);
  DataVector& dv_4073 = temps.at(1652);
  dv_4073 = d[114] * dv_3912;
  DataVector& dv_231 = temps.at(228);
  DataVector& dv_4074 = temps.at(983);
  dv_4074 = dv_4067 + dv_4069 - dv_4070 - dv_4071 - dv_4072 +
            d[2920] * (d[586] * dv_231 + dv_4073);
  DataVector& dv_1587 = temps.at(277);
  DataVector& dv_4075 = temps.at(277);
  dv_4075 = dv_1587 * d[2921];
  DataVector& dv_1588 = temps.at(50);
  DataVector& dv_4076 = temps.at(918);
  dv_4076 = d[1063] * dv_1588;
  DataVector& dv_3794 = temps.at(993);
  DataVector& dv_4077 = temps.at(993);
  dv_4077 = dv_3794 * d[2921];
  DataVector& dv_174 = temps.at(172);
  DataVector& dv_4078 = temps.at(942);
  dv_4078 = d[4] * dv_174;
  DataVector& dv_4079 = temps.at(1180);
  dv_4079 = d[395] * dv_4078;
  DataVector& dv_240 = temps.at(236);
  DataVector& dv_4080 = temps.at(3142);
  dv_4080 = dv_240 * dv_4;
  DataVector& dv_1829 = temps.at(1614);
  DataVector& dv_4081 = temps.at(1614);
  dv_4081 = dv_1829 * d[2921] + dv_4080;
  DataVector& dv_4082 = temps.at(947);
  dv_4082 = d[2426] * dv_4081;
  DataVector& dv_4083 = temps.at(1614);
  dv_4083 = d[873] * dv_4081;
  DataVector& dv_110 = temps.at(108);
  DataVector& dv_3921 = temps.at(3105);
  DataVector& dv_4084 = temps.at(108);
  dv_4084 = dv_110 * dv_3921;
  DataVector& dv_4085 = temps.at(1647);
  DataVector& dv_582 = temps.at(534);
  dv_4085 = Dx * dv_582;
  DataVector& dv_3804 = temps.at(336);
  DataVector& dv_4086 = temps.at(1060);
  dv_4086 = d[205] * (dv_3804 * d[2922] + dv_4085);
  DataVector& dv_1818 = temps.at(1604);
  DataVector& dv_2031 = temps.at(1663);
  DataVector& dv_2965 = temps.at(2493);
  DataVector& dv_4087 = temps.at(1211);
  dv_4087 = (-d[2921]) * (dv_1818 + dv_2031 * d[2923]) + dv_2965;
  DataVector& dv_3340 = temps.at(2824);
  DataVector& dv_4088 = temps.at(3216);
  dv_4088 = d[2820] * dv_3340;
  DataVector& dv_2388 = temps.at(1978);
  DataVector& dv_4089 = temps.at(371);
  dv_4089 = dv_0 * dv_2388;
  DataVector& dv_25 = temps.at(25);
  DataVector& dv_403 = temps.at(383);
  DataVector& dv_4090 = temps.at(383);
  DataVector& dv_594 = temps.at(545);
  dv_4090 = dv_25 + dv_403 + dv_594;
  DataVector& dv_1706 = temps.at(133);
  DataVector& dv_4091 = temps.at(133);
  dv_4091 = -dv_1706;
  DataVector& dv_1679 = temps.at(1513);
  DataVector& dv_1688 = temps.at(1520);
  DataVector& dv_1713 = temps.at(138);
  DataVector& dv_1717 = temps.at(1538);
  DataVector& dv_801 = temps.at(749);
  DataVector& sc_5 = temps.at(3227);
  sc_5 =
      d[2571] * ((-d[1622] - d[9]) * dv_801 + (d[2922] * d[2921]) * dv_1688) +
      d[2590] * (dv_1679 * d[2922] + dv_4085) +
      d[51] * (-dv_2965 + d[2921] * (dv_1713 + dv_1717));
  DataVector& dv_1981 = temps.at(149);
  DataVector& dv_3373 = temps.at(2855);
  DataVector& dv_527 = temps.at(482);
  sc_5 += d[55] * (dv_3373 + dv_4091 * d[2923]) +
          d[63] * ((-d[2928]) * dv_527 + (-d[2921]) * dv_1981);
  DataVector& dv_4092 = temps.at(1287);
  dv_4092 = d[12] * sc_5;
  DataVector& dv_3792 = temps.at(1485);
  DataVector& dv_4093 = temps.at(1485);
  dv_4093 = dv_3792 * dv_4065;
  DataVector& dv_138 = temps.at(136);
  DataVector& dv_14 = temps.at(14);
  DataVector& dv_153 = temps.at(151);
  DataVector& dv_1641 = temps.at(1486);
  DataVector& dv_1651 = temps.at(1496);
  DataVector& dv_210 = temps.at(207);
  DataVector& dv_294 = temps.at(286);
  DataVector& dv_399 = temps.at(380);
  DataVector& dv_556 = temps.at(509);
  DataVector& dv_721 = temps.at(670);
  DataVector& dv_971 = temps.at(877);
  DataVector& sc_7 = temps.at(3229);
  sc_7 = dv_138 * dv_1641 + dv_153 * (68.0 * dv_14 - dv_210 - dv_721) -
         dv_1651 * dv_971 + dv_294 * (-dv_399 - dv_556);
  DataVector& dv_1645 = temps.at(1490);
  DataVector& dv_853 = temps.at(792);
  sc_7 += 15.0 * dv_1645 * dv_853;
  sc_5 = d[1066] * sc_7;
  DataVector& dv_160 = temps.at(158);
  DataVector& dv_171 = temps.at(169);
  DataVector& dv_4094 = temps.at(1647);
  dv_4094 = dv_160 - dv_171 + sc_5;
  DataVector& dv_4095 = temps.at(1271);
  dv_4095 = -dv_4067 - dv_4069 + dv_4070 + dv_4071 + dv_4072 +
            d[2920] * (d[2731] * dv_231 - dv_4073);
  DataVector& dv_10 = temps.at(10);
  DataVector& dv_1781 = temps.at(1435);
  DataVector& dv_4096 = temps.at(1683);
  dv_4096 = (-d[16]) * dv_4075 + d[1063] * dv_1588 + 2.0 * dv_10 * dv_4065 -
            dv_1781 * dv_4063;
  DataVector& dv_3818 = temps.at(1535);
  DataVector& dv_4097 = temps.at(1535);
  dv_4097 = dv_3818 * dv_4096;
  DataVector& dv_1606 = temps.at(191);
  DataVector& dv_229 = temps.at(226);
  DataVector& dv_4098 = temps.at(1652);
  dv_4098 = dv_1606 * dv_229;
  DataVector& dv_1607 = temps.at(189);
  DataVector& dv_3827 = temps.at(249);
  DataVector& dv_4099 = temps.at(249);
  dv_4099 = dv_1607 * dv_3827;
  DataVector& dv_4100 = temps.at(975);
  dv_4100 = dv_4 * d[2921];
  DataVector& dv_167 = temps.at(165);
  DataVector& dv_1678 = temps.at(1512);
  DataVector& dv_1826 = temps.at(1611);
  DataVector& dv_2603 = temps.at(2185);
  DataVector& dv_3469 = temps.at(2944);
  DataVector& dv_490 = temps.at(458);
  sc_7 = (4.0 * d[2574]) * dv_490 +
         d[19] * ((-d[20]) * (dv_0 * dv_3469 + dv_1826 * d[2923]) +
                  d[77] * (dv_167 + dv_1678) + dv_0 * dv_2603);
  DataVector& dv_114 = temps.at(112);
  DataVector& dv_1830 = temps.at(1615);
  DataVector& dv_1833 = temps.at(1618);
  DataVector& dv_241 = temps.at(237);
  DataVector& dv_3353 = temps.at(2837);
  sc_7 += d[20] * ((-d[20]) * (dv_1833 * d[2923] - dv_241) + d[104] * dv_1830 +
                   d[77] * (dv_114 + dv_3353));
  DataVector& dv_163 = temps.at(161);
  DataVector& dv_3225 = temps.at(2712);
  DataVector& dv_731 = temps.at(680);
  sc_7 += d[28] * (d[1] * dv_3225 + d[2922] * ((d[144] + d[91]) * dv_163 +
                                               d[468] * dv_14 + dv_731)) +
          d[55] * (-dv_1818 + dv_4091 * d[2923]);
  sc_5 = (2.0 * d[2928] * d[12]) * sc_7;
  DataVector& dv_1586 = temps.at(1349);
  DataVector& dv_1810 = temps.at(1596);
  DataVector& dv_253 = temps.at(246);
  DataVector& dv_3823 = temps.at(64);
  DataVector& dv_3824 = temps.at(250);
  DataVector& dv_3834 = temps.at(266);
  DataVector& dv_4101 = temps.at(1650);
  dv_4101 = (-d[651]) * dv_3824 +
            (-d[2647] * d[2921] - d[2923] * (d[1] + d[1045])) * dv_3823 +
            (-d[2921]) * dv_3834 + (2.0 * d[22] * d[2923]) * dv_1586 +
            (4.0 * d[48] * d[92] * d[2926] * d[2921]) * dv_253 -
            dv_1810 * dv_4063 + sc_5;
  DataVector& dv_1815 = temps.at(1601);
  dv_4101 +=
      -dv_1815 * ((d[37] + d[2921]) * dv_5 + Dy * d[227] +
                  d[2920] * ((-d[29]) * Dx + (4.0 * d[2928] * d[2922]) * Dy));
  DataVector& dv_119 = temps.at(117);
  DataVector& dv_271 = temps.at(264);
  DataVector& dv_3833 = temps.at(233);
  DataVector& dv_898 = temps.at(830);
  dv_4101 += -dv_3833 *
             (d[1] * ((d[104] + d[20]) * Dy + dv_4100) +
              d[2922] * ((d[48] + d[897]) * dv_231 + d[348] * dv_119 + dv_898) +
              d[2923] * ((d[348] + d[49]) * dv_5 + (-d[112]) * dv_5 + dv_271));
  DataVector& dv_2 = temps.at(2);
  DataVector& dv_6 = temps.at(6);
  dv_4101 += (4.0 * d[22]) * dv_2 * dv_4064 + (4.0 * d[48] * d[92] * d[2926]) *
                                                  dv_6 *
                                                  (Dx * d[87] + d[90] * dv_240);
  DataVector& dv_375 = temps.at(358);
  DataVector& dv_4102 = temps.at(233);
  dv_4102 = (-d[107]) * Dy + (-d[77]) * dv_375 + d[1063] * dv_10;
  DataVector& dv_1874 = temps.at(68);
  DataVector& dv_4103 = temps.at(1618);
  dv_4103 = dv_1874 * dv_4102;
  DataVector& dv_4104 = temps.at(266);
  dv_4104 = dv_4077 + dv_4079 - dv_4082 - dv_4083 - dv_4092;
  DataVector& dv_2181 = temps.at(1788);
  DataVector& dv_3808 = temps.at(1138);
  dv_4104 +=
      d[54] * (d[50] * dv_4087 + d[53] * ((-d[86]) * dv_3808 + dv_4088) +
               d[63] * (dv_2181 + d[2921] * ((-d[2923]) * dv_4090 + dv_4089)) -
               dv_4084 - dv_4086);
  DataVector& dv_1595 = temps.at(135);
  DataVector& dv_4105 = temps.at(680);
  dv_4105 = (-d[16]) * dv_4094 + d[1063] * dv_1595 + dv_10 * dv_4104 + dv_4093;
  DataVector& dv_4106 = temps.at(458);
  dv_4106 = (3.0 * d[142]) * dv_119;
  DataVector& dv_4107 = temps.at(940);
  dv_4107 = 3.0 * dv_4;
  DataVector& dv_34 = temps.at(34);
  DataVector& dv_4108 = temps.at(64);
  dv_4108 = d[2923] * (dv_34 + dv_4107);
  DataVector& dv_1582 = temps.at(1454);
  DataVector& dv_4109 = temps.at(246);
  dv_4109 = dv_1582 + dv_4;
  DataVector& dv_4110 = temps.at(504);
  DataVector& dv_550 = temps.at(504);
  dv_4110 = (-2.0 * d[2920]) * Dy + d[29] * dv_550;
  DataVector& dv_4111 = temps.at(1615);
  dv_4111 = 14.0 * dv_119;
  DataVector& dv_4112 = temps.at(133);
  dv_4112 = (24.0 * d[142]) * dv_119 + (-d[1627] * d[2923]) * dv_4109 +
            d[2923] * (25.0 * dv_4 + 36.0 * dv_5);
  DataVector& dv_289 = temps.at(281);
  DataVector& dv_4113 = temps.at(264);
  dv_4113 = (d[1137] * d[539] * d[839] * d[2926]) * dv_289;
  DataVector& dv_4114 = temps.at(1596);
  dv_4114 = (-5.0 * d[19]) * Dy + Dy * d[348] + d[354] * dv_4;
  DataVector& dv_4115 = temps.at(956);
  dv_4115 = d[697] * dv_5 + d[938] * dv_1531;
  DataVector& dv_4116 = temps.at(824);
  DataVector& dv_892 = temps.at(824);
  dv_4116 = d[2584] * dv_892;
  DataVector& dv_4117 = temps.at(1611);
  DataVector& dv_752 = temps.at(700);
  dv_4117 = 6.0 * dv_752;
  DataVector& dv_4118 = temps.at(830);
  dv_4118 = Dy * d[535];
  DataVector& dv_3485 = temps.at(2960);
  DataVector& dv_4119 = temps.at(1512);
  dv_4119 = (d[1384] + d[1548]) * dv_4118 + dv_3485;
  DataVector& dv_4120 = temps.at(1349);
  dv_4120 = d[306] * dv_231;
  DataVector& dv_1514 = temps.at(1393);
  DataVector& dv_3911 = temps.at(3131);
  DataVector& dv_4121 = temps.at(250);
  dv_4121 = dv_1514 + dv_3911;
  DataVector& dv_3910 = temps.at(1130);
  DataVector& dv_4122 = temps.at(953);
  dv_4122 = (d[1248] + 29.0) * dv_3910;
  DataVector& dv_2725 = temps.at(2278);
  DataVector& dv_4123 = temps.at(877);
  DataVector& dv_751 = temps.at(699);
  dv_4123 = (-d[1786] - d[256]) * dv_3910 + d[1959] * dv_751 + dv_2725;
  DataVector& dv_2667 = temps.at(2230);
  DataVector& dv_4124 = temps.at(169);
  dv_4124 =
      d[2542] * ((d[1421] + d[2100]) * dv_752 + d[1961] * dv_751 + dv_2667);
  DataVector& dv_4125 = temps.at(670);
  DataVector& dv_833 = temps.at(776);
  dv_4125 = d[20] * ((-d[258]) * dv_1115 +
                     (-20.0 * d[1] - 20.0 * d[1386]) * dv_833 + d[250] * dv_0);
  DataVector& dv_3567 = temps.at(3040);
  DataVector& dv_4126 = temps.at(158);
  dv_4126 =
      (d[2928] + 19.0 * d[3]) * dv_833 + (-d[245]) * dv_3567 + d[254] * dv_0;
  DataVector& dv_4127 = temps.at(838);
  DataVector& dv_906 = temps.at(838);
  dv_4127 = d[398] * dv_906;
  DataVector& dv_2008 = temps.at(167);
  DataVector& dv_4128 = temps.at(380);
  dv_4128 = (-26.0 * d[2922]) * Dy + dv_2008;
  DataVector& dv_1709 = temps.at(1534);
  DataVector& dv_4129 = temps.at(1490);
  dv_4129 = (d[2928] + d[1932]) * dv_1709 + (-60.0 * d[2928] * d[6]) * Dy +
            d[211] * dv_0;
  DataVector& dv_4130 = temps.at(138);
  dv_4130 = d[1468] * dv_263;
  DataVector& dv_2234 = temps.at(1835);
  DataVector& dv_4131 = temps.at(149);
  dv_4131 = (-d[2609] + d[320] + d[321]) * dv_2234 + dv_4130;
  DataVector& dv_2208 = temps.at(1811);
  DataVector& dv_4132 = temps.at(1513);
  dv_4132 = -dv_2208 + dv_751;
  DataVector& dv_265 = temps.at(258);
  DataVector& dv_4133 = temps.at(1520);
  dv_4133 = (-2.0 * d[48]) * dv_4132 +
            d[20] * (113.0 * dv_751 + 60.0 * dv_752) +
            d[259] * ((d[2033] + d[2596]) * Dx + 38.0 * dv_265);
  DataVector& dv_2949 = temps.at(2477);
  DataVector& dv_4134 = temps.at(1513);
  dv_4134 =
      d[53] * (d[20] * (221.0 * dv_751 + 80.0 * dv_752) + d[221] * dv_4132 +
               d[77] * ((d[2048] + d[2828]) * Dx + dv_2949));
  DataVector& dv_3637 = temps.at(3107);
  DataVector& dv_4135 = temps.at(1339);
  dv_4135 =
      d[182] * ((d[1298] * d[2923] + d[1511] * (d[1276] - 1.0) + d[563]) * Dy +
                d[217] * dv_0 + 56.0 * dv_3637);
  DataVector& dv_1 = temps.at(1);
  DataVector& dv_4136 = temps.at(935);
  dv_4136 = d[27] * dv_1;
  DataVector& dv_4137 = temps.at(1654);
  dv_4137 = (-d[185] - d[279] + d[48]) * dv_3567;
  DataVector& dv_3528 = temps.at(3002);
  DataVector& dv_4138 = temps.at(1573);
  DataVector& dv_780 = temps.at(728);
  dv_4138 = -dv_3528 + dv_780;
  DataVector& dv_4139 = temps.at(1307);
  dv_4139 = 62.0 * dv_752;
  DataVector& dv_4140 = temps.at(3192);
  dv_4140 =
      d[122] * ((-d[2831]) * dv_3910 + (d[267] * d[2923]) * Dx +
                (-d[1123] * d[2640]) * Dx) +
      d[205] * ((d[2605] * d[6] + d[297]) * Dx +
                (d[1682] + d[259] * (d[1368] + 41.0) - d[2606]) * dv_3910) +
      d[2591] * dv_3912;
  dv_4140 += d[53] *
             (d[101] * (Dx * d[2832] - 164.0 * dv_265) + d[271] * dv_3921 +
              d[273] * ((d[535] + 13.0) * dv_2008 + d[268] * dv_752 + dv_4139) +
              d[50] * (dv_2949 + dv_4138));
  DataVector& dv_1762 = temps.at(1565);
  DataVector& dv_3896 = temps.at(3126);
  dv_4140 +=
      d[57] * ((-d[88] + 19.0 * d[2921] * d[2923]) * dv_1762 + d[2602] * dv_0) +
      d[63] *
          ((d[2671] - 180.0 * d[277] + d[2834] * d[99]) * dv_0 +
           (d[2833] + d[364]) * dv_3896 + (d[1] * (d[1668] + d[2551])) * dv_5);
  DataVector& dv_4141 = temps.at(1661);
  dv_4141 =
      d[1199] * ((d[2928] * d[20] * (d[2145] + 163.0 * d[6]) - d[2146]) * Dx +
                 (d[1768] + d[2130] + d[2611]) * dv_752);
  DataVector& dv_1720 = temps.at(1540);
  dv_4141 += d[122] * (d[104] * (-dv_3911 - 11.0 * dv_751) +
                       d[348] * (-dv_1720 - dv_3911) +
                       d[77] * ((d[1674] + 86.0 * d[6]) * Dx + 70.0 * dv_265));
  dv_4141 +=
      d[2591] * (Dx * d[2123] + d[1088] * dv_752 - 27.0 * dv_231) +
      d[337] *
          ((2.0 * d[2928] * d[20] * (d[2140] + 77.0 * d[6]) - d[2143]) * Dx +
           (-66.0 * d[101] + d[1693] * d[273] + d[2131] - d[544]) * dv_2208);
  DataVector& dv_4142 = temps.at(3131);
  dv_4142 = d[1450] * (dv_751 + 42.0 * dv_752);
  DataVector& dv_4143 = temps.at(3115);
  dv_4143 = d[221] * dv_751;
  DataVector& dv_4144 = temps.at(2947);
  dv_4144 = 338.0 * dv_265;
  DataVector& dv_4145 = temps.at(3202);
  dv_4145 = d[155] * dv_752;
  DataVector& dv_4146 = temps.at(1684);
  dv_4146 = d[116] * ((d[1248] + 11.0) * dv_1514 + dv_4145 + 47.0 * dv_752);
  DataVector& dv_4147 = temps.at(1350);
  dv_4147 = d[2831] * dv_1531;
  DataVector& dv_4148 = temps.at(1090);
  dv_4148 = 8.0 * dv_1115;
  DataVector& dv_1530 = temps.at(1407);
  DataVector& dv_3902 = temps.at(1441);
  DataVector& dv_4149 = temps.at(3204);
  dv_4149 =
      (d[2674] * d[362]) * dv_1530 +
      d[123] * (d[104] * (dv_2008 + dv_3902 + dv_4117) +
                d[116] * ((d[1248] + 23.0) * dv_751 + dv_4145 + 9.0 * dv_752) +
                d[259] * ((-d[2110] + d[288] + 37.0) * Dx - 328.0 * dv_265));
  DataVector& dv_3847 = temps.at(3128);
  dv_4149 += d[124] *
             ((-d[259]) * ((d[2115] + d[2842]) * Dx + 842.0 * dv_265) +
              d[104] * ((d[1596] + 5.0) * dv_751 + 10.0 * dv_752) +
              d[116] * ((d[1248] + 17.0) * dv_1514 + dv_3847 + 56.0 * dv_752));
  DataVector& dv_3374 = temps.at(2856);
  DataVector& dv_3375 = temps.at(2857);
  DataVector& dv_3858 = temps.at(245);
  dv_4149 +=
      d[127] * ((d[2107] + d[2118] * d[259] + d[350] * d[2923]) * dv_240 +
                (-d[1315] + 24.0 * d[2921] * d[2923]) * dv_3858 + dv_3374) +
      d[128] * ((d[2112] - d[2120] * d[259] + d[749]) * dv_240 +
                (d[1628] + d[1637] + d[77]) * dv_4148 + dv_3375);
  dv_4149 += d[299] *
             ((-d[2928] - d[1958]) * dv_3567 + d[1002] * dv_1 + d[2724] * dv_0);
  DataVector& dv_4150 = temps.at(3128);
  dv_4150 = 13.0 * dv_751;
  DataVector& dv_4151 = temps.at(3202);
  dv_4151 = -46.0 * dv_752;
  DataVector& dv_4152 = temps.at(342);
  dv_4152 = 8.0 * dv_752;
  DataVector& dv_3914 = temps.at(1148);
  DataVector& dv_4153 = temps.at(306);
  dv_4153 = d[1333] * dv_3914;
  DataVector& dv_4154 = temps.at(1651);
  dv_4154 = -39.0 * dv_265;
  DataVector& dv_101 = temps.at(99);
  DataVector& dv_165 = temps.at(163);
  DataVector& dv_4155 = temps.at(3195);
  dv_4155 = dv_101 + dv_165;
  DataVector& dv_4156 = temps.at(163);
  DataVector& dv_822 = temps.at(767);
  dv_4156 = dv_165 - dv_822;
  DataVector& dv_1705 = temps.at(107);
  DataVector& dv_4157 = temps.at(107);
  DataVector& dv_601 = temps.at(552);
  dv_4157 = dv_1705 + dv_601;
  DataVector& dv_3038 = temps.at(446);
  DataVector& dv_4158 = temps.at(1648);
  DataVector& dv_602 = temps.at(553);
  dv_4158 = -dv_3038 + dv_602;
  DataVector& dv_4159 = temps.at(1030);
  DataVector& dv_648 = temps.at(598);
  dv_4159 = -dv_648;
  DataVector& dv_168 = temps.at(166);
  DataVector& dv_4160 = temps.at(986);
  dv_4160 = dv_168 + dv_4159;
  DataVector& dv_4161 = temps.at(1344);
  dv_4161 = Dx * d[129];
  DataVector& dv_4162 = temps.at(407);
  DataVector& dv_432 = temps.at(407);
  DataVector& dv_496 = temps.at(461);
  dv_4162 = dv_432 + dv_496;
  DataVector& dv_4163 = temps.at(489);
  dv_4163 = Dy * d[120];
  DataVector& dv_4164 = temps.at(1653);
  dv_4164 = -178.0 * dv_14;
  DataVector& dv_4165 = temps.at(1658);
  dv_4165 = 49.0 * dv_14;
  DataVector& dv_4166 = temps.at(1300);
  dv_4166 = -142.0 * dv_14;
  DataVector& dv_4001 = temps.at(348);
  DataVector& dv_4167 = temps.at(348);
  dv_4167 = dv_1607 * dv_4001;
  DataVector& dv_4168 = temps.at(994);
  dv_4168 = d[2853] * dv_114;
  DataVector& dv_143 = temps.at(141);
  DataVector& dv_4169 = temps.at(988);
  dv_4169 = d[591] * dv_143;
  DataVector& dv_1218 = temps.at(1113);
  DataVector& dv_4170 = temps.at(996);
  dv_4170 = d[2788] * dv_1218;
  DataVector& dv_4171 = temps.at(3132);
  dv_4171 = d[2781] * dv_1218;
  DataVector& dv_1406 = temps.at(1295);
  DataVector& dv_4172 = temps.at(1197);
  dv_4172 = d[1825] * dv_1406;
  DataVector& dv_4173 = temps.at(1017);
  dv_4173 = d[680] * dv_4171;
  DataVector& dv_1400 = temps.at(1289);
  DataVector& dv_4174 = temps.at(1289);
  dv_4174 = d[2855] * dv_1400;
  DataVector& dv_4175 = temps.at(3143);
  dv_4175 = d[151] * dv_4171;
  DataVector& dv_1402 = temps.at(1291);
  DataVector& dv_4176 = temps.at(1291);
  dv_4176 = d[2781] * dv_1402;
  DataVector& dv_4177 = temps.at(1105);
  dv_4177 = d[22] * dv_4174;
  DataVector& dv_1215 = temps.at(1110);
  DataVector& dv_4178 = temps.at(1174);
  dv_4178 = d[2855] * dv_1215;
  DataVector& dv_4179 = temps.at(3148);
  dv_4179 = d[40] * dv_4178;
  DataVector& dv_4180 = temps.at(1293);
  dv_4180 = d[2855] * dv_1404;
  DataVector& dv_3746 = temps.at(3210);
  DataVector& dv_4181 = temps.at(398);
  dv_4181 = d[2855] * dv_3746;
  DataVector& dv_1230 = temps.at(1125);
  DataVector& dv_4182 = temps.at(1125);
  dv_4182 = d[36] * dv_1230;
  DataVector& dv_4183 = temps.at(1657);
  dv_4183 = d[2860] * dv_163;
  DataVector& dv_1950 = temps.at(1666);
  DataVector& dv_4184 = temps.at(443);
  dv_4184 = (d[22] * d[852]) * dv_1950;
  DataVector& dv_4185 = temps.at(3140);
  dv_4185 = d[2855] * dv_1426;
  DataVector& dv_4186 = temps.at(1019);
  dv_4186 = d[619] * dv_4185;
  DataVector& dv_4187 = temps.at(1010);
  dv_4187 = d[2799] * dv_4185;
  DataVector& dv_3028 = temps.at(2182);
  DataVector& dv_4188 = temps.at(1003);
  dv_4188 = d[83] * dv_3028;
  DataVector& dv_4189 = temps.at(3142);
  dv_4189 = d[83] * dv_4080;
  DataVector& dv_4190 = temps.at(928);
  dv_4190 = d[1669] * dv_4176;
  DataVector& dv_4191 = temps.at(367);
  dv_4191 = d[2721] * dv_4171;
  DataVector& dv_4192 = temps.at(1352);
  dv_4192 = d[2746] * dv_4175;
  DataVector& dv_4193 = temps.at(920);
  dv_4193 = dv_4190 * d[2926];
  DataVector& dv_4194 = temps.at(1656);
  dv_4194 = d[31] * dv_4178;
  DataVector& dv_4195 = temps.at(3210);
  dv_4195 = d[2864] * dv_3746;
  DataVector& dv_3757 = temps.at(3183);
  DataVector& dv_4196 = temps.at(3183);
  dv_4196 = d[2864] * dv_3757;
  DataVector& dv_1432 = temps.at(1315);
  DataVector& dv_4197 = temps.at(1053);
  dv_4197 = (d[2865] * d[3]) * dv_1432;
  DataVector& dv_4008 = temps.at(1066);
  DataVector& dv_4198 = temps.at(1066);
  dv_4198 = d[151] * dv_4008;
  DataVector& dv_4199 = temps.at(1061);
  dv_4199 = d[2565] * dv_4185;
  DataVector& dv_4200 = temps.at(982);
  dv_4200 = d[2867] * dv_1432;
  DataVector& dv_4201 = temps.at(3133);
  dv_4201 = d[648] * dv_4185;
  DataVector& dv_1079 = temps.at(977);
  DataVector& dv_4040 = temps.at(1644);
  DataVector& dv_4202 = temps.at(1194);
  dv_4202 = (-d[1155]) * dv_4197 + (-d[1182]) * dv_1079 + (-d[124]) * dv_4040 +
            (-d[2393]) * dv_4168 + (-d[2447]) * dv_4181 + (-d[2490]) * dv_4179 +
            (-d[2494]) * dv_4180 + (-d[2496]) * dv_4178 + (-d[2498]) * dv_4180 +
            (-d[2500]) * dv_4174;
  DataVector& dv_1015 = temps.at(915);
  DataVector& dv_1137 = temps.at(1033);
  dv_4202 += (-d[2501]) * dv_4179 + (-d[2513]) * dv_4185 + (-d[252]) * dv_1015 +
             (-d[252]) * dv_1137 + (-d[2525]) * dv_4185 + (-d[2526]) * dv_4185 +
             (-d[2530]) * dv_4186 + (-d[2532]) * dv_4185 +
             (-d[2535]) * dv_4185 + (-d[2537]) * dv_4187;
  dv_4202 += (-d[2538]) * dv_4187 + (-d[2544]) * dv_4185 +
             (-d[2757]) * dv_4188 + (-d[2758]) * dv_4189 +
             (-d[2763]) * dv_4195 + (-d[2764]) * dv_4193 +
             (-d[2765]) * dv_4195 + (-d[2768]) * dv_4193 +
             (-d[2775]) * dv_4190 + (-d[2775]) * dv_4194;
  DataVector& dv_1109 = temps.at(1006);
  DataVector& dv_3724 = temps.at(3190);
  DataVector& dv_3730 = temps.at(3194);
  DataVector& dv_3733 = temps.at(3197);
  DataVector& dv_3751 = temps.at(3215);
  DataVector& dv_3999 = temps.at(313);
  DataVector& dv_4010 = temps.at(1045);
  DataVector& dv_4034 = temps.at(1238);
  DataVector& dv_994 = temps.at(894);
  dv_4202 += (-d[2816]) * dv_4191 + (-d[2853]) * dv_994 + (-d[2855]) * dv_3724 +
             (-d[2855]) * dv_3730 + (-d[2855]) * dv_3733 +
             (-d[2860]) * dv_1109 + (-d[2861]) * dv_4034 +
             (-d[2863]) * dv_3751 + (-d[2864]) * dv_3999 + (-d[2866]) * dv_4010;
  DataVector& dv_1037 = temps.at(937);
  DataVector& dv_1051 = temps.at(950);
  DataVector& dv_1136 = temps.at(1032);
  DataVector& dv_3988 = temps.at(1097);
  DataVector& dv_4011 = temps.at(1131);
  DataVector& dv_557 = temps.at(510);
  dv_4202 += (-d[2866]) * dv_4011 + (-d[319]) * dv_1051 + (-d[636]) * dv_557 +
             (-d[703]) * dv_4199 + (-d[719]) * dv_1136 + (-d[721]) * dv_4180 +
             (-d[722]) * dv_4183 + (-d[765]) * dv_4177 + (-d[2921]) * dv_1037 +
             (-d[2921]) * dv_3988;
  DataVector& dv_2271 = temps.at(1867);
  DataVector& dv_3995 = temps.at(1195);
  DataVector& dv_4017 = temps.at(1660);
  DataVector& dv_4032 = temps.at(1642);
  dv_4202 += (-d[2921]) * dv_3995 + (-64.0 * d[123]) * dv_4032 +
             (d[106] * d[632]) * dv_4183 + (d[1147] * d[670]) * dv_2271 +
             (d[1182] * d[2415]) * dv_4017 + (d[1183] * d[724]) * dv_4186 +
             (d[119] * d[2768]) * dv_4201 + (d[124] * d[2171]) * dv_4038 +
             (d[128] * d[2795]) * dv_4172 + (d[13] * d[582]) * dv_4200;
  DataVector& dv_1005 = temps.at(905);
  DataVector& dv_1044 = temps.at(944);
  DataVector& dv_131 = temps.at(129);
  DataVector& dv_1465 = temps.at(1347);
  DataVector& dv_15 = temps.at(15);
  DataVector& dv_3679 = temps.at(3149);
  DataVector& dv_836 = temps.at(779);
  dv_4202 += (d[132] * d[2855]) * dv_1465 + (d[157] * d[2400]) * dv_15 +
             (d[180] * d[647]) * dv_3679 + (d[1882] * d[637]) * dv_143 +
             (d[19] * d[2369]) * dv_4170 + (d[19] * d[252]) * dv_1005 +
             (d[19] * d[633]) * dv_836 + (d[2] * d[2923]) * dv_1044 +
             (d[2087] * d[608]) * dv_131 + (d[22] * d[750]) * dv_4174;
  DataVector& dv_1110 = temps.at(1007);
  DataVector& dv_1953 = temps.at(1669);
  DataVector& dv_3665 = temps.at(3135);
  DataVector& dv_4007 = temps.at(1374);
  DataVector& dv_4015 = temps.at(1338);
  DataVector& dv_739 = temps.at(688);
  dv_4202 += (d[2389] * d[601]) * dv_3665 + (d[2415] * d[821]) * dv_739 +
             (d[2430] * d[6]) * dv_1110 + (d[2491] * d[2856]) * dv_4198 +
             (d[271] * d[2855]) * dv_1953 + (d[273] * d[632]) * dv_3730 +
             (d[2734] * d[2853]) * dv_594 + (d[2735] * d[652]) * dv_4169 +
             (d[2782] * d[462]) * dv_4007 + (d[2809] * d[683]) * dv_4015;
  DataVector& dv_1433 = temps.at(1316);
  DataVector& dv_4044 = temps.at(896);
  dv_4202 += (d[2813] * d[53]) * dv_4185 + (d[2814] * d[720]) * dv_4192 +
             (d[2855] * d[291]) * dv_1433 + (d[2855] * d[654]) * dv_1109 +
             (d[2855] * d[679]) * dv_1215 + (d[2858] * d[63]) * dv_4172 +
             (d[3] * d[52]) * dv_4044 + (d[332] * d[770]) * dv_4180 +
             (d[420] * d[675]) * dv_4174 + (d[52] * d[839]) * dv_4032;
  DataVector& dv_1075 = temps.at(973);
  DataVector& dv_1180 = temps.at(1075);
  DataVector& dv_1410 = temps.at(1299);
  DataVector& dv_4048 = temps.at(1334);
  dv_4202 += (d[616] * d[654]) * dv_4169 + (d[632] * d[762]) * dv_4179 +
             (d[83] * d[976]) * dv_1180 + (-d[101] * d[2]) * dv_1005 +
             (-d[101] * d[2384]) * dv_1410 + (-d[101] * d[2401]) * dv_4018 +
             (-d[1017] * d[2861]) * dv_4048 + (-d[19] * d[2854]) * dv_1075 +
             (-d[190] * d[2810]) * dv_4173 + (-d[2435] * d[2860]) * dv_143;
  DataVector& dv_669 = temps.at(619);
  DataVector& dv_760 = temps.at(708);
  dv_4202 += (-d[2442] * d[2817]) * dv_4185 + (-d[2485] * d[598]) * dv_4200 +
             (-d[259] * d[583]) * dv_669 + (-d[259] * d[585]) * dv_760 +
             (-d[2719] * d[2921]) * dv_4170 + (-d[2756] * d[652]) * dv_163 +
             (-d[276] * d[2795]) * dv_4184 + (-d[2789] * d[2804]) * dv_1218 +
             (-d[2801] * d[601]) * dv_4171 + (-d[2811] * d[725]) * dv_4173;
  DataVector& dv_3761 = temps.at(3217);
  DataVector& dv_837 = temps.at(780);
  dv_4202 += (-d[2858] * d[3]) * dv_4184 + (-d[2859] * d[2862]) * dv_1432 +
             (-d[2862] * d[619]) * dv_1953 + (-d[2867] * d[595]) * dv_3761 +
             (-d[35] * d[580]) * dv_163 +
             (d[1017] * d[354] * d[772]) * dv_4174 +
             (d[115] * d[1178] * d[55]) * dv_4175 +
             (d[1171] * d[2717] * d[40]) * dv_4171 +
             (d[131] * d[2737] * d[48]) * dv_837 +
             (d[1339] * d[2859] * d[370]) * dv_4176;
  DataVector& dv_1452 = temps.at(1335);
  DataVector& dv_4039 = temps.at(1236);
  DataVector& dv_4047 = temps.at(1353);
  dv_4202 += (d[1417] * d[2529] * d[598]) * dv_4185 +
             (d[1417] * d[2855] * d[752]) * dv_1950 +
             (d[1803] * d[2863] * d[764]) * dv_1432 +
             (d[190] * d[2812] * d[50]) * dv_4175 +
             (d[2171] * d[2802] * d[752]) * dv_1426 +
             (d[2425] * d[758] * d[764]) * dv_4035 +
             (d[2446] * d[370] * d[52]) * dv_4047 +
             (d[2499] * d[48] * d[529]) * dv_4185 +
             (d[2808] * d[57] * d[595]) * dv_1452 +
             (d[582] * d[605] * d[758]) * dv_4039;
  DataVector& dv_1396 = temps.at(1285);
  DataVector& dv_1442 = temps.at(1325);
  DataVector& dv_48 = temps.at(48);
  dv_4202 += (-d[101] * d[2784] * d[2803]) * dv_131 +
             (-d[1147] * d[2740] * d[669]) * dv_739 +
             (-d[162] * d[599] * d[63]) * dv_48 +
             (-d[162] * d[634] * d[2921]) * dv_48 +
             (-d[1660] * d[319] * d[752]) * dv_4176 +
             (-d[2790] * d[2792] * d[319]) * dv_1396 +
             (-d[2799] * d[2808] * d[50]) * dv_1442 +
             (-d[319] * d[682] * d[775]) * dv_4175 +
             (d[1007] * d[2781] * d[760] * d[2926]) * dv_4198 +
             (d[1339] * d[22] * d[2856] * d[619]) * dv_1396;
  DataVector& dv_1192 = temps.at(1087);
  DataVector& dv_243 = temps.at(66);
  DataVector& dv_4006 = temps.at(324);
  dv_4202 += (d[1390] * d[2793] * d[370] * d[384]) * dv_4171 +
             (d[2760] * d[2815] * d[462] * d[2920]) * dv_4006 +
             (-d[1921] * d[2781] * d[2857] * d[582]) * dv_1406 +
             d[101] * dv_1192 + d[1155] * dv_4190 + d[1155] * dv_4194 +
             d[2589] * dv_4044 + d[259] * dv_3983 + d[2744] * dv_243 +
             d[2746] * dv_4182;
  DataVector& dv_1008 = temps.at(908);
  DataVector& dv_4042 = temps.at(1659);
  dv_4202 += d[2748] * dv_4199 + d[2760] * dv_4195 + d[2764] * dv_4196 +
             d[2768] * dv_4196 + d[2775] * dv_4197 + d[2796] * dv_4042 +
             d[2800] * dv_4176 + d[2817] * dv_4191 + d[2818] * dv_4201 +
             d[2854] * dv_1008;
  DataVector& dv_1143 = temps.at(1039);
  DataVector& dv_3646 = temps.at(3116);
  DataVector& dv_3669 = temps.at(3139);
  DataVector& dv_3985 = temps.at(370);
  dv_4202 += d[2855] * dv_3985 + d[2857] * dv_4192 + d[319] * dv_1143 +
             d[319] * dv_3669 + d[50] * dv_999 + d[584] * dv_557 +
             d[587] * dv_4168 + d[601] * dv_1091 + d[624] * dv_243 +
             d[63] * dv_3646;
  DataVector& dv_2314 = temps.at(1906);
  DataVector& dv_3655 = temps.at(3125);
  dv_4202 += (-d[580]) * dv_3567 * dv_4 + d[63] * dv_999 + d[703] * dv_4182 +
             d[705] * dv_4177 + d[728] * dv_4181 + d[742] * dv_4180 +
             d[847] * dv_4189 + d[856] * dv_4188 + d[909] * dv_2314 +
             dv_3655 * d[2921];
  dv_4202 += d[580] * dv_0 * dv_3485;
  DataVector& dv_4203 = temps.at(908);
  dv_4203 = (-d[17] * d[2868]) * dv_10 + dv_6;
  DataVector& dv_316 = temps.at(217);
  DataVector& dv_4204 = temps.at(217);
  dv_4204 = 30.0 * dv_316;
  DataVector& dv_188 = temps.at(185);
  DataVector& dv_3789 = temps.at(1500);
  DataVector& dv_4205 = temps.at(185);
  dv_4205 =
      d[92] *
      (dv_188 + d[2926] * ((-d[19]) * dv_1641 + d[20] * dv_3789 + dv_4204));
  DataVector& dv_4206 = temps.at(1500);
  dv_4206 = (-d[2868]) * dv_10 + dv_375;
  DataVector& dv_4207 = temps.at(1486);
  dv_4207 = -dv_4206;
  DataVector& dv_4208 = temps.at(1291);
  dv_4208 = d[92] * dv_4207;
  DataVector& dv_1783 = temps.at(1555);
  DataVector& dv_4209 = temps.at(1555);
  dv_4209 = dv_1783 * dv_4208;
  DataVector& dv_1800 = temps.at(1586);
  DataVector& dv_2107 = temps.at(1725);
  DataVector& dv_267 = temps.at(260);
  sc_5 = (d[1394] + d[1697] + d[1757]) * dv_2107 + (-d[92]) * dv_1800 +
         (-d[12] - d[159] - d[571]) * dv_1815 + d[1227] * dv_267;
  DataVector& dv_112 = temps.at(110);
  DataVector& dv_1173 = temps.at(1068);
  DataVector& dv_254 = temps.at(247);
  DataVector& dv_284 = temps.at(276);
  DataVector& dv_3779 = temps.at(1927);
  sc_5 += d[13] * ((-d[2922]) *
                       ((-d[91]) * dv_3779 + (-d[121] * (d[144] + d[49])) * Dx +
                        dv_112 + dv_284) +
                   (-d[2923]) * (d[116] * dv_854 - dv_1173 + dv_138 + dv_294) +
                   (4.0 * d[2928] * d[92]) * dv_6) +
          d[271] * dv_254;
  DataVector& dv_4210 = temps.at(988);
  dv_4210 = d[306] * sc_5;
  DataVector& dv_1741 = temps.at(170);
  DataVector& dv_4211 = temps.at(1586);
  dv_4211 = d[2869] * dv_1741;
  DataVector& dv_1751 = temps.at(241);
  DataVector& dv_4212 = temps.at(1291);
  dv_4212 = dv_1751 * dv_4208;
  DataVector& dv_4213 = temps.at(247);
  dv_4213 = 15.0 * dv_265;
  DataVector& dv_1580 = temps.at(1452);
  DataVector& dv_2189 = temps.at(1795);
  DataVector& dv_558 = temps.at(511);
  sc_7 = d[123] * ((d[1244] + d[1813]) * Dx + dv_265) +
         d[124] * ((d[284] + 19.0 * d[6] - 22.0) * Dx + 5.0 * dv_265) +
         d[125] * ((d[1620] - 8.0) * Dy + dv_2189 + dv_4118) +
         d[127] * ((d[284] - 1.0) * dv_558 + dv_1580);
  DataVector& dv_3913 = temps.at(3186);
  sc_7 += d[128] * ((d[1705] - 22.0) * Dy + d[1338] * dv_0 + dv_3567) +
          d[129] * ((d[152] + 7.0 * d[6] - 8.0) * Dx + dv_3913);
  sc_5 = (-d[2906]) * sc_7;
  DataVector& dv_1630 = temps.at(1475);
  DataVector& dv_3137 = temps.at(2626);
  DataVector& sc_4 = temps.at(3226);
  sc_4 = d[19] * ((d[1701] - d[1813] * d[2921]) * Dy +
                  (d[512] + d[556]) * dv_1630 + dv_3896) +
         d[52] * ((d[1544] - d[2876] + 6.0) * Dx - dv_4213) +
         d[2920] * ((d[104] + d[1612] + d[20] * (d[152] + d[531] + 6.0)) * Dx +
                    (d[1622] + d[306]) * dv_3137);
  DataVector& dv_3916 = temps.at(1104);
  DataVector& dv_3941 = temps.at(1188);
  sc_4 += d[2921] * ((d[104] - d[20] * (d[1580] - 6.0) + d[2559]) * Dy +
                     (d[2509] + d[9]) * dv_3916 + 14.0 * dv_3941);
  sc_7 = (-d[845]) * sc_4;
  DataVector& dv_1031 = temps.at(931);
  DataVector& dv_790 = temps.at(738);
  DataVector& sc_1 = temps.at(3223);
  sc_1 = (-d[19]) * ((-d[403] - d[2921] * (d[1189] + 11.0)) * dv_240 +
                     d[2916] * dv_5 + 196.0 * dv_1031) +
         (-d[206]) *
             ((-d[2921]) * ((d[1248] - d[1728] + 22.0) * Dx - 196.0 * dv_265) +
              d[1083] * dv_790);
  DataVector& dv_1696 = temps.at(1527);
  sc_1 +=
      d[20] * ((d[2509] + d[46]) * dv_1696 +
               (-d[314] + d[2921] * (44.0 - d[1997])) * Dy + d[1909] * dv_5) +
      d[52] * ((d[268] - d[2828] + 44.0) * Dx - dv_2949);
  sc_4 = d[2428] * sc_1;
  DataVector& dv_1564 = temps.at(1437);
  DataVector& dv_69 = temps.at(69);
  DataVector& sc_0 = temps.at(3222);
  sc_0 = (-d[19]) * ((d[2915] + d[88]) * dv_1564 + (-d[1721] - d[9]) * dv_69 +
                     d[288] * dv_5) +
         (-d[134] * d[508]) * Dx +
         d[2920] * ((-d[1299] + d[1504] * d[48] - d[2835]) * dv_3910 +
                    (-d[1358] - d[1478] - d[563]) * dv_751 +
                    (10.0 * d[280] * d[6]) * Dx);
  sc_0 += d[2921] * ((-d[892]) * dv_3567 + (-d[1511] - d[2873]) * dv_1 +
                     (2.0 * d[2928] * d[2922] * (d[1519] + d[354])) * Dx);
  sc_1 = d[43] * sc_0;
  DataVector& dv_538 = temps.at(493);
  DataVector& sc_3 = temps.at(3225);
  sc_3 = (-d[2899]) * dv_4161 +
         d[191] * ((-d[1750] * d[48] - d[1855] + d[2850] * d[6] + d[639]) * Dx +
                   (d[2922] * (d[2848] + d[2917])) * Dy) +
         d[50] * ((d[20] + d[91]) * dv_3567 + (d[29] + d[31]) * dv_4068 +
                  (d[2575] * d[48] + d[2901] + d[3] * d[42]) * dv_538);
  DataVector& dv_2125 = temps.at(1741);
  sc_3 += d[52] * ((d[132] + d[1574] + d[1817] * d[91] + d[2849]) * Dx +
                   (d[2917] + d[321]) * dv_3910) +
          d[55] * ((d[1724] + d[261]) * dv_1630 + (d[3] + d[46]) * dv_69 +
                   dv_2125) +
          d[63] * ((d[2206] * d[3] + d[2849] + d[49] * (d[1680] - 1.0)) * Dy +
                   (d[2928] * d[2122] + d[1512] + d[1638]) * dv_0 +
                   (-d[280]) * dv_4118);
  sc_0 = d[873] * sc_3;
  DataVector& dv_264 = temps.at(257);
  DataVector& sc_2 = temps.at(3224);
  sc_2 = (-d[55]) * ((-173.0 * d[3] + d[622]) * dv_0 +
                     (d[403] + d[2921] * (d[1276] + 11.0)) * dv_240 +
                     (-d[2593]) * dv_5) +
         (-d[59]) * ((-d[27] * (d[1593] - 11.0) + d[36]) * Dy + d[1266] * dv_5 +
                     3.0 * dv_264);
  DataVector& dv_2170 = temps.at(1778);
  DataVector& dv_3922 = temps.at(1139);
  sc_2 += d[129] * ((-d[1593] + d[1868] - 11.0) * dv_2170 +
                    (3.0 * d[2922] * d[2923]) * Dy) +
          d[205] * (Dx * d[1401] + d[1668] * dv_752 +
                    d[2921] * ((d[1371] + d[2916] - 66.0) * Dx + dv_3922));
  sc_2 += d[337] * ((-d[2921]) * ((d[1596] + d[1785] + 11.0) * dv_2170 -
                                  173.0 * dv_265) +
                    Dx * d[1284] + d[1474] * dv_752) +
          d[60] * ((d[1728] - 66.0) * dv_5 + (-d[1867] - d[88]) * dv_0 +
                   d[2905] * dv_5);
  sc_3 = d[899] * sc_2;
  DataVector& dv_3931 = temps.at(52);
  DataVector& dv_4214 = temps.at(260);
  dv_4214 = (-d[2882]) *
                ((d[1801] * d[2921] - d[314]) * Dy + d[2915] * dv_0 + dv_3931 +
                 d[2920] * ((d[2646] + d[288] + 12.0) * Dx + dv_4213)) +
            (d[2908] * d[2909]) * dv_6 +
            (d[2907] * d[2910] * (d[1387] + d[557])) * dv_6 + sc_1 + sc_4 +
            sc_5 + sc_7;
  dv_4214 += d[465] * ((d[2914] * d[2921] - d[314]) * Dy + (-d[88]) * dv_0 +
                       (d[2920] * (d[286] + d[2914])) * Dx +
                       (2.0 * d[6] * d[2921]) * Dy) +
             sc_0 + sc_3;
  DataVector& dv_4215 = temps.at(1927);
  DataVector& dv_869 = temps.at(805);
  DataVector& dv_93 = temps.at(91);
  dv_4215 =
      d[2722] * ((-d[2926]) * dv_854 + d[2918] * dv_4100 + d[2919] * dv_613) +
      d[2723] * ((-d[2926]) * dv_869 + d[2918] * dv_3779 + d[2919] * dv_93) +
      d[2918] * dv_2107;
  DataVector& dv_3982 = temps.at(372);
  DataVector& dv_4216 = temps.at(50);
  dv_4216 = (-d[16]) * dv_4205 + d[2870] * dv_1588 + dv_3982 * dv_4215;
  DataVector& dv_19 = temps.at(19);
  DataVector& dv_4217 = temps.at(975);
  dv_4217 = 4.0 * dv_19;
  DataVector& dv_321 = temps.at(304);
  DataVector& dv_382 = temps.at(365);
  DataVector& dv_395 = temps.at(376);
  DataVector& dv_414 = temps.at(389);
  DataVector& dv_445 = temps.at(418);
  DataVector& dv_475 = temps.at(444);
  DataVector& sc_9 = temps.at(3231);
  sc_9 = d[119] * (d[118] * dv_395 + d[54] * dv_414 + dv_174 * dv_382) +
         d[2922] * (d[100] * dv_321 + d[114] * dv_445 + d[118] * dv_475);
  DataVector& dv_477 = temps.at(445);
  DataVector& dv_500 = temps.at(211);
  DataVector& dv_523 = temps.at(415);
  sc_9 += d[2923] *
          ((-d[114]) * dv_500 + (-d[29]) * dv_477 + (3.0 * d[12]) * dv_523);
  DataVector& sc_6 = temps.at(3228);
  sc_6 = dv_10 * sc_9;
  DataVector& dv_313 = temps.at(125);
  DataVector& dv_342 = temps.at(325);
  DataVector& dv_374 = temps.at(357);
  DataVector& dv_44 = temps.at(44);
  sc_2 = d[1] * dv_313 * dv_44 +
         dv_375 * ((-d[114]) * dv_374 + (3.0 * d[12]) * dv_342 - dv_321) + sc_6;
  DataVector& dv_308 = temps.at(239);
  sc_5 = 2.0 * sc_2 * pow(dv_308, 2.0);
  DataVector& dv_227 = temps.at(224);
  DataVector& dv_228 = temps.at(225);
  DataVector& dv_315 = temps.at(183);
  sc_7 = d[2928] * dv_227 * dv_315 + d[111] * dv_228 + sc_5;
  sc_4 = (-d[108]) * sc_7;
  DataVector& dv_1333 = temps.at(1227);
  DataVector& dv_1335 = temps.at(1229);
  DataVector& dv_1336 = temps.at(1230);
  DataVector& dv_1337 = temps.at(1228);
  DataVector& dv_1339 = temps.at(1232);
  DataVector& dv_1342 = temps.at(1235);
  DataVector& dv_1344 = temps.at(1237);
  DataVector& dv_1346 = temps.at(1239);
  DataVector& dv_1347 = temps.at(1240);
  DataVector& dv_1348 = temps.at(1241);
  sc_5 = -dv_1333 - dv_1335 - dv_1336 - dv_1337 - dv_1339 - dv_1342 - dv_1344 -
         dv_1346 - dv_1347 - dv_1348;
  DataVector& dv_1349 = temps.at(1242);
  DataVector& dv_1350 = temps.at(1226);
  DataVector& dv_1352 = temps.at(1244);
  DataVector& dv_1354 = temps.at(1246);
  DataVector& dv_1355 = temps.at(1247);
  DataVector& dv_1356 = temps.at(1248);
  DataVector& dv_1357 = temps.at(1249);
  DataVector& dv_1359 = temps.at(1251);
  DataVector& dv_1361 = temps.at(1253);
  DataVector& dv_1363 = temps.at(1255);
  sc_5 += -dv_1349 - dv_1350 - dv_1352 - dv_1354 - dv_1355 - dv_1356 - dv_1357 -
          dv_1359 - dv_1361 - dv_1363;
  DataVector& dv_1364 = temps.at(227);
  DataVector& dv_1366 = temps.at(1257);
  DataVector& dv_1368 = temps.at(1259);
  DataVector& dv_1370 = temps.at(1261);
  DataVector& dv_1371 = temps.at(1262);
  DataVector& dv_1372 = temps.at(1263);
  DataVector& dv_1373 = temps.at(1264);
  DataVector& dv_1374 = temps.at(1265);
  DataVector& dv_1375 = temps.at(1266);
  DataVector& dv_1376 = temps.at(1267);
  sc_5 += -dv_1364 - dv_1366 - dv_1368 - dv_1370 - dv_1371 - dv_1372 - dv_1373 -
          dv_1374 - dv_1375 - dv_1376;
  DataVector& dv_1377 = temps.at(1268);
  DataVector& dv_1378 = temps.at(1269);
  DataVector& dv_1379 = temps.at(1270);
  DataVector& dv_1381 = temps.at(1272);
  DataVector& dv_1383 = temps.at(1274);
  DataVector& dv_1385 = temps.at(1276);
  DataVector& dv_1387 = temps.at(1278);
  DataVector& dv_1388 = temps.at(1279);
  DataVector& dv_1390 = temps.at(1281);
  DataVector& dv_1391 = temps.at(1280);
  sc_5 += -dv_1377 - dv_1378 - dv_1379 - dv_1381 - dv_1383 - dv_1385 - dv_1387 -
          dv_1388 - dv_1390 - dv_1391;
  DataVector& dv_1394 = temps.at(1284);
  DataVector& dv_1395 = temps.at(1058);
  DataVector& dv_1397 = temps.at(1286);
  DataVector& dv_1399 = temps.at(1288);
  DataVector& dv_1401 = temps.at(1290);
  DataVector& dv_1403 = temps.at(1292);
  DataVector& dv_1405 = temps.at(1294);
  DataVector& dv_1407 = temps.at(1296);
  DataVector& dv_1408 = temps.at(1297);
  DataVector& dv_1409 = temps.at(1298);
  sc_5 += -dv_1394 - dv_1395 - dv_1397 - dv_1399 - dv_1401 - dv_1403 - dv_1405 -
          dv_1407 - dv_1408 - dv_1409;
  DataVector& dv_1412 = temps.at(1301);
  DataVector& dv_1413 = temps.at(1302);
  DataVector& dv_1415 = temps.at(1303);
  DataVector& dv_1417 = temps.at(1304);
  DataVector& dv_1421 = temps.at(666);
  DataVector& dv_1422 = temps.at(1306);
  DataVector& dv_1424 = temps.at(1308);
  DataVector& dv_1428 = temps.at(1312);
  DataVector& dv_1429 = temps.at(1311);
  DataVector& dv_1431 = temps.at(1314);
  sc_5 += -dv_1412 - dv_1413 - dv_1415 - dv_1417 - dv_1421 - dv_1422 - dv_1424 -
          dv_1428 - dv_1429 - dv_1431;
  DataVector& dv_1434 = temps.at(1317);
  DataVector& dv_1437 = temps.at(1320);
  DataVector& dv_1439 = temps.at(1322);
  DataVector& dv_1441 = temps.at(1324);
  DataVector& dv_1444 = temps.at(1327);
  DataVector& dv_1445 = temps.at(1328);
  DataVector& dv_1446 = temps.at(1329);
  DataVector& dv_1447 = temps.at(1330);
  DataVector& dv_1453 = temps.at(1336);
  DataVector& dv_1457 = temps.at(1340);
  sc_5 += -dv_1434 - dv_1437 - dv_1439 - dv_1441 - dv_1444 - dv_1445 - dv_1446 -
          dv_1447 - dv_1453 - dv_1457;
  DataVector& dv_1458 = temps.at(1326);
  DataVector& dv_1460 = temps.at(1342);
  DataVector& dv_1463 = temps.at(1345);
  DataVector& dv_1464 = temps.at(1346);
  DataVector& dv_1466 = temps.at(1348);
  DataVector& dv_1469 = temps.at(1351);
  DataVector& dv_1472 = temps.at(1354);
  DataVector& dv_1473 = temps.at(1355);
  DataVector& dv_1476 = temps.at(1358);
  DataVector& dv_1477 = temps.at(535);
  sc_5 += -dv_1458 - dv_1460 - dv_1463 - dv_1464 - dv_1466 - dv_1469 - dv_1472 -
          dv_1473 - dv_1476 - dv_1477;
  DataVector& dv_1479 = temps.at(1360);
  DataVector& dv_1481 = temps.at(1362);
  DataVector& dv_1482 = temps.at(1363);
  DataVector& dv_1483 = temps.at(1364);
  DataVector& dv_1485 = temps.at(1366);
  DataVector& dv_1489 = temps.at(1369);
  DataVector& dv_1490 = temps.at(1370);
  DataVector& dv_1491 = temps.at(1371);
  DataVector& dv_1492 = temps.at(1372);
  DataVector& dv_1493 = temps.at(1373);
  sc_5 += -dv_1479 - dv_1481 - dv_1482 - dv_1483 - dv_1485 - dv_1489 - dv_1490 -
          dv_1491 - dv_1492 - dv_1493;
  DataVector& dv_1418 = temps.at(1305);
  DataVector& dv_1420 = temps.at(1282);
  DataVector& dv_1425 = temps.at(1309);
  DataVector& dv_1438 = temps.at(1321);
  DataVector& dv_1450 = temps.at(1333);
  DataVector& dv_1454 = temps.at(1337);
  DataVector& dv_1486 = temps.at(1367);
  DataVector& dv_1487 = temps.at(1368);
  DataVector& dv_1495 = temps.at(1375);
  sc_5 += (-d[1000]) * dv_1420 + (-d[1000]) * dv_1425 + (-d[1003]) * dv_1438 +
          (-d[1011]) * dv_1438 + (-d[1012]) * dv_1454 + (-d[662]) * dv_1418 +
          (-d[932]) * dv_1486 + (-d[945]) * dv_1450 +
          (-d[40] * d[714]) * dv_1487 - dv_1495;
  DataVector& dv_1419 = temps.at(1256);
  DataVector& dv_16 = temps.at(16);
  sc_5 += (-d[6] * d[658]) * dv_1419 + (-d[658] * d[999]) * dv_16 +
          (2.0 * d[48] * d[71] * d[733]) * dv_14 +
          (d[142] * d[579] * d[83] * d[2920]) * dv_15 +
          (d[147] * d[579] * d[83] * d[2921]) * dv_14 +
          (d[71] * d[944] * d[2921] * d[2923]) * dv_14 +
          (2.0 * d[2928] * d[236] * d[6] * d[714]) * dv_15 +
          (2.0 * d[2928] * d[236] * d[6] * d[714]) * dv_16 +
          (2.0 * d[108] * d[20] * d[48] * d[695]) * dv_15 +
          (8.0 * d[151] * d[43] * d[2921] * d[2923]) * dv_1215;
  sc_5 += (8.0 * d[151] * d[43] * d[2920] * d[2922]) * dv_1218 +
          (8.0 * d[151] * d[75] * d[2920] * d[2922]) * dv_16 +
          (8.0 * d[151] * d[75] * d[2921] * d[2923]) * dv_16 +
          (d[236] * d[6] * d[793] * d[2921] * d[2923]) * dv_16 +
          (d[236] * d[7] * d[803] * d[2920] * d[2922]) * dv_16 +
          (d[579] * d[6] * d[83] * d[2921] * d[2923]) * dv_15 +
          (d[579] * d[7] * d[83] * d[2920] * d[2922]) * dv_14 +
          (2.0 * d[2928] * d[108] * d[628] * d[7] * d[733]) * dv_14 +
          (4.0 * d[2928] * d[106] * d[6] * d[610] * d[626]) * dv_16 +
          (4.0 * d[2928] * d[106] * d[610] * d[628] * d[7]) * dv_16;
  sc_5 += (4.0 * d[2928] * d[594] * d[75] * d[2920] * d[2922]) * dv_16 +
          (4.0 * d[2928] * d[594] * d[75] * d[2921] * d[2923]) * dv_16 +
          (4.0 * d[19] * d[20] * d[236] * d[48] * d[586]) * dv_14 +
          (4.0 * d[19] * d[20] * d[236] * d[48] * d[586]) * dv_15 +
          (8.0 * d[106] * d[19] * d[48] * d[582] * d[6]) * dv_16 +
          (8.0 * d[106] * d[20] * d[48] * d[582] * d[7]) * dv_16 +
          (8.0 * d[142] * d[330] * d[48] * d[626] * d[2920]) * dv_15 +
          (8.0 * d[142] * d[330] * d[48] * d[626] * d[2920]) * dv_16 +
          (8.0 * d[147] * d[330] * d[48] * d[628] * d[2921]) * dv_14 +
          (8.0 * d[147] * d[330] * d[48] * d[628] * d[2921]) * dv_16;
  sc_5 +=
      (d[2928] * d[207] * d[582] * d[733] * d[2920] * d[2922]) * dv_14 +
      (d[2928] * d[549] * d[632] * d[695] * d[2921] * d[2923]) * dv_15 +
      (2.0 * d[2928] * d[19] * d[50] * d[71] * d[969] * d[2923]) * dv_14 +
      (2.0 * d[2928] * d[20] * d[52] * d[71] * d[969] * d[2922]) * dv_15 +
      (2.0 * d[2928] * d[207] * d[582] * d[733] * d[2921] * d[2923]) * dv_14 +
      (8.0 * d[142] * d[20] * d[330] * d[48] * d[598] * d[2920]) * dv_15 +
      (8.0 * d[147] * d[19] * d[330] * d[48] * d[598] * d[2921]) * dv_14 +
      (8.0 * d[330] * d[48] * d[6] * d[626] * d[2921] * d[2923]) * dv_16 +
      (8.0 * d[330] * d[48] * d[6] * d[628] * d[2921] * d[2923]) * dv_15 +
      (8.0 * d[330] * d[48] * d[626] * d[7] * d[2920] * d[2922]) * dv_14;
  sc_5 +=
      (8.0 * d[330] * d[48] * d[628] * d[7] * d[2920] * d[2922]) * dv_16 +
      (20.0 * d[2928] * d[131] * d[2920] * d[2922] * d[2921] * d[2923]) *
          dv_14 +
      (20.0 * d[2928] * d[131] * d[2920] * d[2922] * d[2921] * d[2923]) *
          dv_15 +
      (40.0 * d[2928] * d[131] * d[2920] * d[2922] * d[2921] * d[2923]) *
          dv_16 +
      (d[2928] * d[20] * d[549] * d[595] * d[695] * d[2920] * d[2922]) * dv_15 +
      (2.0 * d[2928] * d[19] * d[330] * d[582] * d[586] * d[2921] * d[2923]) *
          dv_14 +
      (2.0 * d[2928] * d[19] * d[586] * d[632] * d[71] * d[2921] * d[2923]) *
          dv_15 +
      (2.0 * d[2928] * d[20] * d[330] * d[582] * d[586] * d[2920] * d[2922]) *
          dv_15 +
      (2.0 * d[2928] * d[20] * d[586] * d[621] * d[71] * d[2920] * d[2922]) *
          dv_14 +
      (4.0 * d[2928] * d[19] * d[20] * d[236] * d[586] * d[598] * d[6]) * dv_15;
  sc_5 +=
      (4.0 * d[2928] * d[19] * d[20] * d[236] * d[586] * d[598] * d[7]) *
          dv_14 +
      (8.0 * d[19] * d[330] * d[48] * d[598] * d[6] * d[2921] * d[2923]) *
          dv_15 +
      (8.0 * d[20] * d[330] * d[48] * d[598] * d[7] * d[2920] * d[2922]) *
          dv_14 +
      (16.0 * d[106] * d[48] * d[582] * d[2920] * d[2922] * d[2921] * d[2923]) *
          dv_16 +
      (16.0 * d[19] * d[330] * d[48] * d[598] * d[6] * d[2921] * d[2923]) *
          dv_16 +
      (16.0 * d[20] * d[330] * d[48] * d[598] * d[7] * d[2920] * d[2922]) *
          dv_16 +
      (2.0 * d[2928] * d[108] * d[598] * d[733] * d[2920] * d[2922] * d[2921] *
       d[2923]) *
          dv_14 +
      (4.0 * d[2928] * d[236] * d[586] * d[626] * d[2920] * d[2922] * d[2921] *
       d[2923]) *
          dv_14 +
      (4.0 * d[2928] * d[236] * d[586] * d[628] * d[2920] * d[2922] * d[2921] *
       d[2923]) *
          dv_15 +
      (8.0 * d[2928] * d[106] * d[598] * d[610] * d[2920] * d[2922] * d[2921] *
       d[2923]) *
          dv_16;
  DataVector& dv_1169 = temps.at(1064);
  DataVector& dv_1174 = temps.at(1069);
  sc_5 +=
      (8.0 * d[151] * d[22] * d[50] * d[2922]) * dv_1169 * dv_1174 +
      (8.0 * d[151] * d[22] * d[52] * d[2923]) * dv_1169 * dv_1174 +
      (d[2928] * d[621] * d[695] * d[2922] * d[2921]) * dv_1169 * dv_1174 +
      (2.0 * d[2928] * d[40] * d[703] * d[2922] * d[2923]) * dv_1169 * dv_1174 +
      (2.0 * d[48] * d[695] * d[2926] * d[2920] * d[2921]) * dv_1169 * dv_1174 +
      (4.0 * d[40] * d[48] * d[50] * d[586] * d[2920]) * dv_1169 * dv_1174 +
      (4.0 * d[40] * d[48] * d[52] * d[586] * d[2921]) * dv_1169 * dv_1174;
  sc_5 +=
      (8.0 * d[130] * d[48] * d[610] * d[2920] * d[2921]) * dv_1169 * dv_1174 +
      (24.0 * d[151] * d[19] * d[22] * d[2922] * d[2921]) * dv_1169 * dv_1174 +
      (24.0 * d[151] * d[20] * d[22] * d[2920] * d[2923]) * dv_1169 * dv_1174 +
      (d[2928] * d[12] * d[582] * d[695] * d[2920] * d[2923]) * dv_1169 *
          dv_1174 +
      (d[2928] * d[20] * d[595] * d[695] * d[2920] * d[2923]) * dv_1169 *
          dv_1174 +
      (2.0 * d[2928] * d[43] * d[582] * d[610] * d[2920] * d[2923]) * dv_1169 *
          dv_1174 +
      (2.0 * d[2928] * d[43] * d[582] * d[610] * d[2922] * d[2921]) * dv_1169 *
          dv_1174;
  sc_5 += (2.0 * d[2928] * d[626] * d[695] * d[2926] * d[2922] * d[2923]) *
              dv_1169 * dv_1174 +
          (4.0 * d[19] * d[40] * d[48] * d[621] * d[2922] * d[2923]) * dv_1169 *
              dv_1174 +
          (4.0 * d[20] * d[40] * d[48] * d[632] * d[2922] * d[2923]) * dv_1169 *
              dv_1174 +
          (4.0 * d[40] * d[48] * d[50] * d[595] * d[6] * d[2920]) * dv_1169 *
              dv_1174 +
          (4.0 * d[40] * d[48] * d[52] * d[595] * d[7] * d[2921]) * dv_1169 *
              dv_1174 +
          (4.0 * d[40] * d[48] * d[6] * d[621] * d[2920] * d[2921]) * dv_1169 *
              dv_1174 +
          (4.0 * d[40] * d[48] * d[632] * d[7] * d[2920] * d[2921]) * dv_1169 *
              dv_1174;
  sc_5 +=
      (2.0 * d[2928] * d[0] * d[19] * d[50] * d[586] * d[595] * d[2922]) *
          dv_1169 * dv_1174 +
      (2.0 * d[2928] * d[0] * d[19] * d[586] * d[621] * d[2922] * d[2921]) *
          dv_1169 * dv_1174 +
      (2.0 * d[2928] * d[0] * d[20] * d[52] * d[586] * d[595] * d[2923]) *
          dv_1169 * dv_1174 +
      (2.0 * d[2928] * d[0] * d[20] * d[586] * d[632] * d[2920] * d[2923]) *
          dv_1169 * dv_1174 +
      (2.0 * d[2928] * d[598] * d[695] * d[7] * d[2926] * d[2920] * d[2921]) *
          dv_1169 * dv_1174 +
      (4.0 * d[2928] * d[22] * d[582] * d[6] * d[621] * d[2920] * d[2923]) *
          dv_1169 * dv_1174 +
      (4.0 * d[2928] * d[22] * d[582] * d[632] * d[7] * d[2922] * d[2921]) *
          dv_1169 * dv_1174;
  sc_5 += (8.0 * d[19] * d[20] * d[40] * d[48] * d[595] * d[2922] * d[2923]) *
              dv_1169 * dv_1174 +
          (4.0 * d[2928] * d[19] * d[22] * d[582] * d[595] * d[7] * d[2922] *
           d[2921]) *
              dv_1169 * dv_1174 +
          (4.0 * d[2928] * d[20] * d[22] * d[582] * d[595] * d[6] * d[2920] *
           d[2923]) *
              dv_1169 * dv_1174;
  DataVector& dv_1498 = temps.at(1378);
  sc_7 = -dv_1498 * sc_5;
  sc_6 = d[131];
  DataVector& dv_106 = temps.at(104);
  DataVector& dv_22 = temps.at(22);
  DataVector& dv_31 = temps.at(31);
  DataVector& dv_533 = temps.at(488);
  DataVector& dv_537 = temps.at(492);
  DataVector& dv_539 = temps.at(63);
  DataVector& dv_540 = temps.at(494);
  DataVector& dv_546 = temps.at(500);
  DataVector& dv_62 = temps.at(62);
  sc_6 *=
      (d[141] * d[143]) * dv_106 + d[146] * dv_540 + dv_546 +
      d[2922] * ((-d[2920]) * dv_533 + dv_537) +
      d[2923] * (d[139] * ((-d[19]) * dv_31 + d[20] * dv_62 - dv_22) + dv_539);
  DataVector& dv_536 = temps.at(491);
  DataVector& dv_562 = temps.at(515);
  DataVector& dv_570 = temps.at(523);
  DataVector& dv_578 = temps.at(530);
  DataVector& dv_584 = temps.at(536);
  DataVector& dv_592 = temps.at(544);
  DataVector& sc_8 = temps.at(3230);
  sc_8 =
      (-d[6]) * (d[19] * dv_578 + dv_584) + dv_592 +
      d[2922] * ((-d[156]) * dv_536 + (-d[2920]) * dv_570 +
                 (86.0 * d[2928] * d[12] * d[2923] - d[157] - d[158] * d[50]) *
                     Dx * Dy +
                 d[52] * dv_562);
  DataVector& dv_383 = temps.at(366);
  DataVector& dv_573 = temps.at(526);
  DataVector& dv_575 = temps.at(528);
  sc_8 += d[2923] * ((-d[161]) * dv_383 + (-d[58]) * dv_573 +
                     (-d[2921]) * dv_575 + (25.0 * d[20] * d[2920]) * Dx * Dy);
  sc_9 = d[168] * sc_8;
  DataVector& dv_675 = temps.at(625);
  DataVector& dv_689 = temps.at(639);
  DataVector& dv_706 = temps.at(656);
  DataVector& dv_710 = temps.at(660);
  DataVector& dv_715 = temps.at(590);
  DataVector& sc_10 = temps.at(3232);
  sc_10 = (-d[6]) * dv_706 + d[167] * dv_710 + dv_715 +
          d[2922] * (d[52] * dv_675 + dv_689);
  DataVector& dv_161 = temps.at(159);
  DataVector& dv_690 = temps.at(640);
  DataVector& dv_692 = temps.at(642);
  DataVector& dv_695 = temps.at(645);
  DataVector& dv_699 = temps.at(649);
  sc_10 += d[2923] * ((-d[120]) * dv_690 + (-d[218]) * dv_161 +
                      d[63] * ((-d[91]) * dv_692 + d[20] * dv_695) + dv_699);
  sc_8 = d[237] * sc_10;
  sc_10 = d[331];
  DataVector& dv_446 = temps.at(419);
  DataVector& dv_644 = temps.at(594);
  DataVector& dv_814 = temps.at(760);
  DataVector& dv_815 = temps.at(761);
  DataVector& dv_825 = temps.at(770);
  DataVector& dv_826 = temps.at(771);
  DataVector& dv_827 = temps.at(764);
  DataVector& dv_832 = temps.at(565);
  DataVector& dv_845 = temps.at(788);
  sc_10 *= d[313] * ((-d[19]) * dv_446 + d[20] * dv_644) + dv_845 +
           d[2922] * ((-d[55]) * dv_815 + d[52] * dv_814 + dv_825) +
           d[2923] * (d[312] * ((-d[19]) * dv_826 + dv_827) + dv_832);
  DataVector& dv_547 = temps.at(501);
  DataVector& dv_548 = temps.at(502);
  DataVector& dv_549 = temps.at(503);
  DataVector& dv_551 = temps.at(505);
  DataVector& dv_555 = temps.at(508);
  DataVector& sc_12 = temps.at(3234);
  sc_12 =
      (-d[2920]) * dv_548 +
      d[6] * ((-d[3]) * dv_547 + (6.0 * d[2920] * d[2923]) * Dx * Dy - dv_549) +
      d[2922] * (dv_551 + dv_555 * d[2920]);
  DataVector& dv_29 = temps.at(29);
  DataVector& dv_468 = temps.at(437);
  DataVector& dv_553 = temps.at(507);
  DataVector& dv_559 = temps.at(512);
  sc_12 += d[2923] * (dv_559 + d[2921] * (-dv_29 - dv_468 - dv_553 - dv_556));
  DataVector& sc_11 = temps.at(3233);
  sc_11 = d[74] * sc_12;
  DataVector& dv_743 = temps.at(691);
  DataVector& dv_806 = temps.at(705);
  DataVector& dv_951 = temps.at(819);
  DataVector& dv_963 = temps.at(869);
  DataVector& dv_965 = temps.at(871);
  DataVector& dv_968 = temps.at(874);
  DataVector& dv_983 = temps.at(883);
  sc_2 = (d[34] * (d[523] * (-d[19] * d[522] - d[521]) + d[525]) + d[517] +
          d[2922] * (-d[36] * d[520] + 4.0 * d[518] * d[92] * d[2921])) *
             dv_965 +
         (-d[205] * d[537] - d[337] * d[533] - d[530] +
          d[55] * d[2922] * d[2923] + d[57] * d[2922] * d[2923]) *
             dv_968 +
         (-d[260]) * dv_743 + (-d[301]) * dv_806 +
         (-d[469] * d[6] * d[2920] * d[2921] + d[472] * d[2922] +
          d[2920] * (-d[180] * d[473] + d[474])) *
             dv_951 +
         (d[495] * d[496] + d[498] * d[499] - d[515]) * dv_963 + dv_983 + sc_6 +
         sc_9;
  DataVector& dv_593 = temps.at(35);
  DataVector& dv_597 = temps.at(548);
  DataVector& dv_599 = temps.at(550);
  DataVector& dv_604 = temps.at(555);
  DataVector& dv_606 = temps.at(557);
  DataVector& dv_607 = temps.at(558);
  DataVector& dv_611 = temps.at(562);
  sc_2 += d[171] *
          ((-d[2923]) * dv_604 + (-d[2922]) * (-dv_607 + dv_611 * d[2920]) +
           (8.0 * d[2928] * d[7]) * dv_16 + d[6] * dv_606 - dv_593 - dv_597 -
           dv_599);
  DataVector& dv_667 = temps.at(617);
  DataVector& dv_846 = temps.at(188);
  DataVector& dv_849 = temps.at(648);
  DataVector& dv_855 = temps.at(794);
  DataVector& dv_856 = temps.at(795);
  DataVector& dv_861 = temps.at(799);
  DataVector& dv_872 = temps.at(808);
  sc_2 +=
      d[208] * dv_667 +
      d[351] * ((-d[189]) * dv_861 + (-d[2922]) * dv_872 +
                (2.0 * d[2928] * d[92]) * dv_846 +
                (d[142] * d[174] * d[52]) * dv_16 +
                (d[147] * d[178] * d[50]) * dv_16 - dv_849 - dv_855 * dv_856) +
      sc_10 + sc_11 + sc_8;
  DataVector& dv_940 = temps.at(814);
  DataVector& dv_941 = temps.at(832);
  DataVector& dv_942 = temps.at(817);
  DataVector& dv_943 = temps.at(811);
  DataVector& dv_944 = temps.at(847);
  DataVector& dv_945 = temps.at(821);
  DataVector& dv_946 = temps.at(859);
  DataVector& dv_947 = temps.at(857);
  DataVector& dv_948 = temps.at(850);
  sc_2 += dv_948 *
          ((-d[357]) * dv_940 +
           d[2920] * ((-d[180]) * dv_944 + d[31] * (d[92] * dv_943 + dv_942) +
                      d[449] * dv_941) +
           d[2922] * ((-d[321]) * dv_945 + d[411] * dv_946 + dv_947));
  DataVector& dv_883 = temps.at(815);
  DataVector& dv_886 = temps.at(818);
  DataVector& dv_891 = temps.at(823);
  DataVector& dv_893 = temps.at(825);
  DataVector& dv_894 = temps.at(826);
  DataVector& dv_895 = temps.at(827);
  DataVector& dv_901 = temps.at(833);
  DataVector& dv_904 = temps.at(836);
  DataVector& dv_905 = temps.at(837);
  sc_2 += -dv_905 *
          ((-d[357]) * dv_883 +
           d[206] * ((-d[356]) * ((-d[92]) * dv_901 + dv_895) +
                     (d[92] * (-d[227] * d[369] + d[372])) * dv_894 + dv_904) +
           d[2922] * (d[36] * ((-d[92]) * dv_891 + dv_886) + dv_893));
  DataVector& dv_908 = temps.at(840);
  DataVector& dv_912 = temps.at(844);
  DataVector& dv_916 = temps.at(813);
  DataVector& dv_917 = temps.at(796);
  DataVector& dv_922 = temps.at(852);
  DataVector& dv_924 = temps.at(845);
  DataVector& dv_925 = temps.at(846);
  sc_2 += -dv_925 * (d[34] * ((-d[259]) * dv_912 + dv_916) + dv_908 +
                     d[2922] * ((-d[36]) * (d[49] * dv_922 + dv_917) + dv_924));
  DataVector& dv_928 = temps.at(854);
  DataVector& dv_930 = temps.at(856);
  DataVector& dv_934 = temps.at(858);
  DataVector& dv_935 = temps.at(855);
  DataVector& dv_937 = temps.at(860);
  DataVector& dv_938 = temps.at(843);
  sc_2 += -dv_935 * ((-d[2922]) * dv_930 - dv_928 + dv_934 * d[2920]) -
          dv_937 * dv_938;
  DataVector& dv_952 = temps.at(791);
  DataVector& dv_953 = temps.at(842);
  DataVector& dv_956 = temps.at(863);
  DataVector& dv_959 = temps.at(866);
  DataVector& dv_961 = temps.at(868);
  sc_2 += -dv_961 *
          ((-d[478]) * dv_953 + (-d[2922]) * dv_959 +
           (2.0 * d[479] * d[6] * d[2920] * d[2921]) * dv_6 - dv_952 - dv_956);
  DataVector& dv_224 = temps.at(221);
  sc_5 = (6.0 * d[130] * d[16]) * dv_10 * dv_224 * sc_2;
  DataVector& dv_1330 = temps.at(1224);
  DataVector& dv_1331 = temps.at(1225);
  DataVector& dv_249 = temps.at(242);
  DataVector& dv_250 = temps.at(243);
  DataVector& dv_304 = temps.at(296);
  DataVector& dv_307 = temps.at(299);
  sc_1 = (-d[73]) * dv_228 - dv_1330 * dv_1331 - dv_1330 * dv_249 -
         dv_250 * dv_304 - dv_304 * dv_307 + sc_4 + sc_7;
  DataVector& dv_11 = temps.at(11);
  DataVector& dv_248 = temps.at(55);
  sc_1 += (-d[109]) * dv_227 * (d[47] * pow(dv_227, 2.0) + dv_315) +
          (4.0 * d[2928] * d[16] * d[74]) * dv_227 * dv_229 * dv_248 +
          (12.0 * d[2928] * d[75]) * dv_11 * dv_19 * dv_227 * dv_248 + sc_5;
  DataVector& dv_237 = temps.at(234);
  DataVector& dv_239 = temps.at(235);
  DataVector& dv_67 = temps.at(67);
  DataVector& dv_72 = temps.at(72);
  DataVector& dv_74 = temps.at(70);
  sc_1 += (8.0 * d[2928] * d[16] * d[74]) * dv_10 * dv_229 * dv_44 *
          ((-d[13]) * dv_74 + dv_237 - dv_72 + d[2926] * (dv_239 + dv_67));
  DataVector& dv_225 = temps.at(222);
  sc_0 = (-d[1027]) * dv_225 * sc_1;
  DataVector& dv_20 = temps.at(20);
  DataVector& dv_223 = temps.at(220);
  DataVector& dv_33 = temps.at(33);
  sc_3 = (d[23] * d[24]) * dv_20 * (d[17] * dv_10 * dv_44 - dv_33 * dv_6) +
         (-d[26] * d[70]) * dv_20 * dv_223 + sc_0 + 1.0;
  get(get<CurvedScalarWave::Tags::Psi>(*result)) = sc_3 / sqrt(dv_19);
  DataVector& dv_1589 = temps.at(49);
  DataVector& dv_1596 = temps.at(176);
  DataVector& dv_1597 = temps.at(145);
  DataVector& dv_1623 = temps.at(1469);
  DataVector& dv_1624 = temps.at(1391);
  DataVector& dv_1665 = temps.at(1506);
  DataVector& dv_180 = temps.at(9);
  sc_7 = (-d[166] * d[67]) * dv_1624 + (-d[1067] * d[166] * d[39]) * dv_1589 +
         (8.0 * d[1051] * d[43]) * dv_1596 * dv_180 +
         d[1068] * dv_1597 * dv_180 + dv_1623 + dv_1665 * d[2926];
  DataVector& dv_1505 = temps.at(1385);
  DataVector& dv_1521 = temps.at(1399);
  DataVector& dv_1522 = temps.at(1400);
  DataVector& dv_1605 = temps.at(205);
  DataVector& dv_1625 = temps.at(1470);
  DataVector& dv_1626 = temps.at(1471);
  sc_7 += d[1071] * dv_1596 *
          ((-d[1069]) * dv_1605 + (-d[2928] * d[1028] * d[1070]) * dv_1626 +
           d[5] * dv_1505 + d[5] * dv_1522 + d[68] * dv_1625 + dv_1521);
  DataVector& dv_12 = temps.at(12);
  DataVector& dv_1666 = temps.at(212);
  DataVector& dv_1704 = temps.at(1533);
  DataVector& dv_1739 = temps.at(3);
  sc_7 += d[1071] * dv_180 *
          (-dv_1739 +
           d[2922] * ((6.0 * d[48] * d[2920]) * dv_12 * dv_1666 - dv_1704));
  sc_5 = d[2926] * sc_7;
  DataVector& dv_1507 = temps.at(1387);
  DataVector& dv_1602 = temps.at(178);
  DataVector& dv_1612 = temps.at(1462);
  DataVector& dv_1613 = temps.at(1463);
  DataVector& dv_1614 = temps.at(1464);
  DataVector& dv_1617 = temps.at(1467);
  DataVector& dv_1618 = temps.at(1396);
  DataVector& dv_1619 = temps.at(32);
  DataVector& dv_1621 = temps.at(1395);
  DataVector& dv_1622 = temps.at(1468);
  DataVector& dv_177 = temps.at(175);
  sc_1 = (-d[510]) * dv_1622 *
             ((-d[17]) * dv_1619 + (-d[1064] * d[39]) * dv_1613 +
              d[1035] * dv_1612 - dv_1507 * dv_1621 + dv_1614 * d[2927] +
              dv_1617 + dv_1618 * dv_177) +
         (4.0 * d[2926] * d[2927]) * dv_1602 + sc_5;
  DataVector& dv_1535 = temps.at(47);
  DataVector& dv_1537 = temps.at(1412);
  DataVector& dv_1538 = temps.at(1413);
  DataVector& dv_1585 = temps.at(1436);
  DataVector& dv_1610 = temps.at(1460);
  DataVector& dv_78 = temps.at(76);
  sc_1 += (-d[1048] * d[40]) * dv_1535 * dv_78 +
          (-d[1050] * d[17]) * dv_1535 * dv_19 +
          (-d[1051] * d[390]) * dv_1538 * dv_19 +
          (-d[17] * d[420]) * dv_1538 * dv_1610 +
          (-d[22] * d[39]) * dv_1537 * dv_1538 +
          (4.0 * d[17] * d[420]) * dv_1585 * dv_19 +
          (4.0 * d[22] * d[39]) * dv_11 * dv_1585;
  sc_1 += (-d[1066]) * dv_1602 * dv_1610 * dv_20 +
          (-16.0 * d[1025] * d[1042] * d[420]) * dv_11 * dv_1535;
  DataVector& dv_1500 = temps.at(1380);
  DataVector& dv_1501 = temps.at(1381);
  DataVector& dv_1509 = temps.at(1388);
  DataVector& dv_1611 = temps.at(1461);
  DataVector& dv_1616 = temps.at(1466);
  sc_1 += (2.0 * d[12]) * dv_1602 * dv_20 *
              ((-d[1028]) * dv_1500 + (-d[18]) * dv_1537 +
               (2.0 * d[2928] * d[5]) * dv_1509 * dv_6 +
               (3.0 * d[1028] * d[17] * d[2927]) * dv_11 - dv_1501 - dv_2) +
          (18.0 * d[2928] * d[2926]) * dv_1610 * dv_1611 * dv_1616;
  sc_0 = (-d[1047] * d[70]) * sc_1;
  sc_12 = (-d[1199]);
  DataVector& dv_195 = temps.at(192);
  DataVector& dv_2047 = temps.at(1692);
  DataVector& dv_2056 = temps.at(1698);
  DataVector& dv_2069 = temps.at(1706);
  DataVector& dv_213 = temps.at(210);
  DataVector& dv_488 = temps.at(456);
  DataVector& dv_489 = temps.at(457);
  DataVector& dv_517 = temps.at(475);
  sc_12 *= -dv_751 * ((-d[2921]) * dv_2069 + dv_240 * dv_489) +
           d[2922] * (-dv_14 * (dv_2047 + dv_2056 + 165.0 * dv_5) +
                      34.0 * dv_195 + dv_213 + dv_488 * dv_5 + dv_517);
  DataVector& dv_2070 = temps.at(411);
  DataVector& sc_14 = temps.at(3236);
  sc_14 = dv_1564 * ((-d[2921]) * dv_2070 + dv_489 * dv_538);
  DataVector& dv_2055 = temps.at(308);
  DataVector& dv_2071 = temps.at(1707);
  DataVector& dv_2073 = temps.at(1709);
  DataVector& dv_2074 = temps.at(1710);
  DataVector& dv_326 = temps.at(309);
  sc_14 += d[2923] * (dv_14 * (-534.0 * dv_15 + dv_2074 + 356.0 * dv_5) -
                      dv_1531 * dv_2073 + 69.0 * dv_195 + dv_2055 + dv_2071 +
                      33.0 * dv_326);
  DataVector& sc_13 = temps.at(3235);
  sc_13 = d[60] * sc_14;
  DataVector& dv_1832 = temps.at(1617);
  DataVector& dv_2060 = temps.at(148);
  DataVector& dv_479 = temps.at(447);
  DataVector& dv_493 = temps.at(399);
  DataVector& dv_495 = temps.at(203);
  DataVector& dv_541 = temps.at(495);
  DataVector& dv_624 = temps.at(575);
  sc_6 = (-d[55]) * (dv_1564 * (dv_2060 * d[2921] + dv_479 * dv_624) +
                     d[2923] * (-dv_14 * (dv_493 - 110.0 * dv_5) -
                                dv_1832 * dv_541 + dv_495)) +
         sc_12;
  DataVector& dv_2068 = temps.at(451);
  DataVector& dv_420 = temps.at(395);
  DataVector& dv_484 = temps.at(452);
  DataVector& dv_498 = temps.at(463);
  DataVector& dv_499 = temps.at(441);
  DataVector& dv_989 = temps.at(889);
  sc_6 += d[1198] * (dv_751 * ((-d[2921]) * dv_2068 + 4.0 * dv_484) +
                     d[2922] * (-dv_14 * (dv_498 + 75.0 * dv_5) -
                                dv_420 * dv_5 + dv_499 + 41.0 * dv_989));
  DataVector& dv_2058 = temps.at(1700);
  DataVector& dv_2059 = temps.at(1701);
  sc_6 += d[126] * (dv_2058 * dv_752 + dv_2059 * dv_751);
  DataVector& dv_202 = temps.at(199);
  DataVector& dv_205 = temps.at(202);
  DataVector& dv_2062 = temps.at(154);
  DataVector& dv_2063 = temps.at(1703);
  DataVector& dv_2064 = temps.at(1687);
  DataVector& dv_2066 = temps.at(1705);
  DataVector& dv_993 = temps.at(893);
  sc_6 += d[57] * ((-d[2923]) *
                       (dv_14 * (dv_2064 + dv_993) + 45.0 * dv_202 +
                        10.0 * dv_205 + dv_2063 + dv_2066 - 36.0 * dv_989) +
                   dv_1564 * (dv_2062 * d[2921] + dv_484)) +
          sc_13;
  sc_9 = (-d[956]) * sc_6;
  DataVector& dv_126 = temps.at(124);
  DataVector& dv_2026 = temps.at(1531);
  DataVector& dv_2084 = temps.at(1713);
  DataVector& dv_336 = temps.at(319);
  DataVector& dv_516 = temps.at(414);
  sc_14 = d[2922] * (dv_126 * dv_2084 + 24.0 * dv_202 + dv_2026 -
                     dv_29 * (57.0 * dv_5 + dv_516) + 32.0 * dv_326 + dv_336);
  DataVector& dv_2083 = temps.at(472);
  DataVector& dv_512 = temps.at(473);
  sc_14 += -dv_751 * ((-d[2921]) * dv_2083 + dv_240 * dv_512);
  sc_12 = (-d[1199]) * sc_14;
  DataVector& dv_2081 = temps.at(1711);
  DataVector& sc_15 = temps.at(3237);
  sc_15 = dv_1564 * ((-d[2921]) * dv_2081 + dv_512 * dv_538);
  DataVector& dv_2015 = temps.at(1537);
  DataVector& dv_2054 = temps.at(1697);
  DataVector& dv_21 = temps.at(21);
  DataVector& dv_363 = temps.at(346);
  DataVector& dv_571 = temps.at(524);
  sc_15 +=
      d[2923] * (dv_14 * (-426.0 * dv_15 + dv_2074 + 284.0 * dv_5) + dv_2015 +
                 63.0 * dv_202 + dv_2054 - dv_21 * (dv_16 + dv_571) + dv_363);
  sc_14 = d[60] * sc_15;
  DataVector& dv_1860 = temps.at(468);
  DataVector& dv_2078 = temps.at(467);
  DataVector& dv_451 = temps.at(424);
  DataVector& dv_520 = temps.at(462);
  DataVector& dv_522 = temps.at(439);
  sc_13 = (-d[1198]) * ((-d[2922]) * (-dv_14 * (93.0 * dv_5 + dv_520) -
                                      dv_451 * dv_5 + dv_522 + 24.0 * dv_989) +
                        dv_751 * (-4.0 * dv_1860 + dv_2078 * d[2921]));
  DataVector& dv_2044 = temps.at(1689);
  DataVector& dv_2079 = temps.at(115);
  DataVector& dv_376 = temps.at(359);
  DataVector& dv_501 = temps.at(353);
  DataVector& dv_518 = temps.at(476);
  sc_13 += (-d[55]) * (dv_1564 * (dv_2044 * dv_501 + dv_2079 * d[2921]) +
                       d[2923] * (-dv_14 * (-98.0 * dv_5 + dv_516) -
                                  dv_1832 * dv_376 + dv_518)) +
           sc_12;
  DataVector& dv_184 = temps.at(181);
  DataVector& dv_2075 = temps.at(142);
  sc_13 += d[129] * (Dy * dv_184 * d[2922] + dv_2075 * dv_751);
  DataVector& dv_2033 = temps.at(1245);
  DataVector& dv_2076 = temps.at(1702);
  DataVector& dv_568 = temps.at(521);
  sc_13 += d[57] * ((-d[2923]) *
                        (dv_14 * (dv_2064 + dv_568) + 20.0 * dv_195 +
                         10.0 * dv_202 + dv_2033 + 20.0 * dv_205 + dv_2066) +
                    dv_1564 * (dv_1860 + dv_2076 * d[2921])) +
           sc_14;
  sc_6 = (-3.0 * d[78]) * sc_13;
  sc_12 = d[63];
  DataVector& dv_17 = temps.at(17);
  DataVector& dv_1980 = temps.at(1250);
  DataVector& dv_1994 = temps.at(109);
  DataVector& dv_1995 = temps.at(103);
  DataVector& dv_2027 = temps.at(90);
  DataVector& dv_2028 = temps.at(349);
  DataVector& dv_2029 = temps.at(1523);
  DataVector& dv_412 = temps.at(354);
  sc_12 *=
      dv_0 * (-33.0 * Dy * dv_17 + dv_14 * dv_2027 + dv_2029 * d[2921]) +
      d[2923] * (-dv_14 * (dv_1995 + dv_2028) - dv_1980 - dv_1994 - dv_412);
  DataVector& dv_1997 = temps.at(114);
  DataVector& dv_1998 = temps.at(1610);
  DataVector& dv_2034 = temps.at(1671);
  DataVector& dv_207 = temps.at(204);
  DataVector& dv_400 = temps.at(381);
  sc_14 =
      (-d[50]) *
      ((-d[2923]) * (dv_14 * (dv_114 + dv_1997) + dv_1998 + dv_2034 + dv_207) +
       dv_0 * (Dy * dv_400 + dv_2031 * d[2921]));
  DataVector& dv_1657 = temps.at(1501);
  DataVector& dv_2025 = temps.at(28);
  DataVector& dv_27 = temps.at(27);
  DataVector& dv_405 = temps.at(385);
  DataVector& dv_567 = temps.at(520);
  DataVector& dv_96 = temps.at(94);
  sc_14 += d[52] * (dv_751 * (-Dy * dv_405 + dv_2025 * d[2921]) +
                    d[2922] * (dv_17 * (dv_1657 + dv_27) + dv_2026 -
                               dv_96 * (dv_143 + 17.0 * dv_5 + dv_567)));
  DataVector& dv_1659 = temps.at(1503);
  DataVector& dv_2035 = temps.at(1668);
  DataVector& dv_411 = temps.at(208);
  sc_14 += d[53] * ((-d[2922]) * (-dv_14 * (dv_411 + 33.0 * dv_5) - dv_1659 +
                                  dv_412 + 17.0 * dv_989) +
                    dv_751 * (dv_2035 * d[2921] - dv_400 * dv_538));
  DataVector& dv_2023 = temps.at(312);
  DataVector& dv_2024 = temps.at(1526);
  sc_14 += d[55] * (dv_0 * dv_2023 - dv_1 * dv_2024) + sc_12;
  sc_13 = (d[12] * d[88]) * sc_14;
  DataVector& dv_1652 = temps.at(1497);
  DataVector& dv_39 = temps.at(39);
  DataVector& dv_391 = temps.at(374);
  DataVector& dv_51 = temps.at(51);
  DataVector& dv_528 = temps.at(483);
  sc_12 = (-d[50]) * ((-d[2923]) * (dv_14 * (dv_126 + dv_25 + dv_528) + dv_391 +
                                    dv_5 * dv_51 - dv_989) +
                      dv_0 * (dv_1652 + dv_1709 * dv_39));
  DataVector& dv_1646 = temps.at(1491);
  DataVector& dv_1700 = temps.at(1530);
  sc_12 += (-d[52]) * ((-d[2922]) * (-dv_115 * (dv_1531 + dv_17) +
                                     dv_17 * (dv_17 + dv_34) + dv_195) +
                       dv_1646 * dv_1700);
  DataVector& dv_2019 = temps.at(95);
  DataVector& dv_2020 = temps.at(1521);
  sc_12 += (-d[55]) * (Dy * dv_2019 * d[2923] - dv_0 * dv_2020);
  DataVector& dv_155 = temps.at(153);
  DataVector& dv_1643 = temps.at(1488);
  DataVector& dv_2021 = temps.at(1543);
  DataVector& dv_393 = temps.at(345);
  DataVector& dv_543 = temps.at(497);
  DataVector& dv_639 = temps.at(589);
  sc_12 += (-d[63]) * (dv_1564 * ((-d[2921]) * dv_2021 + 6.0 * dv_1643) +
                       d[2923] * (-dv_1531 * dv_543 +
                                  dv_29 * (dv_155 + dv_16 + dv_639) + dv_393));
  DataVector& dv_1653 = temps.at(1498);
  DataVector& dv_1979 = temps.at(1313);
  DataVector& dv_392 = temps.at(375);
  sc_12 +=
      (d[20] * d[2920]) *
      ((-d[2922]) * (-dv_1653 * dv_34 + dv_29 * (-dv_21 - dv_392) + dv_393) +
       dv_1700 * dv_1979);
  sc_14 = d[1197] * sc_12;
  DataVector& sc_18 = temps.at(3240);
  sc_18 = (-d[1198]);
  DataVector& dv_2013 = temps.at(1532);
  DataVector& dv_2090 = temps.at(408);
  DataVector& dv_2091 = temps.at(1715);
  DataVector& dv_455 = temps.at(428);
  DataVector& dv_459 = temps.at(432);
  DataVector& dv_798 = temps.at(746);
  sc_18 *= (-d[2923]) * (dv_14 * (dv_2091 + 49.0 * dv_5 - dv_798) + dv_2013 +
                         dv_2034 - 54.0 * dv_326 + dv_336 + dv_455 * dv_5) +
           dv_0 * (dv_2090 * d[2921] + dv_240 * dv_459);
  DataVector& dv_1124 = temps.at(1021);
  DataVector& dv_2057 = temps.at(1699);
  DataVector& dv_469 = temps.at(438);
  DataVector& sc_20 = temps.at(3242);
  sc_20 =
      (-d[2922]) * (-54.0 * dv_1124 + 63.0 * dv_195 + dv_2015 + dv_2057 +
                    dv_2071 - dv_96 * (dv_469 + 84.0 * dv_5) + 76.0 * dv_989);
  DataVector& dv_1727 = temps.at(1546);
  DataVector& dv_2095 = temps.at(635);
  sc_20 += dv_1727 * (-Dy * dv_459 + dv_2095 * d[2921]);
  DataVector& sc_19 = temps.at(3241);
  sc_19 = d[60] * sc_20;
  DataVector& sc_17 = temps.at(3239);
  sc_17 = sc_18;
  DataVector& dv_2052 = temps.at(1696);
  DataVector& dv_2087 = temps.at(1484);
  DataVector& dv_2092 = temps.at(436);
  DataVector& dv_2093 = temps.at(1716);
  DataVector& dv_2094 = temps.at(1717);
  DataVector& dv_473 = temps.at(442);
  DataVector& dv_988 = temps.at(888);
  sc_17 +=
      d[1199] * (-dv_0 * ((-d[2921]) * dv_2093 + 4.0 * dv_2087) +
                 d[2923] * ((54.0 * d[2921]) * dv_988 -
                            dv_14 * (dv_2092 + dv_2094) - dv_2052 - dv_473));
  DataVector& dv_2085 = temps.at(1673);
  DataVector& dv_2086 = temps.at(1674);
  sc_17 += d[129] * (dv_0 * dv_2085 - dv_1 * dv_2086);
  DataVector& dv_2043 = temps.at(1688);
  DataVector& dv_2088 = temps.at(425);
  DataVector& dv_819 = temps.at(765);
  sc_17 += d[55] * (dv_1530 * (-dv_2087 + dv_2088 * d[2921]) +
                    d[2922] * (-dv_14 * (dv_2043 + 144.0 * dv_5) +
                               2.0 * dv_17 * (31.0 * dv_5 + dv_819) + dv_2063));
  DataVector& dv_1011 = temps.at(911);
  DataVector& dv_1598 = temps.at(156);
  DataVector& dv_2089 = temps.at(1714);
  DataVector& dv_447 = temps.at(420);
  DataVector& dv_464 = temps.at(433);
  DataVector& dv_466 = temps.at(435);
  sc_17 += d[57] * (-dv_1530 * ((-d[2921]) * dv_1598 - dv_2044 * dv_447) +
                    d[2922] * (d[27] * dv_988 - dv_1011 -
                               dv_14 * (dv_2089 + dv_464) + dv_466)) +
           sc_19;
  DataVector& sc_16 = temps.at(3238);
  sc_16 = d[12] * sc_17;
  DataVector& dv_1736 = temps.at(261);
  DataVector& dv_1857 = temps.at(269);
  DataVector& dv_1985 = temps.at(1528);
  DataVector& dv_476 = temps.at(303);
  sc_15 = d[1073] * dv_1857 + d[34] * dv_1985 + dv_1736 * dv_476 * d[2920] +
          dv_477 * d[2922] + sc_16;
  sc_12 = d[1200] * sc_15;
  DataVector& dv_324 = temps.at(307);
  DataVector& dv_426 = temps.at(401);
  DataVector& dv_463 = temps.at(193);
  sc_19 = (-d[2923]) *
          (dv_14 * (-dv_2047 + dv_463 + 55.0 * dv_5) + 34.0 * dv_202 + dv_213 +
           dv_324 - 46.0 * dv_326 + dv_426 * dv_5 - 34.0 * dv_989);
  DataVector& dv_2046 = temps.at(1691);
  DataVector& dv_428 = temps.at(403);
  sc_19 += dv_0 * (dv_2046 * d[2921] + dv_240 * dv_428);
  sc_17 = (-d[1198]) * sc_19;
  DataVector& dv_439 = temps.at(412);
  sc_18 =
      (-d[2922]) * (69.0 * dv_202 + dv_2054 + dv_2055 - dv_2056 * dv_5 +
                    dv_2057 - dv_96 * (dv_439 + 82.0 * dv_5) + 110.0 * dv_989);
  DataVector& dv_2053 = temps.at(402);
  sc_18 += dv_1530 * (dv_2053 * d[2921] - dv_428 * dv_538);
  sc_19 = (d[19] * d[20]) * sc_18;
  sc_18 = (2.0 * d[52] * d[2921]);
  DataVector& dv_2040 = temps.at(1686);
  DataVector& dv_2050 = temps.at(1694);
  DataVector& dv_2051 = temps.at(1695);
  DataVector& dv_443 = temps.at(416);
  sc_18 *= -dv_0 * ((-d[2921]) * dv_2051 + 4.0 * dv_2040) +
           d[2923] * ((46.0 * d[2921]) * dv_988 -
                      dv_14 * (dv_2050 + 178.0 * dv_5) - dv_2052 - dv_443);
  DataVector& dv_2037 = temps.at(1341);
  DataVector& dv_2039 = temps.at(1685);
  sc_16 = (-d[129]) * (Dy * dv_2039 * d[2923] - dv_0 * dv_2037) + sc_17;
  DataVector& dv_1061 = temps.at(960);
  DataVector& dv_2032 = temps.at(498);
  DataVector& dv_2045 = temps.at(1690);
  DataVector& dv_415 = temps.at(390);
  DataVector& dv_435 = temps.at(352);
  DataVector& dv_437 = temps.at(410);
  DataVector& dv_94 = temps.at(92);
  sc_16 += (-d[57]) * ((-d[2922]) * (-dv_1061 + dv_14 * (-dv_435 - dv_94) +
                                     dv_2032 + dv_437) +
                       dv_1700 * (-dv_2044 * dv_415 + dv_2045 * d[2921])) +
           sc_19;
  sc_16 += sc_18;
  DataVector& dv_1698 = temps.at(1529);
  DataVector& dv_2041 = temps.at(396);
  sc_16 +=
      d[55] * (dv_1530 * (-dv_2040 + dv_2041 * d[2921]) +
               d[2922] * (-dv_14 * (dv_2043 + 246.0 * dv_5) +
                          10.0 * dv_17 * (dv_1698 + dv_17) + 45.0 * dv_195));
  sc_15 = d[954] * sc_16;
  sc_16 = (8.0 * d[2926]) * dv_174;
  DataVector& dv_1513 = temps.at(1392);
  DataVector& dv_1553 = temps.at(1427);
  DataVector& dv_1984 = temps.at(131);
  DataVector& dv_2018 = temps.at(139);
  DataVector& dv_379 = temps.at(362);
  DataVector& dv_381 = temps.at(364);
  DataVector& dv_728 = temps.at(677);
  sc_16 *=
      (-d[19]) * (dv_1553 + dv_1984) +
      (-d[2920]) * (-dv_1513 * dv_2018 + d[2922] * (dv_379 + dv_541)) +
      d[2921] * ((-d[2923]) * (dv_1698 + dv_381) + dv_0 * (d[29] + dv_728));
  DataVector& dv_1852 = temps.at(265);
  DataVector& dv_1853 = temps.at(41);
  DataVector& dv_1854 = temps.at(1598);
  DataVector& dv_1855 = temps.at(1424);
  DataVector& dv_1983 = temps.at(1583);
  sc_8 = (d[1103] * d[88]) * dv_1854 + (d[1192] * d[12]) * dv_1853 +
         (d[94] * d[2927]) * dv_1855 + (-d[1070] * d[955]) * dv_523 +
         (-d[118] * d[2925]) * dv_523 + d[1068] * dv_1852 + d[1188] * dv_477 +
         d[1189] * dv_477 + d[1190] * dv_1983 + sc_13 + sc_6 + sc_9;
  DataVector& dv_1510 = temps.at(1389);
  DataVector& dv_1851 = temps.at(1616);
  DataVector& dv_1856 = temps.at(18);
  DataVector& dv_1858 = temps.at(1405);
  DataVector& dv_2017 = temps.at(214);
  sc_8 += d[1191] * dv_1985 + d[1194] * dv_500 + d[1195] * dv_500 +
          d[1196] * (dv_1858 + dv_477 * d[2920]) + dv_1510 * dv_1851 * dv_2017 +
          dv_1856 * d[2924] + sc_12 + sc_14 + sc_15 + sc_16;
}
}  // namespace CurvedScalarWave::Worldtube::detail
