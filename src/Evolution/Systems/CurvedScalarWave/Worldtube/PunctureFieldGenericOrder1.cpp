// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/Worldtube/PunctureField.hpp"

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/DynamicBuffer.hpp"
#include "DataStructures/Tags/TempTensor.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Utilities/Gsl.hpp"

namespace CurvedScalarWave::Worldtube {

void puncture_field_generic_1(
    gsl::not_null<Variables<tmpl::list<
        CurvedScalarWave::Tags::Psi, ::Tags::dt<CurvedScalarWave::Tags::Psi>,
        ::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>,
                      Frame::Inertial>>>*>
        result,
    const tnsr::I<DataVector, 3, Frame::Inertial>& centered_coords,
    const tnsr::I<double, 3>& particle_position,
    const tnsr::I<double, 3>& particle_velocity,
    const tnsr::I<double, 3>& particle_acceleration, const double BH_mass) {
  const size_t grid_size = get<0>(centered_coords).size();
  result->initialize(grid_size);
  const double xp = particle_position[0];
  const double yp = particle_position[1];
  const double xpdot = particle_velocity[0];
  const double ypdot = particle_velocity[1];
  const double xpddot = particle_acceleration[0];
  const double ypddot = particle_acceleration[1];
  const double rp = get(magnitude(particle_position));
  const double rpdot = (xp * xpdot + yp * ypdot) / rp;

  const auto& Dx = get<0>(centered_coords);
  const auto& Dy = get<1>(centered_coords);
  const auto& z = get<2>(centered_coords);

  const double M = BH_mass;

  DynamicBuffer<DataVector> temps(314, grid_size);

  const double d_0 = 2 * M;
  const double d_1 = rp * rp * rp;
  const double d_2 = xp * xpdot;
  const double d_3 = yp * ypdot;
  const double d_4 = d_2 + d_3;
  const double d_5 = 1.0 / d_1;
  const double d_6 = rp * rp;
  const double d_7 = 4 * M;
  const double d_8 = d_4 * d_4;
  const double d_9 = ypdot * ypdot;
  const double d_10 = xpdot * xpdot;
  const double d_11 = d_10 - 1;
  const double d_12 =
      d_0 * d_6 + d_0 * d_8 + d_1 * (d_11 + d_9) + d_4 * d_7 * rp;
  const double d_13 = 1.0 / d_12;
  const double d_14 = d_13 * d_5;
  const double d_15 = rp * rp * rp * rp * rp;
  const double d_16 = 1.0 / d_15;
  const double d_17 = xp * xp;
  const double d_18 = yp * yp;
  const double d_19 = 2 * rp;
  const double d_20 = rp * rp * rp * rp * rp * rp * rp * rp * rp * rp;
  const double d_21 = 1.0 / d_20;
  const double d_22 = rp * rp * rp * rp * rp * rp * rp;
  const double d_23 = 1.0 / (d_12 * d_12);
  const double d_24 = rp * rp * rp * rp * rp * rp;
  const double d_25 = 8 * M;
  const double d_26 = rp * rp * rp * rp;
  const double d_27 = xp * yp;
  const double d_28 = 2 * xpdot;
  const double d_29 = rp * rp * rp * rp * rp * rp * rp * rp * rp;
  const double d_30 = M * d_13;
  const double d_31 = 2 * yp;
  const double d_32 = M * d_31;
  const double d_33 = d_6 * ypdot;
  const double d_34 = 3 * d_3;
  const double d_35 = 3 * yp;
  const double d_36 = d_0 * ypdot;
  const double d_37 = d_35 + d_36;
  const double d_38 = 8 * d_30;
  const double d_39 = rp * rp * rp * rp * rp * rp * rp * rp;
  const double d_40 = M * M;
  const double d_41 = 2 * d_40;
  const double d_42 = 3 * d_6;
  const double d_43 = yp * yp * yp;
  const double d_44 = xp * xp * xp;
  const double d_45 = xp * xp * xp * xp;
  const double d_46 = yp * yp * yp * yp;
  const double d_47 = d_17 * d_18;
  const double d_48 = 4 * d_6;
  const double d_49 = d_17 * yp;
  const double d_50 = M * rp;
  const double d_51 = 2 * xp;
  const double d_52 = d_18 * xp;
  const double d_53 = 3 * d_44;
  const double d_54 = 8 * d_6;
  const double d_55 = 6 * d_26;
  const double d_56 = M * d_1;
  const double d_57 = 1.0 / d_6;
  const double d_58 = 1.0 / d_24;
  const double d_59 = 1.0 / d_26;
  const double d_60 = 3 * rpdot;
  const double d_61 = d_13 * d_60;
  const double d_62 = xp * xpddot;
  const double d_63 = yp * ypddot;
  const double d_64 = 2 * d_9;
  const double d_65 = 2 * d_10 + 2 * d_62 + 2 * d_63 + d_64 - rpdot;
  const double d_66 = M * d_57;
  const double d_67 = ypddot * ypdot;
  const double d_68 = 2 * rpdot;
  const double d_69 = d_62 + d_63;
  const double d_70 = d_10 + d_9;
  const double d_71 = d_69 + d_70;
  const double d_72 = -M * d_59 * d_60 * d_8 +
                      d_0 * d_4 * d_5 * (-d_68 + d_71) + d_67 + xpddot * xpdot;
  const double d_73 = d_65 * d_66 + d_72;
  const double d_74 = d_23 * d_73;
  const double d_75 = d_6 * d_60;
  const double d_76 = d_50 * d_51;
  const double d_77 = M * d_4;
  const double d_78 = d_51 * d_77;
  const double d_79 = d_1 * xpdot + d_78;
  const double d_80 = d_32 * d_4;
  const double d_81 = d_1 * ypdot + d_80;
  const double d_82 = 36 * M;
  const double d_83 = -yp;
  const double d_84 = -d_10 + d_69 - d_9;
  const double d_85 = d_13 * rp;
  const double d_86 = pow(rp, 11);
  const double d_87 = 1.0 / d_86;
  const double d_88 = 1.0 / (d_12 * d_12 * d_12);
  const double d_89 = 16 * M;
  const double d_90 = d_18 * ypddot;
  const double d_91 = d_9 + 1;
  const double d_92 = d_68 * rp;
  const double d_93 = 2 * d_18;
  const double d_94 = ypdot * ypdot * ypdot;
  const double d_95 = M * d_63;
  const double d_96 = M * d_9;
  const double d_97 = d_18 * d_67;
  const double d_98 = d_18 * d_9;
  const double d_99 = -d_73;
  const double d_100 = d_23 * d_99;
  const double d_101 = -d_84;
  const double d_102 = 2 * d_6;
  const double d_103 = d_1 * rpdot;
  const double d_104 = 4 * d_13 * d_26;
  const double d_105 = 4 * yp;
  const double d_106 = 2 * ypdot;
  const double d_107 = 3 * d_52;
  const double d_108 = 12 * d_16;
  const double d_109 = 16 * d_40;
  const double d_110 = d_109 * d_14;
  const double d_111 = 8 * d_66;
  const double d_112 = -d_76 + d_79;
  const double d_113 = d_109 * d_58;
  const double d_114 = M * d_16;
  const double d_115 = 1.0 / rp;
  const double d_116 = d_16 * d_82;
  const double d_117 = 8 * d_14 * d_40;
  const double d_118 = 24 * d_40 * d_58;
  const double d_119 = d_21 * d_40;
  const double d_120 = 45 * d_119;
  const double d_121 = 1.0 / d_29;
  const double d_122 = d_0 * d_121;
  const double d_123 = 12 * d_15;
  const double d_124 = d_0 * d_26;
  const double d_125 = 6 * M;
  const double d_126 = 2 * M * rp * yp - d_81;
  const double d_127 = d_17 + d_18;
  const double d_128 = -d_10 * (d_17 - d_93) + 2 * d_2 * (d_0 - d_34) +
                       ypdot * (d_106 * d_17 + yp * (-d_3 + d_7));
  const double d_129 = d_4 - rp;
  const double d_130 = -d_129;
  DataVector& dv_0 = temps.at(0);
  dv_0 = Dx * xp;
  DataVector& dv_1 = temps.at(1);
  dv_1 = Dy * yp;
  DataVector& dv_2 = temps.at(2);
  dv_2 = dv_0 + dv_1;
  DataVector& dv_3 = temps.at(3);
  dv_3 = d_0 * dv_2;
  DataVector& dv_4 = temps.at(4);
  dv_4 = dv_3 * rp;
  DataVector& dv_5 = temps.at(5);
  dv_5 = Dx * xpdot;
  DataVector& dv_6 = temps.at(6);
  dv_6 = Dy * ypdot;
  DataVector& dv_7 = temps.at(7);
  dv_7 = dv_5 + dv_6;
  DataVector& dv_8 = temps.at(8);
  dv_8 = d_1 * dv_7;
  DataVector& dv_9 = temps.at(9);
  dv_9 = d_4 * dv_3;
  DataVector& dv_10 = temps.at(10);
  dv_10 = dv_8 + dv_9;
  DataVector& dv_11 = temps.at(11);
  dv_11 = dv_10 + dv_4;
  DataVector& dv_12 = temps.at(12);
  dv_12 = dv_11 * dv_11;
  DataVector& dv_13 = temps.at(13);
  dv_13 = d_14 * dv_12;
  DataVector& dv_14 = temps.at(14);
  dv_14 = dv_2 * dv_2;
  DataVector& dv_15 = temps.at(15);
  dv_15 = d_0 * dv_14;
  DataVector& dv_16 = temps.at(16);
  dv_16 = Dx * Dx;
  DataVector& dv_17 = temps.at(17);
  dv_17 = Dy * Dy;
  DataVector& dv_18 = temps.at(18);
  dv_18 = z * z;
  DataVector& dv_19 = temps.at(19);
  dv_19 = dv_17 + dv_18;
  DataVector& dv_20 = temps.at(20);
  dv_20 = dv_16 + dv_19;
  DataVector& dv_21 = temps.at(21);
  dv_21 = d_5 * dv_15 + dv_20;
  DataVector& dv_22 = temps.at(22);
  dv_22 = -dv_13 + dv_21;
  DataVector& dv_23 = temps.at(23);
  dv_23 = sqrt(dv_22);
  DataVector& dv_24 = temps.at(24);
  dv_24 = 1.0 / dv_23;
  DataVector& dv_25 = temps.at(25);
  dv_25 = 1.0 / dv_22;
  DataVector& dv_26 = temps.at(26);
  dv_26 = M * dv_25;
  DataVector& dv_27 = temps.at(27);
  dv_27 = 6 * dv_1;
  DataVector& dv_28 = temps.at(28);
  dv_28 = dv_0 * dv_27;
  DataVector& dv_29 = temps.at(29);
  dv_29 = -dv_16;
  DataVector& dv_30 = temps.at(30);
  dv_30 = 2 * dv_17;
  DataVector& dv_31 = temps.at(31);
  dv_31 = 2 * dv_18;
  DataVector& dv_32 = temps.at(32);
  dv_32 = dv_30 + dv_31;
  DataVector& dv_33 = temps.at(33);
  dv_33 = dv_29 + dv_32;
  DataVector& dv_34 = temps.at(34);
  dv_34 = -dv_17;
  DataVector& dv_35 = temps.at(35);
  dv_35 = 2 * dv_16;
  DataVector& dv_36 = temps.at(36);
  dv_36 = dv_31 + dv_35;
  DataVector& dv_37 = temps.at(37);
  dv_37 = dv_34 + dv_36;
  DataVector& dv_38 = temps.at(38);
  dv_38 = d_17 * dv_33 + d_18 * dv_37 - dv_28;
  DataVector& dv_39 = temps.at(39);
  dv_39 = dv_2 * dv_38;
  DataVector& dv_40 = temps.at(10);
  dv_40 = dv_10 - dv_4;
  DataVector& dv_41 = temps.at(40);
  dv_41 = 4 * dv_1;
  DataVector& dv_42 = temps.at(41);
  dv_42 = dv_19 + dv_29;
  DataVector& dv_43 = temps.at(42);
  dv_43 = dv_16 + dv_18;
  DataVector& dv_44 = temps.at(43);
  dv_44 = dv_34 + dv_43;
  DataVector& dv_45 = temps.at(44);
  dv_45 = d_17 * dv_42 + d_18 * dv_44 - dv_0 * dv_41;
  DataVector& dv_46 = temps.at(45);
  dv_46 = d_19 * dv_45;
  DataVector& dv_47 = temps.at(46);
  dv_47 = -dv_38;
  DataVector& dv_48 = temps.at(47);
  dv_48 = d_4 * dv_47 + dv_46;
  DataVector& dv_49 = temps.at(48);
  dv_49 = d_13 * dv_48;
  DataVector& dv_50 = temps.at(49);
  dv_50 = dv_40 * dv_49;
  DataVector& dv_51 = temps.at(50);
  dv_51 = dv_39 + dv_50;
  DataVector& dv_52 = temps.at(51);
  dv_52 = d_16 * dv_51;
  DataVector& dv_53 = temps.at(52);
  dv_53 = 4 * dv_18;
  DataVector& dv_54 = temps.at(53);
  dv_54 = -dv_40;
  DataVector& dv_55 = temps.at(54);
  dv_55 = dv_54 * dv_54;
  DataVector& dv_56 = temps.at(55);
  dv_56 = d_23 * dv_55;
  DataVector& dv_57 = temps.at(56);
  dv_57 = d_22 * dv_56;
  DataVector& dv_58 = temps.at(57);
  dv_58 = 8 * dv_18;
  DataVector& dv_59 = temps.at(58);
  dv_59 = d_24 * dv_56;
  DataVector& dv_60 = temps.at(59);
  dv_60 = d_13 * dv_22;
  DataVector& dv_61 = temps.at(60);
  dv_61 = d_20 * dv_60;
  DataVector& dv_62 = temps.at(61);
  dv_62 = Dx * yp;
  DataVector& dv_63 = temps.at(62);
  dv_63 = Dy * xp;
  DataVector& dv_64 = temps.at(63);
  dv_64 = Dx * Dy;
  DataVector& dv_65 = temps.at(64);
  dv_65 = -d_27 * dv_18 + d_6 * dv_64;
  DataVector& dv_66 = temps.at(65);
  dv_66 = d_17 * dv_19;
  DataVector& dv_67 = temps.at(66);
  dv_67 = d_18 * dv_17 + dv_66;
  DataVector& dv_68 = temps.at(67);
  dv_68 = d_17 * dv_16;
  DataVector& dv_69 = temps.at(68);
  dv_69 = d_10 * dv_67 - d_28 * dv_65 * ypdot + d_9 * (d_18 * dv_43 + dv_68) +
          pow(dv_62 - dv_63, 2);
  DataVector& dv_70 = temps.at(69);
  dv_70 = dv_56 * dv_69;
  DataVector& dv_71 = temps.at(70);
  dv_71 = d_30 * dv_22;
  DataVector& dv_72 = temps.at(71);
  dv_72 = -dv_30;
  DataVector& dv_73 = temps.at(72);
  dv_73 = -d_9 * dv_43 * yp;
  DataVector& dv_74 = temps.at(73);
  dv_74 = d_7 * dv_43 * ypdot + dv_73 + yp * (dv_16 + dv_72);
  DataVector& dv_75 = temps.at(74);
  dv_75 = -dv_31;
  DataVector& dv_76 = temps.at(75);
  dv_76 = dv_16 + dv_75;
  DataVector& dv_77 = temps.at(65);
  dv_77 = d_18 * (dv_17 + dv_75) + dv_66;
  DataVector& dv_78 = temps.at(63);
  dv_78 = dv_64 * (d_32 - d_33) - xp * (d_0 * dv_19 - d_34 * dv_18);
  DataVector& dv_79 = temps.at(76);
  dv_79 = 2 * Dy;
  DataVector& dv_80 = temps.at(77);
  dv_80 = dv_0 * dv_79;
  DataVector& dv_81 = temps.at(78);
  dv_81 = d_10 * dv_77 + d_28 * dv_78 + d_37 * dv_80;
  DataVector& dv_82 = temps.at(79);
  dv_82 = -d_17 * (-d_9 * dv_76 + dv_17 - dv_35) - dv_74 * yp + dv_81;
  DataVector& dv_83 = temps.at(80);
  dv_83 = 4 * dv_82;
  DataVector& dv_84 = temps.at(81);
  dv_84 = d_15 * dv_56;
  DataVector& dv_85 = temps.at(82);
  dv_85 = dv_22 * dv_69;
  DataVector& dv_86 = temps.at(83);
  dv_86 = dv_2 * dv_47;
  DataVector& dv_87 = temps.at(84);
  dv_87 = dv_49 * dv_54;
  DataVector& dv_88 = temps.at(85);
  dv_88 = dv_86 + dv_87;
  DataVector& dv_89 = temps.at(86);
  dv_89 = dv_88 * dv_88;
  DataVector& dv_90 = temps.at(87);
  dv_90 = dv_2 * dv_2 * dv_2 * dv_2;
  DataVector& dv_91 = temps.at(88);
  dv_91 = d_41 * dv_90;
  DataVector& dv_92 = temps.at(89);
  dv_92 = 3 * dv_16;
  DataVector& dv_93 = temps.at(90);
  dv_93 = -dv_92;
  DataVector& dv_94 = temps.at(91);
  dv_94 = 4 * dv_17;
  DataVector& dv_95 = temps.at(92);
  dv_95 = dv_53 + dv_94;
  DataVector& dv_96 = temps.at(93);
  dv_96 = dv_93 + dv_95;
  DataVector& dv_97 = temps.at(94);
  dv_97 = 3 * dv_17;
  DataVector& dv_98 = temps.at(95);
  dv_98 = -dv_97;
  DataVector& dv_99 = temps.at(96);
  dv_99 = 4 * dv_16;
  DataVector& dv_100 = temps.at(97);
  dv_100 = dv_53 + dv_99;
  DataVector& dv_101 = temps.at(98);
  dv_101 = dv_100 + dv_98;
  DataVector& dv_102 = temps.at(99);
  dv_102 = -14 * Dx * Dy * xp * yp + d_17 * dv_96 + d_18 * dv_101;
  DataVector& dv_103 = temps.at(100);
  dv_103 = M * dv_14;
  DataVector& dv_104 = temps.at(101);
  dv_104 = dv_103 * rp;
  DataVector& dv_105 = temps.at(102);
  dv_105 = d_30 * d_42 * (dv_48 * dv_48);
  DataVector& dv_106 = temps.at(103);
  dv_106 = Dy * dv_44;
  DataVector& dv_107 = temps.at(104);
  dv_107 = 30 * d_43 * dv_0 * dv_106;
  DataVector& dv_108 = temps.at(105);
  dv_108 = Dx * d_44;
  DataVector& dv_109 = temps.at(106);
  dv_109 = dv_1 * dv_108;
  DataVector& dv_110 = temps.at(107);
  dv_110 = 30 * dv_109;
  DataVector& dv_111 = temps.at(108);
  dv_111 = dv_110 * dv_42;
  DataVector& dv_112 = temps.at(109);
  dv_112 = Dx * Dx * Dx * Dx;
  DataVector& dv_113 = temps.at(110);
  dv_113 = 2 * dv_112;
  DataVector& dv_114 = temps.at(111);
  dv_114 = 11 * dv_16;
  DataVector& dv_115 = temps.at(112);
  dv_115 = dv_113 - dv_114 * dv_19 + 2 * (dv_19 * dv_19);
  DataVector& dv_116 = temps.at(113);
  dv_116 = 11 * dv_17;
  DataVector& dv_117 = temps.at(114);
  dv_117 = -dv_53;
  DataVector& dv_118 = temps.at(115);
  dv_118 = dv_116 + dv_117;
  DataVector& dv_119 = temps.at(116);
  dv_119 = Dy * Dy * Dy * Dy;
  DataVector& dv_120 = temps.at(117);
  dv_120 = z * z * z * z;
  DataVector& dv_121 = temps.at(110);
  dv_121 = dv_113 - dv_116 * dv_18 + 2 * dv_119 + 2 * dv_120;
  DataVector& dv_122 = temps.at(118);
  dv_122 = 68 * dv_17;
  DataVector& dv_123 = temps.at(119);
  dv_123 = 7 * dv_18;
  DataVector& dv_124 = temps.at(120);
  dv_124 = dv_122 - dv_123;
  DataVector& dv_125 = temps.at(117);
  dv_125 = 4 * dv_120;
  DataVector& dv_126 = temps.at(121);
  dv_126 = 11 * dv_112 + 11 * dv_119 + dv_123 * dv_17 - dv_125;
  DataVector& dv_127 = temps.at(122);
  dv_127 = -dv_124 * dv_16 + dv_126;
  DataVector& dv_128 = temps.at(123);
  dv_128 = -d_45 * dv_115 - d_46 * (-dv_118 * dv_16 + dv_121) + d_47 * dv_127 +
           dv_107 + dv_111;
  DataVector& dv_129 = temps.at(124);
  dv_129 = -dv_33;
  DataVector& dv_130 = temps.at(125);
  dv_130 = Dy * d_45;
  DataVector& dv_131 = temps.at(126);
  dv_131 = 5 * dv_16;
  DataVector& dv_132 = temps.at(127);
  dv_132 = 5 * dv_18;
  DataVector& dv_133 = temps.at(128);
  dv_133 = dv_131 + dv_132;
  DataVector& dv_134 = temps.at(129);
  dv_134 = dv_133 + dv_34;
  DataVector& dv_135 = temps.at(130);
  dv_135 = -dv_134;
  DataVector& dv_136 = temps.at(131);
  dv_136 = d_46 * dv_79;
  DataVector& dv_137 = temps.at(132);
  dv_137 = 6 * dv_17;
  DataVector& dv_138 = temps.at(133);
  dv_138 = dv_137 + dv_53;
  DataVector& dv_139 = temps.at(134);
  dv_139 = dv_138 + dv_29;
  DataVector& dv_140 = temps.at(135);
  dv_140 = d_35 * dv_108;
  DataVector& dv_141 = temps.at(136);
  dv_141 = -34 * dv_16 + dv_58;
  DataVector& dv_142 = temps.at(137);
  dv_142 = dv_116 + dv_141;
  DataVector& dv_143 = temps.at(138);
  dv_143 = Dy * d_18;
  DataVector& dv_144 = temps.at(139);
  dv_144 = d_17 * dv_143;
  DataVector& dv_145 = temps.at(140);
  dv_145 = -9 * dv_17;
  DataVector& dv_146 = temps.at(141);
  dv_146 = dv_100 + dv_145;
  DataVector& dv_147 = temps.at(142);
  dv_147 = 3 * dv_0;
  DataVector& dv_148 = temps.at(143);
  dv_148 = d_43 * dv_146 * dv_147 + dv_142 * dv_144;
  DataVector& dv_149 = temps.at(144);
  dv_149 = dv_129 * dv_130 - dv_135 * dv_136 + dv_139 * dv_140 + dv_148;
  DataVector& dv_150 = temps.at(145);
  dv_150 = dv_2 * dv_2 * dv_2;
  DataVector& dv_151 = temps.at(146);
  dv_151 = d_40 * dv_150;
  DataVector& dv_152 = temps.at(147);
  dv_152 = d_31 * dv_151;
  DataVector& dv_153 = temps.at(148);
  dv_153 = dv_0 * dv_143;
  DataVector& dv_154 = temps.at(149);
  dv_154 = dv_36 + dv_98;
  DataVector& dv_155 = temps.at(150);
  dv_155 = -dv_131 + dv_31 + dv_94;
  DataVector& dv_156 = temps.at(151);
  dv_156 = d_43 * dv_154 + d_49 * dv_155 + dv_108 * dv_79 - 12 * dv_153;
  DataVector& dv_157 = temps.at(152);
  dv_157 = d_50 * dv_2;
  DataVector& dv_158 = temps.at(153);
  dv_158 = dv_156 * dv_157;
  DataVector& dv_159 = temps.at(154);
  dv_159 = -dv_152 + dv_158;
  DataVector& dv_160 = temps.at(155);
  dv_160 = d_6 * dv_149 + dv_159;
  DataVector& dv_161 = temps.at(156);
  dv_161 = d_51 * dv_151;
  DataVector& dv_162 = temps.at(157);
  dv_162 = Dx * d_43;
  DataVector& dv_163 = temps.at(158);
  dv_163 = dv_162 * dv_79;
  DataVector& dv_164 = temps.at(159);
  dv_164 = Dx * d_17;
  DataVector& dv_165 = temps.at(160);
  dv_165 = 12 * dv_1;
  DataVector& dv_166 = temps.at(161);
  dv_166 = dv_164 * dv_165;
  DataVector& dv_167 = temps.at(162);
  dv_167 = dv_32 + dv_93;
  DataVector& dv_168 = temps.at(163);
  dv_168 = 5 * dv_17;
  DataVector& dv_169 = temps.at(164);
  dv_169 = -dv_168 + dv_31 + dv_99;
  DataVector& dv_170 = temps.at(165);
  dv_170 = -d_44 * dv_167 - d_52 * dv_169 - dv_163 + dv_166;
  DataVector& dv_171 = temps.at(163);
  dv_171 = dv_132 + dv_168;
  DataVector& dv_172 = temps.at(127);
  dv_172 = 2 * Dx;
  DataVector& dv_173 = temps.at(166);
  dv_173 = d_45 * dv_172;
  DataVector& dv_174 = temps.at(167);
  dv_174 = Dx * d_46;
  DataVector& dv_175 = temps.at(168);
  dv_175 = 6 * dv_16;
  DataVector& dv_176 = temps.at(169);
  dv_176 = dv_175 + dv_34;
  DataVector& dv_177 = temps.at(170);
  dv_177 = -dv_176 - dv_53;
  DataVector& dv_178 = temps.at(171);
  dv_178 = d_43 * dv_63;
  DataVector& dv_179 = temps.at(172);
  dv_179 = -9 * dv_16;
  DataVector& dv_180 = temps.at(173);
  dv_180 = dv_179 + dv_95;
  DataVector& dv_181 = temps.at(174);
  dv_181 = -dv_180;
  DataVector& dv_182 = temps.at(175);
  dv_182 = -34 * dv_17 + dv_58;
  DataVector& dv_183 = temps.at(176);
  dv_183 = dv_114 + dv_182;
  DataVector& dv_184 = temps.at(177);
  dv_184 = d_18 * dv_164;
  DataVector& dv_185 = temps.at(178);
  dv_185 = d_53 * dv_1 * dv_181 + dv_173 * (-dv_171 - dv_29) + dv_174 * dv_37 +
           3 * dv_177 * dv_178 - dv_183 * dv_184;
  DataVector& dv_186 = temps.at(179);
  dv_186 = d_6 * dv_185 + dv_157 * dv_170 + dv_161;
  DataVector& dv_187 = temps.at(180);
  dv_187 = d_40 * dv_14;
  DataVector& dv_188 = temps.at(15);
  dv_188 = dv_15 * rp;
  DataVector& dv_189 = temps.at(181);
  dv_189 = d_54 * dv_14 - d_55 * dv_20 - d_56 * dv_20 + dv_187 + dv_188;
  DataVector& dv_190 = temps.at(182);
  dv_190 = d_19 * dv_2;
  DataVector& dv_191 = temps.at(183);
  dv_191 = -dv_186 * xpdot + dv_189 * dv_190;
  DataVector& dv_192 = temps.at(184);
  dv_192 = dv_160 * ypdot + dv_191;
  DataVector& dv_193 = temps.at(185);
  dv_193 = d_1 * d_13 * dv_192;
  DataVector& dv_194 = temps.at(9);
  dv_194 = d_5 * dv_9 - d_57 * dv_3 + dv_7;
  DataVector& dv_195 = temps.at(186);
  dv_195 =
      -d_48 * dv_128 - dv_102 * dv_104 - dv_105 + 4 * dv_193 * dv_194 + dv_91;
  DataVector& dv_196 = temps.at(187);
  dv_196 = dv_195 * dv_25;
  DataVector& dv_197 = temps.at(188);
  dv_197 = -d_14 * dv_55 + dv_21;
  DataVector& dv_198 = temps.at(189);
  dv_198 = dv_197 * rp;
  DataVector& dv_199 = temps.at(80);
  dv_199 = M * dv_58 * dv_59 + d_22 * d_38 * dv_85 + d_25 * d_26 * dv_70 +
           d_29 * dv_58 * dv_71 + d_39 * dv_60 * dv_83 + dv_196 * dv_198 -
           9 * dv_26 * dv_89 - dv_53 * dv_57 - dv_53 * dv_61 + dv_83 * dv_84;
  DataVector& dv_200 = temps.at(190);
  dv_200 = d_21 * dv_199;
  DataVector& dv_201 = temps.at(191);
  dv_201 = M * dv_23;
  DataVector& dv_202 = temps.at(192);
  dv_202 = d_58 * dv_201;
  DataVector& dv_203 = temps.at(4);
  dv_203 = dv_4 * dv_7;
  DataVector& dv_204 = temps.at(193);
  dv_204 = d_61 * dv_12;
  DataVector& dv_205 = temps.at(12);
  dv_205 = d_26 * dv_12;
  DataVector& dv_206 = temps.at(194);
  dv_206 = Dy * ypddot;
  DataVector& dv_207 = temps.at(195);
  dv_207 = Dx * xpddot + dv_206;
  DataVector& dv_208 = temps.at(196);
  dv_208 = d_0 * dv_7;
  DataVector& dv_209 = temps.at(197);
  dv_209 = M * dv_2;
  DataVector& dv_210 = temps.at(198);
  dv_210 = d_75 * dv_7;
  DataVector& dv_211 = temps.at(196);
  dv_211 = d_1 * dv_207 + d_4 * dv_208 + d_68 * dv_209 + d_71 * dv_3 +
           dv_208 * rp + dv_210;
  DataVector& dv_212 = temps.at(199);
  dv_212 = d_14 * dv_11;
  DataVector& dv_213 = temps.at(200);
  dv_213 = Dx + d_5 * d_51 * dv_209;
  DataVector& dv_214 = temps.at(201);
  dv_214 = -dv_212 * (d_76 + d_79) + dv_213;
  DataVector& dv_215 = temps.at(202);
  dv_215 = d_5 * dv_2;
  DataVector& dv_216 = temps.at(203);
  dv_216 = Dy + d_32 * dv_215;
  DataVector& dv_217 = temps.at(199);
  dv_217 = -dv_212 * (d_32 * rp + d_81) + dv_216;
  DataVector& dv_218 = temps.at(204);
  dv_218 = dv_214 * xpdot + dv_217 * ypdot;
  DataVector& dv_219 = temps.at(205);
  dv_219 = -d_59 * (3 * M * dv_14 * rpdot + d_13 * dv_11 * dv_211 * rp -
                    d_74 * dv_205 - dv_203 - dv_204) +
           dv_218;
  DataVector& dv_220 = temps.at(206);
  dv_220 = dv_219 * dv_24;
  DataVector& dv_221 = temps.at(207);
  dv_221 = Dx - xp;
  DataVector& dv_222 = temps.at(208);
  dv_222 = Dy + d_83;
  DataVector& dv_223 = temps.at(209);
  dv_223 = dv_221 * xpdot + dv_222 * ypdot;
  DataVector& dv_224 = temps.at(210);
  dv_224 = dv_223 * rp;
  DataVector& dv_225 = temps.at(211);
  dv_225 = 3 * Dy;
  DataVector& dv_226 = temps.at(212);
  dv_226 = 3 * dv_1;
  DataVector& dv_227 = temps.at(213);
  dv_227 = 3 * Dx;
  DataVector& dv_228 = temps.at(214);
  dv_228 = xpdot * (dv_164 - dv_62 * (d_31 + dv_225) + xp * (dv_226 + dv_33)) +
           ypdot * (dv_143 - dv_63 * (d_51 + dv_227) + yp * (dv_147 + dv_37));
  DataVector& dv_229 = temps.at(195);
  dv_229 = d_70 + dv_207;
  DataVector& dv_230 = temps.at(198);
  dv_230 = -2 * M * d_4 * dv_223 - 2 * M * dv_2 * rpdot + d_0 * dv_224 +
           d_1 * dv_229 + dv_210;
  DataVector& dv_231 = temps.at(215);
  dv_231 = d_84 * dv_3 + dv_230;
  DataVector& dv_232 = temps.at(216);
  dv_232 = d_26 * dv_48;
  DataVector& dv_233 = temps.at(217);
  dv_233 = dv_79 + yp;
  DataVector& dv_234 = temps.at(218);
  dv_234 = 2 * dv_1;
  DataVector& dv_235 = temps.at(219);
  dv_235 = xpdot * (dv_164 - dv_233 * dv_62 + xp * (dv_234 + dv_42)) +
           ypdot * (dv_143 - dv_63 * (dv_172 + xp) + yp * (2 * dv_0 + dv_44));
  DataVector& dv_236 = temps.at(220);
  dv_236 = 2 * d_4 * dv_228;
  DataVector& dv_237 = temps.at(221);
  dv_237 = d_68 * dv_45 + d_84 * dv_47 - 4 * dv_235 * rp + dv_236;
  DataVector& dv_238 = temps.at(222);
  dv_238 = 20 * dv_18;
  DataVector& dv_239 = temps.at(223);
  dv_239 = d_73 * d_88 * dv_55;
  DataVector& dv_240 = temps.at(224);
  dv_240 = 32 * M * dv_239;
  DataVector& dv_241 = temps.at(225);
  dv_241 = d_23 * dv_58;
  DataVector& dv_242 = temps.at(226);
  dv_242 = -dv_231;
  DataVector& dv_243 = temps.at(227);
  dv_243 = d_74 * dv_22;
  DataVector& dv_244 = temps.at(228);
  dv_244 = dv_19 * xp * (xpdot * xpdot * xpdot);
  DataVector& dv_245 = temps.at(229);
  dv_245 = ypdot * (Dy * (-d_17 + d_6) - dv_143 + dv_19 * yp);
  DataVector& dv_246 = temps.at(230);
  dv_246 = Dy * dv_0;
  DataVector& dv_247 = temps.at(231);
  dv_247 = Dy * d_17;
  DataVector& dv_248 = temps.at(232);
  dv_248 = d_9 * dv_164;
  DataVector& dv_249 = temps.at(233);
  dv_249 = d_9 * dv_43;
  DataVector& dv_250 = temps.at(234);
  dv_250 = d_6 * (d_9 + dv_206) + d_92 * dv_6;
  DataVector& dv_251 = temps.at(64);
  dv_251 =
      16 * d_10 * dv_245 + 16 * dv_244 +
      16 * xpdot *
          (Dx * (dv_250 - yp * (Dy + d_91 * yp)) - dv_248 - dv_67 * xpddot +
           xp * (Dy * (Dy + yp) - d_63 * dv_18 + dv_249)) +
      16 * ypdot *
          (Dx * xp * yp - d_90 * dv_16 - d_90 * dv_18 + dv_16 * yp - dv_246 -
           dv_247 + dv_65 * xpddot - dv_68 * ypddot - dv_73);
  DataVector& dv_252 = temps.at(66);
  dv_252 = dv_164 * (d_9 + 2);
  DataVector& dv_253 = temps.at(72);
  dv_253 = dv_34 + dv_35;
  DataVector& dv_254 = temps.at(235);
  dv_254 = dv_225 + yp;
  DataVector& dv_255 = temps.at(233);
  dv_255 =
      -8 * M * d_64 * dv_0 - 8 * d_0 * dv_0 * dv_206 -
      8 * d_10 * (d_0 * (dv_1 + dv_19) - dv_245) + 8 * d_17 * d_67 * dv_31 +
      8 * d_17 * dv_6 - 8 * d_3 * dv_147 - 8 * d_3 * dv_16 + 8 * d_3 * dv_30 -
      8 * d_67 * dv_68 - 8 * d_93 * dv_6 + 8 * d_94 * dv_16 * yp +
      8 * d_94 * dv_18 * yp + 8 * d_95 * dv_31 + 8 * d_95 * dv_35 -
      8 * d_96 * dv_31 - 8 * d_96 * dv_35 - 8 * d_97 * dv_16 -
      8 * d_97 * dv_18 + 8 * dv_147 * dv_6 + 8 * dv_244 - 8 * dv_78 * xpddot +
      8 * xpdot *
          (Dx * (d_36 * dv_233 - d_98 + dv_250 + dv_254 * yp) - dv_252 -
           dv_77 * xpddot +
           xp *
               (d_0 * dv_6 - d_35 * dv_18 * ypddot - dv_226 + dv_249 + dv_253));
  DataVector& dv_256 = temps.at(204);
  dv_256 = -d_59 * (d_100 * dv_205 + d_60 * dv_103 + d_85 * dv_11 * dv_211 -
                    dv_203 - dv_204) +
           dv_218;
  DataVector& dv_257 = temps.at(196);
  dv_257 = dv_22 * dv_22;
  DataVector& dv_258 = temps.at(4);
  dv_258 = 1.0 / dv_257;
  DataVector& dv_259 = temps.at(207);
  dv_259 = -dv_221 * xpdot - dv_222 * ypdot;
  DataVector& dv_260 = temps.at(193);
  dv_260 = d_13 * dv_54;
  DataVector& dv_261 = temps.at(12);
  dv_261 = 18 * dv_88;
  DataVector& dv_262 = temps.at(11);
  dv_262 = dv_195 * dv_197;
  DataVector& dv_263 = temps.at(44);
  dv_263 = -dv_45;
  DataVector& dv_264 = temps.at(124);
  dv_264 = d_17 * dv_129 - d_18 * dv_37 + dv_28;
  DataVector& dv_265 = temps.at(37);
  dv_265 = d_19 * dv_263 - d_4 * dv_264;
  DataVector& dv_266 = temps.at(28);
  dv_266 = dv_0 * dv_1;
  DataVector& dv_267 = temps.at(63);
  dv_267 = -d_17 * dv_96 - d_18 * dv_101 + 14 * dv_266;
  DataVector& dv_268 = temps.at(112);
  dv_268 = d_45 * dv_115 + d_46 * (-dv_118 * dv_16 + dv_121) - d_47 * dv_127 -
           dv_107 - dv_111;
  DataVector& dv_269 = temps.at(110);
  dv_269 = -dv_192;
  DataVector& dv_270 = temps.at(104);
  dv_270 = -dv_194;
  DataVector& dv_271 = temps.at(122);
  dv_271 = dv_269 * dv_270;
  DataVector& dv_272 = temps.at(108);
  dv_272 = -dv_99;
  DataVector& dv_273 = temps.at(115);
  dv_273 = 11 * dv_18;
  DataVector& dv_274 = temps.at(194);
  dv_274 = Dx * d_45;
  DataVector& dv_275 = temps.at(67);
  dv_275 = dv_274 * (dv_116 + dv_272 + dv_273);
  DataVector& dv_276 = temps.at(113);
  dv_276 = dv_100 - dv_116;
  DataVector& dv_277 = temps.at(41);
  dv_277 = Dy * dv_42;
  DataVector& dv_278 = temps.at(118);
  dv_278 = -dv_122 + dv_123;
  DataVector& dv_279 = temps.at(217);
  dv_279 = 22 * dv_16 + dv_278;
  DataVector& dv_280 = temps.at(228);
  dv_280 = d_17 * dv_62;
  DataVector& dv_281 = temps.at(109);
  dv_281 = 4 * dv_112;
  DataVector& dv_282 = temps.at(142);
  dv_282 = 45 * dv_1;
  DataVector& dv_283 = temps.at(65);
  dv_283 = 22 * dv_17;
  DataVector& dv_284 = temps.at(234);
  dv_284 = 15 * dv_1;
  DataVector& dv_285 = temps.at(229);
  dv_285 = -dv_114 + dv_95;
  DataVector& dv_286 = temps.at(236);
  dv_286 = dv_18 + dv_29 + dv_97;
  DataVector& dv_287 = temps.at(237);
  dv_287 = dv_43 + dv_98;
  DataVector& dv_288 = temps.at(238);
  dv_288 = 15 * dv_0;
  DataVector& dv_289 = temps.at(239);
  dv_289 = Dy * Dy * Dy;
  DataVector& dv_290 = temps.at(240);
  dv_290 = 68 * dv_1;
  DataVector& dv_291 = temps.at(241);
  dv_291 = dv_2 * dv_223;
  DataVector& dv_292 = temps.at(145);
  dv_292 = d_41 * dv_150;
  DataVector& dv_293 = temps.at(180);
  dv_293 = 6 * dv_187;
  DataVector& dv_294 = temps.at(242);
  dv_294 = dv_259 * dv_293;
  DataVector& dv_295 = temps.at(243);
  dv_295 = dv_209 * rpdot;
  DataVector& dv_296 = temps.at(244);
  dv_296 = M * dv_224;
  DataVector& dv_297 = temps.at(245);
  dv_297 = Dx * d_18;
  DataVector& dv_298 = temps.at(246);
  dv_298 = 12 * dv_17;
  DataVector& dv_299 = temps.at(247);
  dv_299 = dv_179 + dv_298 + dv_53;
  DataVector& dv_300 = temps.at(170);
  dv_300 = dv_177 * dv_225;
  DataVector& dv_301 = temps.at(248);
  dv_301 = yp * (dv_175 + dv_53 + dv_98);
  DataVector& dv_302 = temps.at(163);
  dv_302 = 2 * d_45 * (dv_171 + dv_93);
  DataVector& dv_303 = temps.at(169);
  dv_303 = -dv_176 - dv_31;
  DataVector& dv_304 = temps.at(175);
  dv_304 = 33 * dv_16 + dv_182;
  DataVector& dv_305 = temps.at(249);
  dv_305 = Dx * ypdot;
  DataVector& dv_306 = temps.at(250);
  dv_306 = Dy * xpdot;
  DataVector& dv_307 = temps.at(251);
  dv_307 = -dv_94;
  DataVector& dv_308 = temps.at(252);
  dv_308 = dv_137 + dv_31;
  DataVector& dv_309 = temps.at(29);
  dv_309 = dv_29 + dv_308;
  DataVector& dv_310 = temps.at(253);
  dv_310 = 4 * Dy;
  DataVector& dv_311 = temps.at(133);
  dv_311 = dv_138 + dv_93;
  DataVector& dv_312 = temps.at(254);
  dv_312 = 12 * dv_16;
  DataVector& dv_313 = temps.at(255);
  dv_313 = 12 * dv_18;
  DataVector& dv_314 = temps.at(128);
  dv_314 = dv_133 + dv_98;
  DataVector& dv_315 = temps.at(136);
  dv_315 = dv_141 + 33 * dv_17;
  DataVector& dv_316 = temps.at(137);
  dv_316 = dv_142 * dv_79;
  DataVector& dv_317 = temps.at(95);
  dv_317 = dv_145 + dv_312 + dv_53;
  DataVector& dv_318 = temps.at(256);
  dv_318 = (1.0 / 24.0) * dv_258;
  DataVector& dv_319 = temps.at(212);
  dv_319 = -d_18 * dv_172 + dv_164 + dv_226 * xp;
  DataVector& dv_320 = temps.at(257);
  dv_320 = 24 * dv_23;
  DataVector& dv_321 = temps.at(197);
  dv_321 = d_16 * dv_209;
  DataVector& dv_322 = temps.at(258);
  dv_322 = dv_320 * dv_321;
  DataVector& dv_323 = temps.at(259);
  dv_323 = dv_38 * xp;
  DataVector& dv_324 = temps.at(260);
  dv_324 = d_108 * dv_201;
  DataVector& dv_325 = temps.at(232);
  dv_325 = -d_6 * dv_6 * xpdot + d_91 * dv_297 - dv_1 * xp + dv_248;
  DataVector& dv_326 = temps.at(261);
  dv_326 = pow(dv_22, 3.0 / 2.0);
  DataVector& dv_327 = temps.at(262);
  dv_327 = d_110 * dv_326;
  DataVector& dv_328 = temps.at(66);
  dv_328 = Dx * d_98 + d_37 * dv_63 + dv_252 - yp * (-d_0 * dv_306 + dv_62) -
           ypdot * (d_6 * dv_306 + d_7 * dv_62);
  DataVector& dv_329 = temps.at(263);
  dv_329 = d_111 * d_13 * dv_326;
  DataVector& dv_330 = temps.at(264);
  dv_330 = d_112 * dv_40;
  DataVector& dv_331 = temps.at(265);
  dv_331 = d_5 * dv_201 * dv_241;
  DataVector& dv_332 = temps.at(266);
  dv_332 = d_23 * dv_330;
  DataVector& dv_333 = temps.at(267);
  dv_333 = d_109 * d_59 * dv_18 * dv_23;
  DataVector& dv_334 = temps.at(268);
  dv_334 = dv_40 * dv_40;
  DataVector& dv_335 = temps.at(269);
  dv_335 = d_23 * dv_334;
  DataVector& dv_336 = temps.at(270);
  dv_336 = d_113 * dv_23;
  DataVector& dv_337 = temps.at(271);
  dv_337 = -dv_297;
  DataVector& dv_338 = temps.at(272);
  dv_338 = d_19 * (d_51 * dv_1 + dv_164 + dv_337) - d_4 * dv_319;
  DataVector& dv_339 = temps.at(273);
  dv_339 = d_114 * dv_320;
  DataVector& dv_340 = temps.at(274);
  dv_340 = 8 * d_16 * dv_201;
  DataVector& dv_341 = temps.at(275);
  dv_341 = d_112 * dv_49;
  DataVector& dv_342 = temps.at(276);
  dv_342 = d_13 * dv_214;
  DataVector& dv_343 = temps.at(277);
  dv_343 = dv_201 * dv_53;
  DataVector& dv_344 = temps.at(278);
  dv_344 = d_115 * d_40 * dv_58;
  DataVector& dv_345 = temps.at(279);
  dv_345 = dv_23 * dv_342;
  DataVector& dv_346 = temps.at(280);
  dv_346 = dv_214 * dv_24;
  DataVector& dv_347 = temps.at(281);
  dv_347 = d_116 * dv_346;
  DataVector& dv_348 = temps.at(282);
  dv_348 = dv_23 * dv_69;
  DataVector& dv_349 = temps.at(283);
  dv_349 = d_113 * dv_348;
  DataVector& dv_350 = temps.at(284);
  dv_350 = d_5 * dv_313;
  DataVector& dv_351 = temps.at(280);
  dv_351 = dv_335 * dv_346;
  DataVector& dv_352 = temps.at(285);
  dv_352 = M * dv_351;
  DataVector& dv_353 = temps.at(286);
  dv_353 = 24 * dv_18;
  DataVector& dv_354 = temps.at(287);
  dv_354 = d_40 * dv_353;
  DataVector& dv_355 = temps.at(282);
  dv_355 = d_117 * dv_348;
  DataVector& dv_356 = temps.at(72);
  dv_356 = d_9 * dv_76 + dv_253;
  DataVector& dv_357 = temps.at(73);
  dv_357 = -dv_74 * yp + dv_81;
  DataVector& dv_358 = temps.at(78);
  dv_358 = d_17 * dv_356 + dv_357;
  DataVector& dv_359 = temps.at(75);
  dv_359 = dv_340 * dv_358;
  DataVector& dv_360 = temps.at(288);
  dv_360 = d_57 * d_7 * dv_358;
  DataVector& dv_361 = temps.at(261);
  dv_361 = 1.0 / dv_326;
  DataVector& dv_362 = temps.at(289);
  dv_362 = dv_214 * dv_361;
  DataVector& dv_363 = temps.at(290);
  dv_363 = d_119 * dv_24;
  DataVector& dv_364 = temps.at(87);
  dv_364 = 4 * d_1 * d_13 * dv_194 *
               (dv_191 + ypdot * (d_6 * (-dv_130 * dv_33 + dv_134 * dv_136 +
                                         dv_139 * dv_140 + dv_148) +
                                  dv_159)) +
           2 * d_40 * dv_90 - d_48 * dv_128 - dv_102 * dv_104 - dv_105;
  DataVector& dv_365 = temps.at(268);
  dv_365 = d_121 * (-d_14 * dv_334 + dv_21);
  DataVector& dv_366 = temps.at(21);
  dv_366 = -dv_276;
  DataVector& dv_367 = temps.at(135);
  dv_367 = M * d_42 * dv_49;
  DataVector& dv_368 = temps.at(154);
  dv_368 = 48 * d_1 * dv_14;
  DataVector& dv_369 = temps.at(180);
  dv_369 = d_4 * dv_293;
  DataVector& dv_370 = temps.at(129);
  dv_370 = dv_1 * dv_172 + xp * (dv_19 + dv_92);
  DataVector& dv_371 = temps.at(131);
  dv_371 = dv_174 * dv_6;
  DataVector& dv_372 = temps.at(99);
  dv_372 = dv_164 * (d_125 - 17 * d_3);
  DataVector& dv_373 = temps.at(102);
  dv_373 = M * dv_175;
  DataVector& dv_374 = temps.at(252);
  dv_374 = -15 * dv_16 + dv_308;
  DataVector& dv_375 = temps.at(183);
  dv_375 = M * dv_137;
  DataVector& dv_376 = temps.at(168);
  dv_376 = -15 * dv_17 + dv_175 + dv_31;
  DataVector& dv_377 = temps.at(230);
  dv_377 = d_43 * dv_246;
  DataVector& dv_378 = temps.at(134);
  dv_378 = d_7 * dv_24;
  DataVector& dv_379 = temps.at(143);
  dv_379 = -d_17 * dv_79 + d_35 * dv_0 + dv_143;
  DataVector& dv_380 = temps.at(9);
  dv_380 = dv_47 * yp;
  DataVector& dv_381 = temps.at(123);
  dv_381 = dv_0 * yp;
  DataVector& dv_382 = temps.at(291);
  dv_382 = -d_10 * dv_143 + d_33 * dv_5 - dv_247 * (d_10 + 1) + dv_381;
  DataVector& dv_383 = temps.at(292);
  dv_383 = Dx * d_6 * xpdot * ypdot - d_11 * dv_247 -
           xp * (Dx * d_35 + d_0 * dv_305 - d_7 * dv_306) -
           yp * (d_0 * dv_5 + d_10 * dv_1 + dv_234);
  DataVector& dv_384 = temps.at(293);
  dv_384 = d_126 * dv_54;
  DataVector& dv_385 = temps.at(294);
  dv_385 = d_23 * dv_384;
  DataVector& dv_386 = temps.at(231);
  dv_386 = -dv_247;
  DataVector& dv_387 = temps.at(295);
  dv_387 = d_19 * (d_31 * dv_0 + dv_143 + dv_386) - d_4 * dv_379;
  DataVector& dv_388 = temps.at(296);
  dv_388 = d_126 * dv_49;
  DataVector& dv_389 = temps.at(297);
  dv_389 = d_13 * dv_217;
  DataVector& dv_390 = temps.at(298);
  dv_390 = dv_23 * dv_389;
  DataVector& dv_391 = temps.at(83);
  dv_391 = d_116 * dv_86;
  DataVector& dv_392 = temps.at(299);
  dv_392 = dv_217 * dv_24;
  DataVector& dv_393 = temps.at(300);
  dv_393 = M * dv_56;
  DataVector& dv_394 = temps.at(301);
  dv_394 = dv_350 * dv_393;
  DataVector& dv_395 = temps.at(302);
  dv_395 = d_59 * dv_56;
  DataVector& dv_396 = temps.at(303);
  dv_396 = dv_354 * dv_395;
  DataVector& dv_397 = temps.at(69);
  dv_397 = d_118 * dv_70;
  DataVector& dv_398 = temps.at(304);
  dv_398 = d_116 * dv_87;
  DataVector& dv_399 = temps.at(300);
  dv_399 = d_108 * dv_393;
  DataVector& dv_400 = temps.at(305);
  dv_400 = d_120 * dv_89;
  DataVector& dv_401 = temps.at(261);
  dv_401 = dv_217 * dv_361;
  DataVector& dv_402 = temps.at(306);
  dv_402 = 2 * dv_260;
  DataVector& dv_403 = temps.at(307);
  dv_403 = d_122 * dv_195;
  DataVector& dv_404 = temps.at(308);
  dv_404 = 5 * d_121 * dv_262;
  DataVector& dv_405 = temps.at(309);
  dv_405 = d_44 * dv_62;
  DataVector& dv_406 = temps.at(94);
  dv_406 = dv_80 + yp * (dv_43 + dv_97);
  DataVector& dv_407 = temps.at(42);
  dv_407 = 12 * dv_6;
  DataVector& dv_408 = temps.at(77);
  dv_408 = 4 * Dx;
  DataVector& dv_409 = temps.at(310);
  dv_409 = 48 * d_127 * dv_22;
  DataVector& dv_410 = temps.at(311);
  dv_410 = d_13 * dv_257;
  DataVector& dv_411 = temps.at(312);
  dv_411 = dv_22 * dv_56;
  DataVector& dv_412 = temps.at(313);
  dv_412 = d_25 * dv_411;

  get(get<CurvedScalarWave::Tags::Psi>(*result)) =
      dv_24 * (-1.0 / 24.0 * dv_200 * dv_26 - 1.0 / 2.0 * dv_26 * dv_52 + 1);
  get(get<::Tags::dt<CurvedScalarWave::Tags::Psi>>(*result)) =
      -dv_318 *
      (M * d_87 * dv_23 *
           (16 * M * d_13 * d_20 * dv_18 * dv_256 +
            64 * M * d_13 * d_22 * dv_22 * dv_69 * rpdot +
            80 * M * d_13 * d_29 * dv_18 * dv_22 * rpdot +
            16 * M * d_13 * d_39 * dv_256 * dv_69 +
            16 * M * d_15 * d_23 * dv_242 * dv_54 * dv_69 +
            16 * M * d_22 * d_23 * dv_18 * dv_242 * dv_54 +
            32 * M * d_23 * d_24 * dv_18 * dv_55 * rpdot +
            16 * M * d_23 * d_26 * dv_55 * dv_69 * rpdot - M * dv_251 * dv_84 +
            18 * M * dv_256 * dv_258 * dv_89 * rp -
            d_102 * dv_196 *
                (3 * M * d_59 * dv_14 * rpdot +
                 d_13 * d_5 * dv_54 * (2 * M * d_101 * dv_2 - dv_230) +
                 d_23 * d_99 * dv_55 - d_5 * dv_259 * dv_3 -
                 d_59 * d_61 * dv_55 - dv_7) -
            d_102 * dv_256 * dv_258 * dv_262 +
            8 * d_13 * d_29 * dv_256 * dv_82 +
            36 * d_13 * d_39 * dv_22 * dv_82 * rpdot -
            d_13 * d_86 * dv_256 * dv_58 +
            12 * d_15 * d_23 * dv_55 * dv_82 * rpdot - d_20 * dv_18 * dv_240 +
            8 * d_23 * d_24 * dv_242 * dv_54 * dv_82 +
            8 * d_23 * d_73 * dv_18 * dv_22 * pow(rp, 14) -
            16 * d_29 * dv_239 * dv_82 - d_29 * dv_255 * dv_60 -
            d_39 * dv_240 * dv_69 - d_39 * dv_241 * dv_242 * dv_54 -
            d_39 * dv_251 * dv_71 + 16 * d_73 * d_86 * d_88 * dv_18 * dv_55 -
            d_74 * d_86 * d_89 * dv_85 - d_89 * dv_18 * dv_243 * pow(rp, 13) -
            44 * dv_18 * dv_61 * rpdot +
            4 * dv_195 * dv_197 * dv_25 * rp * rpdot -
            dv_198 * dv_25 *
                (6 * M * d_1 * d_13 * dv_265 *
                     (d_101 * dv_264 + d_68 * dv_263 + 4 * dv_235 * rp -
                      dv_236) +
                 6 * M * d_23 * d_24 * d_99 * (dv_265 * dv_265) -
                 8 * d_100 * d_22 * dv_271 - 4 * d_103 * d_13 * dv_271 -
                 d_104 * dv_269 *
                     (2 * M * d_101 * d_5 * dv_2 +
                      6 * M * d_4 * d_59 * dv_2 * rpdot +
                      2 * M * d_57 * dv_259 - d_0 * d_4 * d_5 * dv_259 -
                      d_7 * dv_215 * rpdot - dv_229) -
                 d_104 * dv_270 *
                     (-d_68 * dv_189 * dv_2 - dv_160 * ypddot +
                      dv_186 * xpddot + 2 * dv_189 * dv_223 * rp -
                      dv_190 *
                          (-M * d_75 * dv_20 + 2 * M * dv_14 * rpdot -
                           d_0 * dv_8 - 24 * d_103 * dv_20 - 12 * d_26 * dv_7 -
                           d_41 * dv_291 - 16 * d_6 * dv_291 -
                           d_7 * dv_2 * dv_224 + 16 * dv_14 * rp * rpdot) -
                      xpdot *
                          (d_6 * (xpdot *
                                      (d_43 * (dv_300 + dv_303 * yp) -
                                       d_44 * dv_172 *
                                           (27 * dv_1 + 20 * dv_17 + dv_238 +
                                            dv_272) +
                                       d_49 * (-9 * Dy * dv_180 + dv_304 * yp) +
                                       d_93 * dv_0 * (18 * Dy * yp - dv_183) +
                                       dv_302) +
                                  ypdot *
                                      (20 * Dx * dv_130 +
                                       d_107 * (dv_300 + dv_301) +
                                       d_43 * dv_172 * (dv_1 + dv_100 + dv_72) +
                                       d_53 * (Dy * dv_181 + dv_299 * yp) -
                                       2 * dv_280 * (34 * dv_1 + dv_183))) -
                           d_92 * dv_185 +
                           dv_157 *
                               (d_106 * (-d_27 * (5 * dv_1 + dv_169) +
                                         d_44 * dv_79 + 6 * dv_164 * dv_222 +
                                         dv_297 * (-d_83 - dv_225)) +
                                xpdot * (8 * Dx * dv_254 * xp * yp -
                                         3 * d_17 * (dv_167 + dv_41) +
                                         d_18 * (2 * Dy * yp - dv_169) -
                                         6 * dv_108)) -
                           dv_170 * dv_295 + dv_170 * dv_296 + dv_292 * xpdot -
                           dv_294 * xp) -
                      ypdot *
                          (d_6 *
                               (d_43 * (-dv_5 * (-20 * dv_1 - 27 * dv_17 +
                                                 dv_312 + dv_313) +
                                        2 * ypdot *
                                            (dv_135 * dv_310 + dv_314 * yp)) +
                                d_44 *
                                    (3 * dv_305 *
                                         (dv_117 - dv_137 + dv_16 + dv_165) +
                                     xpdot * (d_35 * dv_311 + dv_310 * dv_33)) +
                                d_45 * (-dv_309 * ypdot + dv_5 * dv_79) -
                                d_49 * (dv_5 * (54 * dv_17 + dv_179 +
                                                36 * dv_18 + dv_290) +
                                        ypdot * (-dv_315 * yp + dv_316)) -
                                d_52 * (9 * dv_305 * (dv_146 + dv_27) +
                                        xpdot * (-d_35 * dv_317 + dv_316))) +
                           d_92 * dv_149 + dv_156 * dv_295 - dv_156 * dv_296 +
                           dv_157 * (-d_17 * (2 * dv_5 * (dv_225 + 5 * yp) +
                                              ypdot * (-8 * dv_1 - dv_131 -
                                                       dv_307 - dv_75)) +
                                     d_18 * (4 * dv_254 * dv_5 -
                                             3 * ypdot * (dv_154 + dv_234)) -
                                     d_31 * xp *
                                         (6 * dv_305 * (-d_83 - dv_79) +
                                          xpdot * (dv_155 + dv_27)) +
                                     2 * d_44 * (dv_305 + dv_306)) +
                           dv_292 * ypdot - dv_294 * yp)) -
                 rp * (d_54 *
                           (-xpdot *
                                (d_44 *
                                     (-dv_16 * (22 * dv_18 + dv_282 + dv_283) +
                                      dv_19 * (dv_284 + dv_95) + dv_281) +
                                 d_52 * (15 * Dy * yp * (-dv_17 + dv_18) -
                                         dv_126 + dv_16 * (dv_124 + dv_282)) -
                                 dv_162 * (15 * dv_106 + dv_276 * yp) + dv_275 +
                                 dv_280 * (-45 * dv_277 + dv_279 * yp)) +
                            ypdot *
                                (-d_18 * dv_288 *
                                     (-dv_225 * dv_44 + dv_287 * yp) -
                                 d_43 * (-d_105 * dv_289 + dv_1 * dv_273 +
                                         4 * dv_119 + dv_125 +
                                         dv_16 * (11 * dv_1 - dv_283 + dv_58) -
                                         dv_18 * dv_283 + dv_281) +
                                 d_49 * (-dv_1 * dv_123 + dv_126 +
                                         dv_16 * (dv_278 + dv_290) -
                                         22 * dv_289 * yp) +
                                 15 * dv_108 * (dv_277 - dv_286 * yp) +
                                 dv_130 * dv_285)) +
                       dv_103 * dv_267 * rpdot - 8 * dv_151 * dv_223 +
                       dv_188 *
                           (xpdot * (3 * dv_164 - dv_62 * (7 * Dy + d_105) +
                                     xp * (7 * dv_1 + dv_96)) +
                            ypdot * (d_18 * dv_225 - dv_63 * (7 * Dx + 4 * xp) +
                                     yp * (7 * dv_0 + dv_101))) -
                       dv_224 * dv_267 * dv_3 + 8 * dv_268 * rp * rpdot) -
                 rpdot * (d_48 * dv_268 + dv_104 * dv_267 + dv_91)) -
            dv_238 * dv_57 * rpdot - 8 * dv_243 * dv_82 * pow(rp, 12) -
            dv_255 * dv_59 -
            dv_26 * dv_261 *
                (d_13 * dv_242 * dv_48 * rp - 2 * d_74 * dv_232 * dv_54 +
                 2 * dv_2 * dv_228 * rp + 2 * dv_2 * dv_47 * rpdot -
                 dv_224 * dv_47 + dv_237 * dv_260 * rp - dv_87 * rpdot)) -
       3 * M * dv_200 * dv_220 - d_82 * dv_220 * dv_52 -
       14 * d_87 * dv_199 * dv_201 * rpdot - 72 * dv_202 * dv_51 * rpdot -
       12 * dv_202 *
           (2 * d_23 * dv_232 * dv_40 * (d_65 * d_66 + d_72) + d_68 * dv_50 -
            d_85 * dv_237 * dv_40 + dv_190 * dv_228 + dv_224 * dv_38 -
            dv_231 * dv_49 * rp - dv_39 * rpdot) +
       24 * dv_219 * dv_23);
  get<0>(get<::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>,
                           Frame::Inertial>>(*result)) =
      -dv_318 *
      (-5 * M * dv_362 * dv_364 * dv_365 -
       d_108 * dv_352 * (d_17 * dv_356 + dv_357) - d_118 * dv_351 * dv_69 +
       d_120 * dv_362 * (dv_51 * dv_51) +
       d_122 * dv_24 * dv_364 * (-d_14 * dv_330 + dv_213) -
       d_13 * dv_338 * dv_339 * dv_40 - d_59 * dv_351 * dv_354 +
       dv_214 * dv_320 - dv_214 * dv_355 - dv_319 * dv_322 + dv_323 * dv_324 +
       dv_324 * dv_341 + dv_325 * dv_327 + dv_325 * dv_335 * dv_336 +
       dv_328 * dv_329 + dv_328 * dv_335 * dv_340 - dv_330 * dv_331 +
       dv_332 * dv_333 + dv_332 * dv_349 + dv_332 * dv_359 + dv_342 * dv_343 -
       dv_344 * dv_345 - dv_345 * dv_360 - dv_347 * dv_39 - dv_347 * dv_50 +
       dv_350 * dv_352 +
       18 * dv_363 * dv_51 *
           (2 * d_13 * dv_338 * dv_40 + 2 * dv_2 * dv_319 - dv_323 - dv_341) -
       dv_365 * dv_378 *
           (M * dv_2 * rp * (d_44 * dv_167 + d_52 * dv_169 + dv_163 - dv_166) +
            d_1 * d_13 * dv_192 * (2 * M * d_57 * xp - d_5 * d_78 - xpdot) +
            2 * d_6 *
                (d_44 * dv_284 * (dv_19 + dv_93) + dv_174 * dv_366 +
                 15 * dv_178 * (dv_18 + dv_34 + dv_92) + dv_184 * dv_279 +
                 dv_275) -
            dv_161 -
            dv_260 * (d_123 * dv_370 + d_124 * dv_370 -
                      d_50 * (d_44 * (d_3 * dv_374 + dv_373) +
                              d_52 * (d_3 * dv_376 + dv_375) + dv_234 * dv_372 +
                              4 * dv_274 * dv_6 + 4 * dv_371 +
                              xpdot * (d_45 * (dv_179 + dv_32) + d_46 * dv_30 +
                                       d_47 * (-17 * dv_17 + dv_31 + dv_312) -
                                       dv_110 + 12 * dv_377)) -
                      d_6 * (d_107 * (M * dv_94 + d_3 * dv_317) +
                             d_53 * (M * dv_99 + d_3 * dv_311) + dv_173 * dv_6 +
                             20 * dv_371 + dv_372 * dv_41 +
                             xpdot * (d_46 * dv_303 + d_47 * dv_304 -
                                      54 * dv_109 + dv_302 + 36 * dv_377)) -
                      dv_368 * xp + dv_369 * xp) -
            dv_338 * dv_367));
  get<1>(get<::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>,
                           Frame::Inertial>>(*result)) =
      dv_318 *
      (M * dv_401 * dv_404 -
       d_121 * dv_197 * dv_378 *
           (d_102 *
                (Dy * d_45 * dv_285 - Dy * d_46 * (dv_114 + dv_273 + dv_307) -
                 d_43 * dv_287 * dv_288 -
                 dv_144 * (dv_123 - 68 * dv_16 + dv_283) -
                 15 * dv_286 * dv_405) +
            dv_152 - dv_158 -
            dv_193 * (2 * M * d_57 * yp - d_5 * d_80 - ypdot) +
            dv_260 *
                (d_123 * dv_406 + d_124 * dv_406 -
                 d_50 *
                     (d_43 * (d_3 * (dv_145 + dv_36) + dv_375) +
                      d_45 * dv_35 * ypdot +
                      d_49 * (d_3 * (-17 * dv_16 + dv_298 + dv_31) + dv_373) +
                      6 * dv_153 * (d_0 - 5 * d_3) + dv_405 * dv_407 +
                      xpdot * (d_43 * dv_376 * xp + d_44 * dv_374 * yp +
                               dv_130 * dv_408 - 34 * dv_143 * dv_164 +
                               dv_174 * dv_310)) -
                 d_6 *
                     (d_107 * (Dx * dv_79 * (-9 * d_3 + d_7) + dv_301 * xpdot) +
                      d_35 * d_44 * (Dx * dv_407 + dv_299 * xpdot) +
                      2 * d_43 *
                          (Dy * (Dy * d_125 + dv_5 * yp) + d_3 * dv_314) +
                      d_45 * (20 * Dy * dv_5 - dv_309 * ypdot) +
                      d_49 * (d_3 * dv_315 +
                              dv_408 * (M * dv_227 - 17 * dv_1 * xpdot))) -
                 dv_368 * yp + dv_369 * yp) +
            dv_367 * dv_387) -
       dv_217 * dv_320 + dv_217 * dv_355 -
       dv_24 * dv_403 * (-d_126 * d_14 * dv_54 + dv_216) -
       dv_260 * dv_339 * dv_387 +
       dv_261 * dv_363 *
           (2 * dv_2 * dv_379 + dv_380 - dv_387 * dv_402 + dv_388) +
       dv_322 * dv_379 + dv_324 * dv_380 + dv_324 * dv_388 + dv_327 * dv_382 +
       dv_329 * dv_383 + dv_331 * dv_384 - dv_333 * dv_385 +
       dv_336 * dv_382 * dv_56 + dv_340 * dv_383 * dv_56 - dv_343 * dv_389 +
       dv_344 * dv_390 - dv_349 * dv_385 - dv_359 * dv_385 + dv_360 * dv_390 -
       dv_391 * dv_392 - dv_392 * dv_394 + dv_392 * dv_396 + dv_392 * dv_397 -
       dv_392 * dv_398 + dv_392 * dv_399 * dv_82 - dv_400 * dv_401);
  get<2>(get<::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>,
                           Frame::Inertial>>(*result)) =
      -1.0 / 24.0 * z *
      (d_109 * d_115 * dv_410 + d_109 * dv_22 * dv_395 + d_110 * d_8 * dv_257 -
       d_111 * d_128 * dv_410 + d_113 * d_8 * dv_411 +
       d_114 * d_129 * dv_260 * dv_409 - d_117 * dv_85 +
       72 * d_119 * d_127 * dv_88 * (-d_130 * dv_260 + dv_2) -
       d_128 * d_16 * dv_412 -
       d_25 * dv_197 *
           (3 * d_127 * d_13 * d_130 * d_50 * (-d_4 * dv_38 + dv_46) +
            d_127 *
                (dv_103 + rp * (-d_17 * dv_285 + d_18 * dv_366 + 30 * dv_266)) -
            dv_402 * (-d_127 * d_77 * dv_2 -
                      d_127 * rp *
                          (xpdot * (5 * dv_164 + dv_27 * xp + dv_337) +
                           ypdot * (5 * dv_143 + 6 * dv_381 + dv_386)) +
                      d_55 * dv_2 + d_56 * dv_2)) *
           1.0 / d_39 -
       d_38 * dv_257 + 48 * d_5 * dv_103 - d_5 * dv_412 - 24 * dv_13 +
       24 * dv_16 + 24 * dv_17 + dv_25 * dv_400 - dv_26 * dv_404 +
       dv_321 * dv_409 - dv_344 * dv_60 + dv_353 - dv_358 * dv_399 -
       dv_360 * dv_60 + dv_391 + dv_394 - dv_396 - dv_397 + dv_398 + dv_403 +
       dv_53 * dv_71) /
      pow(dv_22, 5.0 / 2.0);
}
}  // namespace CurvedScalarWave::Worldtube
