// Distributed under the MIT License.
// See LICENSE.txt for details.
#pragma once
#include <optional>
#include <string>
#include <tuple>
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/SchwarzschildReference.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Options/Auto.hpp"
#include "Options/Context.hpp"
#include "Options/ParseError.hpp"
#include "Options/String.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/PupStlCpp17.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::worldtube {
/// Diagnostic replay on a static, conforming spherical-shell face at R=2.5.
/// Record the missing inner neighbor, including multistep startup evaluations.
struct FaceReplayParameters {
  struct File { using type=std::string; static constexpr Options::String help{"Face tape filename; must not exist when recording."}; };
  struct Mode { using type=std::string; static constexpr Options::String help{"Record or Replay."}; };
  struct ModelSector { using type=std::string; static constexpr Options::String help{"None (all reference), Gauge, CP, Physical, CPAndPhysical (reference gauge), or All (all model)."}; };
  struct Stride { using type=size_t; static constexpr Options::String help{"Retain every nth normal record; startup is always exact. Stride 1 requires exact sample times; larger strides use local degree-five interpolation."}; static size_t lower_bound(){return 1;} };
  struct SchwarzschildReference { using type=Options::Auto<BoundaryConditions::detail::SchwarzschildReferenceParameters,Options::AutoLabel::None>; static constexpr Options::String help{"Gauge model parameters, or None for Algebraic zero-rate radiation."}; };
  using options=tmpl::list<File,Mode,ModelSector,Stride,SchwarzschildReference>;
  static constexpr Options::String help{"Diagnostic sector replay; fixed static spherical-shell mesh and global AB stepping only."};
  FaceReplayParameters()=default;
  FaceReplayParameters(std::string file_in,std::string mode,std::string sector,size_t stride_in,
      std::optional<BoundaryConditions::detail::SchwarzschildReferenceParameters> reference_in,const Options::Context& context={})
      :file(std::move(file_in)),record(mode=="Record"),model_sector(std::move(sector)),stride(stride_in),reference(std::move(reference_in)) {
    if(file.empty() or (mode!="Record" and mode!="Replay") or stride==0 or
       (model_sector!="None" and model_sector!="Gauge" and model_sector!="CP" and model_sector!="Physical" and model_sector!="CPAndPhysical" and model_sector!="All") or
       (record and (stride!=1 or model_sector!="None"))) {
      PARSE_ERROR(context,"Invalid FaceReplay file, mode, sector or stride. Record requires sector None and stride 1.");
    }
  }
  std::string file{};
  bool record{false};
  std::string model_sector{"None"};
  size_t stride{1};
  std::optional<BoundaryConditions::detail::SchwarzschildReferenceParameters> reference{};
  void pup(PUP::er& p){p|file;p|record;p|model_sector;p|stride;p|reference;}
};
inline bool operator==(const FaceReplayParameters& a,const FaceReplayParameters& b){return std::tie(a.file,a.record,a.model_sector,a.stride,a.reference)==std::tie(b.file,b.record,b.model_sector,b.stride,b.reference);}
inline bool operator!=(const FaceReplayParameters& a,const FaceReplayParameters& b){return not(a==b);}
struct ReferenceReplayGauge {
  struct ReferenceReplay { using type=FaceReplayParameters; static constexpr Options::String help{"Reference interface replay with selected model sectors. Applies to the complete boundary, not just gauge."}; };
  using options=tmpl::list<ReferenceReplay>;
  static constexpr Options::String help{"Diagnostic complete-boundary replay."};
  ReferenceReplayGauge()=default;
  explicit ReferenceReplayGauge(FaceReplayParameters p):parameters(std::move(p)){}
  FaceReplayParameters parameters{};
};

template<size_t Dim> struct FaceReplayData {
  tnsr::aa<DataVector,Dim> metric{},pi{};
  tnsr::iaa<DataVector,Dim> phi{};
  tnsr::i<DataVector,Dim> raw_normal{};
  size_t radial_points{0};
  void pup(PUP::er& p){p|metric;p|pi;p|phi;p|raw_normal;p|radial_points;}
};
template<size_t Dim> bool operator==(const FaceReplayData<Dim>& a,const FaceReplayData<Dim>& b){return std::tie(a.metric,a.pi,a.phi,a.raw_normal,a.radial_points)==std::tie(b.metric,b.pi,b.phi,b.raw_normal,b.radial_points);}
template<size_t Dim> bool operator!=(const FaceReplayData<Dim>& a,const FaceReplayData<Dim>& b){return not(a==b);}

void update_face_replay(gsl::not_null<std::optional<FaceReplayData<3>>*> data,
    const FaceReplayParameters& options,const tnsr::aa<DataVector,3>& metric,
    const tnsr::aa<DataVector,3>& pi,const tnsr::iaa<DataVector,3>& phi,
    const tnsr::I<DataVector,3>& coordinates,const Mesh<3>& mesh,
    const InverseJacobian<DataVector,3,Frame::ElementLogical,Frame::Inertial>& inverse_jacobian,
    const TimeStepId& time_id,double time);
} // namespace gh::worldtube
