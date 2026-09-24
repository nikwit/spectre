// Distributed under the MIT License.
// See LICENSE.txt for details.
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/FaceReplay.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <mutex>
#include <vector>
#include "Utilities/ErrorHandling/Error.hpp"

namespace gh::worldtube {
namespace {
constexpr std::array<char,8> magic{{'G','H','R','P','L','Y','0','1'}};
template<class T> void put(std::ostream& s,const T& value){s.write(reinterpret_cast<const char*>(&value),sizeof(T));if(not s)ERROR("Face replay write failed");}
template<class T> void get_value(std::istream& s,T& value){s.read(reinterpret_cast<char*>(&value),sizeof(T));if(not s)ERROR("Truncated face replay file");}
struct Key { int64_t slab; uint64_t substep; double step_time,time; };
void write_key(std::ostream& s,const Key& k){put(s,k.slab);put(s,k.substep);put(s,k.step_time);put(s,k.time);}
Key read_key(std::istream& s){Key k{};get_value(s,k.slab);get_value(s,k.substep);get_value(s,k.step_time);get_value(s,k.time);return k;}
bool near(double a,double b){return std::abs(a-b)<2.e-11;}
struct Writer {
  std::ofstream stream;
  size_t points;
  Writer(const std::string& file,const std::vector<double>& xyz):points(xyz.size()/3){
    if(std::filesystem::exists(file))ERROR("Refusing to overwrite face replay tape "<<file);
    stream.open(file,std::ios::binary);stream.write(magic.data(),8);put(stream,static_cast<uint64_t>(points));
    stream.write(reinterpret_cast<const char*>(xyz.data()),static_cast<std::streamsize>(xyz.size()*sizeof(double)));
  }
  void append(const Key& k,const std::vector<double>& values){
    if(values.size()!=50*points)ERROR("Face replay resolution changed");
    write_key(stream,k);stream.write(reinterpret_cast<const char*>(values.data()),static_cast<std::streamsize>(values.size()*sizeof(double)));
    if(not stream)ERROR("Face replay recording failed");
    // Flush each record so an interrupted run leaves a detectable complete prefix.
    stream.flush();
  }
};
struct Record { Key key; std::streamoff offset; };
struct Reader {
  std::ifstream stream;
  size_t points;
  std::vector<Record> startup{},main{};
  Reader(const std::string& file,const std::vector<double>& xyz,size_t stride){
    stream.open(file,std::ios::binary);std::array<char,8> version{};stream.read(version.data(),8);
    if(version!=magic)ERROR("Invalid face replay format: "<<file);
    uint64_t nf=0;get_value(stream,nf);points=nf;
    if(3*points!=xyz.size())ERROR("Face replay angular resolution mismatch");
    for(const auto expected:xyz){double observed=0.;get_value(stream,observed);if(std::abs(observed-expected)>2.e-13)ERROR("Face replay coordinate/order mismatch");}
    const auto end=static_cast<std::streamoff>(std::filesystem::file_size(file));
    size_t index=0;Record last{};
    while(static_cast<std::streamoff>(stream.tellg())<end){
      const Key k=read_key(stream);const auto offset=static_cast<std::streamoff>(stream.tellg());
      if(k.substep!=0)ERROR("Face replay supports multistep evaluations without RK substeps only");
      if(k.slab<0)startup.push_back({k,offset});
      else {
        if(index%stride==0)main.push_back({k,offset});
        last={k,offset};++index;
      }
      stream.seekg(static_cast<std::streamoff>(50*points*sizeof(double)),std::ios::cur);
      if(not stream or static_cast<std::streamoff>(stream.tellg())>end)ERROR("Truncated face replay record");
    }
    if(main.empty())ERROR("Face replay has no normal evolution samples");
    if(main.back().offset!=last.offset)main.push_back(last);
    for(size_t i=1;i<main.size();++i)if(main[i].key.time<=main[i-1].key.time)ERROR("Nonmonotonic normal face replay times");
  }
  std::vector<double> read(const Record& r){
    std::vector<double> v(50*points);stream.clear();stream.seekg(r.offset);
    stream.read(reinterpret_cast<char*>(v.data()),static_cast<std::streamsize>(v.size()*sizeof(double)));
    if(not stream)ERROR("Face replay data read failed");return v;
  }
  std::vector<double> evaluate(const Key& k,size_t stride){
    if(k.slab<0){
      for(const auto& r:startup)if(r.key.slab==k.slab and r.key.substep==k.substep and near(r.key.step_time,k.step_time) and near(r.key.time,k.time))return read(r);
      ERROR("Missing face replay self-start stage at "<<k.time<<", slab "<<k.slab);
    }
    auto it=std::lower_bound(main.begin(),main.end(),k.time,[](const Record& r,double t){return r.key.time<t;});
    if(it!=main.end() and near(it->key.time,k.time))return read(*it);
    if(it!=main.begin() and near((it-1)->key.time,k.time))return read(*(it-1));
    if(stride==1)ERROR("Exact face replay has no sample at "<<k.time);
    if(k.time<main.front().key.time or k.time>main.back().key.time or main.size()<6)ERROR("Face replay interpolation would extrapolate or has insufficient data");
    const auto j=static_cast<size_t>(it-main.begin());
    const size_t lo=std::min(j>3?j-3:0,main.size()-6);
    const auto base=read(main[lo]);auto result=base;
    for(size_t i=0;i<6;++i){
      double w=1.;for(size_t m=0;m<6;++m)if(i!=m)w*=(k.time-main[lo+m].key.time)/(main[lo+i].key.time-main[lo+m].key.time);
      const auto v=read(main[lo+i]);for(size_t c=0;c<v.size();++c)result[c]+=w*(v[c]-base[c]);
    }
    return result;
  }
};
std::mutex tape_mutex;
std::map<std::string,std::unique_ptr<Writer>> writers;
std::map<std::pair<std::string,size_t>,std::unique_ptr<Reader>> readers;
template<class T> void append_face(std::vector<double>& v,const T& field,size_t nr,size_t offset,size_t nf){for(const auto& component:field)for(size_t p=0;p<nf;++p)v.push_back(component[p*nr+offset]);}
template<class T> void unpack(T& result,const std::vector<double>& v,size_t& offset,size_t nf){result=T(nf,0.);for(auto& component:result){std::copy_n(v.begin()+static_cast<std::ptrdiff_t>(offset),nf,component.begin());offset+=nf;}}
} // namespace

void update_face_replay(const gsl::not_null<std::optional<FaceReplayData<3>>*> data,
    const FaceReplayParameters& options,const tnsr::aa<DataVector,3>& metric,
    const tnsr::aa<DataVector,3>& pi,const tnsr::iaa<DataVector,3>& phi,
    const tnsr::I<DataVector,3>& coordinates,const Mesh<3>& mesh,
    const InverseJacobian<DataVector,3,Frame::ElementLogical,Frame::Inertial>& inverse_jacobian,
    const TimeStepId& id,const double time){
  const size_t nr=mesh.extents(0),nf=mesh.extents(1)*mesh.extents(2),off=options.record?nr-1:0;
  if(mesh.quadrature(0)!=Spectral::Quadrature::GaussLobatto or mesh.basis(1)!=Spectral::Basis::SphericalHarmonic)ERROR("Face replay needs radial LGL and spherical harmonics");
  std::vector<double> xyz;xyz.reserve(3*nf);append_face(xyz,coordinates,nr,off,nf);
  for(size_t p=0;p<nf;++p){double r2=0.;for(size_t i=0;i<3;++i)r2+=xyz[i*nf+p]*xyz[i*nf+p];if(std::abs(std::sqrt(r2)-2.5)>2.e-12)ERROR("Diagnostic face replay requires coordinate radius 2.5M");}
  const Key key{id.slab_number(),id.substep(),id.step_time().value(),time};
  const std::lock_guard<std::mutex> lock(tape_mutex);
  if(options.record){
    auto& w=writers[options.file];if(not w)w=std::make_unique<Writer>(options.file,xyz);
    std::vector<double> v;v.reserve(50*nf);append_face(v,metric,nr,off,nf);append_face(v,pi,nr,off,nf);append_face(v,phi,nr,off,nf);w->append(key,v);data->reset();return;
  }
  auto& reader=readers[{options.file,options.stride}];if(not reader)reader=std::make_unique<Reader>(options.file,xyz,options.stride);
  const auto v=reader->evaluate(key,options.stride);FaceReplayData<3> face{};size_t pos=0;
  unpack(face.metric,v,pos,nf);unpack(face.pi,v,pos,nf);unpack(face.phi,v,pos,nf);
  face.raw_normal=tnsr::i<DataVector,3>(nf,0.);face.radial_points=nr;
  for(size_t i=0;i<3;++i)for(size_t p=0;p<nf;++p)face.raw_normal.get(i)[p]=-inverse_jacobian.get(0,i)[p*nr];
  *data=std::move(face);
}
} // namespace gh::worldtube
