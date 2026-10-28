# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Generate reference fixtures by executing arithmetic from a local SpEC tree.

Run with python-spectre GenerateReference.py /path/to/SpEC --compiler clang++.
No SpEC installation is needed by the unit tests. Original source fragments
are compiled in a temporary directory, never copied into this repository.
The compatibility classes below supply containers and profiling stubs, not
the numerical algorithms. Even the linear fit is extracted from SpEC.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile


def between(source, start, end):
    return source[source.index(start) : source.index(end, source.index(start))]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("spec_root", type=Path)
    parser.add_argument("--compiler", default="c++")
    args = parser.parse_args()
    root = args.spec_root.resolve()
    files = [
        "Spectral/SpectralFunctionals/ComputeGridDiagnostics.cpp",
        "Utils/Math/PolynomialFit.cpp",
        "Evolution/FoshSystem/DualFrameSystem/NormOfProjectedDerivConstraint.cpp",
        "ComputeItems/NumberOfFilteredModes.cpp",
        "Evolution/AmrDrivers/ChangeSpectralGrid.cpp",
        "Utils/LowLevelUtils/SimpleNorms.cpp",
        "Spectral/BasisFunctions/TensorSphericalHarmonics/YlmStorageIterator.hpp",
        "Spectral/BasisFunctions/TensorSphericalHarmonics/YlmStorageIterator.cpp",
        "Spectral/BasisFunctions/TensorSphericalHarmonics/TensorYlmStorageIterator.hpp",
        "Spectral/BasisFunctions/TensorSphericalHarmonics/TensorYlmStorageIterator.cpp",
        "Spectral/BasisFunctions/TensorSphericalHarmonics/TensorYlmCoefficient.cpp",
    ]
    sources = {name: (root / name).read_text() for name in files}
    spectrum = between(
        sources[files[0]],
        "namespace {\n//===",
        "//========================================================================="
        "\n\nComputeDetailedGridDiagnostics::",
    )
    fit = between(
        sources[files[1]], "  void fit(", "  // Fit to N-1th order polynomial"
    )
    l2_norm = between(
        sources[files[5]],
        "double L2Norm(",
        "//-------------------------------------------------------------------------",
    )
    grouping = between(sources[files[2]], "    const int Ntopol=", "\n  }\n}")
    i1_sum = between(
        sources[files[0]][sources[files[0]].index("::GetPowerMonitor(") :],
        "      case I1:",
        "      case S1:",
    )
    # Execute the actual permission expressions, with neighbor inhibition and
    # cooldown disabled, to isolate the shared local p-controller rules.
    driver = sources[files[4]]
    drop = between(
        driver,
        "      mAllowDropPoints[k][bf] =\n          MaxNext",
        "\n      // If this is first call",
    )
    add = between(
        driver,
        "      mAllowAddPoints[k][bf] =\n          (NeighborPileUpModes",
        "\n      // We (Mark and Nils)",
    )
    forced_add = between(
        driver,
        "        mAllowAddPoints[k][bf] =\n            (mTopErrorTol",
        "\n        // We (Mark and Nils)",
    )
    filter_count = between(
        sources[files[3]],
        "      double nf = 0;",
        "      result.insert({key, nf});",
    )

    # Execute SpEC's actual storage iterator and spin weights. Compatibility
    # stubs only supply its mesh-registration base and IPoint container.
    iterator_source = "\n".join(
        "\n".join(
            line
            for line in sources[name].splitlines()
            if not line.startswith(("#include", "#pragma once"))
        )
        for name in files[6:10]
    )
    spin_weight = between(
        sources[files[10]],
        "int TensorYlmCoefSpinWeight(",
        "IPoint TensorYlmCoefAstar(",
    )

    cases = [
        (
            "smooth",
            [
                10 ** (-0.45 * i) * (1 + 0.04 * math.sin(0.7 * i))
                for i in range(12)
            ],
            1,
        ),
        (
            "unfiltered",
            [10 ** (-0.35 * i) * (1 + 0.07 * math.sin(i)) for i in range(12)],
            0,
        ),
        (
            "three_filtered",
            [10 ** (-0.4 * i) * (1 + 0.05 * math.sin(i)) for i in range(14)],
            3,
        ),
        (
            "pile_up",
            [
                1.0,
                0.3,
                0.09,
                0.027,
                0.0081,
                0.00243,
                0.0008,
                0.0008,
                0.0008,
                0.0008,
                0.0008,
                0.0008,
            ],
            1,
        ),
        ("zero", [0.0] * 12, 1),
        ("flat", [1.0] * 12, 1),
        (
            "parity_gaps",
            [10 ** (-0.4 * i) if i % 2 == 0 else 0.0 for i in range(13)],
            1,
        ),
        ("leading_zero", [0.0] + [10 ** (-0.4 * i) for i in range(1, 12)], 1),
        ("resolved_polynomial", [1.0, 0.2, 0.03] + [0.0] * 9, 1),
        (
            "roundoff_tail",
            [
                1.0,
                0.1,
                0.01,
                0.001,
                0.0001,
                1.0e-5,
                1.0e-6,
                1.0e-7,
                1.0e-8,
                1.0e-16,
                1.0e-17,
                1.0e-18,
            ],
            1,
        ),
    ]
    norm_cases = [
        ("three_I1", [1.0e-4, 4.0e-4, 9.0e-4], [1, 1, 1]),
        ("I1_S2", [1.0e-4, 4.0e-4, 9.0e-4], [1, 2]),
    ]
    permission_cases = [
        # current log error, next log error, constraint/target, pile-up, force
        ("add", -3.7, -3.4, 0.01, 0.0, True),
        ("drop", -5.0, -4.7, 0.01, 0.0, True),
        ("constraint_blocks_drop", -5.0, -4.7, 1.0, 0.0, True),
        ("pile_up_blocks_add", -3.7, -3.4, 0.01, 2.0, True),
        ("unforced_pile_up", -3.7, -3.4, 0.01, 2.0, False),
    ]
    # Two independent tensor components, four transverse nodal points, and
    # four spectral coefficients along the monitored direction.
    coefficient_components = [
        [
            1.0,
            2.0,
            3.0,
            4.0,
            2.0,
            0.0,
            -1.0,
            0.5,
            1.0,
            0.5,
            2.0,
            0.0,
            0.0,
            1.0,
            0.5,
            2.0,
        ],
        [
            0.5,
            1.0,
            0.0,
            2.0,
            -0.5,
            3.0,
            1.0,
            0.0,
            2.0,
            1.0,
            0.5,
            1.0,
            1.0,
            0.0,
            2.0,
            0.5,
        ],
    ]
    preamble = r"""
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <list>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;
#define ASSERT(test, message) do { if (!(test)) throw runtime_error("reference assertion"); } while(false)
#define REQUIRE(test, message) ASSERT(test, message)
namespace MV {
struct Size { explicit Size(int n):n(n){} int n; };
template<class C> struct Seq { const C& c; };
template<class C> Seq<C> seq(const C& c){return {c};}
}
template<class T> struct MyVector : vector<T> {
  using vector<T>::vector;
  MyVector() = default;
  explicit MyVector(MV::Size n):vector<T>(n.n){}
  MyVector(MV::Size n,T value):vector<T>(n.n,value){}
  template<class C> explicit MyVector(MV::Seq<C> s):vector<T>(s.c.begin(),s.c.end()){}
  template<class C> void assign(MV::Seq<C> s){vector<T>::assign(s.c.begin(),s.c.end());}
  void assign(MV::Size n,T value){vector<T>::assign(n.n,value);}
  int Size()const{return this->size();}
};
template<class T> using MyList = list<T>;
template<class T> T sqr(T x){return x*x;}
int IRound(double x){return static_cast<int>(floor(x+0.5));}
struct SimpleProfiler { explicit SimpleProfiler(const char*){} };
struct SaverOfMeshInformation {
 explicit SaverOfMeshInformation(const char*){}
 void Unregister()const{}
};
using IPoint=MyVector<int>;
"""
    classes = r"""
class PolynomialFit {
 public:
  PolynomialFit(const MyVector<double>& x,const MyVector<double>& y,int order)
      :coef(MV::Size(2)),err(MV::Size(2)) {
    ASSERT(order==1 && x.Size()>2, "linear reference only");
    double chi,q;
    fit(x.data(),y.data(),x.Size(),&coef[0],&coef[1],&err[0],&err[1],&chi,&q);
  }
  const auto& Coefs()const{return coef;}
  const auto& ErrorsInCoefs()const{return err;}
 private: MyVector<double> coef,err;
};
class ComputeGridDiagnostics {
 public:
  ComputeGridDiagnostics(const MyVector<double>&,double);
  double operator()(const string& k)const{return mDiag.at(k);}
 private:
  MyVector<double> mAllCoefs;
  int mNrFilteredModes;
  unordered_map<string,double> mDiag;
};
struct Topology {string name; int dim; string MyClassID()const{return name;} int Dim()const{return dim;}};
struct Subdomain {
  MyVector<Topology*> tops;
  const auto& IrreducibleTopologies()const{return tops;}
  vector<int> OffsetIrreducibleTopologies()const{
    vector<int> result;int n=0;for(auto t:tops){result.push_back(n);n+=t->dim;}return result;
  }
};
struct ResponseFilter {
  vector<double> response;
  const ResponseFilter& FilterFunction(int)const{return *this;}
  void Filter(double& value,int j,int)const{value*=response[j];}
};
struct StrideIterator {
  int position=0;
  void Reset(){position=0;}
  explicit operator bool()const{return position<4;}
  void operator++(){++position;}
  int Offset()const{return 4*position;}
  int Stride()const{return 1;}
};
"""

    def numbers(values):
        return ",".join(format(x, ".17g") for x in values)

    main_code = ["int main(){cout<<setprecision(17);"]
    for name, power, filtered in cases:
        main_code.append(
            f"""{{ComputeGridDiagnostics d(MyVector<double>{{{numbers(power)}}},{filtered});
          cout<<"S {name} "<<d("TruncationError")<<" "<<d("PredictedTruncationErrorForOneLessMode")
          <<" "<<d("ConvergenceFactor")<<" "<<d("NumberOfPiledUpModes")<<" "<<d("SpectrumIsDegenerate")<<"\\n";}}"""
        )
    for name, norms, dims in norm_cases:
        top_defs = ";".join(
            f'Topology t{i}{{"{"I1" if d==1 else "S2"}",{d}}}'
            for i, d in enumerate(dims)
        )
        top_ptrs = ",".join(f"&t{i}" for i in range(len(dims)))
        main_code.append(f"""{{{top_defs}; Subdomain sd{{{{{top_ptrs}}}}};
          MyVector<double> norms{{{numbers(norms)}}},mResult;
          for(auto& value:norms){{double nodal[3]={{value,value,value}};value=L2Norm(nodal,3);}}
          {grouping}
          cout<<"C {name}";for(auto x:mResult)cout<<" "<<x;cout<<"\\n";}}""")
    for name, error, next_error, constraint, pile, force in permission_cases:
        main_code.append(
            f"""{{int k=0,bf=0;
          array<array<bool,1>,1> mAllowDropPoints{{}},mAllowAddPoints{{}};
          array<array<double,1>,1> MaxNextLog10TruncError{{{{{{{next_error}}}}}}},
            MaxNextTopError{{{{{{{constraint}}}}}}},MaxTopError{{{{{{{constraint}}}}}}},
            MaxLog10TruncError{{{{{{{error}}}}}}},MinPileUpModes{{{{{{{pile}}}}}}};
          array<double,1> mLog10TruncErrorMax{{{{-4.}}}};
          double cBufferForLog10NextTruncError=.1,cBufferForNextTopError=1.25,
            mTopErrorTol=1.,mMaxAllowedPileUpModes=3.,cBufferForPileUpModes=1.5;
          bool NeighborPileUpModes=false;
          {drop}
          {add}
          if({str(force).lower()}){{{forced_add}}}
          cout<<"P {name} "<<mAllowAddPoints[0][0]<<" "<<mAllowDropPoints[0][0]<<"\\n";}}"""
        )
    main_code.append(f"""{{int bf=0,N=8;
      MyVector<ResponseFilter> sf{{ResponseFilter{{{{1.,1.,1.,1.,.9,.5,0.,0.}}}}}};
      {filter_count}
      cout<<"F filter "<<nf<<" "<<IRound(nf)<<"\\n";}}""")
    main_code.append(f"""{{
      MyVector<double> SumOfSquares(MV::Size(4));
      MyVector<int> Counts(MV::Size(4));
      vector<vector<double>> components{{{{{numbers(coefficient_components[0])}}},{{{numbers(coefficient_components[1])}}}}};
      StrideIterator iterator;auto* mpStrideIter=&iterator;
      int mNumberOfCoeffs=4;constexpr int I1=0;
      for(const auto& Coefficients:components){{switch(I1){{{i1_sum}}}}}
      cout<<"A I1";for(int k=0;k<4;++k)cout<<" "<<sqrt(SumOfSquares[k]/Counts[k]);cout<<"\\n";
    }}""")
    main_code.append(r"""{
      vector<IPoint> metric{{},{0},{1},{2},{0,0},{0,1},{0,2},{1,1},{1,2},{2,2}};
      vector<IPoint> phi;
      for(int i=0;i<3;++i)for(auto component:metric){
        component.insert(component.begin(),i);phi.push_back(component);
      }
      for(const auto& group:vector<vector<IPoint>>{metric,metric,phi}){
        vector<int> counts(9,0);
        for(const auto& component:group){
          TensorYlmStorageIterator it(TensorYlmCoefSpinWeight(component),9,17);
          for(it.Reset();it;++it)++counts[it.l()];
        }
        cout<<"Y angular";for(auto n:counts)cout<<" "<<n;cout<<"\n";
      }
    }""")
    main_code.append("}")
    cpp = (
        preamble
        + fit
        + l2_norm
        + classes
        + iterator_source
        + spin_weight
        + spectrum
        + "\n"
        + "\n".join(main_code)
    )
    with tempfile.TemporaryDirectory(prefix="spec-amr-reference-") as temp:
        source = Path(temp) / "reference.cpp"
        source.write_text(cpp)
        exe = Path(temp) / "reference"
        subprocess.run(
            [args.compiler, "-std=c++17", "-O0", str(source), "-o", str(exe)],
            check=True,
        )
        lines = subprocess.check_output([str(exe)], text=True).splitlines()
    angular_counts = [
        list(map(int, line.split()[2:]))
        for line in lines
        if line.startswith("Y ")
    ]
    results = {line.split()[1]: line.split()[2:] for line in lines}
    revision = subprocess.check_output(
        ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
    ).strip()
    metadata = {
        "spec_revision": revision,
        "source_sha256": {
            n: hashlib.sha256(s.encode()).hexdigest()
            for n, s in sources.items()
        },
        "source_execution": (
            "Unmodified arithmetic fragments; original linear fit;"
            " container/profiler compatibility only."
        ),
        "spectra": [
            {
                "name": n,
                "power": p,
                "filtered_modes": f,
                "log10_error": float(results[n][0]),
                "log10_next_error": float(results[n][1]),
                "convergence": float(results[n][2]),
                "pile_up": float(results[n][3]),
                "degenerate": bool(int(results[n][4])),
            }
            for n, p, f in cases
        ],
        "constraints": [
            {
                "name": n,
                "directional_rms": ns,
                "topology_dimensions": ds,
                "spec_aggregated": list(map(float, results[n])),
            }
            for n, ns, ds in norm_cases
        ],
        "permissions": [
            {
                "name": n,
                "log10_error": e,
                "log10_next_error": ne,
                "constraint_ratio": c,
                "pile_up": p,
                "force": f,
                "allow_add": bool(int(results[n][0])),
                "allow_drop": bool(int(results[n][1])),
            }
            for n, e, ne, c, p, f in permission_cases
        ],
        "filter_example": {
            "response": [1.0, 1.0, 1.0, 1.0, 0.9, 0.5, 0.0, 0.0],
            "effective_filtered_modes": float(results["filter"][0]),
            "rounded_filtered_modes": int(results["filter"][1]),
        },
        "angular_coefficient_counts": {
            "l_max": 8,
            "radial_points": 1,
            "spacetime_metric": angular_counts[0],
            "pi": angular_counts[1],
            "phi": angular_counts[2],
        },
        "power_aggregation": {
            "components": coefficient_components,
            "spec_power": list(map(float, results["I1"])),
        },
    }
    out = Path(__file__).parent
    (out / "ReferenceData.json").write_text(
        json.dumps(metadata, indent=2, allow_nan=False) + "\n"
    )
    # The generated header avoids a JSON dependency in the C++ unit tests.
    header = [
        "// Distributed under the MIT License.",
        "// See LICENSE.txt for details.",
        (
            "// Generated by GenerateReference.py. See ReferenceData.json for"
            " provenance."
        ),
        "#pragma once",
        "#include <array>",
        "#include <cstddef>",
        "#include <vector>",
        "namespace spec_amr_reference {",
        (
            "struct SpectrumCase {const char* name; std::vector<double> power;"
            " size_t filtered_modes; std::array<double,4> expected; bool"
            " degenerate;};"
        ),
        "inline const std::vector<SpectrumCase> spectra{",
    ]
    for x in metadata["spectra"]:
        expected = [
            x[k]
            for k in [
                "log10_error",
                "log10_next_error",
                "convergence",
                "pile_up",
            ]
        ]
        header.append(
            f'{{"{x["name"]}",{{{numbers(x["power"])}}},{x["filtered_modes"]},{{{{{numbers(expected)}}}}},{str(x["degenerate"]).lower()}}},'
        )
    header += [
        "};",
        (
            "struct ConstraintCase {const char* name; std::array<double,3>"
            " directional_rms; std::vector<double> spec_aggregated;};"
        ),
        "inline const std::vector<ConstraintCase> constraints{",
    ]
    for x in metadata["constraints"]:
        header.append(
            f'{{"{x["name"]}",{{{{{numbers(x["directional_rms"])}}}}},{{{numbers(x["spec_aggregated"])}}}}},'
        )
    header += [
        "};",
        (
            "struct PermissionCase {double error, next_error, constraint_ratio,"
            " pile_up; bool force, allow_add, allow_drop;};"
        ),
        "inline const std::vector<PermissionCase> permissions{",
    ]
    for x in metadata["permissions"]:
        header.append(
            f'{{{numbers([10**x["log10_error"],10**x["log10_next_error"],x["constraint_ratio"],x["pile_up"]])},{str(x["force"]).lower()},{str(x["allow_add"]).lower()},{str(x["allow_drop"]).lower()}}},'
        )
    header += [
        "};",
        "inline const std::array<std::array<double,16>,2>"
        " coefficient_components{{"
        + ",".join("{{" + numbers(c) + "}}" for c in coefficient_components)
        + "}};",
        "inline const std::array<double,4> i1_power{{"
        + numbers(metadata["power_aggregation"]["spec_power"])
        + "}};",
        "inline const std::array<std::array<size_t,9>,3> angular_counts{{"
        + ",".join(
            "{{" + ",".join(map(str, counts)) + "}}"
            for counts in angular_counts
        )
        + "}};",
        "}  // namespace spec_amr_reference",
        "",
    ]
    (out / "ReferenceData.hpp").write_text("\n".join(header))
    print(f"Generated fixtures from SpEC {revision}; source hashes recorded.")


if __name__ == "__main__":
    main()
