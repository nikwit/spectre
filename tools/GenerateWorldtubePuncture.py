#!/usr/bin/env python3

# Distributed under the MIT License.
# See LICENSE.txt for details.

"""Generate the equatorial Schwarzschild puncture in an instantaneous frame.

Input: the five trusted SymPy expression caches produced by parse_all.py in
https://github.com/nikwit/puncture-field (cpp-port). These must contain the
already differentiated, generic expressions, including the exact inertial
time derivative. Specializing the field before differentiation is incorrect.

Example:
  python GenerateWorldtubePuncture.py --expressions-dir /path/to/cpp-port \
      --output /path/to/PunctureFieldOrder2.cpp
Then run clang-format on the output. Requires SymPy 1.14.
"""

import argparse
import hashlib
import json
import pickle
from collections import deque
from pathlib import Path

import sympy as sp
from sympy.printing.c import C99CodePrinter


class Printer(C99CodePrinter):
    def _print_Integer(self, expr):
        return f"{int(expr)}.0"


def generate(expressions_dir, output):
    names = ["PhiP2", "dtPhiP2", "dxPhiP2", "dyPhiP2", "dzPhiP2"]
    xp, yp, rp, vr, vx = sp.symbols("xp yp rp rpdot xpdot")
    expressions = []
    hashes = {}
    for name in names:
        path = expressions_dir / (name + ".pkl")
        hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()
        with path.open("rb") as stream:
            expr = pickle.load(stream)
        expressions.append(expr.xreplace({xp: rp, yp: sp.S.Zero, vr: vx}))
    replacements, outputs = sp.cse(expressions, optimizations="basic")
    scalar_assignments, vector_assignments = [], []
    vector_symbols = set(sp.symbols("Dx Dy z"))
    scalar_cache, rename = {}, {}
    printer = Printer()

    def is_vector(expr):
        return bool(expr.free_symbols & vector_symbols)

    def scalar(expr):
        if expr.is_Atom:
            return expr
        if expr not in scalar_cache:
            name = sp.Symbol(f"d_{len(scalar_assignments)}")
            scalar_cache[expr] = name
            scalar_assignments.append((name, expr))
        return scalar_cache[expr]

    def vector(expr):
        name = sp.Symbol(f"v_{len(vector_assignments)}")
        vector_symbols.add(name)
        vector_assignments.append((name, expr))
        return name

    def cost(expr):
        return 1 + sum(cost(a) for a in expr.args)

    def bounded(expr):
        if not is_vector(expr):
            return scalar(expr)
        if expr.is_Atom:
            return expr
        args = [bounded(a) for a in expr.args]
        # Bound Blaze expression-template size without changing mathematical
        # order. Materialize partial sums/products, then reuse their buffers
        # after the last consumer. Hoisted scalar subtrees cost one node.
        if isinstance(expr, (sp.Add, sp.Mul)):
            while 1 + sum(cost(a) for a in args) > 32:
                group, remaining = [], list(args)
                while (
                    remaining
                    and sum(cost(a) for a in group) + cost(remaining[0]) < 31
                ):
                    group.append(remaining.pop(0))
                if len(group) < 2:
                    largest = max(range(len(args)), key=lambda i: cost(args[i]))
                    assert cost(args[largest]) > 1
                    args[largest] = vector(args[largest])
                else:
                    args = [
                        vector(expr.func(*group, evaluate=False))
                    ] + remaining
            return expr.func(*args, evaluate=False)
        args = [vector(a) if cost(a) > 16 else a for a in args]
        return expr.func(*args)

    for name, expr in replacements:
        expr = expr.xreplace(rename)
        rename[name] = (
            vector(bounded(expr)) if is_vector(expr) else scalar(expr)
        )
    output_exprs = [bounded(expr.xreplace(rename)) for expr in outputs]
    records = vector_assignments + list(enumerate(output_exprs))
    last_use = {}
    for i, (_, expr) in enumerate(records):
        for symbol in expr.free_symbols & vector_symbols:
            last_use[symbol] = i
    free, slots, next_slot = deque(), {}, 0
    lines = []
    output_lhs = [
        "get(get<CurvedScalarWave::Tags::Psi>(*result))",
        "get(get<::Tags::dt<CurvedScalarWave::Tags::Psi>>(*result))",
    ] + [
        f"get<{i}>(get<::Tags::deriv<CurvedScalarWave::Tags::Psi, "
        "tmpl::size_t<3>, Frame::Inertial>>(*result))"
        for i in range(3)
    ]
    for i, (name, expr) in enumerate(records):
        dying = sorted(
            (s for s in expr.free_symbols if s in slots and last_use[s] == i),
            key=str,
        )
        # All assignments are elementwise, so a dying input can be the output.
        free.extend(slots[s] for s in dying)
        if isinstance(name, int):
            lhs = output_lhs[name]
        else:
            if not free:
                free.append(next_slot)
                next_slot += 1
            slots[name] = free.popleft()
            lines.append(f"  DataVector& {name} = temps.at({slots[name]});")
            lhs = str(name)
        lines.append(f"  {lhs} = {printer.doprint(expr)};")

    scalar_lines = "\n".join(
        f"  const double {s} = {printer.doprint(e)};"
        for s, e in scalar_assignments
    )
    body = "\n".join(lines)
    text = f"""// Distributed under the MIT License.
// See LICENSE.txt for details.

// Generated by tools/GenerateWorldtubePuncture.py with SymPy {sp.__version__}.
// Specialize AFTER inertial differentiation: xp=rp, yp=0, rpdot=xpdot.
// Velocity and acceleration retain both equatorial components.
#include "Evolution/Systems/CurvedScalarWave/Worldtube/PunctureField.hpp"

#include <cmath>

#include "Evolution/Systems/CurvedScalarWave/Worldtube/PunctureFieldWorkspace.hpp"
#include "Utilities/EqualWithinRoundoff.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"

namespace CurvedScalarWave::Worldtube {{

// NOLINTNEXTLINE(google-readability-function-size, readability-function-size)
void puncture_field_2(
    const gsl::not_null<Variables<tmpl::list<
        CurvedScalarWave::Tags::Psi, ::Tags::dt<CurvedScalarWave::Tags::Psi>,
        ::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>,
                      Frame::Inertial>>>*> result,
    const tnsr::I<DataVector, 3, Frame::Inertial>& centered_coords,
    const tnsr::I<double, 3>& particle_position,
    const tnsr::I<double, 3>& particle_velocity,
    const tnsr::I<double, 3>& particle_acceleration, const double bh_mass,
    PunctureFieldWorkspace* workspace) {{
  const size_t grid_size = get<0>(centered_coords).size();
  result->initialize(grid_size);
  if (grid_size == 0) {{
    return;
  }}
  ASSERT(equal_within_roundoff(get<2>(particle_position), 0.),
         "The second-order puncture requires an equatorial particle position.");
  const double rp = hypot(get<0>(particle_position), get<1>(particle_position));
  const double cos_phi = get<0>(particle_position) / rp;
  const double sin_phi = get<1>(particle_position) / rp;
  const double xpdot = cos_phi * get<0>(particle_velocity) +
                       sin_phi * get<1>(particle_velocity);
  const double ypdot = -sin_phi * get<0>(particle_velocity) +
                       cos_phi * get<1>(particle_velocity);
  const double xpddot = cos_phi * get<0>(particle_acceleration) +
                        sin_phi * get<1>(particle_acceleration);
  const double ypddot = -sin_phi * get<0>(particle_acceleration) +
                        cos_phi * get<1>(particle_acceleration);
  const double M = bh_mass;
  PunctureFieldWorkspace local_workspace{{}};
  auto& temps = (workspace == nullptr ? local_workspace : *workspace)
                    .buffer({next_slot + 3}, grid_size);
  auto& Dx = temps.at({next_slot});
  auto& Dy = temps.at({next_slot + 1});
  const auto& z = get<2>(centered_coords);
  Dx = cos_phi * get<0>(centered_coords) + sin_phi * get<1>(centered_coords);
  Dy = -sin_phi * get<0>(centered_coords) + cos_phi * get<1>(centered_coords);

{scalar_lines}

{body}

  // The rotation is fixed for this evaluation. The inertial time derivative
  // above is unchanged; only the spatial covector is rotated back.
  auto& gradient = get<::Tags::deriv<CurvedScalarWave::Tags::Psi,
                                    tmpl::size_t<3>, Frame::Inertial>>(*result);
  auto& rotated_dx = temps.at({next_slot + 2});
  rotated_dx = cos_phi * get<0>(gradient) - sin_phi * get<1>(gradient);
  get<1>(gradient) = sin_phi * get<0>(gradient) + cos_phi * get<1>(gradient);
  get<0>(gradient) = rotated_dx;
}}
}}  // namespace CurvedScalarWave::Worldtube
"""
    output.write_text(text)
    print(
        json.dumps(
            {
                "input_sha256": hashes,
                "cse_temporaries": len(replacements),
                "scalar_assignments": len(scalar_assignments),
                "vector_assignments": len(vector_assignments),
                "buffers_including_rotation": next_slot + 3,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expressions-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    generate(args.expressions_dir, args.output)
