# Distributed under the MIT License.
# See LICENSE.txt for details.

spectre_add_python_bindings_test(
  "Unit.LinearOperators.Python.Chalmers"
  Test_Chalmers.py
  "Unit;LinearOperators;Python"
  "PyLinearOperators;PyDomainCreators;PyH5"
  TIMEOUT 10)

spectre_add_python_bindings_test(
  "Unit.LinearOperators.Python.DefiniteIntegral"
  Test_DefiniteIntegral.py
  "Unit;LinearOperators;Python"
  PyLinearOperators)


spectre_add_python_bindings_test(
  "Unit.LinearOperators.Python.PartialDerivatives"
  Test_PartialDerivatives.py
  "Unit;LinearOperators;Python"
  PyLinearOperators)

spectre_add_python_bindings_test(
  "Unit.LinearOperators.Python.PowerMonitors"
  Test_PowerMonitors.py
  "Unit;LinearOperators;Python"
  PyLinearOperators)
