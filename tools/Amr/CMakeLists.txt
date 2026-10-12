# Distributed under the MIT License.
# See LICENSE.txt for details.

spectre_python_add_module(
  Amr
  MODULE_PATH "tools"
  PYTHON_FILES
  __init__.py
  Chalmers.py
  ChalmersSymmetry.py
  ChalmersStudy.py
  ChalmersEvolution.py
  ChalmersCube.py
  ChalmersConstraint.py
  ChalmersKerr.py
  ChalmersKerrAnalysis.py
  ChalmersKerrReversals.py
  ChalmersSphere.py
  ChalmersWedges.py
  ChalmersSingleHole.py
  ChalmersGaugeRelaxation.py
  ChalmersGaugeRelaxationAnalysis.py
  ChalmersGhFeedbackAnalysis.py
  ChalmersPileUpAudit.py
  ChalmersCoarseningExperiment.py
  ChalmersSpecCoarsening.py
  ChalmersTiltedKerr.py
  ChalmersGaugeRelaxationReference.py
  ChalmersSphericalReference.py
  ChalmersMeshHistory.py
)
