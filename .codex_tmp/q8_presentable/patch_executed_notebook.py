#!/usr/bin/env python3
"""Apply explanatory-cell updates while preserving executed notebook output."""

from __future__ import annotations

import sys
from pathlib import Path

import nbformat


def heading(cell) -> str:
    if cell.cell_type != "markdown":
        return ""
    first_line = cell.source.lstrip().splitlines()[0]
    return first_line.strip()


def find_heading(notebook, text: str) -> int:
    matches = [
        index for index, cell in enumerate(notebook.cells)
        if heading(cell) == text
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one heading {text!r}, found {len(matches)}")
    return matches[0]


def main(executed_path: Path, generated_path: Path) -> None:
    with executed_path.open() as stream:
        executed = nbformat.read(stream, as_version=4)
    with generated_path.open() as stream:
        generated = nbformat.read(stream, as_version=4)

    generated_by_heading = {
        heading(cell): cell
        for cell in generated.cells
        if heading(cell)
    }

    # Replace the title cell.
    executed.cells[0].source = generated.cells[0].source

    # Insert the navigation cell immediately after the title.
    navigation_heading = "## How to navigate this notebook"
    if not any(
            heading(cell) == navigation_heading
            for cell in executed.cells):
        executed.cells.insert(1, generated_by_heading[navigation_heading])

    # Insert the hierarchy explanation before the first-order response section.
    hierarchy_heading = "## Perturbative hierarchy used below"
    if not any(
            heading(cell) == hierarchy_heading
            for cell in executed.cells):
        first_order_index = find_heading(
            executed, "## 4. First-order response blocks")
        executed.cells.insert(
            first_order_index, generated_by_heading[hierarchy_heading])

    # Replace the old closing caution with the two new final sections.
    final_headings = (
        "## 17. Companion time-derivative analysis",
        "## Current status",
    )
    final_generated_index = find_heading(generated, final_headings[0])
    final_cell = generated.cells[final_generated_index]

    # The generated source places both final headings in one markdown cell.
    old_last = executed.cells[-1]
    if old_last.cell_type != "markdown":
        raise RuntimeError("Expected the executed notebook to end in markdown")
    old_last.source = final_cell.source

    executed.metadata["title"] = (
        "q=8 hierarchical harmonic-worldtube matching")
    executed.metadata.setdefault("codex", {})
    executed.metadata["codex"]["companion_report"] = (
        "q8_harmonic_worldtube_matching_report.pdf")
    nbformat.validate(executed)
    with executed_path.open("w") as stream:
        nbformat.write(executed, stream)

    print(
        f"patched {executed_path}: {len(executed.cells)} cells, "
        "outputs preserved")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(
            "usage: patch_executed_notebook.py EXECUTED GENERATED")
    main(Path(sys.argv[1]), Path(sys.argv[2]))
