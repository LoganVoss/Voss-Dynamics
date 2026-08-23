#!/usr/bin/env python3
"""Assemble the front matter and canonical papers without dropping navigation."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from pypdf import PdfReader, PdfWriter


logging.getLogger("pypdf").setLevel(logging.ERROR)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEFAULT_FRONT = HERE / "main.pdf"
DEFAULT_OUTPUT = HERE / "Voss-Dynamics-Information-Representation-and-Discovery.pdf"
PAPERS = (
    (
        "I. The Principle of Full Invertibility",
        ROOT
        / "(1) The Principle of Full Invertibility"
        / "The_Principle_of_Full_Invertibility.pdf",
    ),
    (
        "II. Emergent Predictive Representation",
        ROOT
        / "(2) Emergent Predictive Representation"
        / "Emergent-Predictive-Representation.pdf",
    ),
    (
        "III. Pressure-Driven Invariant Synthesis",
        ROOT
        / "(3) Pressure-Driven Invariant Synthesis"
        / "Pressure-Driven-Invariant-Synthesis.pdf",
    ),
)


def import_outline(
    writer: PdfWriter,
    reader: PdfReader,
    items: list,
    *,
    page_offset: int,
    parent,
) -> None:
    """Copy an outline while translating every destination into the volume."""

    previous = None
    for item in items:
        if isinstance(item, list):
            if previous is not None:
                import_outline(
                    writer,
                    reader,
                    item,
                    page_offset=page_offset,
                    parent=previous,
                )
            continue

        source_page = reader.get_destination_page_number(item)
        if source_page < 0:
            previous = None
            continue
        previous = writer.add_outline_item(
            item.title,
            page_offset + source_page,
            parent=parent,
            color=getattr(item, "color", None),
            bold=bool(getattr(item, "bold", False)),
            italic=bool(getattr(item, "italic", False)),
        )


def assemble(front_matter: Path, output: Path) -> None:
    front = PdfReader(front_matter)
    if len(front.pages) != 2:
        raise ValueError(f"front matter must contain exactly 2 pages, found {len(front.pages)}")

    writer = PdfWriter()
    writer.append(front, import_outline=False)
    writer.add_outline_item("Cover", 0)
    writer.add_outline_item("One problem, three directions", 1)

    for title, paper in PAPERS:
        if not paper.is_file():
            raise FileNotFoundError(paper)
        reader = PdfReader(paper)
        page_offset = len(writer.pages)
        writer.append(reader, import_outline=False)
        paper_outline = writer.add_outline_item(title, page_offset)
        import_outline(
            writer,
            reader,
            reader.outline,
            page_offset=page_offset,
            parent=paper_outline,
        )

    writer.add_metadata(
        {
            "/Title": "Voss Dynamics: Information, Representation, and Discovery",
            "/Author": "Logan Voss",
            "/Subject": "The complete Voss Dynamics trilogy",
            "/Keywords": "information, representation, observability, invariant synthesis",
        }
    )
    writer.page_mode = "/UseOutlines"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as handle:
        writer.write(handle)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--front-matter", type=Path, default=DEFAULT_FRONT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    assemble(args.front_matter.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
