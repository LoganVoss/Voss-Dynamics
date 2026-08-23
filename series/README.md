# Unified volume

The unified edition places the exact requested two-page front matter before the three canonical papers in order. The assembler imports each paper’s outline and link annotations, so the volume remains navigable instead of becoming a flat print-only concatenation.

## Build

~~~bash
tectonic main.tex
uv run --with pypdf==6.10.0 python build_volume.py
~~~

`main.tex` produces the two-page front matter. `build_volume.py` appends the three canonical PDFs and writes `Voss-Dynamics-Information-Representation-and-Discovery.pdf` with a hierarchical outline, working internal links, and release metadata.
