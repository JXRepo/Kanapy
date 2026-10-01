# Kanapy grain boundary proposal

Open `talk.pdf`. The presentation has a title, 14 main slides and six technical
backup slides. `talk.tex` is the editable source; `speaker_notes_zh.md` explains
every page in Chinese.

All presentation files and the 15 figures used by the slides are in this
directory. All figures are PNG images. `slide_overview.png`
shows all pages together. `logos/` is the only subdirectory and contains the
five RUB/ICAMS branding images used by the theme, also in PNG format.

## Build

From this directory, run:

```bash
make
```

This compiles `talk.tex` and removes the intermediate LaTeX files. The PDF stays
in this directory. To compile while retaining the log for inspection, run:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error talk.tex
```

## Reproduce the figures

Use the repository's configured Python environment:

```bash
/home/users/xuejungs/anaconda3/envs/knpy/bin/python generate_figures.py
```

The script imports Kanapy from the repository's `src` directory, then saves
the PNGs here. It also supports `--only controls`, `--only local`,
`--only context` and `--only periodic`. It creates no figure or data directories.
Required packages are the normal Kanapy dependencies and Matplotlib.

`figure_provenance.json` records the source version, parameters and measurements.
`context_apd.npz` contains the frozen contextual APD and voxel labels.
`context_descriptor.json` and `context_generation.json` record its input and
generation settings. `local_before.npy`, `local_greedy.npy` and `local_joint.npy`
contain the constructed voxel candidates. `context_rve_low_fill.pkl` caches
the packed RVE so redrawing does not repeat packing; deleting it makes the
script regenerate that example with the fixed seed.

## Evidence and limits

The figures use Kanapy 6.5.5 from repository commit
`19713e5612fdaaea2acc4378b2c897be3f54dd2d`. The contextual RVE has 23 grains in
a 24 micrometre box; packing relaxation reached zero remaining overlap contacts.
The controlled APD figures reproduce the existing analytical pinch regression.

In the constructed voxel example, connecting A alone disconnects B. The joint
candidate changes ten of 1331 cells, leaves all three grains face connected and
passes Kanapy's local voxel boundary checks. The figure script asserts these
results and the regular/pinched vertex-link counts.

This is a research proposal. The joint candidate was constructed manually;
the automatic search, continuous-partition repair, periodic repair, FE volume
meshing and mechanical validation remain proposed work. No figure is presented
as Hartmaier's unavailable failing specimen.

The complete PDF was compiled, rendered and checked for layout and page bounds.
PNG files were decoded to verify their integrity. The slides retain ordinary
multi-grain junctions and distinguish volume connectivity, local boundary
topology and geometric validity.
