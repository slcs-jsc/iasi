# AGENTS.md

## Purpose

This repository contains the IASI Code Collection, a small C codebase for working with Eumetsat IASI observations. Most development happens in `src/`, with vendored dependencies in `libs/`, reference data in `data/`, documentation in `docs/`, and a runnable example pipeline in `projects/example/`.

## Repository layout

- `src/`: C sources, headers, and the project `Makefile`. Binaries are built in place.
  `src/jurassic.c` and `src/jurassic.h` are a special case: they are copied into this repository for convenience from the external JURASSIC radiative transfer model at `https://github.com/slcs-jsc/jurassic` and should be treated as vendored upstream code.
- `libs/`: vendored source archives plus `build.sh` to produce local dependency installs under `libs/build/`.
- `data/`: CODA definition files used at runtime.
- `tests/data/`: helper script to download the sample `.nat` test input.
- `projects/example/`: end-to-end example workflow and expected plot outputs.
- `docs/`: Doxygen and MkDocs configuration plus the manual sources.

## Build and validation

Dependencies are GCC, `make`, GSL, netCDF, CODA, and for some workflows MPI/OpenMP. If local libraries are needed, build them from the repository root:

```bash
cd libs
./build.sh
```

Compile from `src/`:

```bash
cd src
make
```

Useful build variants from `src/`:

- `make STATIC=1`: request a static build.
- `make PROF=1`: enable profiling flags.
- `make COV=1`: enable coverage instrumentation.
- `make INFO=1`: emit compiler optimization info.
- `make clean`: remove binaries, objects, and coverage artifacts.

The `Makefile` treats warnings as errors via `-Werror`, so keep the tree warning-free.

There is no dedicated unit test target. Validate changes with the narrowest relevant command first, then the example workflow when behavior changes:

```bash
cd projects/example
./run.sh
```

Notes:

- `projects/example/run.sh` expects compiled binaries in `src/`.
- The example sets `CODA_DEFINITION=../../data/` and `LD_LIBRARY_PATH=../../libs/build/lib/:$LD_LIBRARY_PATH`.
- Example data is downloaded by `tests/data/download.sh` into `tests/data/*.nat`.
- Plot generation also expects `gnuplot` and an external coastline file at `~/wrk/coast/wcl.tab`.

## Code conventions

- Follow the existing C style in `src/*.c` and `src/*.h`: function declarations are split across lines, comments are short and procedural, and naming is mostly snake_case.
- Keep edits minimal and local. This is an old-school C codebase with large static buffers and explicit loops; do not refactor broadly unless the task requires it.
- Preserve compatibility with the current warning set in `src/Makefile`, especially conversion, prototype, shadowing, and declaration warnings.
- If you reformat code, use the existing `make indent` target rather than introducing a new formatter.
- Keep runtime paths relative where the project already does so.
- Be especially careful in `src/jurassic.c` and `src/jurassic.h`: they originate from the external JURASSIC radiative transfer model repository, so prefer minimal patches, keep provenance clear, and avoid style-only edits there.

## Documentation and generated files

- Update `README.md` when build, dependency, or runtime behavior changes.
- Update files in `docs/` when CLI behavior or outputs change in a user-visible way.
- Do not commit generated artifacts that are already ignored, especially:
  - `libs/build/`
  - binaries and objects in `src/`
  - coverage outputs in `src/`
  - downloaded `.nat` files in `tests/data/`
  - generated `.tab` and `.nc` files in `projects/example/`

## Practical agent guidance

- Start by reading `README.md` and `src/Makefile` before changing build logic.
- When changing a specific executable, inspect the corresponding source in `src/` and any shared helpers in `src/libiasi.c`, `src/libiasi.h`, `src/jurassic.c`, and `src/jurassic.h`.
- Prefer targeted validation over full rebuilds when possible, but rebuild affected binaries before running the example workflow.
- Avoid changing vendored archives in `libs/` unless the task is explicitly about dependency updates.
- Treat `src/jurassic.c` and `src/jurassic.h` similarly to vendored code. If they must be changed, note clearly that the edit diverges from the upstream JURASSIC source.
