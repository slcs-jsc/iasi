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

Run `make clean` before rebuilding with changed compiler options; Make does not track command-line flag changes.

The `Makefile` treats warnings as errors via `-Werror`, so keep the tree warning-free.

`make check` currently builds the executables and reports that no automated checks are defined yet. Validate behavior with the narrowest relevant command first, then the example workflow when behavior changes:

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
- Before changing behavior, check existing callers, the example workflow, and documented expectations; question assumptions and explain non-obvious choices. Add abstractions or dependencies only for a concrete need.
- Reuse existing data structures where possible. Avoid new `static` helper functions; use a macro for a small operation when its arguments can be evaluated safely, or a regular function declared and documented in `src/libiasi.h`.
- Preserve compatibility with the current warning set in `src/Makefile`, especially conversion, prototype, shadowing, and declaration warnings.
- Format C sources and headers from `src/` with `make indent` after a successful build. It rewrites all `*.c` and `*.h` files, including copied JURASSIC sources, so review the diff and keep only intended changes.
- Keep runtime paths relative where the project already does so.
- Be especially careful in `src/jurassic.c` and `src/jurassic.h`: they originate from the external JURASSIC radiative transfer model repository, so prefer minimal patches, keep provenance clear, and avoid style-only edits there.

## Documentation and generated files

- Update `README.md` when build, dependency, or runtime behavior changes.
- Update files in `docs/` when CLI behavior or outputs change in a user-visible way.
- Keep the MkDocs manual in `docs/manual/` and Doxygen API documentation current with relevant code, build, and workflow changes. Build affected docs from `src/` with `make mkdocs` and/or `make doxygen` when the tools are available.
- For new or revised API Doxygen comments, use a substantive description beyond `@brief`, document parameters and behavior where relevant, and include `@author` with the responsible author's name. Follow the detailed `/** ... */` examples in `src/jurassic.h` without editing that copied dependency for style alone; do not invent attribution.
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
- Investigate unexpected example results before changing expected outputs; update only affected references after verifying the new results.
- Avoid changing vendored archives in `libs/` unless the task is explicitly about dependency updates.
- Treat `src/jurassic.c` and `src/jurassic.h` similarly to vendored code. If they must be changed, note clearly that the edit diverges from the upstream JURASSIC source.
