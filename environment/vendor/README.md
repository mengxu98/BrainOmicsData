# Plotting source dependency

`scop-8aec27fc-source.tar.gz` contains unmodified package source files from
`mengxu98/scop` commit `8aec27fc594ec8763eb37cde719d4ef5c719118c` (version 0.9.2).
It includes DESCRIPTION, NAMESPACE, LICENSE.md, R, src, man, cleanup and configure
files. Example datasets and application assets are omitted; the archive is the
code needed to inspect the plotting implementation, not a binary environment.
The GPL-3.0-or-later license in the archive applies to this dependency. The
repository's MIT license does not replace that license.

The archive and `FeatureDimPlot.R` hashes are recorded in
`../scop-plotting.lock.tsv`. For Figure S4, extract to a temporary directory and
set `SCOP_SOURCE_PATH` to the extracted source folder. `pkgload::load_all` loads
it using the installed dependency packages, without modifying the installed
scop package. Figure S4 used this plotting build; the separate mapping build in
`analysis/scop-knn.lock.tsv` describes the external-query calculation.

A small `FeatureDimPlot` call on an already normalized Seurat fixture was checked
after loading the archived R source. Without compiling the native code,
`pkgload` reports an unavailable DLL; the tested plotting path works, but this
source-only loading route is not intended for scop normalization or integration
functions that call native code. Figure S4 computes its normalization explicitly
before constructing the plotting object.
