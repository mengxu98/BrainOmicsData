# Plotting source dependency

`scop-8aec27fc-source.tar.gz` contains unmodified package source files from
`mengxu98/scop` commit `8aec27fc594ec8763eb37cde719d4ef5c719118c` (version 0.9.2).
It includes DESCRIPTION, NAMESPACE, LICENSE.md, R, src, man, cleanup and configure
files. Example datasets and application assets are omitted; the archive is the
code needed to inspect the plotting implementation, not a binary environment.
The GPL-3.0-or-later license in the archive applies to this dependency. The
repository's MIT license does not replace that license.

The archive and `FeatureDimPlot.R` hashes are recorded in
`../scop-plotting.lock.tsv`. `../restore_r.R` installs the package and its
dependencies in the locked R 4.5.1 environment.
