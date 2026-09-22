args <- commandArgs(TRUE); run <- normalizePath(args[1])
library(jsonlite)
lib <- file.path(run,'environment/R-library')
dir.create(lib,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(lib,.libPaths()))
source(file.path(run,'pipeline/functions/processed_object.R'))
manifest <- fromJSON(file.path(run,'software_sources/harmony_source.json'))
package <- file.path(run,'software_sources',paste0('harmony_',manifest$version,'.tar.gz'))
stopifnot(manifest$version=='2.0.5',processed_file_sha256(package)==manifest$sha256)
installed <- suppressWarnings(tryCatch(as.character(packageVersion('harmony')),
                                       error=function(e) 'missing'))
if(installed != manifest$version) {
  install.packages(package,repos=NULL,type='source',lib=lib)
}
stopifnot(as.character(packageVersion('harmony'))==manifest$version)
writeLines(capture.output(sessionInfo()),file.path(run,'environment/R_session.txt'))
