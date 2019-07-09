#'Example of an \code{haplohh} object
#'@description The object contains haplotype data for 140 cattle individuals (280 haplotypes) belonging to the
#'Creole breed from Guadeloupe (CGU) and 1424 markers (mapping to chromosome BTA12).
#'@usage data(haplohh_cgu_bta12)
#'@format An object of \code{\link{haplohh-class}}.
#'@references Gautier, M. and Naves, M. (2011). Footprints of selection in the ancestral admixture of a New World Creole cattle breed. \emph{Molecular Ecology}, \strong{20}, 3128-3143.
#'@seealso \code{\link{data2haplohh}}
"haplohh_cgu_bta12"


#'Copy example input files into current working directory
#'@description This function copies the following example files to the current working directory:
#'\itemize{
#'\item \code{example1.hap} "example 1" haplotype file in "standard format"
#'\item \code{example1.map} "example 1" marker information file
#'\item \code{example1.vcf} "example 1" as vcf file
#'\item \code{example2.hap} "example 2" haplotype file in "standard format"
#'\item \code{example2.map} "example 2" marker information file
#'\item \code{example2.vcf} "example 2" as vcf file
#'\item \code{ms.out output} from a small simulation by the program 'ms'
#'\item \code{bta12_cgu.hap} an haplotype file in "standard format"
#'\item \code{bta12_cgu.thap} an haplotype file in "transposed format"
#'\item \code{bta12_hapguess_switch.out} an haplotype file in fastphase output format
#'\item \code{map.inp} a marker information file for all bta_cgu markers
#'}
#'Example 1 was used in (Gautier 2017) to explain the various EHH derived statistics calculated by this package.
#'Example 2 is an extension containing multi-allelic markers and missing values.
#'
#'The bta12 files contain data for 280 haplotypes, originating from 140 individuals belonging to the
#'Creole cattle breed from Guadeloupe, at 1.424 markers mapping to bovine chromosome 12 (BTA12) (Gautier 2011).
#'@references Gautier, M. and Naves, M. (2011). Footprints of selection in the ancestral admixture of a New World Creole cattle breed. \emph{Molecular Ecology}, \strong{20}, 3128-3143.
#'
#'Gautier, M., Klassmann,  A. and Vitalis, R. (2017). rehh 2.0: a reimplementation of the R package rehh to detect positive selection from haplotype structure. \emph{Molecular Ecology Resources}, \strong{17}, 78-90.
#'@seealso \code{\link{data2haplohh}}, \code{\link{remove.example.files}}
#'@export
#'@import rehh.data
#'@importFrom utils unzip
make.example.files <- function() {
  rehh.data_files <- c(
    "bta12_cgu.hap.zip",
    "bta12_cgu.thap.zip",
    "bta12_hapguess_switch.out.zip",
    "map.inp.zip"
  )

  for (file in rehh.data_files) {
    file.copy(system.file(file, package = 'rehh.data'),
              file)
    unzip(file)
    #remove zipped files
    if (file.exists(file)) {
      file.remove(file)
    }
  }

  extdata_files <- c(
    "example1.hap",
    "example1.map",
    "example1.vcf",
    "example2.hap",
    "example2.map",
    "example2.vcf",
    "ms.out",
    "bta12_cgu.vcf.gz"
  )

  for (file in extdata_files) {
    file.copy(system.file('extdata', file, package = "rehh"), file)
  }
}

#'Remove example files from current working directory.
#'@description Remove example files from current working directory.
#'@details Removes the files created by \code{make.example.files()}. 
#'No error is thrown, if files do not exist.
#'@seealso \code{\link{make.example.files}}
#'@export
remove.example.files <- function() {
  rehh.data_files <- c(
    "bta12_cgu.hap.zip",
    "bta12_cgu.thap.zip",
    "bta12_hapguess_switch.out.zip",
    "map.inp.zip"
  )
  for (file in rehh.data_files) {
    #zipped files
    if (file.exists(file))
      file.remove(file)
    #unzipped files
    if (file.exists(substr(file, 1, nchar(file) - 4)))
      file.remove(substr(file, 1, nchar(file) - 4))
  }

  extdata_files <- c(
    "example1.hap",
    "example1.map",
    "example1.vcf",
    "example2.hap",
    "example2.map",
    "example2.vcf",
    "ms.out",
    "bta12_cgu.vcf.gz"
  )

  for (file in extdata_files) {
    if (file.exists(file))
      file.remove(file)
  }
}
