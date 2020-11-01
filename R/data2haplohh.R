#'Convert data from input file to an object of class haplohh
#'@description Convert input data files to an object of \code{\link{haplohh-class}}.
#'@param hap_file file containing haplotype data (see details below).
#'@param map_file file containing map information (see details below).
#'@param min_perc_geno.hap threshold on percentage of missing data for haplotypes
#'(haplotypes with less than \code{min_perc_geno.hap} percent of markers genotyped are discarded). Default is \code{NA},
#'hence no constraint.
#'@param min_perc_geno.mrk threshold on percentage of missing data for markers (markers genotyped on less than
#'\code{min_perc_geno.mrk} percent of haplotypes are discarded). By default, \code{min_perc_geno.mrk=100},
#'hence only fully genotyped markers are retained.
#'This value cannot be set to \code{NA} or zero.
#'@param min_maf threshold on the Minor Allele Frequency. Markers having a MAF lower than or equal to minmaf are discarded.
#'In case of multi-allelic markers the second-most frequent allele is referred to as minor allele.
#'Setting this value to zero eliminates monomorphic sites. Default is \code{NA},
#'hence no constraint.
#'@param chr.name name of the chromosome considered (relevant if data for several chromosomes is
#'contained in the haplotype or map file).
#'@param popsel code of the population considered (relevant for fastPHASE output which
#'can contain haplotypes from various populations).
#'@param recode.allele *Deprecated*. logical. \code{FALSE} by default. \code{TRUE} forces parameter \code{allele_coding} to \code{"map"},
#'\code{FALSE} leaves it unchanged.
#'@param allele_coding the allele coding provided by the user. Either \code{"12"} (default), \code{"01"}, \code{"map"} or \code{"none"}.
#'The option is irrelevant for vcf files and ms output.
#'@param haplotype.in.columns logical. If \code{TRUE}, phased input haplotypes are assumed to be in columns (as produced
#'by the SHAPEIT2 program (O'Connell et al., 2014).
#'@param remove_multiple_markers logical. If \code{FALSE} (default), conversion
#'stops, if multiple markers with the same chromosomal position are encountered.
#'If \code{TRUE}, duplicated markers are removed (all but the first marker with identical positions).
#'@param polarize_vcf logical. Only of relevance for vcf files. If \code{TRUE} (default), tries to polarize
#'variants with help of the AA entry in the INFO field. Unpolarized alleles are discarded.
#'If \code{FALSE}, allele coding of vcf file is used unchanged as internal coding.
#'@param capitalize_AA logical. Only of relevance for vcf files with ancestral allele information.
#'Low confidence ancestral alleles are usually coded by lower-case letters. If \code{TRUE} (default), these are
#'changed to upper case before the alleles of the sample are matched for polarization.
#'@param vcf_reader library used to read vcf. By default, low-level parsing is
#'performed using the generic package \code{data.table}. In order to read compressed files,
#'the package \code{R.utils} must be installed, too.
#'If the specialized package \code{vcfR} is available, set this parameter to \code{"vcfR"}.
#'@param position_scaling_factor intended primarily for output of ms where
#'positions lie in the interval [0,1]. These can be rescaled to sizes
#'of typical markers in real data.
#'@param verbose logical. If \code{TRUE} (default), report verbose progress.
#'@details Five haplotype input formats are supported:
#'\itemize{
#'\item a "standard format" with haplotypes in rows and markers in columns (with no header, but a haplotype ID/name in
#'the first column).
#'\item a "transposed format" similar to the one produced by the phasing program SHAPEIT2
#'(O'Connell et al., 2014) in which haplotypes are in columns and markers in rows
#'(with neither header nor marker IDs nor haplotype IDs).
#'\item output files from the fastPHASE program (Sheet and Stephens, 2006).
#'If haplotypes from several different population were phased simultaneously (-u fastPHASE option
#'was used), it is necessary to specify the population of interest by parameter \code{popsel}
#'(if this parameter is not or wrongly set, the error message will provide a list of
#'the population numbers contained in the file).
#'\item files in variant call format (vcf). No mapfile is needed is this case. If
#'the file contains several chromosomes, it is necessary to choose one by parameter
#'\code{chr.name}.
#'\item output of the simulation program 'ms'. No mapfile is needed in this case. If the file
#'contains several 'runs', a specific number has to be specified by the
#'parameter \code{chr.name}.
#'}
#'The "transposed format" has to be explicitly set while the other formats
#'are recognized automatically.
#'
#'The map file contains marker information in three, or, if it is used for
#'polarization (see below), five columns:
#'\itemize{
#'\item marker name/id
#'\item chromosome
#'\item position (physical or genetic)
#'\item ancestral allele encoding
#'\item derived allele encoding
#'}
#'The markers must be in the same order as in the haplotype file. If
#'several chromosomes are represented in the map file, it is necessary to choose that
#'which corresponds to the haplotype file by parameter \code{chr.name}.
#'
#'Haplotypes can be given either with alleles already coded as numbers (in two possible ways)
#'or with the actual alleles (e.g. nucleotides) which can be translated into numbers
#'either using the fourth and fifth column of the map file or by their alpha-numeric order.
#'Correspondingly, the parameter \code{allele_coding} has to be set to either \code{"12"},
#'\code{"01"}, \code{"map"} or \code{"none"}:
#'\itemize{
#'\item \code{"12"}: 0 represents missing values, 1 the ancestral allele
#'and 2 (or higher integers) derived allele(s).
#'\item \code{"01"}: \code{NA} or '.' (a point) represent missing values, 0 the
#'ancestral and 1 (or higher integers) derived allele(s).
#'\item \code{"map"}: for each marker, the fourth column of the map file
#'defines the ancestral allele and the fifth column derived alleles.
#'In case of multiple derived alleles, they must be separated by commas without space.
#'Alleles in the haplotype file which do not appear in neither of the two columns
#'of the map file are regarded as missing values (\code{NA}).
#'\item \code{"none"}: \code{NA} or '.' (a point) represent missing values, otherwise for each
#'marker the allele that comes first in alpha-numeric
#'order is coded by 0, the next by 1, etc. Evidently, this coding does not convey
#'any information about allele status as ancestral or derived, hence the alleles
#'cannot be regarded as polarized.
#'}
#'The information of allelic ancestry is exploited only in the frequency-bin-wise
#'standardization of iHS (see \code{\link{ihh2ihs}}). However, although ancestry status does
#'not figure in the formulas of the cross populations statistics
#'Rsb and XP-EHH, their values do depend on the assigned status.
#'
#'The arguments \code{min_perc_geno.hap},
#'\code{min_perc_geno.mrk} and \code{min_maf} are evaluated in this order.
#'@return The returned value is an object of \code{\link{haplohh-class}}.
#'@references Scheet P, Stephens M (2006) A fast and flexible statistical model for large-scale population genotype
#'data: applications to inferring missing genotypes and haplotypic phase. \emph{Am J Hum Genet}, \strong{78}, 629-644.
#'
#'O'Connell J, Gurdasani D, Delaneau O, et al (2014) A general approach for haplotype phasing
#'across the full spectrum of relatedness. \emph{PLoS Genet}, \strong{10}, e1004234.
#'@examples
#'#copy example files into the current working directory.
#'make.example.files()
#'#create object using a haplotype file in "standard format"
#'hap <- data2haplohh(hap_file = "bta12_cgu.hap",
#'                    map_file = "map.inp",
#'                    chr.name = 12,
#'                    allele_coding = "map")
#'#create object using fastPHASE output
#'hap <- data2haplohh(hap_file = "bta12_hapguess_switch.out",
#'                    map_file = "map.inp",
#'                    chr.name = 12,
#'                    popsel = 7,
#'                    allele_coding = "map")
#'#clean up demo files
#'remove.example.files()
#'@export
#'@importFrom methods new
#'@importFrom utils read.table
data2haplohh <-
  function(hap_file,
           map_file = NA,
           min_perc_geno.hap = NA,
           min_perc_geno.mrk = 100,
           min_maf = NA,
           chr.name = NA,
           popsel = NA,
           recode.allele = FALSE,
           allele_coding = "12",
           haplotype.in.columns = FALSE,
           remove_multiple_markers = FALSE,
           polarize_vcf = TRUE,
           capitalize_AA = TRUE,
           vcf_reader = "data.table",
           position_scaling_factor = NA,
           verbose = TRUE) {
    ## check parameters
    ### empty haplotype makes no sense, but is of no harm
    if (!is.na(min_perc_geno.hap) &
        (is.na(min_perc_geno.hap) | min_perc_geno.hap < 0 |
         min_perc_geno.hap > 100)) {
      stop("min_perc_geno.hap should lie in the interval [0,100].",
           call. = FALSE)
    }
    ### empty marker will cause trouble -> forbid zero or NA
    if (is.na(min_perc_geno.mrk) | min_perc_geno.mrk <= 0 |
        min_perc_geno.mrk > 100) {
      stop("min_perc_geno.mrk should lie in the interval (0,100].",
           call. = FALSE)
    }
    ### minor frequency can be maximal 0.5
    if (!is.na(min_maf)) {
      if (!is.numeric(min_maf) | min_maf < 0 |
          min_maf > 0.5) {
        stop("min_maf should lie in the interval [0,0.5].", call. = FALSE)
      }
    }
    ### deprecated function
    if (recode.allele) {
      warning("Deprecated option: recode.allele. Use 'allele_coding' instead.")
      allele_coding <- "map"
    }
    
    ### vcf_readere
    if (!vcf_reader %in% c("vcfR", "data.table")) {
      stop("vcf_reader must be either 'data.table' or 'vcfR'.", call. = FALSE)
    }
    
    ### check possible coding options
    if (is.na(allele_coding) |
        !(allele_coding %in% c("12", "01", "map", "none"))) {
      stop("allele_coding has to be either '12', '01', 'map' or 'none'.",
           call. = FALSE)
    }
    
    if (!is.na(position_scaling_factor)) {
      if (!is.numeric(position_scaling_factor) |
          position_scaling_factor <= 0) {
        stop(
          paste0(
            "position_scaling_factor must be a positive real number.\n",
            "Conversion stopped."
          ),
          call. = FALSE
        )
      }
    }
    
    ## perform conversion
    if (verbose)
      cat("* Reading input file(s) *\n")
    
    if (is.na(hap_file)) {
      stop("No haplotype file specified. Conversion stopped.", call. = FALSE)
    }
    
    if (is.vcf(hap_file)) {
      hh <-
        read.vcf(
          hap_file,
          chr.name = chr.name,
          polarize_vcf = polarize_vcf,
          capitalize_AA =  capitalize_AA,
          vcf_reader = vcf_reader,
          verbose = verbose
        )
    } else if (is.ms(hap_file)) {
      hh <- read.ms(hap_file, chr.name = chr.name, verbose = verbose)
    } else {
      #### begin map_file
      if (is.na(map_file)) {
        stop("No map file specified. Conversion stopped.", call. = FALSE)
      }
      
      map <-
        read.table(
          map_file,
          row.names = 1,
          colClasses = "character",
          stringsAsFactors = FALSE
        )
      
      if (allele_coding == "map" & ncol(map) < 4) {
        stop(
          paste0(
            "Wrong format for map file. ",
            map_file,
            " should contain 5 columns without header:\n"
            ,
            "Marker id, Chromosome name, Position, Ancestral Allele and Derived Allele.\n",
            "Conversion stopped."
          ),
          call. = FALSE
        )
      } else{
        if (ncol(map) < 2) {
          stop(
            paste0(
              "Wrong format for map file. ",
              map_file,
              " should contain 3 columns without header:\n",
              "Marker id, Chromosome name and Position.\n",
              "Conversion stopped."
            ),
            call. = FALSE
          )
        }
      }
      
      chr.name <-
        check_chromosome_names(map_file, unique(as.character(map[, 1])), chr.name)
      
      ### subset map data frame to specified chromosome
      map <- map[as.character(map[, 1]) == chr.name, ]
      
      ### set first slots of haplohh
      hh <- new("haplohh")
      
      hh@chr.name <- chr.name
      hh@positions <- as.numeric(map[[2]])
      mrk.names <- row.names(map)
      names(hh@positions) <- mrk.names
      
      if (verbose)
        cat("Map info:",
            nrow(map),
            "markers declared for chromosome",
            hh@chr.name,
            ".\n")
      
      ### Fichier haplo
      
      if (haplotype.in.columns) {
        tmp_haplo <- read.transposed(hap_file, verbose = verbose)
      } else{
        if (is.fastPhase(hap_file)) {
          tmp_haplo <-
            read.fastPhase(hap_file, popsel = popsel, verbose = verbose)
        } else{
          #fichier au format standard
          tmp_haplo <- read.standard(hap_file, verbose = verbose)
        }
      }
      if (ncol(tmp_haplo) != nrow(map)) {
        stop(
          paste0(
            "The number of markers in the haplotype file (",
            ncol(tmp_haplo),
            ") is not equal\nto the number of markers in the map file (",
            nrow(map),
            ").\nConversion stopped."
          ),
          call. = FALSE
        )
      }
      
      
      if (allele_coding == "map") {
        ### recode with help of map file
        if (verbose)
          cat("Alleles are being recoded according to fourth and fifth column of map file.\n")
        
        ### split into list of alleles
        allele_list <-
          strsplit(paste(map[, 3], map[, 4], sep = ","), ",", fixed = TRUE)
        
        ### returns NA if allele is not found in allele_list
        hh@haplo <- t(apply(tmp_haplo, 1, function(x) {
          # x is a haplotype, element-wise matching
          mapply(match, x, allele_list) - 1L
        }))
        
        ### remove big object
        rm(allele_list)
      } else if (allele_coding == "none") {
        if (verbose)
          cat(
            paste0(
              "Alleles are being recoded at each marker in alpha-numeric order.\n",
              "*** Consequently, coding does not provide information on ancestry status. ***\n"
            )
          )
        
        
        # (only) point is equivalent to NA
        tmp_haplo[tmp_haplo == "."] <- NA
        
        hh@haplo <- apply(tmp_haplo, 2, function(x) {
          alleles <- sort(unique(x))
          # x is a marker, vector-wise matching
          match(x, alleles) - 1L
        })
        
      }
      else if (allele_coding == "12") {
        if (verbose) {
          cat("Assume that alleles in the haplotype file are coded as:\n")
          cat("0: missing value\n1: ancestral allele\n2, 3, ...: derived allele(s).\n")
        }
        ## matrix must contain only integers
        hh@haplo <-
          matrix(suppressWarnings(as.integer(tmp_haplo)),
                 nrow(tmp_haplo),
                 ncol(tmp_haplo)) - 1L
        
        if (anyNA(hh@haplo)) {
          stop(
            paste0(
              "Alleles are not coded in format \"12\".\n",
              "Check your data or use another value for argument 'allele_coding'.\n",
              "Conversion stopped."
            ),
            call. = FALSE
          )
        }
        
        ### replace -1 by NA
        hh@haplo[hh@haplo == -1] <- NA
        
      } else if (allele_coding == "01") {
        ### hapfile is not a vcf file, but alleles are coded like in vcf
        if (verbose) {
          cat("Assume that alleles in the haplotype file are coded as:\n")
          cat("NA or '.': missing value\n0: ancestral allele\n1, 2, ...: derived allele(s).\n")
        }
        # (only) point is equivalent to NA
        tmp_haplo[tmp_haplo == "."] <- NA
        
        ## convert to integer, stop if allele is neither NA nor integer
        tryCatch(
          hh@haplo <-
            matrix(
              as.integer(tmp_haplo),
              nrow(tmp_haplo),
              ncol(tmp_haplo)
            ),
          warning = function(w) {
            stop(
              paste0(
                "Alleles are not coded in format \"01\".\n",
                "Check your data or use another value for argument 'allele_coding'.\n",
                "Conversion stopped."
              ),
              call. = FALSE
            )
          }
        )
        
      }
      
      rownames(hh@haplo) <- rownames(tmp_haplo)
      colnames(hh@haplo) <- mrk.names
      
      #remove big object
      rm(map)
      rm(tmp_haplo)
    }
    
    ## assert that coded allele numbers are positive
    if (any(hh@haplo < 0, na.rm = TRUE)) {
      stop("Found alleles coded by negative numbers. Conversion stopped.",
           call. = FALSE)
    }
    
    ## assert that positions are ordered
    if (sum(diff(positions(hh)) < 0) > 0) {
      stop("Markers must be ordered numerically in the map file.",
           call. = FALSE)
    }
    
    ## scale positions
    if (!is.na(position_scaling_factor)) {
      hh@positions <- hh@positions * position_scaling_factor
    }
    # 
    ## check for multiple markers
    multiple_markers <- duplicated(hh@positions)
    if (sum(multiple_markers) > 0) {
      if (remove_multiple_markers) {
        hh@positions <- hh@positions[!multiple_markers]
        hh@haplo <- hh@haplo[, !multiple_markers, drop = FALSE]
        warning(paste(
          "Removed",
          sum(multiple_markers),
          "markers with non-unique positions."
        ))
      } else{
        stop(paste(
          sum(multiple_markers),
          "markers have non-unique positions. Conversion stopped."
        ),
        call. = FALSE)
      }
    }
    
    
    # filtering
    hh <- subset(
      hh,
      min_perc_geno.hap = min_perc_geno.hap,
      min_perc_geno.mrk = min_perc_geno.mrk,
      min_maf = min_maf,
      verbose = verbose
    )
    
    if (min(dim(haplo(hh))) == 0 &
        !is.vcf(hap_file) & !is.ms(hap_file)) {
      if (verbose)
        cat("Check whether allele_coding = \"",
            allele_coding,
            "\" is appropriate!\n",
            sep = "")
      
    }
    
    return(hh)
  }

check_chromosome_names <-
  function(file, chr.names_in_file, chr.name) {
    ### check for multiple chromosomes in map data frame
    if (is.na(chr.name)) {
      if (length(chr.names_in_file) != 1) {
        cat("More than one chromosome name in file:",
            file,
            ".\n")
        cat("Here is a list of chromosome names found:\n",
            chr.names_in_file,
            "\n")
        stop("Please specify a chromosome name. Conversion stopped.",
             call. = FALSE)
      }
      chr.name <- chr.names_in_file
    } else{
      chr.name <- as.character(chr.name)
      if (!(chr.name %in% chr.names_in_file)) {
        cat("No markers mapping to chromosome ",
            chr.name,
            " are found in the file ",
            file,
            " .\n")
        
        cat("Here is a list of chromosome names in the map file:\n",
            chr.names_in_file,
            "\n")
        stop("Please specify one chromosome. Conversion stopped.",
             call. = FALSE)
      }
    }
    return(chr.name)
  }

is.fastPhase <- function(hap_file) {
  out_fphase <-
    scan(
      hap_file,
      what = "character",
      sep = "\n",
      quiet = TRUE,
      nlines = 15
    )
  test_fphase_1 <- grep("fastPHASE", out_fphase)
  test_fphase_2 <- grep("BEGIN COMMAND_LINE", out_fphase)
  
  return(length(test_fphase_1) > 0 &
           length(test_fphase_2) > 0)
}

is.vcf <- function(hap_file) {
  firstline <-
    scan(
      hap_file,
      what = "character",
      sep = "\n",
      quiet = TRUE,
      nlines = 1
    )
  return(length(grep("fileformat=VCF", firstline) > 0))
}

is.ms <- function(hap_file) {
  firstline <-
    scan(
      hap_file,
      what = "character",
      sep = "\n",
      quiet = TRUE,
      nlines = 1
    )
  # search for "ms" in first "word" in first line
  return(length(grep("^\\S*ms", firstline) > 0))
}

read.standard <- function(hap_file, verbose) {
  #fichier au format standard
  if (verbose)
    cat("Haplotype input file in standard format assumed.\n")
  
  #retrocompatibilite: ca marchait avec la tabulation avant
  tmp.ncol <-
    length(unlist(strsplit(readLines(hap_file, n = 1), split =
                             "\t|\\s+")))
  tmp_haplo <-
    matrix(scan(hap_file, what = "character", quiet = TRUE),
           ncol = tmp.ncol,
           byrow = TRUE)
  
  rownames <- tmp_haplo[, 1]
  if (anyDuplicated(rownames)) {
    warning(
      paste0(
        "Haplotype identifiers were not unique in haplotype file.\n",
        "They have been modified to become unique."
      ),
      call. = FALSE
    )
    rownames <- make.unique(rownames)
  }
  rownames(tmp_haplo) <- rownames
  
  return(tmp_haplo[, -1])
}

read.transposed <- function(hap_file, verbose) {
  if (verbose)
    cat("Haplotype input file in transposed format assumed.\n")
  
  tmp.nhap <-
    length(unlist(strsplit(readLines(hap_file, n = 1), split =
                             "\t|\\s+")))
  tmp_haplo <-
    matrix(scan(hap_file, what = "character", quiet = TRUE), nrow =
             tmp.nhap)
  return(tmp_haplo)
}


read.fastPhase <- function(hap_file, popsel, verbose) {
  #fichier fastphase
  if (verbose)
    cat("Haplotype input file in fastPHASE format assumed.\n")
  
  out_fphase <- scan(hap_file,
                     what = "character",
                     sep = "\n",
                     quiet = TRUE)
  BEGIN_GENO <- grep("BEGIN GENOTYPES", out_fphase)[1] + 1
  END_GENO <- grep("END GENOTYPES", out_fphase)[1] - 1
  
  # subset to actual data (omit header)
  out_fphase <- out_fphase[BEGIN_GENO:END_GENO]
  
  test_poplabel <- grep("subpop. label:", out_fphase)
  if (length(test_poplabel) > 0) {
    nom_hap_cplet <- out_fphase[test_poplabel]
    nhap_tot <- length(nom_hap_cplet)
    pop_label <- numeric(nhap_tot)
    tmp_poplab <-
      (strsplit(nom_hap_cplet, split = "subpop. label:"))
    for (i in 1:nhap_tot) {
      pop_label[i] <-
        as.numeric(unlist(strsplit(tmp_poplab[[i]][2], split = "\\(internally"))[1])
    }
    populations <- unique(pop_label)
    if (verbose)
      cat(
        "Haplotypes in the fastPHASE output file originate from",
        length(populations),
        "populations.\n"
      )
    
    if (is.na(popsel) & length(populations) == 1) {
      popsel <- populations[1]
    }
    
    if (is.na(popsel) | !(popsel %in% pop_label)) {
      stop(
        paste0(
          "Please specify by 'popsel' one of the following population numbers:\n",
          paste(populations, collapse = " "),
          "\n",
          "Conversion stopped."
        ),
        call. = FALSE
      )
    }
    hapsel <- (which(pop_label == popsel) - 1) * 3 + 1
    hapsel <-
      sort(as.numeric(cbind(hapsel, hapsel + 1, hapsel + 2)))
    out_fphase <- out_fphase[hapsel]
  } else{
    # no sub-population identifiers found
    if (verbose) {
      cat("No population identifiers found in fastPHASE file.\n")
      if (!is.na(popsel)) {
        cat("Ignoring argument 'popsel'.\n")
      }
    }
  }
  
  ### for each individual there are 3 lines
  nind <- length(out_fphase) / 3
  ### and two haplotypes (assuming diploid individuals)
  nhap <- 2 * nind
  
  first_hap <- unlist(strsplit(out_fphase[2], split = " "))
  tmp_haplo <- matrix(as.character(NA),
                      nrow = nhap,
                      ncol = length(first_hap))
  hapnames <- vector(mode = "character", nhap)
  
  for (i in 1:nind) {
    id_line <- 3 * (i - 1) + 1
    hap1_line <- id_line + 1
    hap2_line <- id_line + 2
    hap1_index <- 2 * (i - 1) + 1
    hap2_index <- hap1_index + 1
    
    if (length(out_fphase[id_line]) > 0) {
      ind_id <- strsplit(out_fphase[id_line], split = "\\s+")[[1]][1]
      hapnames[hap1_index] <- paste(ind_id, "1", sep = "_")
      hapnames[hap2_index] <- paste(ind_id, "2", sep = "_")
    }
    
    hap1 <- unlist(strsplit(out_fphase[hap1_line], split = " "))
    tmp_haplo[hap1_index, ] <- hap1
    
    hap2 <- unlist(strsplit(out_fphase[hap2_line], split = " "))
    tmp_haplo[hap2_index, ] <- hap2
  }
  
  if (!anyDuplicated(hapnames)) {
    rownames(tmp_haplo) <- hapnames
  }
  
  return(tmp_haplo)
}

read.ms <- function(hap_file, chr.name, verbose) {
  if (!requireNamespace("gap", quietly = TRUE)) {
    stop("Package 'gap' needed to read ms output. Conversion stopped.",
         call. = FALSE)
  }
  
  if (verbose)
    cat("Input file in 'ms' output format assumed.\n")
  
  ms <- gap::read.ms.output(hap_file, verbose = verbose)
  
  if (!is.na(chr.name)) {
    chr.nbr <- suppressWarnings(as.integer(chr.name))
    if (is.na(chr.nbr)) {
      stop("For ms output files 'chr.name' has to be an integer number.",
           call. = FALSE)
    }
  } else{
    chr.nbr <- NA
  }
  
  if (ms$nreps > 1) {
    if (is.na(chr.nbr) | chr.nbr < 1 | chr.nbr > ms$nreps) {
      stop(
        paste(
          "Ms output file contains",
          ms$nreps,
          "simulations.\nPlease select one by specifying its number in 'chr.name'."
        ),
        call. = FALSE
      )
    }
  } else{
    chr.nbr <- 1
  }
  
  hh <- new("haplohh")
  hh@chr.name <- as.character(chr.nbr)
  
  if (verbose)
    cat("Extracting positions.\n")
  
  hh@positions <- ms$positions[[chr.nbr]]
  
  if (verbose)
    cat("Extracting haplotypes.\n")
  
  hh@haplo <- t(ms$gametes[[chr.nbr]])
  
  rownames(hh@haplo) <- paste0("H", seq_len(nrow(hh@haplo)))
  colnames(hh@haplo) <- paste0("s", seq_len(ncol(hh@haplo)))
  
  return(hh)
}

read.vcf <-
  function(vcf_file,
           chr.name,
           polarize_vcf,
           capitalize_AA,
           vcf_reader,
           verbose) {
    if (vcf_reader == "vcfR" &
        !requireNamespace("vcfR", quietly = TRUE)) {
      stop("Package 'vcfR' needed to read vcf files. Conversion stopped.",
           call. = FALSE)
    }
    if (vcf_reader == "data.table" &
        !requireNamespace("data.table", quietly = TRUE)) {
      stop("Package 'data.table' needed to read vcf files. Conversion stopped.",
           call. = FALSE)
    }
    
    if (verbose)
      cat("Using package '", vcf_reader, "' to read vcf.\n", sep = "")
    
    if (verbose)
      cat("Extracting map information.\n")
    
    if (vcf_reader == "vcfR") {
      vcf <- vcfR::read.vcfR(vcf_file, verbose = verbose)
      
      map <- data.frame(
        CHROM = vcfR::getCHROM(vcf),
        POS = vcfR::getPOS(vcf),
        REF = vcfR::getREF(vcf),
        ALT = vcfR::getALT(vcf),
        stringsAsFactors = FALSE
      )
    } else{
      map <- data.table::fread(
        vcf_file,
        skip = "#CHROM",
        select = c("#CHROM", "POS", "ID", "REF", "ALT", "INFO"),
        stringsAsFactors = FALSE,
        showProgress = FALSE
      )
    }
    
    chr.name <-
      check_chromosome_names(vcf_file, as.character(unique(map[, 1])), chr.name)
    selected <- map[, 1] == chr.name
    map <- map[selected, ]
    
    if (polarize_vcf) {
      if (verbose)
        cat("Extracting ancestral allele from info field of vcf file.\n")
      
      if (vcf_reader == "vcfR") {
        if (!("AA" %in% vcfR::vcf_field_names(vcf, tag = "INFO")$ID)) {
          stop("Key 'AA' not found in INFO field of vcf file. Conversion stopped.",
               call. = FALSE)
        }
        #get ancestral allele as nucleotide
        AA <- vcfR::extract.info(vcf, "AA")[selected]
      } else{
        # if key 'AA' is absent or empty, set whole string to NA
        map$INFO[!grepl("AA=[^;]+", map$INFO)] <- NA
        # extract value for key 'AA'
        AA <- sub(";.*$", "", sub("^.*AA=", "", map$INFO))
        # remove column to free memory
        map$INFO <- NULL
      }
      if (sum(is.na(AA)) > 0) {
        warning(paste(
          "Ancestral allele info field is empty for",
          sum(is.na(AA)),
          "markers."
        ))
      }
      if (capitalize_AA) {
        AA <- toupper(AA)
      }
    }
    
    if (verbose)
      cat("Extracting haplotypes.\n")
    
    if (vcf_reader == "vcfR") {
      gt <- vcfR::extract.gt(vcf)[selected, , drop = FALSE]
    } else{
      # extract GT field (always the first) from sample columns
      gt <- sub(":.*$", "", as.matrix(
        data.table::fread(
          vcf_file,
          skip = "#CHROM",
          drop = c(
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT"
          ),
          stringsAsFactors = FALSE,
          showProgress = FALSE,
        )[selected, ]
      ))
    }
    
    # check that ploidy of each individual is the same at all markers
    ind_ploidy <- apply(gt, MARGIN = 2,
                        function(x) {
                          # replace NAs by empty string
                          x[is.na(x)] <- ""
                          l <- strsplit(x, split = "[/|]")
                          t <- tabulate(lengths(l))
                          # if no information at all
                          if (all(t == 0)) {
                            stop("Cannot determine ploidy for at least one individual. Conversion stopped.",
                                 call. = FALSE)
                          }
                          # if multiple ploidies, return NA
                          ifelse(length(t[t != 0]) != 1, NA, which.max(t))
                        })
    
    if (anyNA(ind_ploidy)) {
      stop(paste(
        sum(is.na(ind_ploidy)),
        "individuals have different ploidy at different markers."
      ),
      call. = FALSE)
    }
    
    ### vcfR translates all absent genotypes, irrespective of ploidy, to single NA
    ### However, we need the correct ploidy for absent genotypes, too!
    ### Work-around: replace NA by "." , ".|." , etc. using the ploidy
    ### of the non-absent genotypes at other markers of the same individual,
    ### hence (hopefully) reconstructing the original information in the vcf file
    if (vcf_reader == "vcfR") {
      for (i in seq_len(ncol(gt))) {
        gt[is.na(gt[, i]), i] <- paste0(rep(".|", ind_ploidy[i] - 1), ".")
      }
    }
    ### end of work-around
    
    # report sample statistics
    ploidy <- tabulate(ind_ploidy)
    names(ploidy) <- seq_along(ploidy)
    if (verbose) {
      cat("Number of individuals which are \n")
      cat("Haploid Diploid Triploid, ... : \n")
      cat(names(ploidy), "\n")
      cat(ploidy, "\n")
    }
    
    # parse vcf genotypes into integers
    tmp_haplo <- matrix(apply(gt, MARGIN = 1,
                              function(x) {
                                suppressWarnings(as.integer(unlist(strsplit(x, split = "[/|]"))))
                              }),
                        ncol = nrow(gt))
    
    if (nrow(gt) > 0) {
      # set haplotype names as individual + underscore + 1:ploidy
      rownames(tmp_haplo) <-
        as.vector(unlist(
          mapply(function(x, y) {
            if (y == 1) {
              # haploid -> return un-changed
              x
            } else{
              # return vector c("HG1_1",...,"HG1_y")
              paste(x, 1:y, sep = "_")
            }
          }, colnames(gt), ind_ploidy, USE.NAMES = FALSE),
          use.names = FALSE
        ))
    }
    
    hh <- new("haplohh")
    
    hh@chr.name <- chr.name
    hh@haplo <- tmp_haplo
    hh@positions <- as.numeric(map[, 2])
    
    if (vcf_reader == "vcfR") {
      mrk.names <- vcfR::getID(vcf)[selected]
    } else{
      # replace point by NA
      mrk.names <- sub("\\.", NA, map$ID)
      # check on duplicates (done automatically by vcfR)
      if (anyDuplicated(na.omit(mrk.names))) {
        stop("ID column contains non-unique names.", call. = FALSE)
      }
    }
    
    if (!anyNA(mrk.names)) {
      colnames(hh@haplo) <- mrk.names
      names(hh@positions) <- mrk.names
    } else{
      if (verbose) {
        if (sum(is.na(mrk.names)) == ncol(hh@haplo)) {
          cat("No marker identifiers found in vcf file.\n")
        } else{
          cat(
            sum(is.na(mrk.names)),
            "out of",
            ncol(hh@haplo),
            "markers have no identifier in vcf file.\n"
          )
        }
      }
    }
    
    # polarize
    if (exists("AA")) {
      if (verbose)
        cat("Polarizing variants.\n")
      
      allele_list <-
        strsplit(paste(map[, "REF"], map[, "ALT"], sep = ","), ",", fixed = TRUE)
      
      #if ancestral allele among REF or ALT, get number, otherwise zero
      aan <-
        mapply(match, AA, allele_list, USE.NAMES = FALSE) - 1L
      
      if (sum(is.na(aan)) == ncol(hh@haplo)) {
        stop("No marker could be polarized. Conversion stopped.", call. = FALSE)
      }
      
      #switch allele coding 0 with aan, if aan is not zero.
      #if the ancestral allele is not known/does not match, this yield NA
      hh@haplo <- t(apply(hh@haplo, 1, function(x) {
        aan * (x == 0L) + x * (aan == 0L | (aan > 0L & x != aan))
      }))
      
      # if only one marker then matrix must be explicitly coerxed
      if (nrow(hh@haplo) == 1) {
        hh@haplo <- as.matrix(hh@haplo[, !is.na(aan)])
      } else{
        hh@haplo <- hh@haplo[, !is.na(aan)]
      }
      hh@positions <- hh@positions[!is.na(aan)]
      
      if (verbose) {
        cat(sum(!is.na(aan)), "markers have been polarized.\n")
        cat(sum(is.na(aan)),
            "unpolarized markers have been removed.\n")
      }
    }
    
    return(hh)
  }
