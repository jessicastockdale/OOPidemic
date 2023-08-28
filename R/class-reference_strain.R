#' ReferenceStrain
#' 
#' @description An R6 class representing a reference strain for a pathogen. The purpose of this class is to act as a reference for a Strain's genome. It also is where the mutation rate for Strains is defined.
#' 
#' @examples
#' ReferenceStrain$new(name = "ref_strain", g_len = 1000, mut_rate = 0.00001)
#' ReferenceStrain$new(name = "ref_strain", genome = c("A", "C", "G", "T"))
#' ReferenceStrain$new(name = "ref_strain", genome = c("A", "C", "G", "U"), dna = FALSE)
#' @export
#' 
ReferenceStrain <- R6::R6Class("ReferenceStrain",

    lock_class = TRUE,

    ###########################################################################
    # public members

    public = list(

        #######################################################################
        # public variables
        
        #' @field name Name of this reference strain
        name = NULL, 

        #######################################################################
        # public functions

        #' @description
        #' Create a new reference strain object.
        #' 
        #' @param name Name of this reference strain.
        #' @param genome A character vector composed of "A", "C", "G", & "T" or "U". 
        #'  Can be excluded if a randomly generated genome is acceptable.
        #' @param g_len Length of randomly generate genome. Not required if a genoms is provided
        #' @param mut_rate The probability of a single nucleotide polymorphism per unit time. 
        #'  Mutations follow a Jukes-Cantor model.
        #' @param dna Bool indicating if the genome is DNA or RNA based. 
        #'  'FALSE' indicates it's RNA based.
        #' 
        #' @return A new `ReferenceStrain` object
        initialize = function(
            name,
            genome = NULL,
            g_len = 1000L,
            mut_rate = 0.000008,
            dna = TRUE
        ) {

            # parameter validations
            stopifnot(is.character(name), length(name) == 1)
            stopifnot(any(
                is.null(genome),
                all(is.character(genome), length(genome) > 0)
            ))
            stopifnot(
                is.numeric(g_len), g_len %% 1 == 0,
                g_len > 0, length(g_len) == 1
            )
            stopifnot(
                is.numeric(mut_rate), length(mut_rate) == 1,
                0 <= mut_rate, mut_rate <= 1    
            )
            stopifnot(is.logical(dna), length(dna) == 1)

            if (is.null(genome)) {
                # no genome supplied. Generate a random one
                private$g_len_ <- g_len
                private$genome_ <- private$randomise_genome(g_len)

            } else {
                # genome provided:

                # check if genome is DNA or RNA
                if (dna & (sum(genome == "U") > 0)) {
                    stop("`genome` appears to be RNA. Set `dna` to FALSE.")
                } else if (!dna & (sum(genome == "T") > 0)) {
                    stop("`genome` appears to be DNA. Set `dna` to TRUE.")
                }

                # convert genome to integers and save genome length
                private$genome_ <- genome_to_integers(genome, dna)
                private$g_len_ <- length(genome)

            }

            self$name <- name

            private$dna_ <- dna
            private$mut_rate_ <- mut_rate
        },

        #' @description
        #' Get reference strain's genome. 
        #' 
        #' @param as_character Logical indicating if the genome should be returned as a character or integer vector.
        #'  Default is TRUE. 
        #'  'FALSE' returns an integer vector.
        #' 
        #' @return A character or integer vector. If an integer vector is specified then 1 = "A", 2 = "C", 3 = "G", and 4 = "T" if genome is DNA based or "U" if genome is RNA based.
        genome = function(as_character = TRUE) {

            stopifnot(is.logical(as_character), length(as_character) == 1)

            genome <- private$genome_

            if (as_character) {
                # convert to character if required
                genome <- genome_to_characters(genome, private$dna_)
            }

            return(invisible(genome))
        },

        #' @description 
        #' Get the nucleotides at the specified loci
        #' 
        #' @param loci A integer vector specifiying the loci in the genome that you want the nucleotides for.
        #' @param as_character Logical indicating if the nucleotides should be returned as a character or integer vector.
        #'  Default is FALSE. 
        #'  'FALSE' returns an integer vector.
        #' 
        #' @return An integer vector representation of the desired nucleotides
        #' 
        #' @examples 
        #' ref_strain <- ReferenceStrain$new(name = "ref_strain")
        #' ref_strain$nucleotides_at_loci(c(3, 10, 15))
        #' ref_strain$nucleotides_at_loci(c(3, 10, 15), TRUE)
        nucleotides_at_loci = function(loci, as_character = FALSE) {

            stopifnot(
                is.numeric(loci), all(loci %% 1 == 0),
                length(loci) >= 1
            )
            stopifnot(1 <= min(loci), max(loci) <= private$g_len_)
            stopifnot(is.logical(as_character), length(as_character) == 1)

            nucs <- private$genome_[loci]

            if (as_character) {
                nucs <- genome_to_characters(nucs, private$dna_)
            }

            return(nucs)
        },

        #' @description
        #' Print a description of the ReferenceStrain object
        #' 
        #' @param ... arguments will be ignored
        #' 
        #' @return No return value. Description is printed to the console
        print = function(...) {
            cat("ReferenceStrain: \n")
            cat("   Name:   ", self$name, "\n", sep = "")
            cat("   Length: ", private$g_len_, "\n", sep = "")
            cat("   Genome: ", self$genome(TRUE), sep = "", fill = TRUE)
            cat("   DNA:    ", private$dna_, "\n", sep = "")
        }
    ),

    ###########################################################################
    # active bindings

    active = list(

        #' @field g_len The length of the reference strain's genome
        g_len = function(v) {
            # return length of genome
            if (missing(v)) private$g_len_
            else stop("Can't set `$g_len`", call. = FALSE)
        },

        #' @field is_dna Bool indicating if the genome is DNA based ('FALSE' indicates it is RNA based)
        is_dna = function(v) {
            # return a bool identifying if this reference strain is dna
            if (missing(v)) private$dna_
            else stop("Can't set `$is_dna`.", call. = FALSE)
        },

        #' @field mut_rate The probability of a single nucleotide polymorphism per unit time. 
        #'  Mutations follow a Jukes-Cantor model.
        mut_rate = function(v) {
            # return a bool identifying if this reference strain is dna
            if (missing(v)) private$mut_rate_
            else stop("Can't set `$mut_rate`.", call. = FALSE)
        }
    ),

    ###########################################################################
    # private members

    private = list(

        #######################################################################
        # private variables
        genome_ = NULL,
        g_len_ = NA,
        dna_ = NULL,
        mut_rate_ = NA,

        #######################################################################
        # private functions

        # @description
        # Generate a random genome
        # 
        # @param g_len Integer indicating the desired genome length
        # 
        # @return An integer vector representation of a random genome
        randomise_genome = function(g_len) {

            sample(1:4, g_len, replace = TRUE)
        }
    ),
)

