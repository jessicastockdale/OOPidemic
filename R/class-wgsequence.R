#' WGSequence
#' 
#' @description An R6 class representing a whole genome sequence test. Users are not expected to interact directly with this class but it is made available so that a user can subclass it if desired.
#' 
#' @export
#' 
#' @examples 
#' ref_strain <- ReferenceStrain$new("ref_strain")
#' grp <- Group$new(1, ref_strain)
#' host <- grp$infectious_hosts()[[1]]
#' strain <- host$strains[[1]]
#' wgs <- WGSequence$new(1, strain, host)
#' 
WGSequence <- R6::R6Class("WGSequence",

    cloneable = FALSE, 

    ###########################################################################
    # public members

    public = list(

        #######################################################################
        # public variables

        #######################################################################
        # public functions

        #' @description
        #' Create a new wgsequence object. Generally, a user won't be required to do this explicitly. It will be handled by a `Lab` object
        #' 
        #' @param sample_num Integer representing this sample's number for this `host` at this time (time is retrieved from `strain`)
        #' @param strain The `Strain` object that is being sampled (or a R6Class that inherits from Strain)
        #' @param host The `Host` object that is carrying the `strain` (or a R6Class that inherits from Host)
        #' 
        #' @return A new `WGSequence` object
        #' 
        initialize = function(
            sample_num,
            strain,
            host
        ) {

            # parameter validations
            stopifnot(
                is.numeric(sample_num), sample_num %% 1 == 0, 
                length(sample_num) == 1, sample_num > 0
            )
            stopifnot("Strain" %in% class(strain))
            stopifnot("Host" %in% class(host))


            # note that strains have already realised themselves by this point

            # store variables
            private$host_ <- host
            private$loci_ <- strain$loci 
            private$nucleotides_ <- strain$nucleotides
            private$ref_strain_ <- strain$ref_strain
            private$sample_num_ <- sample_num
            private$strain_ <- strain
            private$time_ <- strain$realisation_time

            invisible(self)
        },

        #' @description 
        #' Get genome stored in this wgsequence
        #' 
        #' @param as_character Logical indicating if the genome should be returned as a character or integer vector.
        #'  Default is TRUE. 
        #'  'FALSE' returns an integer vector.
        #' 
        #' @return A character or integer vector. If an integer vector is specified then 1 = "A", 2 = "C", 3 = "G", and 4 = "T" if genome is DNA based or "U" if genome is RNA based.
        #' 
        genome = function(as_character = TRUE) {
            # returns genome as an integer or character vector

            genome <- private$ref_strain_$genome(FALSE)
            genome[private$loci_] <- private$nucleotides_

            # convert to characters if wanted
            if (as_character) {
                genome <- genome_to_characters(genome, self$is_dna)
            }
            
            return(genome)
        },

        #' @description
        #' Print a description of the `WGSequence` object
        #' 
        #' @param ... arguments will be ignored
        #' 
        #' @return No return value. Description is printed to the console
        #' 
        print = function(...) {
            cat("WGSequence: \n")
            cat("   Host:                   ", private$host_$id, "\n", sep = "")
            cat("   Reference Strain:       ", private$ref_strain_$name, "\n", sep = "")
            cat("   Sample Time :           ", private$time_, "\n", sep = "")
            cat("   Mutation Loci:          ", paste0(private$loci_, collapse = ", "), sep = "", fill = TRUE)
            cat("   Mutation Nucleotides:   ", paste0(private$nucleotides_, collapse = ", "), "\n", sep = "", fill = TRUE)

        }

    ),

    ###########################################################################
    # active bindings

    active = list(

        #' @field host The `Host` object carrying the `Strain` object that this sequence was taken from
        host = function(v) {
            if (missing(v)) private$host_
            else stop("Can't set `$host`", call. = FALSE)
        },

        #' @field loci A vector of loci in `strain`'s `ReferenceStrain` genome where single nucleotide polymorphisms had occured at time of sampling (unsorted) 
        loci = function(v) {
            if (missing(v)) private$loci_
            else stop("Can't set `$loci`", call. = FALSE)
        },

        #' @field is_dna A logical indicating if this sequence is DNA or RNA based (TRUE and FALSE respectively)
        is_dna = function(v) {
            if (missing(v)) private$ref_strain_$is_dna
            else stop("Can't set `$is_dna`", call. = FALSE)
        },

        #' @field metadata Returns a named list of data about this sequence: `name`, `host`'s group id, `hosts`'s id, `ref_strain`'s name, `collection_time`, and `sample_number` 
        metadata = function(v) {
            if (missing(v)) {
                list(
                    name = self$name,
                    group = private$host_$group$id,
                    host = private$host_$id,
                    ref_strain = private$ref_strain_$name,
                    collection_time = private$time_,
                    sample_number = private$sample_num_
                )
            } else stop("Can't set `$metadata`", call. = FALSE)
        },

        #' @field name String in the form `<group_id>_<host_id>_<time>_<sample_number>` padded to have up to 1, 3, 3, and 1 leading zeros respectively
        name = function(v) {
            if (missing(v)) {
                paste0(
                    sprintf("%02d", private$host_$group$id), "_",
                    sprintf("%04d", private$host_$id), "_",
                    sprintf("%04d", private$time_), "_",
                    sprintf("%02d", private$sample_num_)
                )
            } else stop("Can't set `$name`", call. = FALSE)
        },

        #' @field nucleotides The nucleotides that are present at the `loci` with single nucleotide polymorphisms
        nucleotides = function(v) {
            if (missing(v)) private$nucleotides_
            else stop("Can't set `$nucleotides`", call. = FALSE)
        },

        #' @field ref_strain The `ReferenceStrain` object for `strain`
        ref_strain = function(v) {
            if (missing(v)) private$ref_strain_
            else stop("Can't set `$ref_strain`", call. = FALSE)
        },
        
        #' @field sample_num Integer representing this sample's number for `host` at `time`
        sample_num = function(v) {
            # return the sample number of this test for it's host
            if (missing(v)) private$sample_num_
            else stop("Can't set `$sample_num`", call. = FALSE)
        },
        
        #' @field strain `Strain` object that was sequenced
        strain = function(v) {
            if (missing(v)) private$strain_
            else stop("Can't set `$strain`", call. = FALSE)
        },
        
        #' @field time When this sequence was recorded
        time = function(v) {
            if (missing(v)) private$time_
            else stop("Can't set `$time`", call. = FALSE)
        }
        
    ),

    ###########################################################################
    # private members

    private = list(

        #######################################################################
        # private variables 
        host_ = NULL,
        loci_ = c(),
        nucleotides_ = c(),
        ref_strain_ = NULL,
        sample_num_ = NA,
        strain_ = NULL,
        time_ = NA

        #######################################################################
        # private functions

    ),    
)
