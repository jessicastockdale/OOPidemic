#' Strain
#' 
#' @description An R6 class representing a strain of a pathogen.
#' 
#' @examples
#' ref_strain <- ReferenceStrain$new(name = "ref_strain")
#' Strain$new(ref_strain)
#' Strain$new(ref_strain, time = 3)
#' Strain$new(ref_strain, distance = 4)
#' @export
#' 
Strain <- R6::R6Class("Strain",

    ###########################################################################
    # public members

    public = list(

        #######################################################################
        # public variables

        #######################################################################
        # public functions

        #' @description
        #' Create a new strain object. Generally only used when creating strains for index cases, because of this it's unlikely that a user will create these explicitly.
        #' 
        #' @param ref_strain A ReferenceStrain object (or a R6Class object that inherits from ReferenceStrain)
        #' @param time The time that this strain is 'spawned'
        #' @param distance This strains genomic distance from 'ref_strain'. Genomic distance is defined as the number of single nucleotide polymorphisms of this strain compared to it's reference strain.
        #' 
        #' @return A new 'Strain' object
        #' 
        initialize = function(
            ref_strain,
            time = 0L, 
            distance = 0 
        ) {

            # parameter validations
            stopifnot("ReferenceStrain" %in% class(ref_strain))
            stopifnot(
                is.numeric(time), time %% 1 == 0, 
                length(time) == 1, time >= 0
            )
            stopifnot(
                is.numeric(distance), distance %% 1 == 0, 
                length(distance) == 1, distance >= 0
            )
            # store variables
            private$ref_strain_ <- ref_strain
            private$spawn_time_ <- time 
            private$realisation_time_ <- time 

            # mutate distance number of nucleotides for variation
            if (distance > 0) {
                private$mutate(distance)
            }
            
            return(self)
        },

        #' @description 
        #' Get strain's genome
        #' 
        #' @param as_character Logical indicating if the genome should be returned as a character or integer vector.
        #'  Default is TRUE. 
        #'  'FALSE' returns an integer vector.
        #' 
        #' @return A character or integer vector. If an integer vector is specified then 1 = "A", 2 = "C", 3 = "G", and 4 = "T" if genome is DNA based or "U" if genome is RNA based.
        #' 
        #' @examples
        #' ref_strain <- ReferenceStrain$new(name = "ref_strain")
        #' strain <- Strain$new(ref_strain)
        #' strain$genome()
        #' strain$genome(as_character = FALSE)
        genome = function(as_character = TRUE) {

            stopifnot(is.logical(as_character), length(as_character) == 1)

            genome <- private$ref_strain_$genome(FALSE)
            genome[private$loci_] <- private$nucleotides_

            if (as_character) {
                genome <- genome_to_characters(genome, self$is_dna)
            }
            
            return(genome)
        },

        #' @description
        #' Print a description of the `ReferenceStrain` object
        #' 
        #' @param ... arguments will be ignored
        #' 
        #' @return No return value. Description is printed to the console
        #' 
        print = function(...) {
            cat("Strain: \n")
            cat("   Reference Strain:       ", private$ref_strain_$name, "\n", sep = "")
            cat("   Spawn Time :            ", private$spawn_time_, "\n", sep = "")
            cat("   Realisation Time:       ", private$realisation_time_, sep = "", fill = TRUE)
            cat("   Mutation Loci:          ", paste0(private$loci_, collapse = ", "), sep = "", fill = TRUE)
            cat("   Mutation Nucleotides:   ", paste0(private$nucleotides_, collapse = ", "), "\n", sep = "", fill = TRUE)
        },

        #' @description
        #' Realise changes that have occured (ie mutations) that have occured since the last realisation time
        #' 
        #' @param time The current time
        #' 
        #' @return TRUE: The motivation for this is that we can call this function for many strains with vapply if we know what the return value will be, so we've made it something small
        #' 
        realise = function(time) {

            stopifnot(
                length(time) == 1,
                is.numeric(time), time %% 1 == 0  
            )

            # has any time elapsed since the last realisation_time_?
            if (time != private$realisation_time_) {
                num_subs <- private$substitution_number(time)

                # did any mutations happen?
                if (num_subs > 0) private$mutate(num_subs)

                private$realisation_time_ <- time
            } 

            invisible(TRUE)
        },

        #' @description
        #' Trigger this strain to reproduce and spawn a descendent
        #' 
        #' @param time The current time
        #' 
        #' @return This object's descendent
        #' 
        #' @examples
        #' ref_strain <- ReferenceStrain$new(name = "ref_strain")
        #' strain <- Strain$new(ref_strain)
        #' strain$reproduce(time = 12)
        reproduce = function(time) {

            stopifnot(
                length(time) == 1,
                is.numeric(time), time %% 1 == 0  
            )

            # first realise self
            self$realise(time)

            # clone self
            descendent <- self$clone(deep = TRUE)

            # record descendent
            private$descendents_ <- c(private$descendents_, descendent)

            invisible(descendent)
        }
    ),

    ###########################################################################
    # active bindings

    active = list(

        #' @field ancestor This strain's immediate ancestor
        ancestor = function(v) {
            if (missing(v)) private$ancestor_
            else stop("Can't set `$ancestor`", call. = FALSE)
        },

        #' @field descendents A vector of this strain's immediate descendents
        descendents = function(v) {
            if (missing(v)) private$descendents_
            else stop("Can't set `$descendents`", call. = FALSE)
        },

        #' @field is_dna A logical indicating if this strain has a DNA or RNA based genome (TRUE and FALSE respectively)
        is_dna = function(v) {
            if (missing(v)) private$ref_strain_$is_dna
            else stop("Can't set `$nucleotides`", call. = FALSE)
        },

        #' @field loci A vector of loci where single nucleotide polymorphisms have occured (unsorted)
        loci = function(v) {
            if (missing(v)) private$loci_
            else stop("Can't set `$loci`", call. = FALSE)
        },

        #' @field nucleotides The nucleotides that are present at the `loci` with single nucleotide polymorphisms
        nucleotides = function(v) {
            if (missing(v)) private$nucleotides_
            else stop("Can't set `$nucleotides`", call. = FALSE)
        },

        #' @field realisation_time The most recent realisation time for this strain
        realisation_time = function(v) {
            if (missing(v)) private$realisation_time_
            else stop("Can't set `$realisation_time`", call. = FALSE)
        },

        #' @field ref_strain The ReferenceStrain object for this strain
        ref_strain = function(v) {
            if (missing(v)) private$ref_strain_
            else stop("Can't set `$ref_strain`", call. = FALSE)
        },
        
        #' @field spawn_time The time that this strain was spawned
        spawn_time = function(v) {
            if (missing(v)) private$spawn_time_
            else stop("Can't set `$spawn_time`", call. = FALSE)
        }
        
    ),

    ###########################################################################
    # private members

    private = list(

        #######################################################################
        # private variables 
        ancestor_ = NULL,
        descendents_ = c(),
        loci_ = c(), 
        nucleotides_ = c(),
        realisation_time_ = 0, 
        ref_strain_ = NULL,
        spawn_time_ = NA,

        #######################################################################
        # private functions

        # @description 
        # The deep_clone function for this Class. When Strain$clone(deep = True) is called, deep_clone will be called once per field. We want to catch the spawn_time_ and ancestor fields as these will be different for the descendent
        deep_clone = function(name, value) {
            # offsprings spawn time is the current realisation time
            if (name == "spawn_time_") {private$realisation_time_}
            # offsprings ancestor is this strain
            else if (name == "ancestor_") {self}
            # insert additional modifications to descendent here if needed
            else {value}
        },

        # @description
        # Have this strain mutate the specified distance
        # 
        # @param distance An integer specifying the distance the genome should be from it's current state (distance = new SNPs)
        mutate = function(distance) {

            # sample distance number of loci 
            loci <- sample(private$ref_strain_$g_len, distance) 

            # loci that are different from reference strain
            prev_mut_loci_index <- which(private$loci_ %in% loci)

            # mutate and store previously mutated loci (if any)
            if (length(prev_mut_loci_index) != 0) {
                # get the nucleotides for those loci
                prev_nucs <- private$nucleotides_[prev_mut_loci_index]
                # get new nucleotides
                replacement_nucs <- private$mutate_nucleotides(prev_nucs)
                # store new nucleotides
                private$nucleotides_[prev_mut_loci_index] <- replacement_nucs

            }

            # loci that are still the same as the ref strain
            new_mut_loci <- loci[!(loci %in% private$loci_)]

            # mutate and store new mutated loci. if any
            if (length(new_mut_loci) != 0) {
                # get the nucleotides for those loci
                ref_nucs <- private$ref_strain_$nucleotides_at_loci(new_mut_loci)
                
                # get new nucleotides
                new_nucs <- private$mutate_nucleotides(ref_nucs)

                # store new loci and nucleotides
                private$loci_ <- c(private$loci_, new_mut_loci)
                private$nucleotides_ <- c(private$nucleotides_, new_nucs)

            }

            invisible(self)
        },

        # @description
        # Mutate the specified nucleotides
        # 
        # @param old_nucs An integer vector specifying the current nucleotides to be mutated
        # 
        # @return An integer vector of nucleotides that will replace 'old_nucs'
        # 
        mutate_nucleotides = function(old_nucs) {

            if (length(old_nucs) == 0) {
                return(c())
            }

            # generate a 4xlength(old_nucs) array where each column is the integers 1 to 4
            # then add an extra row with the old nucleotides
            potential_nucs <- rbind(array(1:4, c(4, length(old_nucs))), old_nucs)

            # remove elements in each column that match the last entry (old nucleotide)
            potential_nucs <- apply(potential_nucs, 2, function(c) c[c != c[length(c)]])

            # select one entry from each column
            new_nucs <- apply(potential_nucs, 2, function(c) c[sample(1:3, 1)])

            return(new_nucs)
        },

        # @description
        # Get number of substitutions that occured between time and private$realisation_time_
        # 
        # @param time The current time
        # 
        # @return An integer specifying the number of substitutions
        # 
        substitution_number = function(time) {
            # scale mutation rate by number of days and genome length

            mean <- (
                private$ref_strain_$g_len * 
                private$ref_strain_$mut_rate * 
                (time - private$realisation_time_)
            )
            # the number of substitutions is binomially(n, p) distributed 
            # because we expect n to be large (n > 20) and p to be small we can the poison to approximate the number of substitutions
            num_subs <- stats::rpois(1, mean)

            return(num_subs)
        }
    )    
)
