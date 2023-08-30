#' Lab
#' 
#' @description An R6 class representing a lab. A lab handles taking WGSequences from hosts that are due to be sampled. A user will need to create this object but will usually delegate it's control to either a `Group` or `Population` object
#' 
#' @examples
#' #Lab$new()
#' 
#' @export
#' 
Lab <- R6::R6Class("Lab",

    cloneable = FALSE,

    ###########################################################################
    # public members

    public = list(

        #######################################################################
        # public variables

        #######################################################################
        # public functions

        #' @description
        #' Create a new `Lab` object. 
        #' 
        #' @param wgs_class The class used to create whole genome sequences. Default is `WGSequence` though any R6Class that inherits from `WGSequence` can be used here
        #' 
        #' @return A `Lab` object
        #' 
        initialize = function(
            wgs_class = WGSequence
        ) {
            stopifnot(is_ancestor_r6class(WGSequence, wgs_class))

            private$wgs_class_ <- wgs_class
            return(invisible(self))
        },

        #' @description
        #' Print a description of this `Lab` object
        #' 
        #' @param ... arguments will be ignored
        #' 
        #' @return No return value. Description is printed to the console
        #' 
        print = function(...) {
            cat("Lab: \n")
            cat("   Id:             ", private$id_, "\n", sep = "")
            cat("   WG Sequences:   ", self$num_wgs, "\n", sep = "")

        },

        #' @description
        #' Build a dataframe of sequence names and whole genome sequences for saving to a fasta file
        #' 
        #' @return A dataframe
        #' 
        fasta_df = function() {

            fasta_data <- lapply(
                private$wg_sequences_,
                function(wgs) {
                    list(
                        name = wgs$name, 
                        genome = paste(wgs$genome(), collapse = "")
                    )
                }
            )
            fasta_df <- do.call(rbind, lapply(fasta_data, data.frame))

            return(fasta_df)
        },

        #' @description
        #' Build a dataframe with details about hosts that have been sampled
        #' 
        #' @return A dataframe
        #' 
        hostdata_df = function() {
            hostdatas <- lapply(private$wg_sequences_, function(wgs) wgs$host$data)
            hostdata_df <- do.call(rbind, lapply(hostdatas, data.frame))

            # remove duplicates - can happen if we are doing multiple samples
            hostdata_df <- hostdata_df[!duplicated(hostdata_df), ]

            return(hostdata_df)
        },

        #' @description
        #' Collect samples from `host`s
        #' 
        #' @param hosts A vector of `Host` objects to be sampled (or R6Classes that inherit from `Host`)
        #' @param time The time that samples are being collected
        #' 
        #' @return This `Lab` object
        #' 
        #' @examples
        #' ref_strain <- ReferenceStrain$new("ref_strain", g_len = 100)
        #' grp <- Group$new(
        #'     1, ref_strain,
        #'     init_inf = 6,
        #'     inc_shape = 0,
        #'     sample_schedule = "calendar", sample_freq = 1
        #' )
        #' lab <- Lab$new()
        #' lab$sample_hosts(grp$infectious_hosts(1), 1)
        #' 
        sample_hosts = function(hosts, time) {

            stopifnot(is.numeric(time), time %% 1 == 0, length(time) == 1)

            if (length(hosts) > 0) {
                vapply(
                    hosts, 
                    private$sample_host, 
                    logical(1),
                    time = time
                    )
            }

            return(invisible(self))
        },

        #' @description
        #' Build a dataframe with metadata for the whole genome sequences
        #' 
        #' @return A dataframe
        #' 
        metadata_df = function() {
            metadatas <- lapply(private$wg_sequences_, function(wgs) wgs$metadata)
            metadata_df <- do.call(rbind, lapply(metadatas, data.frame))

            return(metadata_df)
        },

        #' @description
        #' Saves host data and sample metadata as csvs and genomes as a fasta file
        #' 
        #' @details
        #' This function creates three files: 
        #' * `<directory>/<stem>_hostdata.csv` and 
        #' * `<directory>/<stem>_sample_metadata.csv` and 
        #' * `<directory>/<stem>_samples.fasta`
        #' 
        #' @param directory The target directory (the directory will be created if it doesn't exist, no warnings will be shown)
        #' @param stem the extensionless form of the file name 
        #' 
        #' @return A dataframe
        #' 
        save = function(directory, stem) {
            
            # check that the directory exists (create it if it doesn't)
            dir.create(directory, recursive = TRUE, showWarnings = FALSE)

            # generate filepaths
            hostdata_filepath <- file.path(directory, paste0(stem, "_hostdata.csv"))
            metadata_filepath <- file.path(directory, paste0(stem, "_wgs_metadata.csv"))
            fasta_filepath <- file.path(directory, paste0(stem, "_wgs.fasta"))
            
            # save genomes
            fasta_df <- self$fasta_df()
            fastafile <- write.fasta(
                as.list(fasta_df$genome),
                fasta_df$name,
                nbchar = 80,
                file.out = fasta_filepath,
                open = "w"
            )

            # save metadata and host data
            write.csv(self$metadata_df(), metadata_filepath)
            write.csv(self$hostdata_df(), hostdata_filepath)
           
            return(invisible(self))
        }

    ),

    ###########################################################################
    # active bindings

    active = list(

        #' @field num_wgs The number of whole genome sequences collected so far
        num_wgs = function(v) {
            if (missing(v)) length(private$wg_sequences_)
            else stop("Can't set `$num_wgs`", call. = FALSE)
        },
        
        #' @field wg_sequences A vector of all whole gnome sequences collected so far
        wg_sequences = function(v) {
            if (missing(v)) private$wg_sequences_
            else stop("Can't set `$wg_sequences`", call. = FALSE)
        }
    ),

    ###########################################################################
    # private members

    private = list(

        #######################################################################
        # private variables 
        wg_sequences_ = c(),
        wgs_class_ = NULL,

        #######################################################################
        # private functions

        # @description
        # Decide which strains carried by a host are to be sampled.
        # 
        # @param hosts A vector of `Host` objects to be sampled (or R6Classes that inherit from `Host`)
        # @param time The time that samples are being collected
        # 
        # @return This `Lab` object
        # 
        pick_strains = function(host) {
            # decide which strains are to be sampled

            # if multiple strains in a host is implemented then insert code to decide which strains get sampled here

            return(host$strains)
        },

        # @description
        # Collect samples from `host`
        # 
        # @param host A `Host` object to be sampled (or a R6Class that inherits from `Host`)
        # @param time The time that samples are being collected
        # 
        # @return TRUE: The motivation for this is that we can call this function for many strains with vapply if we know what the return value will be, so we've made it something small
        # 
        sample_host = function(host, time) {
            
            stopifnot(
                "Host" %in% class(host),
                host$is_sampling_due(time)
            )

            # first have host realise any changes in themselves and update their sample time if needed
            host$realise(time)
            host$update_sample_time()
            
            # get strains to be sampled
            strains <- private$pick_strains(host)

            # trigger strains to realise themselves
            # use unique because if we implement multiple strains infecting a host there be could be duplicates
            vapply(
                unique(strains),
                function(strain) strain$realise(time),
                logical(1)
            )

            new_wg_sequences <- mapply(
                private$wgs_class_$new,
                seq_along(strains),
                strains,
                MoreArgs = list(host = host)
            )

            private$wg_sequences_ <- c(private$wg_sequences_, new_wg_sequences)

            return(invisible(TRUE)) # return something small so we can use vapply
        }
    )
)

