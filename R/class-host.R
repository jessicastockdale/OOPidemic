#' Host
#' 
#' @description An R6 class representing a host
#' 
#' @examples
#' ref_strain <- ReferenceStrain$new(name = 'ref_strain')
#' pop <- Population$new(id = 1, ref_strain = ref_strain)
#' Host$new(id = 2, population = pop)
#' @export
#' 
Host <- R6::R6Class("Host",

    cloneable = FALSE,
    
    ###########################################################################
    # public members

    public = list(

        #######################################################################
        # public variables

        #######################################################################
        # public functions

        #' @description
        #' Create a new host object. Generally, a user won't be required to do this explicitly. It will be handled by a `Population` object
        #' 
        #' @param id The id of this host. Must be unique within each population, otherwise, this will lead to conflicting lab results.
        #' @param population The Population object that this host belongs to (or a R6Class object that inherits from Population). 
        #' 
        #' @return A new `Host` object
        #' 
        initialize = function(
            id,
            population
        ) {
            
            # parameter validations
            stopifnot(
                is.numeric(id), id %% 1 == 0,
                id > 0, length(id) == 1
            )
            stopifnot("Population" %in% class(population))

            private$id_ <- id
            private$population_ <- population

            invisible(self)
        },

        #' @description
        #' Have host contract strains of a pathogen from an infector. The infector chooses which strains are passed on. Note that the S(E)IR times for this host should have already beeen decided for this host by calling prepare_for_infection()
        #' 
        #' @param infector The `Host` object that is infecting this `host`
        #' @param strains A list of `Strain` Objects
        #' @param freq The frequencies of each `Strain` object in `strains`
        #' 
        #' @return Returns this host
        #' 
        contract_strains = function(
            infector, 
            strains, 
            freq
        ) {
            if (!is.null(infector)) {
                stopifnot(
                    "Host" %in% class(infector),
                    infector$is_infectious(private$population_$time)
                )
            }
            stopifnot(all(vapply(
                strains, 
                function(strain) "Strain" %in% class(strain), 
                logical(1L)
            )))
            stopifnot(is.numeric(freq), all(freq %% 1 == 0))
            stopifnot(length(strains) == length(freq))

            # make sure host has had prepare_for_infection run
            if (is.na(private$infectious_time_)) {
                stop("Attempt to give strain to a host that hasn't been prepared for infection.")
            } 

            # make sure host hasn't already been infected
            if (length(private$strains_) != 0) {
                stop("Attempt to infect a non-susceptible host")
            }

            private$infector_ <- infector
            private$strains_ <- strains 
            private$freq_ <- freq
            
            return(invisible(self))
        },


        #' @description
        #' Have host infect an infectee
        #' 
        #' @param infectee The `Host` object that is being infected by this `host`. The infectee has to have been prepared for infection prior to this function being called
        #' @param time The time of the infection
        #' 
        #' @return Returns this `Host` object
        #' 
        infect = function(infectee, time) {
            
            stopifnot(self$is_infectious(private$population_$time))
            stopifnot(
                "Host" %in% class(infectee)
            )
            if (is.na(infectee$exposure_time)) {
                stop("Infectee has not been prepared for infection")
            } else if (length(infectee$strains) != 0) {
                stop("Infectee is not susceptible...")
            }
            stopifnot(is.numeric(time), time %% 1 == 0, length(time) == 1)

            # realise changes in self
            self$realise(time)

            # clone strains to pass
            strains_to_pass <- as.vector(lapply(
                private$strains_, 
                function(strain) strain$reproduce(time)
            ))

            # perform infection
            infectee$contract_strains(self, strains_to_pass, c(1))

            # record infectee
            private$infectees_ <- c(private$infectees_, infectee)

            return(invisible(self))
        },

        #' @description
        #' Test if host is exposed at a given time
        #' 
        #' @param time The time being queried
        #' 
        #' @return Logical indicating if the host was exposed at `time`
        #' 
        is_exposed = function(time) {
            stopifnot(is.numeric(time), length(time) == 1)

            if (length(private$strains_) == 0) FALSE 
            else all(
                private$exposure_time_ <= time,
                time < private$infectious_time_ 
            )
        },

        #' @description
        #' Test if host is infectious at a given time
        #' 
        #' @param time The time being queried
        #' 
        #' @return Logical indicating if the host was infectious at `time`
        #' 
        is_infectious = function(time) {
            stopifnot(is.numeric(time), length(time) == 1)

            if (length(private$strains_) == 0) FALSE
            else all(
                private$infectious_time_ <= time,
                time < private$recovery_time_ 
            )
        },

        #' @description
        #' Test if host is recovered at a given time
        #' 
        #' @param time The time being queried
        #' 
        #' @return Logical 
        #' 
        is_recovered = function(time) {
            stopifnot(is.numeric(time), length(time) == 1)

            # return a bool indicating if hosts is infected
            if (length(private$strains_) == 0) FALSE 
            else private$recovery_time_ <= time
        },

        #' @description
        #' Test if host is due sampling at a given time
        #' 
        #' @param time The time being queried
        #' 
        #' @return Logical 
        #' 
        is_sampling_due = function(time) {
            stopifnot(is.numeric(time), length(time) == 1)

            if (length(private$strains_) == 0) FALSE
            else private$sample_time_ == time
        },

        #' @description
        #' Test if host is susceptible at a given time
        #' 
        #' @param time The time being queried
        #' 
        #' @return Logical indicating if the host was susceptible at `time`
        #' 
        is_susceptible = function(time) {
            stopifnot(is.numeric(time), length(time) == 1)

            if (length(private$strains_) == 0) TRUE
            else time < private$exposure_time_
        },

        #' @description
        #' Prepare this host for an infection by deciding the times at which they enter the infectious and recovered compartments and when they are due to be sampled. Finally the host returns an interval which we use to find a their infector.
        #' 
        #' @param time The time of infection
        #' @param initial If this is an index case (index cases are placed straight in the infectious compartment)
        #' 
        #' @return Time interval that will be used to find a potential infector
        #' 
        prepare_for_infection = function(time, initial = FALSE) {

            stopifnot(is.numeric(time), time %% 1 == 0, length(time) == 1)
            stopifnot(is.logical(initial), length(initial) == 1)

            # check if this host has already been infected
            if (!is.na(private$infectious_time_)) {
                stop("Attempt to infect a non-susceptible host")
            }

            private$exposure_time_ <- time
            private$realisation_time_ <- time
            
            # calculate infectious, recovery and sampling times
            if (initial) {
                # this is an initial infectious case. Put it straight in the infectious compartment
                private$infectious_time_ <- time
                infector_interval <- NA

            } else {
                intervals <- private$population_$pop_interval_stack
                private$infectious_time_ <- time + ceiling(intervals$inc)
                infector_interval <- intervals$infector_interval

            }
            private$recovery_time_ <- private$generate_recovery_time()
            private$sample_time_ <- private$generate_sample_time()

            return(infector_interval)
        },

        #' @description
        #' Print a description of the `Host` object
        #' 
        #' @param ... arguments will be ignored
        #' 
        #' @return No return value. Description is printed to the console
        #' 
        print = function(...) {
            cat("Host: \n")
            cat("   Id:             ", private$id_, "\n", sep = "")
            cat("   Population:     ", private$population_$id, "\n", sep = "")
            cat("   Exposure Time:  ", private$exposure_time_, "\n", sep = "")
            cat("   Infectious Time: ", private$infectious_time_, "\n", sep = "")
            cat("   Recovery Time:  ", private$recovery_time_, "\n", sep = "")
            cat("   Sample Time:    ", private$sample_time_, "\n", sep = "")
            cat("   Realision Time: ", private$realisation_time_, "\n", sep = "")
            cat("   Num Strains:    ", length(private$strains_), "\n", sep = "")
            cat("   Frequencies:    ", paste0(private$freq_, collapse = ", "), "\n", sep = "")

            invisible(self)
        },

        #' @description 
        #' Realise changes in this `Host` object: This involves updating the  realisation time and deciding growth of pathogen population (if we implement it)
        #' 
        #' @param time The time of infection
        #' @param initial If this is an index case (index cases are placed straight in the infectious compartment)
        #' 
        #' @return Returns this `Host` object
        #' 
        realise = function(time) {
            stopifnot(is.numeric(time), time %% 1 == 0, length(time) == 1)
            
            # if growth of pathogen populations is implemented:
            # new strains caused by mutations are decided and temporarily stored
            # Growth of existing strains is performed
            # then newly mutated strains are stored permanently

            private$realisation_time_ <- time

            return(invisible(self))
        },

        #' @description 
        #' Set a new sampling time in the future. It will be set to the current sampling time plus the `sample_freq`. If this new time is after the `Host` recovers then it is set to infinity. If `sample_schedule` is random nothing will happen. 
        #' 
        #' @return Returns this `Host` object
        #' 
        update_sample_time = function() {
            
            if (self$sample_schedule != "random") {
                t <- private$sample_time_ + self$sample_freq

                # if sample times are at least recovery time set sample time to infinity so that no further sampling is attempted
                if (t >= private$recovery_time_) t <- Inf

                private$sample_time_ <- t
            }

            return(invisible(self))
        }
    ),

    ###########################################################################
    # active bindings

    active = list(

        #' @field data returns a list of data: host id and population id, infector id and population id, exposure time, infectious time, and recovery time
        data = function(v) {
            if (missing(v)) {

                if (!is.null(private$infector_)) {
                    infector_id <- private$infector_$id
                    infector_population <- private$infector_$population
                } else {
                    infector_id <- NA
                    infector_population <- NA
                }

                list(
                    host_id = private$id_,
                    population = private$population_$id,
                    infector_id = infector_id,
                    infector_population = infector_population,
                    exposure_time = private$exposure_time_,
                    infectious_time = private$infectious_time_,
                    recovery_time = private$recovery_time_
                )
            } else stop("Can't set `$data`", call. = FALSE)
        },

        #' @field exposure_time The time at which this host was exposed (Returns NA if still susceptible)
        exposure_time = function(v) {
            if (missing(v)) private$exposure_time_
            else stop("Can't set `$exposure_time`", call. = FALSE)
        },

        #' @field freq The frequencies of each `Strain` object in `strains`
        freq = function(v) {
            if (missing(v)) private$freq_
            else stop("Can't set `$freq`", call. = FALSE)
        },

        #' @field id The id of this host
        id = function(v) {
            if (missing(v)) private$id_
            else stop("Can't set `$id`", call. = FALSE)
        },

        #' @field infectious_time The time at which this host became infectious (Returns NA if still susceptible)
        infectious_time = function(v) {
            if (missing(v)) private$infectious_time_
            else stop("Can't set `$infectious_time`", call. = FALSE)
        },

        #' @field infectees A vector of `Host` objects that this `host` infected
        infectees = function(v) {
            if (missing(v)) private$infectees_
            else stop("Can't set `$infectees`", call. = FALSE)
        },

        #' @field infector The `Host` object that infected this `host`
        infector = function(v) {
            if (missing(v)) {
                if (! is.null(private$infector_)) private$infector_
                else NULL
            }
            else stop("Can't set `$infector`", call. = FALSE)
        },

        #' @field is_index Logical indicating if this `Host` was an index case
        is_index = function(v) {
            if (missing(v)) {
                if (is.na(private$recovery_time_)) FALSE
                else is.null(private$infector_)
            } else {
                stop("Can't set `$is_index`.", call. = FALSE)
            }
        },

        #' @field population The `Population` object that this host belongs to
        population = function(v) {
            if (missing(v)) {
                if (! is.null(private$population_)) private$population_
                else NULL
            }else stop("Can't set `$population`", call. = FALSE)
        },

        #' @field realisation_time The time at which this `Host` last realised itself
        realisation_time = function(v) {
            if (missing(v)) private$realisation_time_
            else stop("Can't set `$realisation_time`", call. = FALSE)
        },

        #' @field recovery_time The time at which this host recovered (Returns NA if still susceptible)
        recovery_time = function(v) {
            if (missing(v)) private$recovery_time_
            else stop("Can't set `$recovery_time`", call. = FALSE)
        },

        #' @field sample_freq The `sample_freq` specified in this `Host`'s `Population` object
        sample_freq = function(v) {
            if (missing(v)) private$population_$sample_freq
            else stop("Can't set `$sample_freq`", call. = FALSE)
        },

        #' @field sample_time The next time at which this host is due sampling (Returns NA if still susceptible)
        sample_time = function(v) {
            if (missing(v)) private$sample_time_
            else stop("Can't set `$sample_time`", call. = FALSE)
        },

        #' @field sample_schedule The `sample_schedule` specified in this `Host`'s `Population` object
        sample_schedule = function(v) {
            if (missing(v)) private$population_$sample_schedule
            else stop("Can't set `$sample_schedule`", call. = FALSE)
        },

        #' @field strains A vector of `Strain` objects that this `Host` is carrying
        strains = function(v) {
            if (missing(v)) private$strains_
            else stop("Can't set `$strains`", call. = FALSE)
        }

    ),

    ###########################################################################
    # private members

    private = list(

        #######################################################################
        # private variables 
        exposure_time_ = NA,
        id_ = NA,
        infectees_ = c(),
        infectious_time_ = NA,
        infector_ = NULL,
        population_ = NULL,
        realisation_time_ = NA,
        recovery_time_ = NA,
        sample_time_ = c(),
        strains_ = c(),
        freq_ = c(),

        #######################################################################
        # private functions

        # @description
        # Generate a calendar sample time
        # 
        # @return The sample time
        # 
        calendar_sample_time = function() {
            # check if this is an initial infection (we're starting sim)
            if (private$infectious_time_ == 0) {
                # starting simulation - just use individual method
                sample <- self$sample_freq

            } else {
                sample <- ceiling(private$infectious_time_ / self$sample_freq) * self$sample_freq

            }

            return(sample)
        },

        # @description
        # Generate a recovery time
        # 
        # @return An integer
        # 
        generate_recovery_time = function() {
            # generate a recovery time
            private$infectious_time_ + ceiling(stats::rgamma(
                1, 
                private$population_$rec_shape, 
                private$population_$rec_rate
            ))
        },

        # @description
        # Generate a sample time determined by `samp_schedule` specified in this `host`'s `population` object
        # 
        # @return The sample time
        # 
        generate_sample_time = function() {
            # generates sample times determined by samp_schedule

            if (self$sample_schedule == "individual") {
                sample <- private$individual_sample_time()

            } else if (self$sample_schedule == "calendar") {
                sample <- private$calendar_sample_time()

            } else { # random
                sample <- private$random_sample_time()

            }
            # if sample times are at least recovery time set sample time to infinity so that no further sampling is attempted
            if (sample >= private$recovery_time_) sample <- Inf

            return(sample)
        },

        # @description
        # Generate an individual sample time
        # 
        # @return The sample time
        # 
        individual_sample_time = function() {
            return(private$infectious_time_ + self$sample_freq)
        },

        # @description
        # Generate a random sample time
        # 
        # @return The sample time
        # 
        random_sample_time = function() {
            
            # because R's sample interprets a vector with a single number in it as just a single number and not a vector we need to test if recovery is one day after infection

            if (private$infectious_time_ == private$recovery_time_ - 1) {
                sample <- private$infectious_time_
            } else {
                sample <- sample(private$infectious_time_:(private$recovery_time_ - 1), 1)
            }

            return(sample)
        }

    ),    
)

