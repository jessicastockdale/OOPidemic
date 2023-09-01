#' Group
#' 
#' @description An R6 class representing a group of individuals. The Default model for an disease outbreak is SEIR but can be switched to an SIR model by setting `inc_shape` to zero.
#' 
#' @examples
#' ref_strain <- ReferenceStrain$new(name = "ref_strain")
#' Group$new(id = 1, ref_strain = ref_strain)
#' Group$new(
#'  id = 1, 
#'  ref_strain = ref_strain, 
#'  find_infector_method = "serial",
#'  si_shape = 6,
#'  si_rate = 2 
#' )
#' Group$new(
#'  id = 1, 
#'  ref_strain = ref_strain, 
#'  find_infector_method = "random",
#'  max_init_dist = 20
#' )
#' @export
#' 
Group <- R6::R6Class("Group",

    cloneable = FALSE,

    ###########################################################################
    # public members

    public = list(

        #######################################################################
        # public variables

        #' @field id The id of this group ** this must be unique if multiple groups are being used in an `Population` object
        id = NA, 
        
        #' @field init_sus The number of initial susceptible individuals
        init_sus = NA, 
        #' @field init_inf The number of initial infected individuals
        init_inf = NA, 

        #' @field inf_rate The infectious pressure exerted by one infectious `host` on the rest of this `group`
        inf_rate = NA, 

        #' @field find_infector_method The method with which an infectees infector is chosen (see definition of new for valid choices)
        find_infector_method = NA, 

        #' @field si_shape,si_rate The Gamma shape and rate parameters for the serial interval 
        si_shape = NA, si_rate = NA, 

        #' @field gt_shape,gt_rate The Gamma shape and rate parameters for the generation time 
        gt_shape = NA, gt_rate = NA, 

        #' @field trans_int_shape,trans_int_rate The Gamma shape and rate parameters for the transmission itnerval
        trans_int_shape = NA, trans_int_rate = NA, 
        
        #' @field inc_shape,inc_rate The Gamma shape and rate parameters for the incubation time
        inc_shape = NA, inc_rate = NA, 

        #' @field rec_shape,rec_rate The Gamma shape and rate parameters for the recovery time. 
        rec_shape = NA, rec_rate = NA, 
        
        #' @field ref_strain The `ReferenceStrain` object used to create `Strain` objects for the index cases if there was any
        ref_strain = NULL, 

        #' @field min_init_dist,max_init_dist The minimum and maximum distances that index `Strain` objects should be from the `ref_strain`
        max_init_dist = NA, min_init_dist = NA, 

        #' @field sample_schedule,sample_freq The method used to schedule sampling times and the time interval between sampling (see new definition for valid choices):
        sample_schedule = NULL, sample_freq = NA, 

        #' @field time The current time for this group
        time = 0L,

        #######################################################################
        # public functions

        #' @description
        #' Create a new `group` object
        #' 
        #' @param id The id of this group ** this must be unique if multiple groups are being used in an `Population` object
        #' @param ref_strain The time that this strain is 'spawned'
        #' @param init_sus The number of initial susceptible individuals (default 99)
        #' @param init_inf The number of initial infected individuals (default 1)
        #' @param inf_rate The infectious pressure exerted by one infectious `host` on the rest of this `group`
        #' @param find_infector_method The method with which an infectees infector is chosen (see definition of new for valid choices)
        #' * 'transmission' (Default) Uses the transition interval to pick the infector. A transition interval is the time elapsed between an infector becoming infectious and the infectee being exposed
        #' * 'serial': Uses the serial interval to pick the infector. A serial interval is the time elapsed between an infector and their infectee becoming infectious
        #' * 'generation': Uses the generation time to pick the infector. A generation time is the time elapsed between an infector and their infectee being exposed
        #' * 'random': An infector is randomly selected from all currently infectious hosts
        #' @param si_shape,si_rate The Gamma shape and rate parameters for the serial interval (Defaults are `6` and `2`)
        #' @param gt_shape,gt_rate The Gamma shape and rate parameters for the generation time 
        #' @param trans_int_shape,trans_int_rate The Gamma shape and rate parameters for the transmission interval
        #' @param inc_shape,inc_rate The Gamma shape and rate parameters for the incubation time. If running an SIR simulation then set `inc_shape` to zero so that a host spends 0 time in the exposed compartment. (Defaults are `3` and `2` respectively)
        #' @param rec_shape,rec_rate The Gamma shape and rate parameters for the recovery time. (Defaults are `5` and `1` respectively)
        #' @param min_init_dist,max_init_dist The minimum and maximum distances that index `Strain` objects should be from the `ref_strain` (Default for both is `0`)
        #' @param sample_schedule,sample_freq The method used to schedule sampling times and the time interval between sampling
        #' * 'random': (Default) A random time is chosen between a `host`'s infectious time and their recovery time
        #' * 'calendar': Samples are taken from every infectious host at `sample_freq` intervals (starting from time = `sample_freq`)
        #' * 'individual': A `host` is sampled every `sample_freq` time intervals up until they recover
        #' @param host_class The class used to create hosts. Default is `Host` a subclassed class can be used here instead
        #' @param strain_class The class used to create strains. Default is `Strain` a subclassed class can be used here instead
        #' 
        #' @return A new 'Group' object
        #' 
        initialize = function(
            id,
            ref_strain = NULL,
            init_sus = 99,
            init_inf = 1,
            inf_rate = 1,
            find_infector_method = "transmission",
            si_shape = NA,
            si_rate = NA,
            gt_shape = NA,
            gt_rate = NA,
            trans_int_shape = 3, 
            trans_int_rate = 2, 
            inc_shape = 3, 
            inc_rate = 2,
            rec_shape = 5,
            rec_rate = 1,
            max_init_dist = 0,
            min_init_dist = 0,
            sample_schedule = "random",
            sample_freq = NA,
            host_class = Host,
            strain_class = Strain
        ) {
            
            # parameter validations
            stopifnot(is.numeric(id), id %% 1 == 0, length(id) == 1)
            stopifnot(any(
                is.null(ref_strain),
                "ReferenceStrain" %in% class(ref_strain)
            ))

            stopifnot(init_sus %% 1 == 0, init_sus >= 0, length(init_sus) == 1)
            stopifnot(init_inf %% 1 == 0, init_inf >= 0, length(init_inf) == 1)

            stopifnot(is.numeric(inf_rate), 0 < inf_rate, length(inf_rate) == 1)

            stopifnot(find_infector_method %in% c(
                "serial", 
                "transmission", 
                "random",
                "generation"
            ))
            if (find_infector_method == "transmission") {
                stopifnot(
                    is.numeric(trans_int_shape), 
                    trans_int_shape > 0, 
                    length(trans_int_shape) == 1
                )
                stopifnot(
                    is.numeric(trans_int_rate), 
                    trans_int_rate > 0, 
                    length(trans_int_rate) == 1
                )
            } else if (find_infector_method == "serial") {
                stopifnot(
                    is.numeric(si_shape), 
                    0 < si_shape, 
                    length(si_shape) == 1
                )
                stopifnot(
                    is.numeric(si_rate), 
                    0 < si_rate, 
                    length(si_rate) == 1
                )
            } else if (find_infector_method == "generation") {
                stopifnot(
                    is.numeric(gt_shape), 
                    gt_shape > 0, 
                    length(gt_shape) == 1
                )
                stopifnot(
                    is.numeric(gt_rate), 
                    gt_rate > 0, 
                    length(gt_rate) == 1
                )
            }

            stopifnot(
                is.numeric(inc_shape), 
                0 <= inc_shape, 
                length(inc_shape) == 1
            )
            stopifnot(
                is.numeric(inc_rate), 
                0 < inc_rate, 
                length(inc_rate) == 1
            )
            stopifnot(
                is.numeric(rec_shape), 
                0 < rec_shape, 
                length(rec_shape) == 1
            )
            stopifnot(
                is.numeric(rec_rate), 
                0 < rec_rate, 
                length(rec_rate) == 1
            )

            stopifnot(
                max_init_dist %% 1 == 0, 
                max_init_dist >= 0, 
                length(max_init_dist) == 1
            )
            stopifnot(
                min_init_dist %% 1 == 0, 
                min_init_dist >= 0, 
                length(min_init_dist) == 1
            )
            stopifnot(max_init_dist >= min_init_dist)

            stopifnot(sample_schedule %in% c("calendar", "individual", "random"))
            if (sample_schedule %in% c("calendar", "individual")) {
                stopifnot(
                    sample_freq %% 1 == 0, 
                    sample_freq >= 1, 
                    length(sample_freq) == 1
                )
            }
            
            stopifnot(is_ancestor_r6class(Host, host_class))
            stopifnot(is_ancestor_r6class(Strain, strain_class))

            # store variables
            self$find_infector_method <- find_infector_method
            self$gt_shape <- gt_shape
            self$gt_rate <- gt_rate
            self$id <- id
            self$inc_rate <- inc_rate
            self$inc_shape <- inc_shape
            self$inf_rate <- inf_rate
            self$init_inf <- init_inf
            self$init_sus <- init_sus
            self$max_init_dist <- max_init_dist
            self$min_init_dist <- min_init_dist
            self$trans_int_rate <- trans_int_rate
            self$trans_int_shape <- trans_int_shape
            self$rec_rate <- rec_rate
            self$rec_shape <- rec_shape
            self$ref_strain <- ref_strain
            self$sample_freq <- sample_freq
            self$sample_schedule <- sample_schedule
            self$si_rate <- si_rate
            self$si_shape <- si_shape
            
            private$host_class_ <- host_class
            private$strain_class_ <- strain_class

            # set up stack of SI, incubation, and transfer interval times
            private$build_interval_stack()

            # initialise outbreak
            private$initialise_hosts()

        },

        #' @description
        #' Get all hosts that are in the exposed compartment at `time`
        #' 
        #' @param time The time being querried (Defaults to this group's time)
        #' 
        #' @return A vector of `Host` objects
        #' 
        exposed_hosts = function(time = self$time) {

            stopifnot(is.numeric(time), time %% 1 == 0, length(time) == 1)

            private$hosts_[vapply(
                private$hosts_, 
                function(host) host$is_exposed(time),
                logical(1)
            )]
        },


        #' @description
        #' Select a `num_infectees` `host`s from the susceptible hosts in this `group`. All susceptible hosts are equally likely to be selected
        #' 
        #' @param num_infectees An integer indicating the number of susceptible hosts to select
        #' 
        #' @return A vector of `Host` objects
        #' 
        #' @examples
        #' ref_strain <- ReferenceStrain$new(name = "ref_strain")
        #' group <- Group$new(id = 1, ref_strain = ref_strain)
        #' group$get_infectees(3)
        #' 
        get_infectees = function(num_infectees) {

            stopifnot(
                is.numeric(num_infectees), 
                num_infectees %% 1 == 0, 
                length(num_infectees) == 1
            )

            infectees <- sample(self$susceptible_hosts(), num_infectees)
            return(infectees)
        },

        #' @description
        #' Select infectors. `find_infector_method` determines how `infector_intervals` is used to find the infectors. The length of `infector_intervals` determines how many infectors are chosen.
        #' 
        #' @param infector_intervals A vector of integers indicating how long ago an infector was infected / infectious (depends on `find_infector_method`). If `find_infector_method="random"` then pass a vector of zeros. 
        #' 
        #' @return A vector of `Host` objects
        #' 
        #' @examples
        #' ref_strain <- ReferenceStrain$new(name = "ref_strain")
        #' group <- Group$new(id = 1, ref_strain = ref_strain, init_inf = 10)
        #' group$get_infectors(c(1, 2))
        #' group <- Group$new(
        #'  id = 1, ref_strain = ref_strain, 
        #'  init_inf = 10, find_infector_method = "random"
        #' )
        #' group$get_infectors(c(0, 0, 0))
        #' 
        get_infectors = function(infector_intervals) {

            stopifnot(
                is.numeric(infector_intervals), 
                length(infector_intervals) > 0
            )

            num <- length(infector_intervals)

            # if there's only one infectious host we have no choice
            if (self$infectious_size == 1) {
                return(invisible(rep(self$infectious_hosts(), num)))
            }

            # which method are we using to choose an infector?
            if (self$find_infector_method == "random") {
                infectee_ind <- sample(self$infectious_size, num, replace = TRUE)

            } else {
                # generation, serial or transmission methods
                # use transmission interval to find infector

                # get time since _____ for all infectious hosts
                if (self$find_infector_method == "generation") {
                    time_since <- self$time_since_infectious_hosts_infected
                } else {
                    # serial or transmission method
                    time_since <- self$time_since_infectious_hosts_infectious
                }

                # build a matrix with num columns where every column is time_since
                ts_mat <- replicate(num, time_since)
                # subtract element i of infector_intervals from col i of ts_mat
                closest_times <- abs(sweep(ts_mat, 2, infector_intervals))
                
                # find the index of the minimum element of each column
                # randomly sample if more then one element qualifies
                # these indexes will corespond to the indexes of the current infectives 
                infectee_ind <- apply(
                    closest_times,
                    2,
                    function(col) {
                        inds <- which(col == min(col))
                        inds[sample.int(length(inds), 1)] # catches inds of length 1
                    }
                )
            }
            
            return(invisible(self$infectious_hosts()[infectee_ind]))
        },

        #' @description
        #' Determine number of infectees
        #' 
        #' @param foi The force of infection exerted on this group
        #' 
        #' @return An integer
        #' 
        get_num_infectees = function(foi) {
            stopifnot(
                is.numeric(foi), 
                length(foi) == 1
            )
            n <- rbinom(1, self$susceptible_size, foi)

            return(n)
        },

        #' @description
        #' Trigger an infection cycle in this `group`. All new infectees (if any) are selected, their infectors are selected and infections of hosts are carried out.
        #' 
        #' @return This group object
        #' 
        infect = function() {
            
            #increment time
            self$time <- self$time + 1

            # get infectees 
            num_infectees <- self$get_num_infectees(private$force_of_infection())
            infectees <- self$get_infectees(num_infectees)
            
            # if no one gets infected stop
            if (length(infectees) == 0) return(invisible(self))

            # prepare_for_infection infectees
            trans_ints <- vapply(
                infectees,
                function(infectee) infectee$prepare_for_infection(self$time),
                numeric(1L)
            )

            # get infectors
            infectors <- self$get_infectors(trans_ints)

            # infect hosts
            mapply(
                function(infector, infectee) {
                    infector$infect(infectee, self$time)
                },
                infectors, 
                infectees
            )

            invisible(self)
        },

        #' @description
        #' Get all hosts that are in the infectious compartment at `time`
        #' 
        #' @param time The time being querried (Defaults to this group's time)
        #' 
        #' @return A vector of `Host` objects
        #' 
        infectious_hosts = function(time = self$time) {

            stopifnot(is.numeric(time), time %% 1 == 0, length(time) == 1)

            private$hosts_[vapply(
                private$hosts_, 
                function(host) host$is_infectious(time),
                logical(1)
            )]
        },

        #' @description
        #' Print a description of the `Group` object
        #' 
        #' @param ... arguments will be ignored
        #' 
        #' @return No return value. Description is printed to the console
        #' 
        print = function(...) {
            cat("Group: \n")
            cat("   Id:                 ", self$id, "\n", sep = "")
            cat("   Time:               ", self$time, "\n", sep = "")
            cat("   Group Size:         ", self$size, "\n", sep = "")
            cat("   Susceptible Hosts:  ", self$susceptible_size, "\n", sep = "")
            if (self$inc_shape != 0) {
                # only print exposed compartment size if we are running a SEIR outbreak
                cat("   Exposed Hosts:      ", self$exposed_size, "\n", sep = "")
            }
            cat("   Infectious Hosts:   ", self$infectious_size, "\n", sep = "")
            cat("   Recovered Hosts:    ", self$recovered_size, "\n", sep = "")
        },

        #' @description
        #' Get all hosts that are in the recovered compartment at `time`
        #' 
        #' @param time The time being querried (Defaults to this group's time)
        #' 
        #' @return A vector of `Host` objects
        #' 
        recovered_hosts = function(time = self$time) {

            stopifnot(is.numeric(time), time %% 1 == 0, length(time) == 1)

            private$hosts_[vapply(
                private$hosts_, 
                function(host) host$is_recovered(time),
                logical(1)
            )]
        },

        #' @description
        #' Simulate an outbreak in this group. This is intended to be used if you are siming an outbreak in a single homogenous group. It's recommended that you use an `Population` class if your simulation requires a heterogenous population with more then one groups of hosts 
        #' 
        #' @param lab A lab object to collect samples from this group
        #' @param feedback How frequently to print group summary to the console (turn off by setting to zero)
        #' 
        #' @return This `group` object 
        #' 
        run_simulation = function(lab, feedback = 10) {
            stopifnot("Lab" %in% class(lab))
            stopifnot(
                is.numeric(feedback), length(feedback) == 1,
                feedback %% 1 == 0, feedback >= 0
            )

            # have lab take any samples due at t = 0
            lab$sample_hosts(self$hosts_due_for_sampling, self$time)

            while (self$is_outbreak_active) {
                # trigger a round of infections 
                self$infect()
                # sample hosts 
                lab$sample_hosts(self$hosts_due_for_sampling, self$time)
                
                # feedback 
                if (all(
                    feedback != 0,
                    self$time %% feedback == 0
                )) {
                    print(self)
                }
            }

            if (feedback != 0) {
                print(self)
                print(lab)
            }

            return(invisible(self))
        },

        #' @description
        #' Get all hosts that are in the susceptible compartment at `time`
        #' 
        #' @param time The time being querried (Defaults to this group's time)
        #' 
        #' @return A vector of `Host` objects
        #' 
        susceptible_hosts = function(time = self$time) {

            stopifnot(is.numeric(time), time %% 1 == 0, length(time) == 1)

            private$hosts_[vapply(
                private$hosts_, 
                function(host) host$is_susceptible(time),
                logical(1)
            )]
        }

    ),

    ###########################################################################
    # active bindings

    active = list(

        #' @field exposed_size The number of hosts in the exposed compartment at the groups current time
        exposed_size = function(v) {
            if (missing(v)) length(self$exposed_hosts())
            else stop("Can't set `$exposed_size`", call. = FALSE)
        },

        #' @field hosts All host objects in a vector
        hosts = function(v) {
            if (missing(v)) private$hosts_
            else stop("Can't set `$hosts`", call. = FALSE)
        },

        #' @field hosts_due_for_sampling A vector of hosts that are due for sampling at the group's current time 
        hosts_due_for_sampling = function(v) {

            if (missing(v)) {
                private$hosts_[vapply(
                    private$hosts_, 
                    function(host) host$is_sampling_due(self$time),
                    logical(1)
                )]
            
            } else {stop("Can't set `$get_hosts_due_for_sampling`", call. = FALSE)}
        },
        
        #' @field is_outbreak_active Logical indicating if there are still exposed or infectious individuals at the group's current time
        is_outbreak_active = function(v) {
            if (missing(v)) {!all(
                self$exposed_size == 0,
                self$infectious_size == 0
            )}
            else stop("Can't set `$is_outbreak_active`", call. = FALSE)
        },

        #' @field infectious_size The number of hosts in the infectious compartment at the groups current time
        infectious_size = function(v) {
            if (missing(v)) length(self$infectious_hosts())
            else stop("Can't set `$infectious_size`", call. = FALSE)
            
        },

        #' @field group_interval_stack Pop a set of intervals from the interval stack. The interval stack is a preprepared set of at least incubation periods and the interval being used to find an infector (serial interval and generation time are also present if that method is being used)
        group_interval_stack = function(v) {
            if (missing(v)) {
                if (nrow(private$interval_stack_) == 0) {
                    stop("interval_stack is empty")
                }

                ints <- private$interval_stack_[1, ]
                private$interval_stack_ <- private$interval_stack_[-1, ]

                return(ints)
            } else stop("Can't set `$recovered_size`", call. = FALSE)
        },

        #' @field recovered_size The number of hosts in the recovered compartment at the groups current time
        recovered_size = function(v) {
            if (missing(v)) length(self$recovered_hosts())
            else stop("Can't set `$recovered_size`", call. = FALSE)
            
        },

        #' @field susceptible_size The number of hosts in the susceptible compartment at the groups current time
        susceptible_size = function(v) {
            if (missing(v)) length(self$susceptible_hosts())
            else stop("Can't set `$susceptible_size`", call. = FALSE)
            
        },

        #' @field outbreak_size The number of infected so far
        outbreak_size = function(v) {
            if (missing(v)) {
                self$exposed_size + self$infectious_size + self$recovered_size
            }
            else stop("Can't set `$outbreak_size`", call. = FALSE)
            
        },

        #' @field size The group size
        size = function(v) {
            # get group size
            if (missing(v)) {
                length(private$hosts_)
            }
            else stop("Can't set `$size`", call. = FALSE)
        },

        #' @field time_since_infectious_hosts_infected A vector containing the time elapsed wince all currently infectious hosts were infected
        time_since_infectious_hosts_infected = function(v) {
            if (missing(v)) {
                vapply(
                    self$infectious_hosts(),
                    function(host) {
                        self$time - host$exposure_time
                    },
                    numeric(1)
                )
            }
        },
        
        #' @field time_since_infectious_hosts_infectious A vector containing the time elapsed wince all currently infectious hosts became infectious
        time_since_infectious_hosts_infectious = function(v) {
            if (missing(v)) {
                vapply(
                    self$infectious_hosts(),
                    function(host) {
                        self$time - host$infectious_time
                    },
                    numeric(1)
                )
            }
        }
        
    ),

    ###########################################################################
    # private members

    private = list(

        #######################################################################
        # private variables 
        host_class_ = NULL,
        hosts_ = NULL,
        strain_class_ = NULL,
        interval_stack_ = NULL,

        #######################################################################
        # private functions

        # @description
        # Generate a stack of intervals for hosts to use. At a minimum it will contain an incubation period and an "infector interval". The infector interval will be quantity used to find infectors for an infectee. Depending on the value of `find_infector_method` this stack may also contain a serial interval or generation time
        # 
        # 
        # @return No return value
        # 
        build_interval_stack = function() {

            n <- self$init_inf + self$init_sus
            factor <- 20

            if (self$find_infector_method == "serial") {
                interval_stack <- data.frame(
                    si = stats::rgamma(
                        factor * n, 
                        shape = self$si_shape, 
                        rate = self$si_rate
                    ),
                    inc = stats::rgamma(
                        factor * n, 
                        shape = self$inc_shape, 
                        rate = self$inc_rate
                    )
                )
                interval_stack$infector_interval <- (
                    interval_stack$si  - interval_stack$inc
                ) 

            } else if (self$find_infector_method == "transmission") {
                interval_stack <- data.frame(
                    inc = stats::rgamma(
                        n, 
                        shape = self$inc_shape, 
                        rate = self$inc_rate
                    ),
                    infector_interval = stats::rgamma(
                        n, 
                        shape = self$trans_int_shape, 
                        rate = self$trans_int_rate
                    )
                )
            } else if (self$find_infector_method == "random") {
                interval_stack <- data.frame(
                    inc = stats::rgamma(
                        n, 
                        shape = self$inc_shape, 
                        rate = self$inc_rate
                    ),
                    infector_interval = rep.int(0, n)
                )
            } else if (self$find_infector_method == "generation") {
                interval_stack <- data.frame(
                    inc = stats::rgamma(
                        n, 
                        shape = self$inc_shape, 
                        rate = self$inc_rate
                    ),
                    infector_interval = stats::rgamma(
                        n, 
                        shape = self$gt_shape,
                        rate = self$gt_rate
                    )
                )
            }

            neg_bools <- interval_stack$infector_interval < 0

            # how many negative transmission intervals are there?
            if (sum(neg_bools) > factor * n / 2) {
                warning(
                    "Parameters provided for serial interval and incubation period produce a negative transmission interval more then 50 % of the time. This may alter the expected distribution of the transmission interval. Please consider changing the parameters provided.",
                    immediate. = TRUE
                )
            }
            # drop the negative values
            interval_stack <- interval_stack[! neg_bools, ]

            # only keep one row for each host
            private$interval_stack_ <- utils::head(interval_stack, n)

            return(NULL)
        },

        # @description
        # Initialise infections in index hosts
        # 
        # @param host The `Host` object being initialised
        # @param variation Logical indicating if the strain infecting this host can have variations in it's genome compared to `ref_strain`
        # 
        # @return TRUE: The motivation for this is that we can call this function for many hosts with vapply if we know what the return value will be, so we've made it something small 
        # 
        initialise_index_host = function(host, variation = TRUE) {
            
            # decide strain's genomic distance from reference strain
            if (variation) {
                distance <- private$initial_genomic_distance()
            } else {
                distance <- 0
            }

            # prepare_for_infection 
            host$prepare_for_infection(time = 0, initial = TRUE)
            # now pass on strain
            host$contract_strains(
                infector = NULL,
                # time = 0,
                strains = c(private$strain_class_$new(self$ref_strain, 0, distance)),
                freq = c(1)#,
                # initial = TRUE
            )        

            return(TRUE)
        },

        # @description
        # Initialise all host objects and then initialise infections in the index cases
        # 
        # @return TRUE: The motivation for this is that we can call this function for many hosts with vapply if we know what the return value will be, so we've made it something small 
        # 
        initialise_hosts = function() {

            # can't use the self$size function here as that function requires the hosts to have been initialised
            pop_size <- self$init_sus + self$init_inf
            
            # create all host objects
            private$hosts_ <- lapply(
                c(1:pop_size), 
                private$host_class_$new,
                group = self
            ) 

            # randomly select indexes of initial infectives
            index_hosts <- sample(private$hosts_, self$init_inf, replace = FALSE)

            # set up strains in infectives to start outbreak
            # infectives is a vector containing the indices of the initial infectives
            if (self$init_inf == 1) {
                # only one strain needed - don't create variation
                private$initialise_index_host(index_hosts[[1]], variation = FALSE)
            
            } else {
                # multiple initial hosts
                vapply(index_hosts, private$initialise_index_host, logical(1L))
            }

            return(NULL)
        },

        # @description
        # Calculate the force of infection experienced by each individual in this group
        # 
        # @return The force of infection (numeric)
        # 
        force_of_infection = function() {
            self$inf_rate * self$infectious_size / self$size
        },

        # @description
        # Calculate the initial genomic distance of a strain from the reference strain
        # 
        # @return An integer
        # 
        initial_genomic_distance = function() {
            self$min_init_dist + sample(self$max_init_dist - self$min_init_dist, 1)
        }
    )
)
