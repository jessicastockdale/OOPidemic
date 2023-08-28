#' Outbreak
#' 
#' @description An R6 class representing an outbreak. Used to simulate an outbreak in several populations (or groups)
#' 
#' @export
#' 
Outbreak <- R6::R6Class("Outbreak", 

    cloneable = FALSE,

    ###########################################################################
    # public members

    public = list(

        #' @field inf_rates An matrix with the same number of rows and columns as there is populations. The entry in the `i`th row and `j`the column is how much infectious pressure an individual in population `i` exerts on an individual in population `j`
        inf_rates = NULL, 
        #' @field lab The `Lab` object (or an R6Class that inherits from `Lab`) that takes samples from `populations`
        lab = NULL, 
        #' @field populations A vector of `Population` objects (or R6Classes that inherit from `Population`) 
        populations = NULL, 
        #' @field time An integer representing the current time in the outbreak
        time = 0L,

        #' @description
        #' Create a new `Outbreak` object. 
        #' 
        #' @param populations A vector of `Population` objects (or R6Classes that inherit from `Population`) 
        #' @param inf_rates An matrix with the same number of rows and columns as there is populations. The entry in the `i`th row and `j`the column is how much infectious pressure an individual in population `i` exerts on an individual in population `j`
        #' @param lab A `Lab` object (or R6Class that inherits from `Lab`)
        #' 
        #' @return A `Outbreak` object
        #' 
        initialize = function(
            populations,
            inf_rates,
            lab
        ) {

            stopifnot(all(vapply(
                populations,
                function(population) {"Population" %in% class(population)},
                logical(1)
            )))
            pop_ids <- vapply(populations, function(p) p$id, numeric(1L))
            if (any(
                length(unique(pop_ids)) != length(populations),
                any(pop_ids <= 0)
            )) {
                stop("ID of each population must be unique and at least 1.")
            }

            stopifnot(
                nrow(inf_rates) == length(populations),
                ncol(inf_rates) == length(populations),
                is.numeric(inf_rates),
                all(inf_rates >= 0)
            )
            stopifnot("Lab" %in% class(lab))
            
            self$inf_rates <- inf_rates
            self$lab <- lab
            self$populations <- populations
            private$pop_ids_ <- pop_ids
        },

        #' @description
        #' Select `host`s from the susceptible hosts in each `populations` by calling `get_infectees()` for each
        #' 
        #' @param num_infectees An vector of integers indicating how many infectees to select from each population
        #' 
        #' @return A vector of vectors of `Host` objects
        #' 
        get_infectees = function(num_infectees) {

            stopifnot(
                is.numeric(num_infectees), 
                all(num_infectees %% 1 == 0),
                length(num_infectees) == length(self$populations)
            )

            infectees <- mapply(
                function(population, new_I) {
                    population$get_infectees(new_I)
                },
                self$populations,
                num_infectees
            )

            return(infectees)
        },

        #' @description
        #' Select infectors from the infectious hosts in each population in `populations` by calling `get_infectors()` for each 
        #' 
        #' @param population_of_infectors A vector of  the population id of each infector. 
        #' @param infector_intervals A vector of integers indicating how long ago each an infector was infected / infectious (depends on the `find_infector_method` in each population). If `find_infector_method="random"` then it's a vector of zeros. 
        #' 
        #' @return A vector of `Host` objects
        #' 
        get_infectors = function(population_of_infectors, infector_intervals) {

            stopifnot(
                is.numeric(population_of_infectors),
                all(population_of_infectors %% 1 == 0),
                all(population_of_infectors %in% private$pop_ids_)
            )
            stopifnot(
                is.numeric(infector_intervals),
                all(infector_intervals >= 0)
            )
            stopifnot(length(population_of_infectors) == length(infector_intervals))

            # loop through each population and select infectors
            # the method for selecting infectors is defined in the Population class
            infectors <- list()
            for (population in self$populations) {
                ind <- which(population_of_infectors == population$id)

                # check if any infectors are needed from this population
                if (length(ind) != 0) {
                    # currently Population$get_infectors() gets how many infectors to select from the length of the vector passed to it. It doesn't need any otehr  info

                    infectors[ind] <- population$get_infectors(infector_intervals[ind])
                }
            }

            return(infectors)
        },

        #' @description
        #' Select the population that each infector comes from
        #' 
        #' @param num_infectees An vector of integers indicating how many infectees to select from each population
        #' @param fois_mat A force of infection matrix of the same shape and form as inf_rates
        #' 
        #' @return A vector of vectors of population ids
        #' 
        get_population_of_infectors = function(num_infectees, fois_mat) {

            stopifnot(
                is.numeric(num_infectees), 
                all(num_infectees %% 1 == 0),
                length(num_infectees) == length(self$populations)
            )
            stopifnot(
                nrow(fois_mat) == length(self$populations),
                ncol(fois_mat) == length(self$populations),
                is.numeric(fois_mat),
                all(fois_mat >= 0),
                all(fois_mat <= 1)
            )

            # store the number of populations (we'll use this a few times)
            n <- length(self$populations)

            # get the population that the infector for each infectee will come from
            # for each population use the elements of it's corresponding column in fois_mat for the probabilities of the population of the infectors
            population_for_infectors <- lapply(
                seq_along(num_infectees),
                function(i) {
                    sample(
                        n, 
                        num_infectees[i], 
                        replace = TRUE, 
                        fois_mat[i, ]
                    )
                }
            )

            return(population_for_infectors)
        },

        #' @description
        #' Determine number of infectees in each population
        #' 
        #' @param fois A vector of the force of infection experienced by each population
        #' 
        #' @return An vector of integers
        #' 
        get_num_infectees = function(fois) {
            stopifnot(
                length(fois) == length(self$populations),
                is.numeric(fois),
                all(fois >= 0)
            )

            num_infectes <- mapply(
                function(population, foi) population$get_num_infectees(foi),
                self$populations,
                fois
            )

            return(num_infectes)
        },

        #' @description
        #' Trigger an infection cycle in all `populations`. All new infectees (if any) are selected, their infectors are selected and infections of hosts are carried out.
        #' 
        #' @return This `Outbreak` object
        #' 
        infect = function() {

            # get number of infectious and susceptible for each population
            inf_sizes <- self$infectious_sizes

            # check if there are any infectives and stop if there isn't
            if (sum(inf_sizes) == 0) {return(invisible(self))}

            # get FOI matrix
            fois_mat <- private$force_of_infection(inf_sizes)
            
            # get new infectees for each population
            num_infectees_in_populations <- self$get_num_infectees(rowSums(fois_mat))
            infectees <- self$get_infectees(num_infectees_in_populations)


            # get the population that the infector for each infectee comes from
            population_for_infectors <- self$get_population_of_infectors(
                num_infectees_in_populations, 
                fois_mat
            )

            # population_for_infectors and infectees are both lists of lists of hosts
            # we can unlist them now
            infectees <- unlist(infectees, recursive = FALSE)
            population_for_infectors <- unlist(population_for_infectors, recursive = FALSE)

            # if there are no infectees then stop
            if (length(infectees) == 0) return(invisible(self))

            # prepare infectees for infection
            infector_intervals <- vapply(
                infectees,
                function(infectee) infectee$prepare_for_infection(self$time),
                numeric(1L)
            )

            # get infectors 
            infectors <- self$get_infectors(
                population_for_infectors, 
                infector_intervals
            )
           
            # have infectors infect infectees
            mapply(
                function(infector, infectee) {
                    infector$infect(infectee, self$time)
                },
                infectors,
                infectees
            )

            return(invisible(self))
        },

        #' @description
        #' Print a description of this `Outbreak` object
        #' 
        #' @param ... arguments will be ignored
        #' 
        #' @return No return value. Description is printed to the console
        #' 
        print = function() {

            outbreak_data <- data.frame(
                Size = self$population_sizes,
                Susceptibles = self$susceptible_sizes,
                Exposed = self$exposed_sizes,
                Infectious = self$infectious_sizes,
                Recovered = self$recovered_sizes
            )

            
            cat("Outreak: \n")
            cat("   Time:                   ", self$time, "\n", sep = "")
            cat("   Samples Taken:          ", self$num_wgs, "\n", sep = "")
            cat("   Population Details: \n")
            writeLines(paste0("         ", capture.output(print(t(outbreak_data)))))
        },

        #' @description
        #' Run this simulation.
        #' 
        #' @param feedback How frequently to print population summary to the console (default 10)
        #' 
        #' @return This `Outbreak` object 
        #' 
        run_simulation = function(feedback = 5) {

            stopifnot(is.numeric(feedback), length(feedback) == 1, feedback %% 1 == 0)

            # first have lab perform samples due at t = 0
            self$sample_populations()
            
            while (self$is_outbreak_active) {

                # increment time for self and for each population
                self$time <- self$time + 1
                private$align_population_times()

                # trigger a round of infections
                self$infect()

                # have lab take samples
                self$sample_populations()

                # give feedback if wanted
                if (feedback != 0) {
                    if (self$time %% feedback == 0) {
                        self$print()
                    } 
                }
            }

            if (feedback != 0) self$print()
    
            return(invisible(self))
        },

        #' @description
        #' Collect samples from `host`s in every population that are due for sampling
        #' 
        #' @return This `Outbreak` object
        #' 
        sample_populations = function() {

            # get hosts that are due sampling

            hosts <- self$hosts_due_for_sampling
            self$lab$sample_hosts(hosts, self$time)

            return(invisible(self))
        }#,

    ),

    ###########################################################################
    # active bindings
    active = list(

        #' @field exposed_sizes The number of hosts in the exposed compartment in each population
        exposed_sizes = function(v) {

            if (missing(v)) {
                vapply(
                    self$populations,
                    function(population) {population$exposed_size},
                    integer(1)
                )
            } else {stop("Can't set `$ exposed_sizes`.", call. = FALSE)}
        },

        #' @field hosts_due_for_sampling A vector of hosts from all populations that are due sampling at the current time
        hosts_due_for_sampling = function(v) {
            
            if (missing(v)) {
                hosts <- lapply(
                    self$populations, 
                    function(population) population$hosts_due_for_sampling
                )
                hosts <- unlist(hosts, recursive = FALSE)

                return(hosts)
            } else {stop("Can't set `$get_hosts_due_for_sampling`", call. = FALSE)}
        },

        #' @field is_outbreak_active Logical indicating if there is any host in any population that is still in the exposed or infectious compartment
        is_outbreak_active = function(v) {
            if (missing(v)) {
                any(vapply(
                    self$populations, 
                    function(population) {population$is_outbreak_active},
                    logical(1L)
                ))
            } else {stop("Can't set `$is_outbreak_active`.", call. = FALSE)}
        },

        #' @field infectious_sizes The number of hosts in the infectious compartment in each population
        infectious_sizes = function(v) {
            # get number of infectious in each population
            if (missing(v)) {
                vapply(
                    self$populations,
                    function(population) {population$infectious_size},
                    integer(1)
                )
            } else {stop("Can't set `$ infectious_sizes`.", call. = FALSE)}
        },

        #' @field recovered_sizes The number of hosts in the recovered compartment in each population
        recovered_sizes = function(v) {
            # get number of recovered in each population
            if (missing(v)) {
                vapply(
                    self$populations,
                    function(population) {population$recovered_size},
                    integer(1)
                )
            } else {stop("Can't set `$recovered_sizes`.", call. = FALSE)}
        },

        #' @field susceptible_sizes The number of hosts in the susceptible compartment in each population
        susceptible_sizes = function(v) {
            # get number of susceptible in each population
            if (missing(v)) {
                vapply(
                    self$populations,
                    function(population) {population$susceptible_size},
                    integer(1)
                )
            } else {stop("Can't set `$susceptible_sizes`.", call. = FALSE)}
        },

        #' @field num_wgs The number of whole genome sequences collected so far
        num_wgs = function(v) {
            # get number of tests run by lab
            if (missing(v)) {
                self$lab$num_wgs
            } else {stop("Can't set `$num_wgs`.", call. = FALSE)}
        },

        #' @field outbreak_size The number of non_susceptible hosts across all populations
        outbreak_size = function(v) {
            # get the outbreak size across all populations
            if (missing(v)) {
                sum(self$population_outbreak_sizes)
            } else {stop("Can't set `$outbreak_size`.", call. = FALSE)}
        },

        #' @field population_outbreak_sizes The number of non_susceptible hosts in each population
        population_outbreak_sizes = function(v) {
            # get the outbreak size across all populations
            if (missing(v)) {
                vapply(
                    self$populations,
                    function(population) {population$outbreak_size},
                    integer(1L)
                )
            } else {stop("Can't set `$population_outbreak_sizes`.", call. = FALSE)}
        },

        #' @field population_sizes The size of each population
        population_sizes = function(v) {
            # get size of each population
            if (missing(v)) {
                vapply(
                    self$populations,
                    function(population) {population$size},
                    double(1)
                )
            } else {stop("Can't set `$population_sizes`.", call. = FALSE)}
        }

    ),

    ###########################################################################
    # private members
    private = list( 

        pop_ids_ = c(),

        # @description
        # align all population's time to this outbreaks time
        # 
        # @return This `outbreak` object 
        # 
        align_population_times = function() {
            # align each Population's time to Outbreak's time
            vapply(
                self$populations,
                function(population) {
                    population$time <- self$time
                    TRUE
                },
                logical(1L)
            )

            return(invisible(self))
        },

        # @description
        # Calculate the force_of_infectious matrix
        # 
        # @param inf_sizes The number of infectious hosts in each population
        # 
        # @return A matrix
        # 
        force_of_infection = function(inf_sizes) {
            # calculate the force of infection matrix. 
            # Each element will be of the form beta_{ij}*I_j/N 

            beta_I <- sweep(self$inf_rates, MARGIN = 2, inf_sizes, "*") 
            foi_mat <- beta_I / sum(self$population_sizes)

            return(foi_mat)
        }
    )
)
