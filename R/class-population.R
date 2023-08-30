#' Population
#' 
#' @description An R6 class representing an population. Used to simulate an disease outbreak in one or more
#' 
#' @examples 
#' # An SEIR model outbreak ina population with two groups
#' set.seed(1)
#' ref_strain <- ReferenceStrain$new("ref_strain")
#' population <- Population$new(
#'     c(
#'         Group$new(1, ref_strain, init_inf = 5, init_sus = 95),
#'         Group$new(2, ref_strain, init_inf = 0, init_sus = 100)
#'     ),
#'     matrix(c(0.75, 0.5, 0.25, 0.1), ncol = 2),
#'     Lab$new()
#' )
#' population$run_simulation()
#' # plot each group separately
#' plot_compartment_dynamics(population, group_labels = c("Senior", "Adult"))
#' # plot groups combined
#' plot_compartment_dynamics(population, combine = TRUE)
#' 
#' @export
#' 
Population <- R6::R6Class("Population", 

    cloneable = FALSE,

    ###########################################################################
    # public members

    public = list(

        #' @field inf_rates A matrix with the same number of rows and columns as there is groups. The entry in the `i`th row and `j`the column is how much infectious pressure an individual in group `i` exerts on an individual in group `j`
        inf_rates = NULL, 
        #' @field lab The `Lab` object (or an R6Class that inherits from `Lab`) that takes samples from `groups`
        lab = NULL, 
        #' @field groups A vector of `Group` objects (or R6Classes that inherit from `Group`) 
        groups = NULL, 
        #' @field time An integer representing the current time in the population
        time = 0L,

        #' @description
        #' Create a new `Population` object. 
        #' 
        #' @param groups A vector of `Group` objects (or R6Classes that inherit from `Group`) 
        #' @param inf_rates A matrix with the same number of rows and columns as there is groups. The entry in the `i`th row and `j`the column is how much infectious pressure an individual in group `i` exerts on an individual in group `j`
        #' @param lab A `Lab` object (or R6Class that inherits from `Lab`)
        #' 
        #' @return A `Population` object
        #' 
        initialize = function(
            groups,
            inf_rates,
            lab
        ) {

            stopifnot(all(vapply(
                groups,
                function(group) {"Group" %in% class(group)},
                logical(1)
            )))
            group_ids <- vapply(groups, function(p) p$id, numeric(1L))
            if (any(
                length(unique(group_ids)) != length(groups),
                any(group_ids <= 0)
            )) {
                stop("ID of each group must be unique and at least 1.")
            }

            stopifnot(
                nrow(inf_rates) == length(groups),
                ncol(inf_rates) == length(groups),
                is.numeric(inf_rates),
                all(inf_rates >= 0)
            )
            stopifnot("Lab" %in% class(lab))
            
            self$inf_rates <- inf_rates
            self$lab <- lab
            self$groups <- groups
            private$group_ids_ <- group_ids
        },

        #' @description
        #' Select `host`s from the susceptible hosts in each `groups` by calling `get_infectees()` for each
        #' 
        #' @param num_infectees A vector of integers indicating how many infectees to select from each group
        #' 
        #' @return A vector of vectors of `Host` objects
        #' 
        get_infectees = function(num_infectees) {

            stopifnot(
                is.numeric(num_infectees), 
                all(num_infectees %% 1 == 0),
                length(num_infectees) == length(self$groups)
            )

            infectees <- mapply(
                function(group, new_I) {
                    group$get_infectees(new_I)
                },
                self$groups,
                num_infectees
            )

            return(infectees)
        },

        #' @description
        #' Select infectors from the infectious hosts in each group in `groups` by calling `get_infectors()` for each 
        #' 
        #' @param group_of_infectors A vector of  the group id of each infector. 
        #' @param infector_intervals A vector of integers indicating how long ago each infector was infected / infectious (depends on the `find_infector_method` in each group). If `find_infector_method="random"` then it's a vector of zeros. 
        #' 
        #' @return A vector of `Host` objects
        #' 
        get_infectors = function(group_of_infectors, infector_intervals) {

            stopifnot(
                is.numeric(group_of_infectors),
                all(group_of_infectors %% 1 == 0),
                all(group_of_infectors %in% private$group_ids_)
            )
            stopifnot(
                is.numeric(infector_intervals),
                all(infector_intervals >= 0)
            )
            stopifnot(length(group_of_infectors) == length(infector_intervals))

            # loop through each group and select infectors
            # the method for selecting infectors is defined in the Group class
            infectors <- list()
            for (group in self$groups) {
                ind <- which(group_of_infectors == group$id)

                # check if any infectors are needed from this group
                if (length(ind) != 0) {
                    # currently Group$get_infectors() gets how many infectors to select from the length of the vector passed to it. It doesn't need any otehr  info

                    infectors[ind] <- group$get_infectors(infector_intervals[ind])
                }
            }

            return(infectors)
        },

        #' @description
        #' Select the group that each infector comes from
        #' 
        #' @param num_infectees An vector of integers indicating how many infectees to select from each group
        #' @param fois_mat A force of infection matrix of the same shape and form as inf_rates
        #' 
        #' @return A vector of vectors of group ids
        #' 
        get_group_of_infectors = function(num_infectees, fois_mat) {

            # store the number of groups (we'll use this a few times)
            n <- length(self$groups)

            stopifnot(
                is.numeric(num_infectees), 
                all(num_infectees %% 1 == 0),
                length(num_infectees) == n
            )
            stopifnot(
                nrow(fois_mat) == n,
                ncol(fois_mat) == n,
                is.numeric(fois_mat),
                all(fois_mat >= 0),
                all(fois_mat <= 1)
            )


            # get the group that the infector for each infectee will come from
            # for each group use the elements of it's corresponding column in fois_mat for the probabilities of the group of the infectors
            group_for_infectors <- lapply(
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

            return(group_for_infectors)
        },

        #' @description
        #' Determine number of infectees in each group
        #' 
        #' @param fois A vector of the force of infection experienced by each group
        #' 
        #' @return An vector of integers
        #' 
        get_num_infectees = function(fois) {
            stopifnot(
                length(fois) == length(self$groups),
                is.numeric(fois),
                all(fois >= 0)
            )

            num_infectes <- mapply(
                function(group, foi) group$get_num_infectees(foi),
                self$groups,
                fois
            )

            return(num_infectes)
        },

        #' @description
        #' Trigger an infection cycle in all `groups`. All new infectees (if any) are selected, their infectors are selected and infections of hosts are carried out.
        #' 
        #' @return This `Population` object
        #' 
        infect = function() {

            # get number of infectious for each group
            inf_sizes <- self$infectious_sizes

            # check if there are any infectives and stop if there isn't
            if (sum(inf_sizes) == 0) {return(invisible(self))}

            # get FOI matrix
            fois_mat <- private$force_of_infection(inf_sizes)
            
            # get new infectees for each group
            num_infectees_in_groups <- self$get_num_infectees(rowSums(fois_mat))
            infectees <- self$get_infectees(num_infectees_in_groups)


            # get the group that the infector for each infectee comes from
            group_for_infectors <- self$get_group_of_infectors(
                num_infectees_in_groups, 
                fois_mat
            )

            # group_for_infectors and infectees are both lists of lists of hosts
            # we can unlist them now
            infectees <- unlist(infectees, recursive = FALSE)
            group_for_infectors <- unlist(group_for_infectors, recursive = FALSE)

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
                group_for_infectors, 
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
        #' Print a description of this `Population` object
        #' 
        #' @param ... arguments will be ignored
        #' 
        #' @return No return value. Description is printed to the console
        #' 
        print = function() {

            population_data <- data.frame(
                Size = self$group_sizes,
                Susceptibles = self$susceptible_sizes,
                Exposed = self$exposed_sizes,
                Infectious = self$infectious_sizes,
                Recovered = self$recovered_sizes
            )

            
            cat("Outreak: \n")
            cat("   Time:                   ", self$time, "\n", sep = "")
            cat("   Samples Taken:          ", self$num_wgs, "\n", sep = "")
            cat("   Group Details: \n")
            writeLines(paste0("         ", capture.output(print(t(population_data)))))
        },

        #' @description
        #' Run this simulation.
        #' 
        #' @param feedback How frequently to print group summary to the console (default 10)
        #' 
        #' @return This `Population` object 
        #' 
        run_simulation = function(feedback = 5) {

            stopifnot(is.numeric(feedback), length(feedback) == 1, feedback %% 1 == 0)

            # first have lab perform samples due at t = 0
            self$sample_groups()
            
            while (self$is_outbreak_active) {

                # increment time for self and for each group
                self$time <- self$time + 1
                private$align_group_times()

                # trigger a round of infections
                self$infect()

                # have lab take samples
                self$sample_groups()

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
        #' Collect samples from `host`s in every group that are due for sampling
        #' 
        #' @return This `Population` object
        #' 
        sample_groups = function() {

            # get hosts that are due sampling

            hosts <- self$hosts_due_for_sampling
            self$lab$sample_hosts(hosts, self$time)

            return(invisible(self))
        }#,

    ),

    ###########################################################################
    # active bindings
    active = list(

        #' @field exposed_sizes The number of hosts in the exposed compartment in each group
        exposed_sizes = function(v) {

            if (missing(v)) {
                vapply(
                    self$groups,
                    function(group) {group$exposed_size},
                    integer(1)
                )
            } else {stop("Can't set `$ exposed_sizes`.", call. = FALSE)}
        },

        #' @field hosts_due_for_sampling A vector of hosts from all groups that are due sampling at the current time
        hosts_due_for_sampling = function(v) {
            
            if (missing(v)) {
                hosts <- lapply(
                    self$groups, 
                    function(group) group$hosts_due_for_sampling
                )
                hosts <- unlist(hosts, recursive = FALSE)

                return(hosts)
            } else {stop("Can't set `$get_hosts_due_for_sampling`", call. = FALSE)}
        },

        #' @field is_outbreak_active Logical indicating if there is any host in any group that is still in the exposed or infectious compartment
        is_outbreak_active = function(v) {
            if (missing(v)) {
                any(vapply(
                    self$groups, 
                    function(group) {group$is_outbreak_active},
                    logical(1L)
                ))
            } else {stop("Can't set `$is_outbreak_active`.", call. = FALSE)}
        },

        #' @field infectious_sizes The number of hosts in the infectious compartment in each group
        infectious_sizes = function(v) {
            # get number of infectious in each group
            if (missing(v)) {
                vapply(
                    self$groups,
                    function(group) {group$infectious_size},
                    integer(1)
                )
            } else {stop("Can't set `$ infectious_sizes`.", call. = FALSE)}
        },

        #' @field recovered_sizes The number of hosts in the recovered compartment in each group
        recovered_sizes = function(v) {
            # get number of recovered in each group
            if (missing(v)) {
                vapply(
                    self$groups,
                    function(group) {group$recovered_size},
                    integer(1)
                )
            } else {stop("Can't set `$recovered_sizes`.", call. = FALSE)}
        },

        #' @field susceptible_sizes The number of hosts in the susceptible compartment in each group
        susceptible_sizes = function(v) {
            # get number of susceptible in each group
            if (missing(v)) {
                vapply(
                    self$groups,
                    function(group) {group$susceptible_size},
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

        #' @field outbreak_size The number of non_susceptible hosts across all groups
        outbreak_size = function(v) {
            # get the outbreak size across all groups
            if (missing(v)) {
                sum(self$group_outbreak_sizes)
            } else {stop("Can't set `$outbreak_size`.", call. = FALSE)}
        },

        #' @field group_outbreak_sizes The number of non_susceptible hosts in each group
        group_outbreak_sizes = function(v) {
            # get the outbreak size across all groups
            if (missing(v)) {
                vapply(
                    self$groups,
                    function(group) {group$outbreak_size},
                    integer(1L)
                )
            } else {stop("Can't set `$group_outbreak_sizes`.", call. = FALSE)}
        },

        #' @field group_sizes The size of each group
        group_sizes = function(v) {
            # get size of each group
            if (missing(v)) {
                vapply(
                    self$groups,
                    function(group) {group$size},
                    double(1)
                )
            } else {stop("Can't set `$group_sizes`.", call. = FALSE)}
        }

    ),

    ###########################################################################
    # private members
    private = list( 

        group_ids_ = c(),

        #' @description
        #' align all group's time to this populations time
        #' 
        #' @return This `population` object 
        #'
        #' @noRd 
        align_group_times = function() {
            # align each Group's time to Population's time
            vapply(
                self$groups,
                function(group) {
                    group$time <- self$time
                    TRUE
                },
                logical(1L)
            )

            return(invisible(self))
        },

        #' @description
        #' Calculate the force_of_infectious matrix
        #' 
        #' @param inf_sizes The number of infectious hosts in each group
        #' 
        #' @return A matrix
        #' 
        #' @noRd 
        force_of_infection = function(inf_sizes) {
            # calculate the force of infection matrix. 
            # Each element will be of the form beta_{ij}*I_j/N 

            beta_I <- sweep(self$inf_rates, MARGIN = 2, inf_sizes, "*") 
            foi_mat <- beta_I / sum(self$group_sizes)

            return(foi_mat)
        }
    )
)
