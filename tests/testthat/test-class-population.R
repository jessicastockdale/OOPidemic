###############################################################################
# Method: Population$new()
###############################################################################

test_that("Population$new() errors when invalid id given", {
    
    expect_error(Population$new(id = "a"))
    expect_error(Population$new(id = 1.5))
    expect_error(Population$new(id = c(1, 2)))
})

test_that("Population$new() errors when invalid ref_strain object given", {

    Class <- R6::R6Class("Class")
    class <- Class$new()

    expect_error(Population$new(1, ref_strain = class))
    expect_error(Population$new(1, ref_strain = 3))
})

test_that("Population$new() errors when invalid init_sus given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Population$new(1, ref_strain, init_sus = "a"))
    expect_error(Population$new(1, ref_strain, init_sus = 1.5))
    expect_error(Population$new(1, ref_strain, init_sus = c(1, 2)))
})

test_that("Population$new() errors when invalid init_inf given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Population$new(1, ref_strain, init_inf = "a"))
    expect_error(Population$new(1, ref_strain, init_inf = 1.5))
    expect_error(Population$new(1, ref_strain, init_inf = c(1, 2)))
})

test_that("Population$new() errors when invalid inf_rate given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Population$new(1, ref_strain, inf_rate = "a"))
    expect_error(Population$new(1, ref_strain, inf_rate = 0))
    expect_error(Population$new(1, ref_strain, inf_rate = c(1, 2)))
})

test_that("Population$new() errors when invalid find_infector_method given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Population$new(1, ref_strain, find_infector_method = "a"))
})

test_that("Population$new() errors when invalid transmission interval gamma parameters given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Population$new(1, ref_strain, trans_int_shape = "a"))
    expect_error(Population$new(1, ref_strain, trans_int_shape = 0))
    expect_error(Population$new(1, ref_strain, trans_int_shape = c(1, 2)))
    expect_error(Population$new(1, ref_strain, trans_int_rate = "a"))
    expect_error(Population$new(1, ref_strain, trans_int_rate = 0))
    expect_error(Population$new(1, ref_strain, trans_int_rate = c(1, 2)))
})

test_that("Population$new() errors when invalid serial interval gamma parameters given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Population$new(
        1, ref_strain, 
        find_infector_method = "serial",
        si_shape = "a"
    ))
    expect_error(Population$new(
        1, ref_strain, 
        find_infector_method = "serial",
        si_shape = 0
    ))
    expect_error(Population$new(
        1, ref_strain, 
        find_infector_method = "serial",
        si_shape = c(1, 2)))
    expect_error(Population$new(
        1, ref_strain, 
        find_infector_method = "serial",
        si_shape = 6, si_rate = "a"
    ))
    expect_error(Population$new(
        1, ref_strain, 
        find_infector_method = "serial",
        si_shape = 6, si_rate = 0
    ))
    expect_error(Population$new(
        1, ref_strain, 
        find_infector_method = "serial",
        si_shape = 6, si_rate = c(1, 2)))
})

test_that("Population$new() errors when invalid generation_interval gamma parameters given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Population$new(
        1, ref_strain, 
        find_infector_method = "generation",
        gt_shape = "a"
    ))
    expect_error(Population$new(
        1, ref_strain, 
        find_infector_method = "generation",
        gt_shape = 0
    ))
    expect_error(Population$new(
        1, ref_strain, 
        find_infector_method = "generation",
        gt_shape = c(1, 2)))
    expect_error(Population$new(
        1, ref_strain, 
        find_infector_method = "generation",
        gt_shape = 6, gt_rate = "a"
    ))
    expect_error(Population$new(
        1, ref_strain, 
        find_infector_method = "generation",
        gt_shape = 6, gt_rate = 0
    ))
    expect_error(Population$new(
        1, ref_strain, 
        find_infector_method = "generation",
        gt_shape = 6, gt_rate = c(1, 2)))
})

test_that("Population$new() errors when invalid incubation gamma parameters given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Population$new(
        1, ref_strain, 
        inc_shape = "a"
    ))
    expect_error(Population$new(
        1, ref_strain, 
        inc_shape = -1
    ))
    expect_error(Population$new(
        1, ref_strain, 
        inc_shape = c(1, 2)
    ))
    expect_error(Population$new(
        1, ref_strain, 
        inc_shape = 6, inc_rate = "a"
    ))
    expect_error(Population$new(
        1, ref_strain, 
        inc_shape = 6, inc_rate = 0
    ))
    expect_error(Population$new(
        1, ref_strain, 
        inc_shape = 6, inc_rate = c(1, 2)
    ))
})

test_that("Population$new() errors when invalid recovery gamma parameters given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Population$new(
        1, ref_strain, 
        rec_shape = "a"
    ))
    expect_error(Population$new(
        1, ref_strain, 
        rec_shape = 0
    ))
    expect_error(Population$new(
        1, ref_strain, 
        rec_shape = c(1, 2)
    ))
    expect_error(Population$new(
        1, ref_strain, 
        rec_shape = 6, rec_rate = "a"
    ))
    expect_error(Population$new(
        1, ref_strain, 
        rec_shape = 6, rec_rate = 0
    ))
    expect_error(Population$new(
        1, ref_strain, 
        rec_shape = 6, rec_rate = c(1, 2)
    ))
})

test_that("Population$new() errors when invalid initial distances given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Population$new(
        1, ref_strain, 
        max_init_dist = "a"
    ))
    expect_error(Population$new(
        1, ref_strain, 
        max_init_dist = 1.5
    ))
    expect_error(Population$new(
        1, ref_strain, 
        max_init_dist = c(1, 2)
    ))

    expect_error(Population$new(
        1, ref_strain, 
        max_init_dist = 200, min_init_dist = "a"
    ))
    expect_error(Population$new(
        1, ref_strain, 
        max_init_dist = 200, min_init_dist = 1.5
    ))
    expect_error(Population$new(
        1, ref_strain, 
        max_init_dist = 200, min_init_dist = c(1, 2)
    ))

    expect_error(Population$new(
        1, ref_strain, 
        max_init_dist = 20, min_init_dist = 21
    ))

})

test_that("Population$new() errors when invalid sample schedule and frequency given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Population$new(
        1, ref_strain, 
        sample_schedule = "a"
    ))

    expect_error(Population$new(
        1, ref_strain, 
        sample_schedule = "individual",
        sample_frequency = "a"
    ))
    expect_error(Population$new(
        1, ref_strain, 
        sample_schedule = "individual",
        sample_frequency = 1.5
    ))
    expect_error(Population$new(
        1, ref_strain, 
        sample_schedule = "individual",
        sample_frequency = c(1, 2)
    ))

})

test_that("Population$new() errors when invalid Host class given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    Class <- R6::R6Class("Class")

    expect_error(Population$new(
        1, ref_strain, 
        host_class = Class
    ))

})

test_that("Population$new() errors when invalid Strain class given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    Class <- R6::R6Class("Class")

    expect_error(Population$new(
        1, ref_strain, 
        strain_class = Class
    ))

})

test_that("Population$new() successfully completes", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    TestP <- R6::R6Class("TestP", inherit = Population,
        active = list(interval_stack = function() private$interval_stack_)
    )
    TestH <- R6::R6Class("TestH", inherit = Host,
        public = list(
            calendar_sample_time_ = function() return(private$calendar_sample_time())
        )
    )
    TestS <- R6::R6Class("TestS", inherit = Strain,
        public = list(
            calendar_sample_time_ = function() return(private$calendar_sample_time())
        )
    )


    init_sus <- 95
    init_inf <- 5
    pop <- TestP$new(
        id = 1,
        ref_strain,
        host_class = TestH,
        strain_class = TestS
    )

    expect_false(is.null(pop$interval_stack))
    expect_length(pop$hosts, init_sus + init_inf)
})

###############################################################################
# Method: Population$exposed_hosts()
###############################################################################

test_that("Population$exposed_hosts() errors when invalid time given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)


    expect_error(pop$exposed_hosts("a"))
    expect_error(pop$exposed_hosts(1.5))
    expect_error(pop$exposed_hosts(c(1, 2)))
})

test_that("Population$exposed_hosts() successfully completes", {
    
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        pop <- Population$new(1, ref_strain, inf_rate = 2)
        # trigger an infection so there are exposed hosts
        pop$infect()

        expected <- c()
        for (host in pop$hosts) {
            if (host$is_exposed(pop$time)) {
                expected <- c(expected, host)
            }
        }

        expect_setequal(pop$exposed_hosts(), expected)
    })
})

###############################################################################
# Method: Population$get_infectees()
###############################################################################

test_that("Population$get_infectees() errors when invalid num_infectees given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)


    expect_error(pop$get_infectees("a"))
    expect_error(pop$get_infectees(1.5))
    expect_error(pop$get_infectees(c(1, 2)))
})

test_that("Population$get_infectees() successfully completes", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)

    num_infectees <- 5
    infectees <- pop$get_infectees(num_infectees)

    expect_length(infectees, num_infectees)
    expect_length(unique(infectees), num_infectees) # check there's no duplicates
    # check they are all susceptible
    expect_true(all(
        vapply(
            infectees, 
            function(infectee) infectee$is_susceptible(pop$time), 
            logical(1L)
        )
    ))
})

###############################################################################
# Method: Population$get_infectors()
###############################################################################

test_that("Population$get_infectors() errors when invalid num_infectees given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)


    expect_error(pop$get_infectors("a"))
    expect_error(pop$get_infectors(c()))
})

test_that("Population$get_infectors() successfully completes with a single infective", {
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        pop <- Population$new(1, ref_strain)

        intervals <- c(1:5)
        infectors <- pop$get_infectors(intervals)
        expect_length(unique(infectors), 1)
        expect_true(unique(infectors)[[1]]$is_infectious(pop$time))
    })
})

test_that("Population$get_infectors() successfully completes using random infector method", {
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        init_inf <- 2
        pop <- Population$new(
            1, ref_strain, 
            init_inf = init_inf,
            find_infector_method = "random"
        )

        intervals <- c(1:5)
        infectors <- pop$get_infectors(intervals)
        expect_length(infectors, length(intervals))
        expect_length(unique(infectors), init_inf)
        expect_true(all(
            vapply(
                infectors,
                function(infector) infector$is_infectious(pop$time),
                logical(1L)
            )
        ))
    })
})

test_that("Population$get_infectors() successfully completes using default transmission infector method", {
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        init_inf <- 2
        pop <- Population$new(
            1, ref_strain, 
            init_inf = init_inf
        )

        intervals <- c(5, 6)

        expected_infectors <- c()
        inf_times <- c()
        rec_times <- c()

        # prepare two infectious hosts that aren't infected at time 0
        for (i in seq_along(intervals)) {
            infectee <- pop$susceptible_hosts()[[i]]
            infectee$prepare_for_infection(intervals[i])
            pop$infectious_hosts()[[1]]$infect(infectee, intervals[i])

            expected_infectors <- c(expected_infectors, infectee)
            inf_times <- c(inf_times, infectee$infectious_time)
            rec_times <- c(rec_times, infectee$recovery_time)
        }

        pop$time <- min(rec_times - 1)
        infectors <- pop$get_infectors(pop$time - inf_times)

        expect_setequal(infectors, expected_infectors)
    })
})

###############################################################################
# Method: Population$get_num_infectees()
###############################################################################

test_that("Population$get_num_infectees() errors when invalid foi given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)


    expect_error(pop$get_num_infectees("a"))
    expect_error(pop$get_num_infectees(c(1, 2)))
})

test_that("Population$get_num_infectees() successfully completes", {
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        init_sus <- 100
        pop <- Population$new(1, ref_strain, init_sus = init_sus)

        num_infectees <- pop$get_num_infectees(0.2)

        expect_true(is.integer(num_infectees))
        expect_length(num_infectees, 1)
        expect_true(0 <= num_infectees)
        expect_true(num_infectees <= init_sus)
    })
})

###############################################################################
# Method: Population$infect()
###############################################################################

test_that("Population$infect() successfully completes with a single infective", {
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        init_sus <- 100
        pop <- Population$new(1, ref_strain, init_sus = init_sus)

        expected_host_ids <- c(73)
        pop$infect()

        expect_length(pop$exposed_hosts(), length(expected_host_ids))
        ids <- vapply(pop$exposed_hosts(), function(host) host$id, integer(1L))
        expect_setequal(
            ids,
            expected_host_ids
        )
        expect_equal(pop$time, 1)
    })
})

###############################################################################
# Method: Population$infectious_hosts()
###############################################################################

test_that("Population$infectious_hosts() errors when invalid time given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)


    expect_error(pop$infectious_hosts("a"))
    expect_error(pop$infectious_hosts(1.5))
    expect_error(pop$infectious_hosts(c(1, 2)))
})

test_that("Population$infectious_hosts() successfully completes", {
    
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        pop <- Population$new(1, ref_strain, inc_shape = 0)
        # trigger an infection so there are exposed hosts
        pop$infect()

        expected <- c()
        for (host in pop$hosts) {
            if (host$is_infectious(pop$time)) {
                expected <- c(expected, host)
            }
        }

        expect_setequal(pop$infectious_hosts(), expected)
    })
})

###############################################################################
# Method: Population$recovered_hosts()
###############################################################################

test_that("Population$recovered_hosts() errors when invalid time given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)


    expect_error(pop$recovered_hosts("a"))
    expect_error(pop$recovered_hosts(1.5))
    expect_error(pop$recovered_hosts(c(1, 2)))
})

test_that("Population$recovered_hosts() successfully completes", {
    
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        init_inf <- 5
        pop <- Population$new(1, ref_strain, init_inf = init_inf)

        index_hosts <- pop$infectious_hosts()
        recovery_times <- vapply(
            index_hosts,
            function(host) host$recovery_time,
            numeric(1L)
        )

        expect_setequal(pop$recovered_hosts(max(recovery_times)), index_hosts)
    })
})

###############################################################################
# Method: Population$run_simulation()
###############################################################################

test_that("Population$run_simluation() errors when invalid lab given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)
    Class <- R6::R6Class("Class")
    class <- Class$new()

    expect_error(pop$run_simluation(class))
})

test_that("Population$run_simluation() errors when invalid time given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)
    lab <- Lab$new()

    expect_error(pop$run_simluation(lab, feedback = "a"))
    expect_error(pop$run_simluation(lab, feedback = c(1, 2)))
    expect_error(pop$run_simluation(lab, feedback = 1.5))
    expect_error(pop$run_simluation(lab, feedback = -1.5))
})

test_that("Population$run_simluation() errors when invalid time given", {
    
    withr::with_seed(1234, {
        ref_strain <- ReferenceStrain$new("ref_strain")
        pop <- Population$new(1, ref_strain)
        lab <- Lab$new()

        pop$run_simulation(lab, 0)
        
        expect_gt(pop$time, 0)
        expect_gt(lab$num_wgs, 0)
        expect_gt(pop$recovered_size, 0)
    })
})

###############################################################################
# Method: Population$susceptible_hosts()
###############################################################################

test_that("Population$susceptible_hosts() errors when invalid time given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)


    expect_error(pop$susceptible_hosts("a"))
    expect_error(pop$susceptible_hosts(1.5))
    expect_error(pop$susceptible_hosts(c(1, 2)))
})

test_that("Population$susceptible_hosts() successfully completes", {
    
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        init_sus <- 5
        pop <- Population$new(1, ref_strain, init_sus = init_sus)

        expect_length(pop$susceptible_hosts(), init_sus)
        expect_true(all(vapply(
            pop$susceptible_hosts(),
            function(host) host$is_susceptible(pop$time),
            logical(1L)
        )))
    })
})

###############################################################################
# Active Binding: Population$hosts_due_for_sampling
###############################################################################

test_that("Population$hosts_due_for_sampling successfully completes", {
    
    withr::with_seed(1234, {
        ref_strain <- ReferenceStrain$new("ref_strain")
        init_inf <- 5
        sample_freq <- 5
        pop <- Population$new(
            1, ref_strain, 
            init_inf = init_inf,
            inc_shape = 0, # turn off exposure compartment
            rec_shape = 20, # force a very long recovery time 
            sample_schedule = "calendar", sample_freq = sample_freq
        )

        pop$time <- sample_freq

        expect_equal(pop$infectious_hosts(), pop$hosts_due_for_sampling)
    })
})

###############################################################################
# Active Binding: Population$is_outbreak_active
###############################################################################

test_that("Population$is_outbreak_active successfully completes", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    init_inf <- 5
    pop <- Population$new(
        1, ref_strain, 
        init_inf = init_inf
    )

    time <- 200
    infector <- pop$infectious_hosts()[[1]]
    infectee <- pop$susceptible_hosts()[[1]]
    
    # perform infection
    infectee$prepare_for_infection(time)
    infector$infect(infectee, time)

    pop$time <- time - 1
    expect_false(pop$is_outbreak_active)
    pop$time <- time
    expect_true(pop$is_outbreak_active)
    pop$time <- infectee$recovery_time - 1
    expect_true(pop$is_outbreak_active)
    pop$time <- infectee$recovery_time
    expect_false(pop$is_outbreak_active)
})

###############################################################################
# Active Binding: Population$pop_interval_stack
###############################################################################

test_that("Population$pop_interval_stack successfully completes and then errors when empty", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    init_sus <- 9
    init_inf <- 1
    pop <- Population$new(
        1, ref_strain,
        init_sus = init_sus, init_inf = init_inf
    )

    for (i in seq(1:init_inf + init_sus)) {
        intervals <- pop$pop_interval_stack
    }
    expect_length(intervals, 2)
    expect_setequal(names(intervals), c("inc", "infector_interval"))
    expect_error(pop$pop_interval_stack)    
})

###############################################################################
# Active Binding: Population$time_since_infectious_hosts_infected
###############################################################################

test_that("Population$time_since_infectious_hosts_infected successfully completes", {

    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        init_inf <- 5
        pop <- Population$new(
            1, ref_strain,
            init_inf = init_inf
        )

        # all should be zero since all infectious hosts are index cases
        expect_equal(pop$time_since_infectious_hosts_infected, rep(0, init_inf))

        # create new infections
        pop$infect()
        # get latest time that a new case becomes infectious
        inf_times <- vapply(
            pop$exposed_hosts(),
            function(host) host$infectious_time,
            numeric(1L)
        )
        pop$time <- max(inf_times)
        expected <- vapply(
            pop$infectious_hosts(),
            function(host) max(inf_times) - host$exposure_time,
            numeric(1L)
        )

        expect_equal(pop$time_since_infectious_hosts_infected, expected)

    })
})

###############################################################################
# Active Binding: Population$time_since_infectious_hosts_infectious
###############################################################################

test_that("Population$time_since_infectious_hosts_infectious successfully completes", {

    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        init_inf <- 5
        pop <- Population$new(
            1, ref_strain,
            init_inf = init_inf
        )

        # all should be zero since all infectious hosts are index cases
        expect_equal(pop$time_since_infectious_hosts_infectious, rep(0, init_inf))

        # create new infections
        pop$infect()
        # get latest time that a new case becomes infectious
        inf_times <- vapply(
            pop$exposed_hosts(),
            function(host) host$infectious_time,
            numeric(1L)
        )
        pop$time <- max(inf_times)
        expected <- vapply(
            pop$infectious_hosts(),
            function(host) max(inf_times) - host$infectious_time,
            numeric(1L)
        )

        expect_equal(pop$time_since_infectious_hosts_infectious, expected)

    })
})

###############################################################################
# Private Member: Population$build_interval_stack()
###############################################################################

test_that("Population$build_interval_stack() successfully completes with serial find_infector_method", {

    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        TestP <- R6::R6Class("TestP", inherit = Population,
            active = list(
                interval_stack = function() return(private$interval_stack_)
            )
        )
        pop <- TestP$new(
            1, ref_strain, 
            find_infector_method = "serial",
            si_shape = 8, si_rate = 2    
        )
        stack <- pop$interval_stack

        expect_equal(nrow(stack), pop$size)
        expect_equal(colnames(stack), c("si", "inc", "infector_interval"))
        expect_true(all(stack$infector_interval == stack$si - stack$inc))
    })
})

test_that("Population$build_interval_stack() successfully completes with transmission find_infector_method", {

    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        TestP <- R6::R6Class("TestP", inherit = Population,
            active = list(
                interval_stack = function() return(private$interval_stack_)
            )
        )
        pop <- TestP$new(
            1, ref_strain, 
            find_infector_method = "transmission",
            trans_int_shape = 3, trans_int_rate = 1    
        )
        stack <- pop$interval_stack

        expect_equal(nrow(stack), pop$size)
        expect_equal(colnames(stack), c("inc", "infector_interval"))
    })
})

test_that("Population$build_interval_stack() successfully completes with random find_infector_method", {

    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        TestP <- R6::R6Class("TestP", inherit = Population,
            active = list(
                interval_stack = function() return(private$interval_stack_)
            )
        )
        pop <- TestP$new(1, ref_strain, find_infector_method = "random")
        stack <- pop$interval_stack

        expect_equal(nrow(stack), pop$size)
        expect_equal(colnames(stack), c("inc", "infector_interval"))
        expect_true(all(stack$infector_interval == 0))
    })
})

test_that("Population$build_interval_stack() successfully completes with generation find_infector_method", {

    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        TestP <- R6::R6Class("TestP", inherit = Population,
            active = list(
                interval_stack = function() return(private$interval_stack_)
            )
        )
        pop <- TestP$new(
            1, ref_strain, 
            find_infector_method = "generation",
            gt_shape = 8, gt_rate = 3    
        )
        stack <- pop$interval_stack

        expect_equal(nrow(stack), pop$size)
        expect_equal(colnames(stack), c("inc", "infector_interval"))
    })
})

###############################################################################
# Private Member: Population$initalise_index_host()
###############################################################################

test_that("Population$initalise_index_host() successfully completes with variation TRUE", {


    ref_strain <- ReferenceStrain$new("ref_strain")
    TestP <- R6::R6Class("TestP", inherit = Population,
        public = list(
            initialise_index_host_ = function(...) private$initialise_index_host(...)
        )
    )
    distance <- 5
    pop <- TestP$new(
        1, ref_strain,
        min_init_dist = distance, max_init_dist = distance
    )
    host <- pop$susceptible_hosts()[[1]]
    pop$initialise_index_host_(host)

    expect_false(host$is_susceptible(pop$time))
    expect_false(all(ref_strain$genome() == host$strains[[1]]$genome()))
})

test_that("Population$initalise_index_host() successfully completes with variation FALSE", {


    ref_strain <- ReferenceStrain$new("ref_strain")
    TestP <- R6::R6Class("TestP", inherit = Population,
        public = list(
            initialise_index_host_ = function(...) private$initialise_index_host(...)
        )
    )
    distance <- 5
    pop <- TestP$new(
        1, ref_strain,
        min_init_dist = distance, max_init_dist = distance
    )
    host <- pop$susceptible_hosts()[[1]]
    pop$initialise_index_host_(host, FALSE)

    expect_false(host$is_susceptible(pop$time))
    expect_true(all(ref_strain$genome() == host$strains[[1]]$genome()))
})

###############################################################################
# Private Member: Population$initialise_hosts()
###############################################################################

test_that("Population$initialise_hosts() successfully completes with one initial infective", {


    ref_strain <- ReferenceStrain$new("ref_strain")
    init_inf <- 1
    init_sus <- 5
    pop <- Population$new(
        1, ref_strain,
        init_sus = init_sus, init_inf = init_inf
    )
    
    expect_length(pop$hosts, init_inf + init_sus)
    expect_length(pop$infectious_hosts(), init_inf)
    expect_true(all(
        ref_strain$genome() == pop$infectious_hosts()[[1]]$strains[[1]]$genome()
    ))

})

test_that("Population$initialise_hosts() successfully completes with multiple initial infective", {


    ref_strain <- ReferenceStrain$new("ref_strain")
    init_inf <- 5
    init_sus <- 5
    pop <- Population$new(
        1, ref_strain,
        init_sus = init_sus, init_inf = init_inf
    )
    
    expect_length(pop$hosts, init_inf + init_sus)
    expect_length(pop$infectious_hosts(), init_inf)
})

###############################################################################
# Private Member: Population$force_of_infection()
###############################################################################

test_that("Population$force_of_infection() successfully completes", {


    ref_strain <- ReferenceStrain$new("ref_strain")
    TestP <- R6::R6Class("TestP", inherit = Population,
        public = list(
            force_of_infection_ = function(...) private$force_of_infection(...)
        )
    )
    inf_rate <- 0.1
    init_inf <- 5
    init_sus <- 95
   
    pop <- TestP$new(
        1, ref_strain,
        init_sus = init_sus, init_inf = init_inf,
        inf_rate = inf_rate
    )

    expected <- inf_rate * init_inf / (init_inf + init_sus)
    expect_equal(pop$force_of_infection_(), expected)
})

###############################################################################
# Private Member: Population$initial_genomic_distance()
###############################################################################

test_that("Population$initial_genomic_distance() successfully completes", {


    ref_strain <- ReferenceStrain$new("ref_strain")
    TestP <- R6::R6Class("TestP", inherit = Population,
        public = list(
            initial_genomic_distance_ = function(...) private$initial_genomic_distance(...)
        )
    )
    min_init_dist <- 0
    max_init_dist <- 20
   
    pop <- TestP$new(
        1, ref_strain,
        max_init_dist = max_init_dist,
        min_init_dist = min_init_dist
    )

    distances <- vapply(
        seq(50),
        function(i) pop$initial_genomic_distance_(),
        numeric(1L)
    )

    expect_true(all(distances >= min_init_dist))
    expect_true(all(distances <= max_init_dist))
})

