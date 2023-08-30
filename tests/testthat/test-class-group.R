###############################################################################
# Method: Group$new()
###############################################################################

test_that("Group$new() errors when invalid id given", {
    
    expect_error(Group$new(id = "a"))
    expect_error(Group$new(id = 1.5))
    expect_error(Group$new(id = c(1, 2)))
})

test_that("Group$new() errors when invalid ref_strain object given", {

    Class <- R6::R6Class("Class")
    class <- Class$new()

    expect_error(Group$new(1, ref_strain = class))
    expect_error(Group$new(1, ref_strain = 3))
})

test_that("Group$new() errors when invalid init_sus given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Group$new(1, ref_strain, init_sus = "a"))
    expect_error(Group$new(1, ref_strain, init_sus = 1.5))
    expect_error(Group$new(1, ref_strain, init_sus = c(1, 2)))
})

test_that("Group$new() errors when invalid init_inf given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Group$new(1, ref_strain, init_inf = "a"))
    expect_error(Group$new(1, ref_strain, init_inf = 1.5))
    expect_error(Group$new(1, ref_strain, init_inf = c(1, 2)))
})

test_that("Group$new() errors when invalid inf_rate given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Group$new(1, ref_strain, inf_rate = "a"))
    expect_error(Group$new(1, ref_strain, inf_rate = 0))
    expect_error(Group$new(1, ref_strain, inf_rate = c(1, 2)))
})

test_that("Group$new() errors when invalid find_infector_method given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Group$new(1, ref_strain, find_infector_method = "a"))
})

test_that("Group$new() errors when invalid transmission interval gamma parameters given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Group$new(1, ref_strain, trans_int_shape = "a"))
    expect_error(Group$new(1, ref_strain, trans_int_shape = 0))
    expect_error(Group$new(1, ref_strain, trans_int_shape = c(1, 2)))
    expect_error(Group$new(1, ref_strain, trans_int_rate = "a"))
    expect_error(Group$new(1, ref_strain, trans_int_rate = 0))
    expect_error(Group$new(1, ref_strain, trans_int_rate = c(1, 2)))
})

test_that("Group$new() errors when invalid serial interval gamma parameters given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Group$new(
        1, ref_strain, 
        find_infector_method = "serial",
        si_shape = "a"
    ))
    expect_error(Group$new(
        1, ref_strain, 
        find_infector_method = "serial",
        si_shape = 0
    ))
    expect_error(Group$new(
        1, ref_strain, 
        find_infector_method = "serial",
        si_shape = c(1, 2)))
    expect_error(Group$new(
        1, ref_strain, 
        find_infector_method = "serial",
        si_shape = 6, si_rate = "a"
    ))
    expect_error(Group$new(
        1, ref_strain, 
        find_infector_method = "serial",
        si_shape = 6, si_rate = 0
    ))
    expect_error(Group$new(
        1, ref_strain, 
        find_infector_method = "serial",
        si_shape = 6, si_rate = c(1, 2)))
})

test_that("Group$new() errors when invalid generation_interval gamma parameters given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Group$new(
        1, ref_strain, 
        find_infector_method = "generation",
        gt_shape = "a"
    ))
    expect_error(Group$new(
        1, ref_strain, 
        find_infector_method = "generation",
        gt_shape = 0
    ))
    expect_error(Group$new(
        1, ref_strain, 
        find_infector_method = "generation",
        gt_shape = c(1, 2)))
    expect_error(Group$new(
        1, ref_strain, 
        find_infector_method = "generation",
        gt_shape = 6, gt_rate = "a"
    ))
    expect_error(Group$new(
        1, ref_strain, 
        find_infector_method = "generation",
        gt_shape = 6, gt_rate = 0
    ))
    expect_error(Group$new(
        1, ref_strain, 
        find_infector_method = "generation",
        gt_shape = 6, gt_rate = c(1, 2)))
})

test_that("Group$new() errors when invalid incubation gamma parameters given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Group$new(
        1, ref_strain, 
        inc_shape = "a"
    ))
    expect_error(Group$new(
        1, ref_strain, 
        inc_shape = -1
    ))
    expect_error(Group$new(
        1, ref_strain, 
        inc_shape = c(1, 2)
    ))
    expect_error(Group$new(
        1, ref_strain, 
        inc_shape = 6, inc_rate = "a"
    ))
    expect_error(Group$new(
        1, ref_strain, 
        inc_shape = 6, inc_rate = 0
    ))
    expect_error(Group$new(
        1, ref_strain, 
        inc_shape = 6, inc_rate = c(1, 2)
    ))
})

test_that("Group$new() errors when invalid recovery gamma parameters given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Group$new(
        1, ref_strain, 
        rec_shape = "a"
    ))
    expect_error(Group$new(
        1, ref_strain, 
        rec_shape = 0
    ))
    expect_error(Group$new(
        1, ref_strain, 
        rec_shape = c(1, 2)
    ))
    expect_error(Group$new(
        1, ref_strain, 
        rec_shape = 6, rec_rate = "a"
    ))
    expect_error(Group$new(
        1, ref_strain, 
        rec_shape = 6, rec_rate = 0
    ))
    expect_error(Group$new(
        1, ref_strain, 
        rec_shape = 6, rec_rate = c(1, 2)
    ))
})

test_that("Group$new() errors when invalid initial distances given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Group$new(
        1, ref_strain, 
        max_init_dist = "a"
    ))
    expect_error(Group$new(
        1, ref_strain, 
        max_init_dist = 1.5
    ))
    expect_error(Group$new(
        1, ref_strain, 
        max_init_dist = c(1, 2)
    ))

    expect_error(Group$new(
        1, ref_strain, 
        max_init_dist = 200, min_init_dist = "a"
    ))
    expect_error(Group$new(
        1, ref_strain, 
        max_init_dist = 200, min_init_dist = 1.5
    ))
    expect_error(Group$new(
        1, ref_strain, 
        max_init_dist = 200, min_init_dist = c(1, 2)
    ))

    expect_error(Group$new(
        1, ref_strain, 
        max_init_dist = 20, min_init_dist = 21
    ))

})

test_that("Group$new() errors when invalid sample schedule and frequency given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Group$new(
        1, ref_strain, 
        sample_schedule = "a"
    ))

    expect_error(Group$new(
        1, ref_strain, 
        sample_schedule = "individual",
        sample_frequency = "a"
    ))
    expect_error(Group$new(
        1, ref_strain, 
        sample_schedule = "individual",
        sample_frequency = 1.5
    ))
    expect_error(Group$new(
        1, ref_strain, 
        sample_schedule = "individual",
        sample_frequency = c(1, 2)
    ))

})

test_that("Group$new() errors when invalid Host class given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    Class <- R6::R6Class("Class")

    expect_error(Group$new(
        1, ref_strain, 
        host_class = Class
    ))

})

test_that("Group$new() errors when invalid Strain class given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    Class <- R6::R6Class("Class")

    expect_error(Group$new(
        1, ref_strain, 
        strain_class = Class
    ))

})

test_that("Group$new() successfully completes", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    TestP <- R6::R6Class("TestP", inherit = Group,
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
    grp <- TestP$new(
        id = 1,
        ref_strain,
        host_class = TestH,
        strain_class = TestS
    )

    expect_false(is.null(grp$interval_stack))
    expect_length(grp$hosts, init_sus + init_inf)
})

###############################################################################
# Method: Group$exposed_hosts()
###############################################################################

test_that("Group$exposed_hosts() errors when invalid time given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    grp <- Group$new(1, ref_strain)


    expect_error(grp$exposed_hosts("a"))
    expect_error(grp$exposed_hosts(1.5))
    expect_error(grp$exposed_hosts(c(1, 2)))
})

test_that("Group$exposed_hosts() successfully completes", {
    
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        grp <- Group$new(1, ref_strain, inf_rate = 2)
        # trigger an infection so there are exposed hosts
        grp$infect()

        expected <- c()
        for (host in grp$hosts) {
            if (host$is_exposed(grp$time)) {
                expected <- c(expected, host)
            }
        }

        expect_setequal(grp$exposed_hosts(), expected)
    })
})

###############################################################################
# Method: Group$get_infectees()
###############################################################################

test_that("Group$get_infectees() errors when invalid num_infectees given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    grp <- Group$new(1, ref_strain)


    expect_error(grp$get_infectees("a"))
    expect_error(grp$get_infectees(1.5))
    expect_error(grp$get_infectees(c(1, 2)))
})

test_that("Group$get_infectees() successfully completes", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    grp <- Group$new(1, ref_strain)

    num_infectees <- 5
    infectees <- grp$get_infectees(num_infectees)

    expect_length(infectees, num_infectees)
    expect_length(unique(infectees), num_infectees) # check there's no duplicates
    # check they are all susceptible
    expect_true(all(
        vapply(
            infectees, 
            function(infectee) infectee$is_susceptible(grp$time), 
            logical(1L)
        )
    ))
})

###############################################################################
# Method: Group$get_infectors()
###############################################################################

test_that("Group$get_infectors() errors when invalid num_infectees given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    grp <- Group$new(1, ref_strain)


    expect_error(grp$get_infectors("a"))
    expect_error(grp$get_infectors(c()))
})

test_that("Group$get_infectors() successfully completes with a single infective", {
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        grp <- Group$new(1, ref_strain)

        intervals <- c(1:5)
        infectors <- grp$get_infectors(intervals)
        expect_length(unique(infectors), 1)
        expect_true(unique(infectors)[[1]]$is_infectious(grp$time))
    })
})

test_that("Group$get_infectors() successfully completes using random infector method", {
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        init_inf <- 2
        grp <- Group$new(
            1, ref_strain, 
            init_inf = init_inf,
            find_infector_method = "random"
        )

        intervals <- c(1:5)
        infectors <- grp$get_infectors(intervals)
        expect_length(infectors, length(intervals))
        expect_length(unique(infectors), init_inf)
        expect_true(all(
            vapply(
                infectors,
                function(infector) infector$is_infectious(grp$time),
                logical(1L)
            )
        ))
    })
})

test_that("Group$get_infectors() successfully completes using default transmission infector method", {
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        init_inf <- 2
        grp <- Group$new(
            1, ref_strain, 
            init_inf = init_inf
        )

        intervals <- c(5, 6)

        expected_infectors <- c()
        inf_times <- c()
        rec_times <- c()

        # prepare two infectious hosts that aren't infected at time 0
        for (i in seq_along(intervals)) {
            infectee <- grp$susceptible_hosts()[[i]]
            infectee$prepare_for_infection(intervals[i])
            grp$infectious_hosts()[[1]]$infect(infectee, intervals[i])

            expected_infectors <- c(expected_infectors, infectee)
            inf_times <- c(inf_times, infectee$infectious_time)
            rec_times <- c(rec_times, infectee$recovery_time)
        }

        grp$time <- min(rec_times - 1)
        infectors <- grp$get_infectors(grp$time - inf_times)

        expect_setequal(infectors, expected_infectors)
    })
})

###############################################################################
# Method: Group$get_num_infectees()
###############################################################################

test_that("Group$get_num_infectees() errors when invalid foi given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    grp <- Group$new(1, ref_strain)


    expect_error(grp$get_num_infectees("a"))
    expect_error(grp$get_num_infectees(c(1, 2)))
})

test_that("Group$get_num_infectees() successfully completes", {
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        init_sus <- 100
        grp <- Group$new(1, ref_strain, init_sus = init_sus)

        num_infectees <- grp$get_num_infectees(0.2)

        expect_true(is.integer(num_infectees))
        expect_length(num_infectees, 1)
        expect_true(0 <= num_infectees)
        expect_true(num_infectees <= init_sus)
    })
})

###############################################################################
# Method: Group$infect()
###############################################################################

test_that("Group$infect() successfully completes with a single infective", {
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        init_sus <- 100
        grp <- Group$new(1, ref_strain, init_sus = init_sus)

        expected_host_ids <- c(73)
        grp$infect()

        expect_length(grp$exposed_hosts(), length(expected_host_ids))
        ids <- vapply(grp$exposed_hosts(), function(host) host$id, integer(1L))
        expect_setequal(
            ids,
            expected_host_ids
        )
        expect_equal(grp$time, 1)
    })
})

###############################################################################
# Method: Group$infectious_hosts()
###############################################################################

test_that("Group$infectious_hosts() errors when invalid time given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    grp <- Group$new(1, ref_strain)


    expect_error(grp$infectious_hosts("a"))
    expect_error(grp$infectious_hosts(1.5))
    expect_error(grp$infectious_hosts(c(1, 2)))
})

test_that("Group$infectious_hosts() successfully completes", {
    
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        grp <- Group$new(1, ref_strain, inc_shape = 0)
        # trigger an infection so there are exposed hosts
        grp$infect()

        expected <- c()
        for (host in grp$hosts) {
            if (host$is_infectious(grp$time)) {
                expected <- c(expected, host)
            }
        }

        expect_setequal(grp$infectious_hosts(), expected)
    })
})

###############################################################################
# Method: Group$recovered_hosts()
###############################################################################

test_that("Group$recovered_hosts() errors when invalid time given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    grp <- Group$new(1, ref_strain)


    expect_error(grp$recovered_hosts("a"))
    expect_error(grp$recovered_hosts(1.5))
    expect_error(grp$recovered_hosts(c(1, 2)))
})

test_that("Group$recovered_hosts() successfully completes", {
    
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        init_inf <- 5
        grp <- Group$new(1, ref_strain, init_inf = init_inf)

        index_hosts <- grp$infectious_hosts()
        recovery_times <- vapply(
            index_hosts,
            function(host) host$recovery_time,
            numeric(1L)
        )

        expect_setequal(grp$recovered_hosts(max(recovery_times)), index_hosts)
    })
})

###############################################################################
# Method: Group$run_simulation()
###############################################################################

test_that("Group$run_simluation() errors when invalid lab given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    grp <- Group$new(1, ref_strain)
    Class <- R6::R6Class("Class")
    class <- Class$new()

    expect_error(grp$run_simluation(class))
})

test_that("Group$run_simluation() errors when invalid time given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    grp <- Group$new(1, ref_strain)
    lab <- Lab$new()

    expect_error(grp$run_simluation(lab, feedback = "a"))
    expect_error(grp$run_simluation(lab, feedback = c(1, 2)))
    expect_error(grp$run_simluation(lab, feedback = 1.5))
    expect_error(grp$run_simluation(lab, feedback = -1.5))
})

test_that("Group$run_simluation() errors when invalid time given", {
    
    withr::with_seed(1234, {
        ref_strain <- ReferenceStrain$new("ref_strain")
        grp <- Group$new(1, ref_strain)
        lab <- Lab$new()

        grp$run_simulation(lab, 0)
        
        expect_gt(grp$time, 0)
        expect_gt(lab$num_wgs, 0)
        expect_gt(grp$recovered_size, 0)
    })
})

###############################################################################
# Method: Group$susceptible_hosts()
###############################################################################

test_that("Group$susceptible_hosts() errors when invalid time given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    grp <- Group$new(1, ref_strain)


    expect_error(grp$susceptible_hosts("a"))
    expect_error(grp$susceptible_hosts(1.5))
    expect_error(grp$susceptible_hosts(c(1, 2)))
})

test_that("Group$susceptible_hosts() successfully completes", {
    
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        init_sus <- 5
        grp <- Group$new(1, ref_strain, init_sus = init_sus)

        expect_length(grp$susceptible_hosts(), init_sus)
        expect_true(all(vapply(
            grp$susceptible_hosts(),
            function(host) host$is_susceptible(grp$time),
            logical(1L)
        )))
    })
})

###############################################################################
# Active Binding: Group$hosts_due_for_sampling
###############################################################################

test_that("Group$hosts_due_for_sampling successfully completes", {
    
    withr::with_seed(1234, {
        ref_strain <- ReferenceStrain$new("ref_strain")
        init_inf <- 5
        sample_freq <- 5
        grp <- Group$new(
            1, ref_strain, 
            init_inf = init_inf,
            inc_shape = 0, # turn off exposure compartment
            rec_shape = 20, # force a very long recovery time 
            sample_schedule = "calendar", sample_freq = sample_freq
        )

        grp$time <- sample_freq

        expect_equal(grp$infectious_hosts(), grp$hosts_due_for_sampling)
    })
})

###############################################################################
# Active Binding: Group$is_outbreak_active
###############################################################################

test_that("Group$is_outbreak_active successfully completes", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    init_inf <- 5
    grp <- Group$new(
        1, ref_strain, 
        init_inf = init_inf
    )

    time <- 200
    infector <- grp$infectious_hosts()[[1]]
    infectee <- grp$susceptible_hosts()[[1]]
    
    # perform infection
    infectee$prepare_for_infection(time)
    infector$infect(infectee, time)

    grp$time <- time - 1
    expect_false(grp$is_outbreak_active)
    grp$time <- time
    expect_true(grp$is_outbreak_active)
    grp$time <- infectee$recovery_time - 1
    expect_true(grp$is_outbreak_active)
    grp$time <- infectee$recovery_time
    expect_false(grp$is_outbreak_active)
})

###############################################################################
# Active Binding: Group$group_interval_stack
###############################################################################

test_that("Group$group_interval_stack successfully completes and then errors when empty", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    init_sus <- 9
    init_inf <- 1
    grp <- Group$new(
        1, ref_strain,
        init_sus = init_sus, init_inf = init_inf
    )

    for (i in seq(1:init_inf + init_sus)) {
        intervals <- grp$group_interval_stack
    }
    expect_length(intervals, 2)
    expect_setequal(names(intervals), c("inc", "infector_interval"))
    expect_error(grp$group_interval_stack)    
})

###############################################################################
# Active Binding: Group$time_since_infectious_hosts_infected
###############################################################################

test_that("Group$time_since_infectious_hosts_infected successfully completes", {

    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        init_inf <- 5
        grp <- Group$new(
            1, ref_strain,
            init_inf = init_inf
        )

        # all should be zero since all infectious hosts are index cases
        expect_equal(grp$time_since_infectious_hosts_infected, rep(0, init_inf))

        # create new infections
        grp$infect()
        # get latest time that a new case becomes infectious
        inf_times <- vapply(
            grp$exposed_hosts(),
            function(host) host$infectious_time,
            numeric(1L)
        )
        grp$time <- max(inf_times)
        expected <- vapply(
            grp$infectious_hosts(),
            function(host) max(inf_times) - host$exposure_time,
            numeric(1L)
        )

        expect_equal(grp$time_since_infectious_hosts_infected, expected)

    })
})

###############################################################################
# Active Binding: Group$time_since_infectious_hosts_infectious
###############################################################################

test_that("Group$time_since_infectious_hosts_infectious successfully completes", {

    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        init_inf <- 5
        grp <- Group$new(
            1, ref_strain,
            init_inf = init_inf
        )

        # all should be zero since all infectious hosts are index cases
        expect_equal(grp$time_since_infectious_hosts_infectious, rep(0, init_inf))

        # create new infections
        grp$infect()
        # get latest time that a new case becomes infectious
        inf_times <- vapply(
            grp$exposed_hosts(),
            function(host) host$infectious_time,
            numeric(1L)
        )
        grp$time <- max(inf_times)
        expected <- vapply(
            grp$infectious_hosts(),
            function(host) max(inf_times) - host$infectious_time,
            numeric(1L)
        )

        expect_equal(grp$time_since_infectious_hosts_infectious, expected)

    })
})

###############################################################################
# Private Member: Group$build_interval_stack()
###############################################################################

test_that("Group$build_interval_stack() successfully completes with serial find_infector_method", {

    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        TestP <- R6::R6Class("TestP", inherit = Group,
            active = list(
                interval_stack = function() return(private$interval_stack_)
            )
        )
        grp <- TestP$new(
            1, ref_strain, 
            find_infector_method = "serial",
            si_shape = 8, si_rate = 2    
        )
        stack <- grp$interval_stack

        expect_equal(nrow(stack), grp$size)
        expect_equal(colnames(stack), c("si", "inc", "infector_interval"))
        expect_true(all(stack$infector_interval == stack$si - stack$inc))
    })
})

test_that("Group$build_interval_stack() successfully completes with transmission find_infector_method", {

    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        TestP <- R6::R6Class("TestP", inherit = Group,
            active = list(
                interval_stack = function() return(private$interval_stack_)
            )
        )
        grp <- TestP$new(
            1, ref_strain, 
            find_infector_method = "transmission",
            trans_int_shape = 3, trans_int_rate = 1    
        )
        stack <- grp$interval_stack

        expect_equal(nrow(stack), grp$size)
        expect_equal(colnames(stack), c("inc", "infector_interval"))
    })
})

test_that("Group$build_interval_stack() successfully completes with random find_infector_method", {

    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        TestP <- R6::R6Class("TestP", inherit = Group,
            active = list(
                interval_stack = function() return(private$interval_stack_)
            )
        )
        grp <- TestP$new(1, ref_strain, find_infector_method = "random")
        stack <- grp$interval_stack

        expect_equal(nrow(stack), grp$size)
        expect_equal(colnames(stack), c("inc", "infector_interval"))
        expect_true(all(stack$infector_interval == 0))
    })
})

test_that("Group$build_interval_stack() successfully completes with generation find_infector_method", {

    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        TestP <- R6::R6Class("TestP", inherit = Group,
            active = list(
                interval_stack = function() return(private$interval_stack_)
            )
        )
        grp <- TestP$new(
            1, ref_strain, 
            find_infector_method = "generation",
            gt_shape = 8, gt_rate = 3    
        )
        stack <- grp$interval_stack

        expect_equal(nrow(stack), grp$size)
        expect_equal(colnames(stack), c("inc", "infector_interval"))
    })
})

###############################################################################
# Private Member: Group$initalise_index_host()
###############################################################################

test_that("Group$initalise_index_host() successfully completes with variation TRUE", {


    ref_strain <- ReferenceStrain$new("ref_strain")
    TestP <- R6::R6Class("TestP", inherit = Group,
        public = list(
            initialise_index_host_ = function(...) private$initialise_index_host(...)
        )
    )
    distance <- 5
    grp <- TestP$new(
        1, ref_strain,
        min_init_dist = distance, max_init_dist = distance
    )
    host <- grp$susceptible_hosts()[[1]]
    grp$initialise_index_host_(host)

    expect_false(host$is_susceptible(grp$time))
    expect_false(all(ref_strain$genome() == host$strains[[1]]$genome()))
})

test_that("Group$initalise_index_host() successfully completes with variation FALSE", {


    ref_strain <- ReferenceStrain$new("ref_strain")
    TestP <- R6::R6Class("TestP", inherit = Group,
        public = list(
            initialise_index_host_ = function(...) private$initialise_index_host(...)
        )
    )
    distance <- 5
    grp <- TestP$new(
        1, ref_strain,
        min_init_dist = distance, max_init_dist = distance
    )
    host <- grp$susceptible_hosts()[[1]]
    grp$initialise_index_host_(host, FALSE)

    expect_false(host$is_susceptible(grp$time))
    expect_true(all(ref_strain$genome() == host$strains[[1]]$genome()))
})

###############################################################################
# Private Member: Group$initialise_hosts()
###############################################################################

test_that("Group$initialise_hosts() successfully completes with one initial infective", {


    ref_strain <- ReferenceStrain$new("ref_strain")
    init_inf <- 1
    init_sus <- 5
    grp <- Group$new(
        1, ref_strain,
        init_sus = init_sus, init_inf = init_inf
    )
    
    expect_length(grp$hosts, init_inf + init_sus)
    expect_length(grp$infectious_hosts(), init_inf)
    expect_true(all(
        ref_strain$genome() == grp$infectious_hosts()[[1]]$strains[[1]]$genome()
    ))

})

test_that("Group$initialise_hosts() successfully completes with multiple initial infective", {


    ref_strain <- ReferenceStrain$new("ref_strain")
    init_inf <- 5
    init_sus <- 5
    grp <- Group$new(
        1, ref_strain,
        init_sus = init_sus, init_inf = init_inf
    )
    
    expect_length(grp$hosts, init_inf + init_sus)
    expect_length(grp$infectious_hosts(), init_inf)
})

###############################################################################
# Private Member: Group$force_of_infection()
###############################################################################

test_that("Group$force_of_infection() successfully completes", {


    ref_strain <- ReferenceStrain$new("ref_strain")
    TestP <- R6::R6Class("TestP", inherit = Group,
        public = list(
            force_of_infection_ = function(...) private$force_of_infection(...)
        )
    )
    inf_rate <- 0.1
    init_inf <- 5
    init_sus <- 95
   
    grp <- TestP$new(
        1, ref_strain,
        init_sus = init_sus, init_inf = init_inf,
        inf_rate = inf_rate
    )

    expected <- inf_rate * init_inf / (init_inf + init_sus)
    expect_equal(grp$force_of_infection_(), expected)
})

###############################################################################
# Private Member: Group$initial_genomic_distance()
###############################################################################

test_that("Group$initial_genomic_distance() successfully completes", {


    ref_strain <- ReferenceStrain$new("ref_strain")
    TestP <- R6::R6Class("TestP", inherit = Group,
        public = list(
            initial_genomic_distance_ = function(...) private$initial_genomic_distance(...)
        )
    )
    min_init_dist <- 0
    max_init_dist <- 20
   
    grp <- TestP$new(
        1, ref_strain,
        max_init_dist = max_init_dist,
        min_init_dist = min_init_dist
    )

    distances <- vapply(
        seq(50),
        function(i) grp$initial_genomic_distance_(),
        numeric(1L)
    )

    expect_true(all(distances >= min_init_dist))
    expect_true(all(distances <= max_init_dist))
})

