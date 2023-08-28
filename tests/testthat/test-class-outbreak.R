###############################################################################
# Method: Outbreak$new()
###############################################################################

test_that("Outbreak$new() errors when an invalid population is given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    Class <- R6::R6Class("Class")
    pops <- c(
        Population$new(2, ref_strain),
        class <- Class$new()
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()

    expect_error(Outbreak$new(pops, inf_rates, lab))
})

test_that("Outbreak$new() errors when an invalid population ids are used", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pops_1 <- c(
        Population$new(1, ref_strain),
        Population$new(1, ref_strain)
    )
    pops_2 <- c(
        Population$new(1, ref_strain),
        Population$new(0, ref_strain)
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()

    expect_error(Outbreak$new(pops_1, inf_rates, lab))
    expect_error(Outbreak$new(pops_2, inf_rates, lab))
})

test_that("Outbreak$new() errors when an invalid inf_rate is given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pops <- c(
        Population$new(1, ref_strain),
        Population$new(2, ref_strain)
    )
    lab <- Lab$new()

    expect_error(Outbreak$new(pops, matrix(rep(1, 4), ncol = 1), lab))
    expect_error(Outbreak$new(pops, matrix(rep(1, 4), ncol = 4), lab))
    expect_error(Outbreak$new(pops, matrix(rep("a", 2), ncol = 4), lab))
    expect_error(Outbreak$new(pops, matrix(rep(-1, 2), ncol = 4), lab))
})

test_that("Outbreak$new() errors when an invalid lab is given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    Class <- R6::R6Class("Class")
    pops <- c(
        Population$new(1, ref_strain),
        Population$new(2, ref_strain)
    )
    class <- Class$new()
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    
    expect_error(Outbreak$new(pops, inf_rates, class))
})

test_that("Outbreak$new() successfully completes", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pops <- c(
        Population$new(1, ref_strain),
        Population$new(2, ref_strain)
    )
    lab <- Lab$new()
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    
    outbreak <- Outbreak$new(pops, inf_rates, lab)

    expect_identical(outbreak$lab, lab)
})

###############################################################################
# Method: Outbreak$get_infectees()
###############################################################################

test_that("Outbreak$get_infectees() errors when an invalid num_infectees is given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pops <- c(
        Population$new(1, ref_strain),
        Population$new(2, ref_strain)
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()
    outbreak <- Outbreak$new(pops, inf_rates, lab)
    
    expect_error(outbreak$get_infectees(c("a", "b")))
    expect_error(outbreak$get_infectees(c(1.2, 3)))
    expect_error(outbreak$get_infectees(c(1)))
})

test_that("Outbreak$get_infectees() successfully completes", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pops <- c(
        Population$new(1, ref_strain),
        Population$new(2, ref_strain)
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()
    outbreak <- Outbreak$new(pops, inf_rates, lab)

    num_infectees <- c(2, 3)
    infectees <- outbreak$get_infectees(num_infectees)

    expect_length(infectees, length(num_infectees))
    expect_true(all(vapply(
        seq_along(num_infectees),
        function(i) length(infectees[[i]]) == num_infectees[[i]],
        logical(1L)
    )))
    expect_true(all(vapply(
        infectees,
        function(infs) all(vapply(
            infs, 
            function(inf) inf$is_susceptible(outbreak$time),
            logical(1L)
        )),
        logical(1L)
    )))
})

###############################################################################
# Method: Outbreak$get_infectors()
###############################################################################

test_that("Outbreak$get_infectors() errors when an invalid num_infectees is given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pops <- c(
        Population$new(1, ref_strain),
        Population$new(2, ref_strain)
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()
    outbreak <- Outbreak$new(pops, inf_rates, lab)
    
    expect_error(outbreak$get_infectors(c("a", "b"), c(1, 2)))
    expect_error(outbreak$get_infectors(c(1.2, 3), c(1, 2)))
    expect_error(outbreak$get_infectors(c(3), c(1)))

})

test_that("Outbreak$get_infectors() errors when an invalid infector_intervals is given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pops <- c(
        Population$new(1, ref_strain),
        Population$new(2, ref_strain)
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()
    outbreak <- Outbreak$new(pops, inf_rates, lab)
    
    expect_error(outbreak$get_infectors(c(1), c("A")))
    expect_error(outbreak$get_infectors(c(1), c(-1)))

})

test_that("Outbreak$get_infectors() errors when an incompatible population_of_infectors and infector_intervals is given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pops <- c(
        Population$new(1, ref_strain),
        Population$new(2, ref_strain)
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()
    outbreak <- Outbreak$new(pops, inf_rates, lab)
    
    expect_error(outbreak$get_infectors(c(1, 2), c(-1)))

})

test_that("Outbreak$get_infectors() completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pops <- c(
        Population$new(1, ref_strain, find_infector_method = "random"),
        Population$new(2, ref_strain, find_infector_method = "random")
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()
    outbreak <- Outbreak$new(pops, inf_rates, lab)
    
    pop_of_infs <- c(1, 2, 1)
    inf_ints <- c(0, 0, 0)
    infectors <- outbreak$get_infectors(pop_of_infs, inf_ints)

    expect_true(all(vapply(
        seq_along(pop_of_infs),
        function(i) {all(
            infectors[[i]]$population$id == pop_of_infs[[i]],
            infectors[[i]]$is_infectious(outbreak$time)
        )},
        logical(1L)
    )))
})

###############################################################################
# Method: Outbreak$get_population_of_infectors()
###############################################################################

test_that("Outbreak$get_population_of_infectors() errors when an invalid num_infectees is given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pops <- c(
        Population$new(1, ref_strain),
        Population$new(2, ref_strain)
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()
    outbreak <- Outbreak$new(pops, inf_rates, lab)
 
    fois_mat <- matrix(rep(0.5, 4), ncol = 2)
    
    expect_error(outbreak$get_population_of_infectors(c("a", "b"), fois_mat))
    expect_error(outbreak$get_population_of_infectors(c(1.2, 3), fois_mat))
    expect_error(outbreak$get_population_of_infectors(c(3), fois_mat))

})

test_that("Outbreak$get_population_of_infectors() errors when an invalid fois_mat is given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pops <- c(
        Population$new(1, ref_strain),
        Population$new(2, ref_strain)
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()
    outbreak <- Outbreak$new(pops, inf_rates, lab)

    
    expect_error(
        outbreak$get_population_of_infectors(
            c(1, 2), matrix(rep(0.5, 4), ncol = 4)
    ))
    expect_error(
        outbreak$get_population_of_infectors(
            c(1, 2), matrix(rep(0.5, 4), ncol = 1)
    ))
    expect_error(
        outbreak$get_population_of_infectors(
            c(1, 2), matrix(rep("a", 4), ncol = 2)
    ))
    expect_error(
        outbreak$get_population_of_infectors(
            c(1, 2), matrix(rep(-1, 4), ncol = 2)
    ))
    expect_error(
        outbreak$get_population_of_infectors(
            c(1, 2), matrix(rep(2, 4), ncol = 2)
    ))

})

test_that("Outbreak$get_population_of_infectors() errors when an invalid num_infectees is given", {

    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        pops <- c(
            Population$new(1, ref_strain),
            Population$new(2, ref_strain)
        )
        inf_rates <- matrix(rep(1, 4), ncol = 2)
        lab <- Lab$new()
        TestO <- R6::R6Class("TestO", inherit = Outbreak,
            public = list(pop_ids = function() private$pop_ids_)
        )
        outbreak <- TestO$new(pops, inf_rates, lab)
    
        num_infectees <- c(2, 1)
        fois_mat <- matrix(rep(0.25, 4), ncol = 2)
        pop_infs <- outbreak$get_population_of_infectors(num_infectees, fois_mat)
        

        expect_length(pop_infs, length(num_infectees))
        expect_true(all(vapply(
            seq_along(num_infectees),
            function(i) {all(
                length(pop_infs[[i]]) == num_infectees[[i]],
                all(pop_infs[[i]] %in% outbreak$pop_ids())
            )},
            logical(1L)
        )))
    })
})

###############################################################################
# Method: Outbreak$get_num_infectees()
###############################################################################

test_that("Outbreak$get_num_infectees() errors when an invalid num_infectees is given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pops <- c(
        Population$new(1, ref_strain),
        Population$new(2, ref_strain)
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()
    outbreak <- Outbreak$new(pops, inf_rates, lab)
 
    fois_mat <- matrix(rep(0.5, 4), ncol = 2)
    
    expect_error(outbreak$get_num_infectees(c(1)))
    expect_error(outbreak$get_num_infectees(c("a", 2)))
    expect_error(outbreak$get_num_infectees(c(-1, 2)))

})

test_that("Outbreak$get_num_infectees() completes successfully", {

    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        pops <- c(
            Population$new(1, ref_strain),
            Population$new(2, ref_strain)
        )
        inf_rates <- matrix(rep(1, 4), ncol = 2)
        lab <- Lab$new()
        outbreak <- Outbreak$new(pops, inf_rates, lab)
    
        fois_mat <- matrix(rep(0.025, 4), ncol = 2)
        num_infectees <- outbreak$get_num_infectees(rowSums(fois_mat))
        

        expect_length(num_infectees, length(pops))
        expect_true(all(num_infectees %% 1 == 0))
        expect_true(all(num_infectees >= 0))

    })
})

###############################################################################
# Method: Outbreak$infect()
###############################################################################


test_that("Outbreak$infect() completes successfully", {

    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        init_sus <- 49
        pops <- c(
            Population$new(1, ref_strain, init_sus = init_sus),
            Population$new(2, ref_strain, init_sus = init_sus)
        )
        inf_rates <- matrix(rep(1, 4), ncol = 2)
        lab <- Lab$new()
        TestO <- R6::R6Class("TestO", inherit = Outbreak,
            public = list(align_population_times_ = function() private$align_population_times())
        )
        outbreak <- TestO$new(pops, inf_rates, lab)
    
        fois_mat <- matrix(rep(0.025, 4), ncol = 2)
        outbreak$time <- 1
        outbreak$align_population_times_()
        outbreak$infect()

        expect_true(all(outbreak$susceptible_sizes < init_sus))
    })
})

###############################################################################
# Method: Outbreak$run_simulation()
###############################################################################


test_that("Outbreak$run_simulation() errors when invalid feedback provided", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pops <- c(
        Population$new(1, ref_strain),
        Population$new(2, ref_strain)
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()
    outbreak <- Outbreak$new(pops, inf_rates, lab)

    expect_error(outbreak$run_simulation("a"))
    expect_error(outbreak$run_simulation(c(2, 1)))
    expect_error(outbreak$run_simulation(1.2))

})

test_that("Outbreak$run_simulation() completes successfully", {

    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        pops <- c(
            Population$new(1, ref_strain),
            Population$new(2, ref_strain)
        )
        inf_rates <- matrix(rep(0.025, 4), ncol = 2)
        lab <- Lab$new()
        outbreak <- Outbreak$new(pops, inf_rates, lab)
        outbreak$run_simulation(0)
    
        expect_gt(outbreak$time, 0)
        expect_gt(outbreak$num_wgs, 0)
        expect_gt(sum(outbreak$recovered_sizes), 0)
    })
})

###############################################################################
# Method: Outbreak$sample_populations()
###############################################################################

test_that("Outbreak$sample_populations() completes successfully", {

    withr::with_seed(1234, {
        ref_strain <- ReferenceStrain$new("ref_strain")
        sample_schedule <- "calendar"
        sample_freq <- 2
        init_inf <- 4
        pops <- c(
            Population$new(
                1, ref_strain,
                init_inf = init_inf,
                inc_shape = 0,
                rec_shape = 11,
                sample_schedule = sample_schedule,
                sample_freq = sample_freq    
            ),
            Population$new(
                2, ref_strain,
                init_inf = init_inf,
                inc_shape = 0,
                rec_shape = 11,
                sample_schedule = sample_schedule,
                sample_freq = sample_freq    
            )
        )
        inf_rates <- matrix(rep(0.025, 4), ncol = 2)
        lab <- Lab$new()
        TestO <- R6::R6Class("TestO", inherit = Outbreak,
            public = list(align_population_times_ = function() private$align_population_times())
        )
        outbreak <- TestO$new(pops, inf_rates, lab)
    
        outbreak$time <- sample_freq
        outbreak$align_population_times_()

        outbreak$sample_populations()

        expect_equal(lab$num_wgs, init_inf * 2)
    })
})

###############################################################################
# Active Bindings: Outbreak$exposed_sizes
###############################################################################


test_that("Outbreak$exposed_sizes completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pops <- c(
        Population$new(1, ref_strain),
        Population$new(2, ref_strain)
    )
    inf_rates <- matrix(rep(0.025, 4), ncol = 2)
    lab <- Lab$new()
    TestO <- R6::R6Class("TestO", inherit = Outbreak,
        public = list(align_population_times_ = function() private$align_population_times())
    )
    outbreak <- TestO$new(pops, inf_rates, lab)
    outbreak$time <- 1
    outbreak$align_population_times_()
    outbreak$infect()

    expected <- vapply(pops, function(pop) pop$exposed_size, numeric(1L))
    expect_equal(outbreak$exposed_sizes, expected)

})

###############################################################################
# Method: Outbreak$hosts_due_for_sampling
###############################################################################

test_that("Outbreak$hosts_due_for_sampling() completes successfully", {

    withr::with_seed(1234, {
        ref_strain <- ReferenceStrain$new("ref_strain")
        sample_schedule <- "calendar"
        sample_freq <- 1
        init_inf <- 4
        pops <- c(
            Population$new(
                1, ref_strain,
                init_inf = init_inf,
                inc_shape = 0,
                rec_shape = 11,
                sample_schedule = sample_schedule,
                sample_freq = sample_freq    
            ),
            Population$new(
                2, ref_strain,
                init_inf = init_inf,
                inc_shape = 0,
                rec_shape = 11,
                sample_schedule = sample_schedule,
                sample_freq = sample_freq    
            )
        )
        inf_rates <- matrix(rep(0.025, 4), ncol = 2)
        lab <- Lab$new()
        TestO <- R6::R6Class("TestO", inherit = Outbreak,
            public = list(align_population_times_ = function() private$align_population_times())
        )
        outbreak <- TestO$new(pops, inf_rates, lab)
        outbreak$time <- 1
        outbreak$align_population_times_()
    
        hosts <- outbreak$hosts_due_for_sampling
        expect_length(hosts, init_inf * 2)
        expect_true(all(vapply(
            hosts,
            function(host) {all(
                host$is_infectious(outbreak$time),
                host$is_sampling_due(outbreak$time)
            )},
            logical(1L)
        )))
    })
})

###############################################################################
# Active Bindings: Outbreak$is_outbreak_active
###############################################################################


test_that("Outbreak$is_outbreak_active completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pops <- c(
        Population$new(1, ref_strain),
        Population$new(2, ref_strain)
    )
    inf_rates <- matrix(rep(0.025, 4), ncol = 2)
    lab <- Lab$new()
    TestO <- R6::R6Class("TestO", inherit = Outbreak,
        public = list(align_population_times_ = function() private$align_population_times())
    )
    outbreak <- TestO$new(pops, inf_rates, lab)

    expect_true(outbreak$is_outbreak_active)

    outbreak$time <- 200
    outbreak$align_population_times_()
    outbreak$infect()

    expect_false(outbreak$is_outbreak_active)
    
})

###############################################################################
# Active Bindings: Outbreak$infectious_sizes
###############################################################################


test_that("Outbreak$infectious_sizes completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    init_inf <- 4
    pops <- c(
        Population$new(1, ref_strain, init_inf = init_inf),
        Population$new(2, ref_strain, init_inf = init_inf)
    )
    inf_rates <- matrix(rep(0.025, 4), ncol = 2)
    lab <- Lab$new()
    outbreak <- Outbreak$new(pops, inf_rates, lab)

    expected <- vapply(pops, function(pop) pop$infectious_size, numeric(1L))
    expect_equal(outbreak$infectious_sizes, expected)
})

###############################################################################
# Active Bindings: Outbreak$recovered_sizes
###############################################################################


test_that("Outbreak$recovered_sizes completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    init_inf_1 <- 4
    init_inf_2 <- 6
    pops <- c(
        Population$new(1, ref_strain, init_inf = init_inf_1),
        Population$new(2, ref_strain, init_inf = init_inf_2)
    )
    inf_rates <- matrix(rep(0.025, 4), ncol = 2)
    lab <- Lab$new()
    TestO <- R6::R6Class("TestO", inherit = Outbreak,
        public = list(align_population_times_ = function() private$align_population_times())
    )
    outbreak <- TestO$new(pops, inf_rates, lab)
    outbreak$time <- 200
    outbreak$align_population_times_()

    expect_equal(outbreak$recovered_sizes, c(init_inf_1, init_inf_2))
})

###############################################################################
# Active Bindings: Outbreak$susceptible_sizes
###############################################################################


test_that("Outbreak$susceptible_sizes completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    init_sus_1 <- 4
    init_sus_2 <- 6
    pops <- c(
        Population$new(1, ref_strain, init_sus = init_sus_1),
        Population$new(2, ref_strain, init_sus = init_sus_2)
    )
    inf_rates <- matrix(rep(0.025, 4), ncol = 2)
    lab <- Lab$new()
    outbreak <- Outbreak$new(pops, inf_rates, lab)

    expect_equal(outbreak$susceptible_sizes, c(init_sus_1, init_sus_2))
})

###############################################################################
# Active Bindings: Outbreak$outbreak_sizes
###############################################################################


test_that("Outbreak$outbreak_sizes completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    init_inf <- 4
    pops <- c(
        Population$new(1, ref_strain, init_inf = init_inf),
        Population$new(2, ref_strain, init_inf = init_inf)
    )
    inf_rates <- matrix(rep(0.025, 4), ncol = 2)
    lab <- Lab$new()
    outbreak <- Outbreak$new(pops, inf_rates, lab)

    expect_equal(outbreak$outbreak_size, init_inf * 2)
})

###############################################################################
# Active Bindings: Outbreak$population_outbreak_sizes
###############################################################################


test_that("Outbreak$population_outbreak_sizes completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    init_inf <- 4
    pops <- c(
        Population$new(1, ref_strain, init_inf = init_inf),
        Population$new(2, ref_strain, init_inf = init_inf)
    )
    inf_rates <- matrix(rep(0.025, 4), ncol = 2)
    lab <- Lab$new()
    outbreak <- Outbreak$new(pops, inf_rates, lab)

    expect_equal(outbreak$population_outbreak_sizes, rep(init_inf, 2))
})

###############################################################################
# Active Bindings: Outbreak$population_sizes
###############################################################################


test_that("Outbreak$population_sizes completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    init_inf <- 4
    init_sus <- 6
    pops <- c(
        Population$new(
            1, ref_strain, 
            init_inf = init_inf,
            init_sus = init_sus
        ),
        Population$new(
            2, ref_strain, 
            init_inf = init_inf,
            init_sus = init_sus
        )
    )
    inf_rates <- matrix(rep(0.025, 4), ncol = 2)
    lab <- Lab$new()
    outbreak <- Outbreak$new(pops, inf_rates, lab)

    expect_equal(outbreak$population_sizes, rep(init_inf + init_sus, 2))
})

###############################################################################
# Active Bindings: Outbreak$align_population_times()
###############################################################################


test_that("Outbreak$align_population_times() completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    init_inf_1 <- 4
    init_inf_2 <- 6
    pops <- c(
        Population$new(1, ref_strain, init_inf = init_inf_1),
        Population$new(2, ref_strain, init_inf = init_inf_2)
    )
    inf_rates <- matrix(rep(1, 4), nrow = 2)
    lab <- Lab$new()
    TestO <- R6::R6Class("TestO", inherit = Outbreak,
        public = list(
            align_population_times_ = function() private$align_population_times()
        )
    )
    outbreak <- TestO$new(pops, inf_rates, lab)
   
    outbreak$time
    expect_true(all(vapply(pops, function(p) p$time, numeric(1L)) == 0))
    outbreak$align_population_times_()
    expect_true(all(vapply(pops, function(p) p$time, numeric(1L)) == outbreak$time))
})

###############################################################################
# Active Bindings: Outbreak$force_of_infection()
###############################################################################


test_that("Outbreak$force_of_infection() completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    init_inf_1 <- 1
    init_inf_2 <- 2
    init_inf_3 <- 5
    init_inf <- c(init_inf_1, init_inf_2, init_inf_3)
    init_sus_1 <- 10
    init_sus_2 <- 20
    init_sus_3 <- 30
    init_sus <- c(init_sus_1, init_sus_2, init_sus_3)
    pops <- c(
        Population$new(
            1, ref_strain, 
            init_inf = init_inf_1,
            init_sus = init_sus_1
        ),
        Population$new(
            2, ref_strain, 
            init_inf = init_inf_2,
            init_sus = init_sus_2
        ),
        Population$new(
            3, ref_strain, 
            init_inf = init_inf_3,
            init_sus = init_sus_3
        )
    )
    inf_rates <- matrix(rep(1, 9), nrow = 3)
    lab <- Lab$new()
    TestO <- R6::R6Class("TestO", inherit = Outbreak,
        public = list(
            force_of_infection_ = function(...) private$force_of_infection(...)
        )
    )
    outbreak <- TestO$new(pops, inf_rates, lab)

    expected_output <- sweep(inf_rates, MARGIN = 2, init_inf, "*") / sum(init_sus + init_inf)
    output <- outbreak$force_of_infection_(init_inf)
    
    expect_equal(output, expected_output)
   
})

