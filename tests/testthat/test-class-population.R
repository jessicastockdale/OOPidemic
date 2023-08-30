###############################################################################
# Method: Population$new()
###############################################################################

test_that("Population$new() errors when an invalid group is given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    Class <- R6::R6Class("Class")
    grps <- c(
        Group$new(2, ref_strain),
        class <- Class$new()
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()

    expect_error(Population$new(grps, inf_rates, lab))
})

test_that("Population$new() errors when an invalid group ids are used", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    groups_1 <- c(
        Group$new(1, ref_strain),
        Group$new(1, ref_strain)
    )
    groups_2 <- c(
        Group$new(1, ref_strain),
        Group$new(0, ref_strain)
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()

    expect_error(Population$new(groups_1, inf_rates, lab))
    expect_error(Population$new(groups_2, inf_rates, lab))
})

test_that("Population$new() errors when an invalid inf_rate is given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    grps <- c(
        Group$new(1, ref_strain),
        Group$new(2, ref_strain)
    )
    lab <- Lab$new()

    expect_error(Population$new(grps, matrix(rep(1, 4), ncol = 1), lab))
    expect_error(Population$new(grps, matrix(rep(1, 4), ncol = 4), lab))
    expect_error(Population$new(grps, matrix(rep("a", 2), ncol = 4), lab))
    expect_error(Population$new(grps, matrix(rep(-1, 2), ncol = 4), lab))
})

test_that("Population$new() errors when an invalid lab is given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    Class <- R6::R6Class("Class")
    grps <- c(
        Group$new(1, ref_strain),
        Group$new(2, ref_strain)
    )
    class <- Class$new()
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    
    expect_error(Population$new(grps, inf_rates, class))
})

test_that("Population$new() successfully completes", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    grps <- c(
        Group$new(1, ref_strain),
        Group$new(2, ref_strain)
    )
    lab <- Lab$new()
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    
    population <- Population$new(grps, inf_rates, lab)

    expect_identical(population$lab, lab)
})

###############################################################################
# Method: Population$get_infectees()
###############################################################################

test_that("Population$get_infectees() errors when an invalid num_infectees is given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    grps <- c(
        Group$new(1, ref_strain),
        Group$new(2, ref_strain)
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()
    population <- Population$new(grps, inf_rates, lab)
    
    expect_error(population$get_infectees(c("a", "b")))
    expect_error(population$get_infectees(c(1.2, 3)))
    expect_error(population$get_infectees(c(1)))
})

test_that("Population$get_infectees() successfully completes", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    grps <- c(
        Group$new(1, ref_strain),
        Group$new(2, ref_strain)
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()
    population <- Population$new(grps, inf_rates, lab)

    num_infectees <- c(2, 3)
    infectees <- population$get_infectees(num_infectees)

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
            function(inf) inf$is_susceptible(population$time),
            logical(1L)
        )),
        logical(1L)
    )))
})

###############################################################################
# Method: Population$get_infectors()
###############################################################################

test_that("Population$get_infectors() errors when an invalid num_infectees is given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    grps <- c(
        Group$new(1, ref_strain),
        Group$new(2, ref_strain)
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()
    population <- Population$new(grps, inf_rates, lab)
    
    expect_error(population$get_infectors(c("a", "b"), c(1, 2)))
    expect_error(population$get_infectors(c(1.2, 3), c(1, 2)))
    expect_error(population$get_infectors(c(3), c(1)))

})

test_that("Population$get_infectors() errors when an invalid infector_intervals is given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    grps <- c(
        Group$new(1, ref_strain),
        Group$new(2, ref_strain)
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()
    population <- Population$new(grps, inf_rates, lab)
    
    expect_error(population$get_infectors(c(1), c("A")))
    expect_error(population$get_infectors(c(1), c(-1)))

})

test_that("Population$get_infectors() errors when an incompatible group_of_infectors and infector_intervals is given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    grps <- c(
        Group$new(1, ref_strain),
        Group$new(2, ref_strain)
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()
    population <- Population$new(grps, inf_rates, lab)
    
    expect_error(population$get_infectors(c(1, 2), c(-1)))

})

test_that("Population$get_infectors() completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    grps <- c(
        Group$new(1, ref_strain, find_infector_method = "random"),
        Group$new(2, ref_strain, find_infector_method = "random")
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()
    population <- Population$new(grps, inf_rates, lab)
    
    group_of_infs <- c(1, 2, 1)
    inf_ints <- c(0, 0, 0)
    infectors <- population$get_infectors(group_of_infs, inf_ints)

    expect_true(all(vapply(
        seq_along(group_of_infs),
        function(i) {all(
            infectors[[i]]$group$id == group_of_infs[[i]],
            infectors[[i]]$is_infectious(population$time)
        )},
        logical(1L)
    )))
})

###############################################################################
# Method: Population$get_group_of_infectors()
###############################################################################

test_that("Population$get_group_of_infectors() errors when an invalid num_infectees is given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    grps <- c(
        Group$new(1, ref_strain),
        Group$new(2, ref_strain)
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()
    population <- Population$new(grps, inf_rates, lab)
 
    fois_mat <- matrix(rep(0.5, 4), ncol = 2)
    
    expect_error(population$get_group_of_infectors(c("a", "b"), fois_mat))
    expect_error(population$get_group_of_infectors(c(1.2, 3), fois_mat))
    expect_error(population$get_group_of_infectors(c(3), fois_mat))

})

test_that("Population$get_group_of_infectors() errors when an invalid fois_mat is given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    grps <- c(
        Group$new(1, ref_strain),
        Group$new(2, ref_strain)
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()
    population <- Population$new(grps, inf_rates, lab)

    
    expect_error(
        population$get_group_of_infectors(
            c(1, 2), matrix(rep(0.5, 4), ncol = 4)
    ))
    expect_error(
        population$get_group_of_infectors(
            c(1, 2), matrix(rep(0.5, 4), ncol = 1)
    ))
    expect_error(
        population$get_group_of_infectors(
            c(1, 2), matrix(rep("a", 4), ncol = 2)
    ))
    expect_error(
        population$get_group_of_infectors(
            c(1, 2), matrix(rep(-1, 4), ncol = 2)
    ))
    expect_error(
        population$get_group_of_infectors(
            c(1, 2), matrix(rep(2, 4), ncol = 2)
    ))

})

test_that("Population$get_group_of_infectors() errors when an invalid num_infectees is given", {

    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        grps <- c(
            Group$new(1, ref_strain),
            Group$new(2, ref_strain)
        )
        inf_rates <- matrix(rep(1, 4), ncol = 2)
        lab <- Lab$new()
        TestO <- R6::R6Class("TestO", inherit = Population,
            public = list(group_ids = function() private$group_ids_)
        )
        population <- TestO$new(grps, inf_rates, lab)
    
        num_infectees <- c(2, 1)
        fois_mat <- matrix(rep(0.25, 4), ncol = 2)
        group_infs <- population$get_group_of_infectors(num_infectees, fois_mat)
        

        expect_length(group_infs, length(num_infectees))
        expect_true(all(vapply(
            seq_along(num_infectees),
            function(i) {all(
                length(group_infs[[i]]) == num_infectees[[i]],
                all(group_infs[[i]] %in% population$group_ids())
            )},
            logical(1L)
        )))
    })
})

###############################################################################
# Method: Population$get_num_infectees()
###############################################################################

test_that("Population$get_num_infectees() errors when an invalid num_infectees is given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    grps <- c(
        Group$new(1, ref_strain),
        Group$new(2, ref_strain)
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()
    population <- Population$new(grps, inf_rates, lab)
 
    fois_mat <- matrix(rep(0.5, 4), ncol = 2)
    
    expect_error(population$get_num_infectees(c(1)))
    expect_error(population$get_num_infectees(c("a", 2)))
    expect_error(population$get_num_infectees(c(-1, 2)))

})

test_that("Population$get_num_infectees() completes successfully", {

    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        grps <- c(
            Group$new(1, ref_strain),
            Group$new(2, ref_strain)
        )
        inf_rates <- matrix(rep(1, 4), ncol = 2)
        lab <- Lab$new()
        population <- Population$new(grps, inf_rates, lab)
    
        fois_mat <- matrix(rep(0.025, 4), ncol = 2)
        num_infectees <- population$get_num_infectees(rowSums(fois_mat))
        

        expect_length(num_infectees, length(grps))
        expect_true(all(num_infectees %% 1 == 0))
        expect_true(all(num_infectees >= 0))

    })
})

###############################################################################
# Method: Population$infect()
###############################################################################


test_that("Population$infect() completes successfully", {

    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        init_sus <- 49
        grps <- c(
            Group$new(1, ref_strain, init_sus = init_sus),
            Group$new(2, ref_strain, init_sus = init_sus)
        )
        inf_rates <- matrix(rep(1, 4), ncol = 2)
        lab <- Lab$new()
        TestO <- R6::R6Class("TestO", inherit = Population,
            public = list(align_group_times_ = function() private$align_group_times())
        )
        population <- TestO$new(grps, inf_rates, lab)
    
        fois_mat <- matrix(rep(0.025, 4), ncol = 2)
        population$time <- 1
        population$align_group_times_()
        population$infect()

        expect_true(all(population$susceptible_sizes < init_sus))
    })
})

###############################################################################
# Method: Population$run_simulation()
###############################################################################


test_that("Population$run_simulation() errors when invalid feedback provided", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    grps <- c(
        Group$new(1, ref_strain),
        Group$new(2, ref_strain)
    )
    inf_rates <- matrix(rep(1, 4), ncol = 2)
    lab <- Lab$new()
    population <- Population$new(grps, inf_rates, lab)

    expect_error(population$run_simulation("a"))
    expect_error(population$run_simulation(c(2, 1)))
    expect_error(population$run_simulation(1.2))

})

test_that("Population$run_simulation() completes successfully", {

    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        grps <- c(
            Group$new(1, ref_strain),
            Group$new(2, ref_strain)
        )
        inf_rates <- matrix(rep(0.025, 4), ncol = 2)
        lab <- Lab$new()
        population <- Population$new(grps, inf_rates, lab)
        population$run_simulation(0)
    
        expect_gt(population$time, 0)
        expect_gt(population$num_wgs, 0)
        expect_gt(sum(population$recovered_sizes), 0)
    })
})

###############################################################################
# Method: Population$sample_groups()
###############################################################################

test_that("Population$sample_groups() completes successfully", {

    withr::with_seed(1234, {
        ref_strain <- ReferenceStrain$new("ref_strain")
        sample_schedule <- "calendar"
        sample_freq <- 2
        init_inf <- 4
        grps <- c(
            Group$new(
                1, ref_strain,
                init_inf = init_inf,
                inc_shape = 0,
                rec_shape = 11,
                sample_schedule = sample_schedule,
                sample_freq = sample_freq    
            ),
            Group$new(
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
        TestO <- R6::R6Class("TestO", inherit = Population,
            public = list(align_group_times_ = function() private$align_group_times())
        )
        population <- TestO$new(grps, inf_rates, lab)
    
        population$time <- sample_freq
        population$align_group_times_()

        population$sample_groups()

        expect_equal(lab$num_wgs, init_inf * 2)
    })
})

###############################################################################
# Active Bindings: Population$exposed_sizes
###############################################################################


test_that("Population$exposed_sizes completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    grps <- c(
        Group$new(1, ref_strain),
        Group$new(2, ref_strain)
    )
    inf_rates <- matrix(rep(0.025, 4), ncol = 2)
    lab <- Lab$new()
    TestO <- R6::R6Class("TestO", inherit = Population,
        public = list(align_group_times_ = function() private$align_group_times())
    )
    population <- TestO$new(grps, inf_rates, lab)
    population$time <- 1
    population$align_group_times_()
    population$infect()

    expected <- vapply(grps, function(grp) grp$exposed_size, numeric(1L))
    expect_equal(population$exposed_sizes, expected)

})

###############################################################################
# Method: Population$hosts_due_for_sampling
###############################################################################

test_that("Population$hosts_due_for_sampling() completes successfully", {

    withr::with_seed(1234, {
        ref_strain <- ReferenceStrain$new("ref_strain")
        sample_schedule <- "calendar"
        sample_freq <- 1
        init_inf <- 4
        grps <- c(
            Group$new(
                1, ref_strain,
                init_inf = init_inf,
                inc_shape = 0,
                rec_shape = 11,
                sample_schedule = sample_schedule,
                sample_freq = sample_freq    
            ),
            Group$new(
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
        TestO <- R6::R6Class("TestO", inherit = Population,
            public = list(align_group_times_ = function() private$align_group_times())
        )
        population <- TestO$new(grps, inf_rates, lab)
        population$time <- 1
        population$align_group_times_()
    
        hosts <- population$hosts_due_for_sampling
        expect_length(hosts, init_inf * 2)
        expect_true(all(vapply(
            hosts,
            function(host) {all(
                host$is_infectious(population$time),
                host$is_sampling_due(population$time)
            )},
            logical(1L)
        )))
    })
})

###############################################################################
# Active Bindings: Population$is_outbreak_active
###############################################################################


test_that("Population$is_outbreak_active completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    grps <- c(
        Group$new(1, ref_strain),
        Group$new(2, ref_strain)
    )
    inf_rates <- matrix(rep(0.025, 4), ncol = 2)
    lab <- Lab$new()
    TestO <- R6::R6Class("TestO", inherit = Population,
        public = list(align_group_times_ = function() private$align_group_times())
    )
    population <- TestO$new(grps, inf_rates, lab)

    expect_true(population$is_outbreak_active)

    population$time <- 200
    population$align_group_times_()
    population$infect()

    expect_false(population$is_outbreak_active)
    
})

###############################################################################
# Active Bindings: Population$infectious_sizes
###############################################################################


test_that("Population$infectious_sizes completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    init_inf <- 4
    grps <- c(
        Group$new(1, ref_strain, init_inf = init_inf),
        Group$new(2, ref_strain, init_inf = init_inf)
    )
    inf_rates <- matrix(rep(0.025, 4), ncol = 2)
    lab <- Lab$new()
    population <- Population$new(grps, inf_rates, lab)

    expected <- vapply(grps, function(grp) grp$infectious_size, numeric(1L))
    expect_equal(population$infectious_sizes, expected)
})

###############################################################################
# Active Bindings: Population$recovered_sizes
###############################################################################


test_that("Population$recovered_sizes completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    init_inf_1 <- 4
    init_inf_2 <- 6
    grps <- c(
        Group$new(1, ref_strain, init_inf = init_inf_1),
        Group$new(2, ref_strain, init_inf = init_inf_2)
    )
    inf_rates <- matrix(rep(0.025, 4), ncol = 2)
    lab <- Lab$new()
    TestO <- R6::R6Class("TestO", inherit = Population,
        public = list(align_group_times_ = function() private$align_group_times())
    )
    population <- TestO$new(grps, inf_rates, lab)
    population$time <- 200
    population$align_group_times_()

    expect_equal(population$recovered_sizes, c(init_inf_1, init_inf_2))
})

###############################################################################
# Active Bindings: Population$susceptible_sizes
###############################################################################


test_that("Population$susceptible_sizes completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    init_sus_1 <- 4
    init_sus_2 <- 6
    grps <- c(
        Group$new(1, ref_strain, init_sus = init_sus_1),
        Group$new(2, ref_strain, init_sus = init_sus_2)
    )
    inf_rates <- matrix(rep(0.025, 4), ncol = 2)
    lab <- Lab$new()
    population <- Population$new(grps, inf_rates, lab)

    expect_equal(population$susceptible_sizes, c(init_sus_1, init_sus_2))
})

###############################################################################
# Active Bindings: Population$outbreak_sizes
###############################################################################


test_that("Population$outbreak_sizes completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    init_inf <- 4
    grps <- c(
        Group$new(1, ref_strain, init_inf = init_inf),
        Group$new(2, ref_strain, init_inf = init_inf)
    )
    inf_rates <- matrix(rep(0.025, 4), ncol = 2)
    lab <- Lab$new()
    population <- Population$new(grps, inf_rates, lab)

    expect_equal(population$outbreak_size, init_inf * 2)
})

###############################################################################
# Active Bindings: Population$group_outbreak_sizes
###############################################################################


test_that("Population$group_outbreak_sizes completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    init_inf <- 4
    grps <- c(
        Group$new(1, ref_strain, init_inf = init_inf),
        Group$new(2, ref_strain, init_inf = init_inf)
    )
    inf_rates <- matrix(rep(0.025, 4), ncol = 2)
    lab <- Lab$new()
    population <- Population$new(grps, inf_rates, lab)

    expect_equal(population$group_outbreak_sizes, rep(init_inf, 2))
})

###############################################################################
# Active Bindings: Population$group_sizes
###############################################################################


test_that("Population$group_sizes completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    init_inf <- 4
    init_sus <- 6
    grps <- c(
        Group$new(
            1, ref_strain, 
            init_inf = init_inf,
            init_sus = init_sus
        ),
        Group$new(
            2, ref_strain, 
            init_inf = init_inf,
            init_sus = init_sus
        )
    )
    inf_rates <- matrix(rep(0.025, 4), ncol = 2)
    lab <- Lab$new()
    population <- Population$new(grps, inf_rates, lab)

    expect_equal(population$group_sizes, rep(init_inf + init_sus, 2))
})

###############################################################################
# Active Bindings: Population$align_group_times()
###############################################################################


test_that("Population$align_group_times() completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    init_inf_1 <- 4
    init_inf_2 <- 6
    grps <- c(
        Group$new(1, ref_strain, init_inf = init_inf_1),
        Group$new(2, ref_strain, init_inf = init_inf_2)
    )
    inf_rates <- matrix(rep(1, 4), nrow = 2)
    lab <- Lab$new()
    TestO <- R6::R6Class("TestO", inherit = Population,
        public = list(
            align_group_times_ = function() private$align_group_times()
        )
    )
    population <- TestO$new(grps, inf_rates, lab)
   
    population$time
    expect_true(all(vapply(grps, function(p) p$time, numeric(1L)) == 0))
    population$align_group_times_()
    expect_true(all(vapply(grps, function(p) p$time, numeric(1L)) == population$time))
})

###############################################################################
# Active Bindings: Population$force_of_infection()
###############################################################################


test_that("Population$force_of_infection() completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    init_inf_1 <- 1
    init_inf_2 <- 2
    init_inf_3 <- 5
    init_inf <- c(init_inf_1, init_inf_2, init_inf_3)
    init_sus_1 <- 10
    init_sus_2 <- 20
    init_sus_3 <- 30
    init_sus <- c(init_sus_1, init_sus_2, init_sus_3)
    grps <- c(
        Group$new(
            1, ref_strain, 
            init_inf = init_inf_1,
            init_sus = init_sus_1
        ),
        Group$new(
            2, ref_strain, 
            init_inf = init_inf_2,
            init_sus = init_sus_2
        ),
        Group$new(
            3, ref_strain, 
            init_inf = init_inf_3,
            init_sus = init_sus_3
        )
    )
    inf_rates <- matrix(rep(1, 9), nrow = 3)
    lab <- Lab$new()
    TestO <- R6::R6Class("TestO", inherit = Population,
        public = list(
            force_of_infection_ = function(...) private$force_of_infection(...)
        )
    )
    population <- TestO$new(grps, inf_rates, lab)

    expected_output <- sweep(inf_rates, MARGIN = 2, init_inf, "*") / sum(init_sus + init_inf)
    output <- population$force_of_infection_(init_inf)
    
    expect_equal(output, expected_output)
   
})

