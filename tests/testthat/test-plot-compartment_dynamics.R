###############################################################################
# function: plot_compartment_dynamics()
###############################################################################

test_that("plot_compartment_dynamics() errors when invalid x given", {
    
    Class <- R6::R6Class("Class")
    class <- Class$new()

    expect_error(plot_compartment_dynamics(3))
    expect_error(plot_compartment_dynamics(class))
})

test_that("plot_compartment_dynamics() errors when invalid frequency given", {
    
    group <- Group$new(1, init_inf = 0)
    expect_error(plot_compartment_dynamics(group, frequency = "A"))
    expect_error(plot_compartment_dynamics(group, frequency = 1.2))
    expect_error(plot_compartment_dynamics(group, frequency = c(1, 2)))
})

test_that("plot_compartment_dynamics() errors when invalid combine given", {
    
    group <- Group$new(1, init_inf = 0)
    expect_error(plot_compartment_dynamics(group, combine = "A"))
    expect_error(plot_compartment_dynamics(group, combine = c(TRUE, FALSE)))
})

test_that("plot_compartment_dynamics() errors when invalid group_labels given", {

    ref_strain <- ReferenceStrain$new("ref_strain")    
    population <- Population$new(
        c(
            Group$new(1, ref_strain, init_inf = 5, init_sus = 95),
            Group$new(2, ref_strain, init_inf = 0, init_sus = 100)
        ),
        matrix(c(0.75, 0.5, 0.25, 0.1), ncol = 2),
        Lab$new()
    )
    expect_error(plot_compartment_dynamics(population, group_labels = 3))
    expect_error(plot_compartment_dynamics(population, group_labels = c("A")))
})

###############################################################################
# function: compartment_sizes_for_group()
###############################################################################

test_that("compartment_sizes_for_group() errors when invalid group given", {
    
    expect_error(compartment_sizes_for_group(
        class <- R6::R6Class("Class")$new(), 2, 1
    ))
})

test_that("compartment_sizes_for_group() errors when invalid max_time given", {
    
    group <- Group$new(1, init_inf = 0)
    expect_error(compartment_sizes_for_group(group, max_time = "A"))
    expect_error(compartment_sizes_for_group(group, max_time = 1.2))
    expect_error(compartment_sizes_for_group(group, max_time = c(1, 4)))
})

test_that("compartment_sizes_for_group() errors when invalid max_time given", {
    
    group <- Group$new(1, init_inf = 0)
    expect_error(compartment_sizes_for_group(group, 100, frequency = "A"))
    expect_error(compartment_sizes_for_group(group, 100, frequency = 1.2))
    expect_error(compartment_sizes_for_group(group, 100, frequency = c(1, 4)))
})

test_that("compartment_sizes_for_group() completes successfully", {
    
    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        init_inf <- 5
        group <- Group$new(
            1, ref_strain,
            init_inf = init_inf
        )
        group$run_simulation(Lab$new(), 0)

        time <- floor(group$time / 2)
        row <- compartment_sizes(group, time)
        df <- compartment_sizes_for_group(group, group$time)

        expect_equal(nrow(df), group$time + 1)
        expect_equal(as.list(df[time + 1, ]), row)
    })
})

###############################################################################
# function: compartment_sizes_for_groups()
###############################################################################

test_that("compartment_sizes_for_groups() errors when invalid groups given", {
    
    expect_error(compartment_sizes_for_groups(c(
        class <- R6::R6Class("Class")$new(), 2, 1,
        Group$new(1, init_inf = 0)
    )))
})

test_that("compartment_sizes_for_groups() errors when invalid max_time given", {
    
    group <- Group$new(1, init_inf = 0)
    expect_error(compartment_sizes_for_groups(group, max_time = "A"))
    expect_error(compartment_sizes_for_groups(group, max_time = 1.2))
    expect_error(compartment_sizes_for_groups(group, max_time = c(1, 4)))
})

test_that("compartment_sizes_for_groups() errors when invalid max_time given", {
    
    group <- Group$new(1, init_inf = 0)
    expect_error(compartment_sizes_for_groups(group, 100, frequency = "A"))
    expect_error(compartment_sizes_for_groups(group, 100, frequency = 1.2))
    expect_error(compartment_sizes_for_groups(group, 100, frequency = c(1, 4)))
})

test_that("compartment_sizes_for_groups() completes successfully", {
    
    withr::with_seed(124, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        grps <- c(
            Group$new(1, ref_strain),
            Group$new(2, ref_strain)
        )
        inf_rates <- matrix(rep(0.025, 4), ncol = 2)
        lab <- Lab$new()
        population <- Population$new(grps, inf_rates, lab)
        population$run_simulation(0)

        df <- compartment_sizes_for_groups(grps, population$time)
        expect_equal(nrow(df), length(grps) * (population$time + 1))
        expect_equal(
            df[seq(1, population$time + 1), ],
            compartment_sizes_for_group(grps[[1]], population$time)
        )
    })

})

###############################################################################
# function: compartment_sizes()
###############################################################################

test_that("compartment_sizes() completes successfully", {
    
    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        init_inf <- 5
        init_sus <- 95
        group <- Group$new(
            1, ref_strain,
            init_inf = init_inf
        )
        group$run_simulation(Lab$new(), 0)

        time <- floor(group$time / 2)
        c_sizes <- compartment_sizes(group, time)

        expect_equal(c_sizes$t, time)    
        expect_equal(c_sizes$id, as.character(group$id))    
        expect_equal(c_sizes$S, length(group$susceptible_hosts(time)))    
        expect_equal(c_sizes$E, length(group$exposed_hosts(time)))    
        expect_equal(c_sizes$I, length(group$infectious_hosts(time)))    
        expect_equal(c_sizes$R, length(group$recovered_hosts(time)))    
    })
})

