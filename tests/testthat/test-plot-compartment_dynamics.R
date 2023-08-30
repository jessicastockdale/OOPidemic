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

