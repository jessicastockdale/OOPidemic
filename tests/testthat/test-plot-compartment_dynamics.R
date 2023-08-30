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
    
    population <- Population$new(1, init_inf = 0)
    expect_error(plot_compartment_dynamics(population, frequency = "A"))
    expect_error(plot_compartment_dynamics(population, frequency = 1.2))
    expect_error(plot_compartment_dynamics(population, frequency = c(1, 2)))
})

test_that("plot_compartment_dynamics() errors when invalid combine given", {
    
    population <- Population$new(1, init_inf = 0)
    expect_error(plot_compartment_dynamics(population, combine = "A"))
    expect_error(plot_compartment_dynamics(population, combine = c(TRUE, FALSE)))
})

test_that("plot_compartment_dynamics() errors when invalid population_labels given", {

    ref_strain <- ReferenceStrain$new("ref_strain")    
    outbreak <- Outbreak$new(
        c(
            Population$new(1, ref_strain, init_inf = 5, init_sus = 95),
            Population$new(2, ref_strain, init_inf = 0, init_sus = 100)
        ),
        matrix(c(0.75, 0.5, 0.25, 0.1), ncol = 2),
        Lab$new()
    )
    expect_error(plot_compartment_dynamics(outbreak, population_labels = 3))
    expect_error(plot_compartment_dynamics(outbreak, population_labels = c("A")))
})

