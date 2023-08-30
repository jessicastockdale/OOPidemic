###############################################################################
# function: incubation_time()
###############################################################################

test_that("incubation_time() runs correctly with a susceptible host", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    group <- Group$new(1, init_inf = 0, init_sus = 1)
    host <- group$hosts[[1]]
    expect_identical(incubation_time(host), NA)
})

test_that("incubation_time() runs correctly with an infected case", {
    
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        group <- Group$new(1, ref_strain, init_inf = 1, init_sus = 10)
        index_host <- group$infectious_hosts()[[1]]
        group$run_simulation(Lab$new(), 0)

        # loop through recovered hosts until we get a non index case
        i <- 1
        host <- group$recovered_hosts()[[i]]
        while (host$is_index) {
            i <- i + 1
            host <- group$recovered_hosts()[[i]]
        }

        expect_identical(incubation_time(host), host$infectious_time - host$exposure_time)
    })

})

###############################################################################
# function: infectiousness_duration()
###############################################################################

test_that("infectiousness_duration() runs correctly with a susceptible host", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    group <- Group$new(1, init_inf = 0, init_sus = 1)
    host <- group$hosts[[1]]
    expect_identical(infectiousness_duration(host), NA)
})

test_that("infectiousness_duration() runs correctly with an infected case", {

    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        group <- Group$new(1, ref_strain, init_inf = 1, init_sus = 10)
        index_host <- group$infectious_hosts()[[1]]
        group$run_simulation(Lab$new(), 0)

        # loop through recovered hosts until we get a non index case
        i <- 1
        host <- group$recovered_hosts()[[i]]
        while (host$is_index) {
            i <- i + 1
            host <- group$recovered_hosts()[[i]]
        }

        expect_identical(infectiousness_duration(host), host$recovery_time - host$infectious_time)
    })

})

###############################################################################
# function: generation_time()
###############################################################################

test_that("generation_time() runs correctly with a susceptible host", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    group <- Group$new(1, init_inf = 0, init_sus = 1)
    host <- group$hosts[[1]]

    expect_identical(generation_time(host), NA)
})

test_that("generation_time() runs correctly with an index case", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    group <- Group$new(1, ref_strain, init_inf = 1, init_sus = 0)
    host <- group$hosts[[1]]

    expect_identical(generation_time(host), NA)
})

test_that("generation_time() runs correctly with an infected case", {

    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        group <- Group$new(1, ref_strain, init_inf = 1, init_sus = 10)
        index_host <- group$infectious_hosts()[[1]]
        group$run_simulation(Lab$new(), 0)

        # loop through recovered hosts until we get a non index case
        i <- 1
        host <- group$recovered_hosts()[[i]]
        while (host$is_index) {
            i <- i + 1
            host <- group$recovered_hosts()[[i]]
        }

        expect_identical(generation_time(host), host$exposure_time - host$infector$exposure_time)
    })

})

###############################################################################
# function: serial_interval()
###############################################################################

test_that("serial_interval() runs correctly with a susceptible host", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    group <- Group$new(1, init_inf = 0, init_sus = 1)
    host <- group$hosts[[1]]

    expect_identical(serial_interval(host), NA)
})

test_that("serial_interval() runs correctly with an index case", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    group <- Group$new(1, ref_strain, init_inf = 1, init_sus = 0)
    host <- group$hosts[[1]]

    expect_identical(serial_interval(host), NA)
})

test_that("serial_interval() runs correctly with an infected case", {

    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        group <- Group$new(1, ref_strain, init_inf = 1, init_sus = 10)
        index_host <- group$infectious_hosts()[[1]]
        group$run_simulation(Lab$new(), 0)

        # loop through recovered hosts until we get a non index case
        i <- 1
        host <- group$recovered_hosts()[[i]]
        while (host$is_index) {
            i <- i + 1
            host <- group$recovered_hosts()[[i]]
        }

        expect_identical(serial_interval(host), host$infectious_time - host$infector$infectious_time)
    })

})

###############################################################################
# function: transmission_interval()
###############################################################################

test_that("transmission_interval() runs correctly with a susceptible host", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    group <- Group$new(1, init_inf = 0, init_sus = 1)
    host <- group$hosts[[1]]

    expect_identical(transmission_interval(host), NA)
})

test_that("transmission_interval() runs correctly with an index case", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    group <- Group$new(1, ref_strain, init_inf = 1, init_sus = 0)
    host <- group$hosts[[1]]

    expect_identical(transmission_interval(host), NA)
})

test_that("transmission_interval() runs correctly with an infected case", {

    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        group <- Group$new(1, ref_strain, init_inf = 1, init_sus = 10)
        index_host <- group$infectious_hosts()[[1]]
        group$run_simulation(Lab$new(), 0)

        # loop through recovered hosts until we get a non index case
        i <- 1
        host <- group$recovered_hosts()[[i]]
        while (host$is_index) {
            i <- i + 1
            host <- group$recovered_hosts()[[i]]
        }

        expect_identical(transmission_interval(host), host$exposure_time - host$infector$infectious_time)
    })

})

###############################################################################
# function: plot_intervals()
###############################################################################

test_that("plot_intervals() errors when invalid x given", {
    
    Class <- R6::R6Class("Class")
    class <- Class$new()

    expect_error(plot_intervals(3))
    expect_error(plot_intervals(class))
})

test_that("plot_intervals() errors when invalid interval given", {
    
    group <- Group$new(1, init_inf = 0)
    expect_error(plot_intervals(group, interval = "A"))
})

test_that("plot_intervals() errors when invalid show_distribution given", {
    
    group <- Group$new(1, init_inf = 0)
    expect_error(plot_intervals(group, show_distribution = "A"))
    expect_error(plot_intervals(group, show_distribution = c(TRUE, FALSE)))
})

test_that("plot_intervals() errors when invalid combine given", {
    
    group <- Group$new(1, init_inf = 0)
    expect_error(plot_intervals(group, combine = "A"))
    expect_error(plot_intervals(group, combine = c(TRUE, FALSE)))
})

test_that("plot_intervals() errors when invalid group_labels given", {

    ref_strain <- ReferenceStrain$new("ref_strain")    
    population <- Population$new(
        c(
            Group$new(1, ref_strain, init_inf = 5, init_sus = 95),
            Group$new(2, ref_strain, init_inf = 0, init_sus = 100)
        ),
        matrix(c(0.75, 0.5, 0.25, 0.1), ncol = 2),
        Lab$new()
    )
    expect_error(plot_intervals(population, group_labels = 3))
    expect_error(plot_intervals(population, group_labels = c("A")))
})
