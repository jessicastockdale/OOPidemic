###############################################################################
# function: node_level()
###############################################################################

test_that("node_level() runs correctly with an index host", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    population <- Population$new(1, ref_strain, init_inf = 1, init_sus = 10)
    index_host <- population$infectious_hosts()[[1]]

    expect_equal(node_level(index_host), 1)
})

test_that("node_level() runs correctly with an non-index case", {
    
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        population <- Population$new(1, ref_strain, init_inf = 1, init_sus = 10)
        population$run_simulation(Lab$new(), 0)

        # loop through recovered hosts until we get a non index case
        i <- 1
        host <- population$recovered_hosts()[[i]]
        while (host$is_index) {
            i <- i + 1
            host <- population$recovered_hosts()[[i]]
        }

        expected <- 2
        infector <- host$infector
        while (!infector$is_index) {
            expected <- expected + 1
            infector <- infector$infector
        }

        expect_equal(node_level(host), expected)
    })

})

###############################################################################
# function: node_shape()
###############################################################################

test_that("node_shape() runs correctly with an index host", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    population <- Population$new(1, ref_strain, init_inf = 1, init_sus = 10)
    index_host <- population$infectious_hosts()[[1]]

    expect_equal(node_shape(index_host), "square")
})

test_that("node_shape() runs correctly with an non-index case", {
    
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        population <- Population$new(1, ref_strain, init_inf = 1, init_sus = 10)
        population$run_simulation(Lab$new(), 0)

        # loop through recovered hosts until we get a non index case
        i <- 1
        host <- population$recovered_hosts()[[i]]
        while (host$is_index) {
            i <- i + 1
            host <- population$recovered_hosts()[[i]]
        }

        expect_equal(node_shape(host), "circle")
    })

})

###############################################################################
# function: node_title()
###############################################################################

###############################################################################
# function: interval()
###############################################################################

test_that("interval() runs correctly with an index host", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    population <- Population$new(1, ref_strain, init_inf = 1, init_sus = 10)
    index_host <- population$infectious_hosts()[[1]]

    expect_equal(interval(index_host), NA)
})

test_that("interval() runs correctly with an non-index case", {
    
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        population <- Population$new(1, ref_strain, init_inf = 1, init_sus = 10)
        population$run_simulation(Lab$new(), 0)

        # loop through recovered hosts until we get a non index case
        i <- 1
        host <- population$recovered_hosts()[[i]]
        while (host$is_index) {
            i <- i + 1
            host <- population$recovered_hosts()[[i]]
        }
        infector <- host$infector

        expect_equal(interval(host, "serial"), host$infectious_time - infector$infectious_time)
        expect_equal(interval(host, "generation"), host$exposure_time - infector$exposure_time)
        expect_equal(interval(host, "transmission"), host$exposure_time - infector$infectious_time)
    })

})

###############################################################################
# function: population_nodes()
###############################################################################

test_that("population_nodes() runs correctly", {
    
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        population <- Population$new(1, ref_strain, init_inf = 1, init_sus = 99)
        population$run_simulation(Lab$new(), 0)
        
        hosts <- population$hosts[vapply(
            population$hosts,
            function(host) host$is_recovered(population$time),
            logical(1L)
        )]

        nodes <- population_nodes(population)
        host_ids <- vapply(hosts, function(host) host$id, integer(1L))

        expect_equal(nodes$group, rep(population$id, length(hosts)))
        expect_equal(nodes$id, paste0(population$id, "-", as.character(host_ids)))
        expect_equal(nodes$label, as.character(host_ids))
    })
})

###############################################################################
# function: population_edges()
###############################################################################

test_that("population_edges() runs correctly", {
    
    withr::with_seed(1234, {
    
        ref_strain <- ReferenceStrain$new("ref_strain")
        population <- Population$new(1, ref_strain, init_inf = 1, init_sus = 99)
        population$run_simulation(Lab$new(), 0)
        
        hosts <- population$hosts[vapply(
            population$hosts,
            function(host) host$is_recovered(population$time),
            logical(1L)
        )]

        edges <- population_edges(population)
        host_ids <- vapply(hosts, function(host) host$id, integer(1L))

        expect_equal(edges$from, vapply(
            hosts,
            function(host) {
                if (!host$is_index) {
                    paste0(host$infector$population$id, "-", host$infector$id)
                } else ""
            },
            character(1L)
        ))
        expect_equal(edges$to, paste0(population$id, "-", as.character(host_ids)))
    })
})

