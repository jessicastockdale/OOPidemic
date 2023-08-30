###############################################################################
# Method: Lab$new()
###############################################################################

test_that("Lab$new() errors when invalid wgs_class given", {
    
    Class <- R6::R6Class("Class")
    class <- Class$new()

    expect_error(Lab$new(class))
})

test_that("Lab$new() comletes successfully", {
    
    lab <- Lab$new()

    expect_true("Lab" %in% class(lab))    
})

###############################################################################
# Method: Lab$fasta_df()
###############################################################################

test_that("Lab$fasta_df() completes successfully", {
    
    withr::with_seed(1234, {
        ref_strain <- ReferenceStrain$new("ref_strain", g_len = 100)
        sample_freq <- 1
        init_inf <- 6
        grp <- Group$new(
            1, ref_strain,
            init_inf = init_inf,
            inc_shape = 0,
            sample_schedule = "calendar", sample_freq = sample_freq
        )
        host <- grp$infectious_hosts()[[1]]
        lab <- Lab$new()

        lab$sample_hosts(grp$infectious_hosts(sample_freq), sample_freq)
        fasta_df <- lab$fasta_df()

        expect_setequal(names(fasta_df), c("name", "genome"))
        expect_equal(nrow(fasta_df), init_inf)

    })
})

###############################################################################
# Method: Lab$hostdata_df()
###############################################################################

test_that("Lab$hostdata_df() completes successfully", {
    
    withr::with_seed(1234, {
        ref_strain <- ReferenceStrain$new("ref_strain", g_len = 100)
        sample_freq <- 1
        init_inf <- 6
        grp <- Group$new(
            1, ref_strain,
            init_inf = init_inf,
            inc_shape = 0,
            sample_schedule = "calendar", sample_freq = sample_freq
        )
        host <- grp$infectious_hosts()[[1]]
        lab <- Lab$new()

        lab$sample_hosts(grp$infectious_hosts(sample_freq), sample_freq)
        hostdata_df <- lab$hostdata_df()

        host_ids <- vapply(
            grp$infectious_hosts(sample_freq),
            function(host) host$id,
            numeric(1L)
        )

        expect_equal(nrow(hostdata_df), init_inf)
        expect_setequal(hostdata_df$host_id, host_ids)
    })
})

###############################################################################
# Method: Lab$sample_hosts()
###############################################################################

test_that("Lab$sample_hosts() errors when invalid hosts given", {
    
    lab <- Lab$new()
    Class <- R6::R6Class("Class")
    class <- Class$new()
    time <- 2

    expect_error(lab$sample_hosts(class, time))
})

test_that("Lab$sample_hosts() errors when invalid time given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain", g_len = 100)
    sample_freq <- 1
    init_inf <- 6
    grp <- Group$new(
        1, ref_strain,
        init_inf = 5
    )
    lab <- Lab$new()

    expect_error(lab$sample_hosts(grp$infectious_hosts(), "a"))
    expect_error(lab$sample_hosts(grp$infectious_hosts(), 1.2))
    expect_error(lab$sample_hosts(grp$infectious_hosts(), c(1, 2)))
})

test_that("Lab$sample_hosts() completes successfully when no hosts are passed", {
    
    lab <- Lab$new()
    num_wgs <- lab$num_wgs
    lab$sample_hosts(c(), 1)

    expect_equal(lab$num_wgs, num_wgs)
})

test_that("Lab$sample_hosts() completes successfully when hosts are passed", {
    
    withr::with_seed(1234, {

        ref_strain <- ReferenceStrain$new("ref_strain", g_len = 100)
        sample_freq <- 1
        init_inf <- 6
        grp <- Group$new(
            1, ref_strain,
            init_inf = init_inf,
            inc_shape = 0,
            sample_schedule = "calendar", sample_freq = sample_freq
        )
        lab <- Lab$new()
        lab$sample_hosts(grp$infectious_hosts(sample_freq), sample_freq)

        expect_equal(lab$num_wgs, init_inf)
    })
})

###############################################################################
# Method: Lab$metadata_df()
###############################################################################

test_that("Lab$metadata_df() completes successfully", {
    
    withr::with_seed(1234, {
        ref_strain <- ReferenceStrain$new("ref_strain", g_len = 100)
        sample_freq <- 1
        init_inf <- 6
        grp <- Group$new(
            1, ref_strain,
            init_inf = init_inf,
            inc_shape = 0,
            sample_schedule = "calendar", sample_freq = sample_freq
        )
        host <- grp$infectious_hosts()[[1]]
        lab <- Lab$new()

        lab$sample_hosts(grp$infectious_hosts(sample_freq), sample_freq)
        metadata_df <- lab$metadata_df()

        wgs_names <- vapply(
            lab$wg_sequences,
            function(wgs) wgs$name,
            character(1L)
        )

        expect_equal(nrow(metadata_df), init_inf)
        expect_setequal(metadata_df$name, wgs_names)
    })
})

###############################################################################
# Private Method: Lab$pick_strains()
###############################################################################

test_that("Lab$pick_strains() runs successfullly", {

    TestL <- R6::R6Class("TestL", inherit = Lab,
        public = list(pick_strains_ = function(...) private$pick_strains(...))
    )    

    ref_strain <- ReferenceStrain$new("ref_strain", g_len = 100)
    sample_freq <- 1
    grp <- Group$new(
        1, ref_strain,
        inc_shape = 0,
        sample_schedule = "calendar", sample_freq = sample_freq
    )
    host <- grp$infectious_hosts()[[1]]
    lab <- TestL$new()

    expect_setequal(lab$pick_strains_(host), host$strains)
})

###############################################################################
# Private Method: Lab$sample_host()
###############################################################################

test_that("Lab$sample_host() errors if host is not due for sampling", {

    TestL <- R6::R6Class("TestL", inherit = Lab,
        public = list(sample_host_ = function(...) private$sample_host(...))
    )    

    ref_strain <- ReferenceStrain$new("ref_strain", g_len = 100)
    grp <- Group$new(1, ref_strain)
    lab <- TestL$new()

    expect_error(lab$sample_host_(grp$susceptible_hosts()[[1]], 3))
})

test_that("Lab$sample_host() errors if incorrect host object passed", {

    TestL <- R6::R6Class("TestL", inherit = Lab,
        public = list(sample_host_ = function(...) private$sample_host(...))
    )    

    Class <- R6::R6Class("Class")
    class <- Class$new()
    lab <- TestL$new()

    expect_error(lab$sample_host_(class, 3))
})

test_that("Lab$sample_host() runs successfully", {

    withr::with_seed(1234, {
        TestL <- R6::R6Class("TestL", inherit = Lab,
            public = list(sample_host_ = function(...) private$sample_host(...))
        )    

        ref_strain <- ReferenceStrain$new("ref_strain", g_len = 100)
        sample_freq <- 1
        grp <- Group$new(
            1, ref_strain,
            inc_shape = 0,
            sample_schedule = "calendar", sample_freq = sample_freq
        )
        host <- grp$infectious_hosts()[[1]]
        strain <- host$strains[[1]]
        lab <- TestL$new()

        lab$sample_host_(host, sample_freq)

        expect_equal(host$realisation_time, sample_freq)
        expect_equal(host$sample_time, sample_freq * 2)
        expect_equal(strain$realisation_time, sample_freq)
        expect_length(lab$wg_sequences, 1)
        expect_identical(lab$wg_sequences[[1]]$host, host)
        expect_identical(lab$wg_sequences[[1]]$strain, strain)

    })
})

