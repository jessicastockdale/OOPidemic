
###############################################################################
# Method: Host$new()
###############################################################################

test_that("Host$new() errors when invalid id given", {
    
    expect_error(Host$new(id = "a"))
    expect_error(Host$new(id = 1.5))
    expect_error(Host$new(id = 0))
    expect_error(Host$new(id = c(1, 2)))
})

test_that("Host$new() errors when invalid population given", {
    
    Class <- R6::R6Class("Class")
    class <- Class$new()

    expect_error(Host$new(1, class))
})

###############################################################################
# Method: Host$contract_strains()
###############################################################################

test_that("Host$contract_strains() errors when invalid infector given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    strain <- Strain$new(ref_strain, 0)
    pop <- Population$new(1, ref_strain)
    infectee <- pop$susceptible_hosts()[[1]]
    infector <- pop$susceptible_hosts()[[2]]

    Class <- R6::R6Class("Class")
    class <- Class$new()

    expect_error(infectee$contract_strains(class, c(strain), c(1)))
    expect_error(infectee$contract_strains(infector, c(strain), c(1)))
})

test_that("Host$contract_strains() errors when invalid strains given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)
    infectee <- pop$susceptible_hosts()[[1]]
    infector <- pop$infectious_hosts()[[1]]

    Class <- R6::R6Class("Class")
    class <- Class$new()

    expect_error(infectee$contract_strains(infector, c(class), c(1)))
})

test_that("Host$contract_strains() errors when invalid freq given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    strain <- Strain$new(ref_strain, 0)
    pop <- Population$new(1, ref_strain)
    infectee <- pop$susceptible_hosts()[[1]]
    infector <- pop$infectious_hosts()[[1]]

    expect_error(infectee$contract_strains(infector, c(strain), c("a")))
    expect_error(infectee$contract_strains(infector, c(strain), c(1.5)))
})

test_that("Host$contract_strains() errors when strains and freq aren't the same length given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    strain <- Strain$new(ref_strain, 0)
    pop <- Population$new(1, ref_strain)
    infectee <- pop$susceptible_hosts()[[1]]
    infector <- pop$infectious_hosts()[[1]]

    expect_error(infectee$contract_strains(infector, c(strain), c(1, 2)))
})

test_that("Host$contract_strains() errors when infectee is not prepared for infection", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    strain <- Strain$new(ref_strain, 0)
    pop <- Population$new(1, ref_strain)
    infectee <- pop$susceptible_hosts()[[1]]
    infector <- pop$infectious_hosts()[[1]]

    expect_error(infectee$contract_strains(infector, c(strain), c(1)))
})

test_that("Host$contract_strains() errors when infectee is already infected", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    strain <- Strain$new(ref_strain, 0)
    pop <- Population$new(1, ref_strain, init_inf = 2)
    infectee <- pop$infectious_hosts()[[1]]
    infector <- pop$infectious_hosts()[[2]]

    expect_error(infectee$contract_strains(infector, c(strain), c(1)))
})

test_that("Host$contract_strains() runs successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    strain <- Strain$new(ref_strain, 0)
    pop <- Population$new(1, ref_strain)
    infectee <- pop$susceptible_hosts()[[1]]
    infector <- pop$infectious_hosts()[[1]]
    
    infectee$prepare_for_infection(1)
    strains <- c(strain)
    freq <- c(1)
    infectee$contract_strains(infector, strains, freq)


    expect_identical(infectee$infector, infector)
    expect_setequal(infectee$strains, strains)
    expect_setequal(infectee$freq, freq)
})

###############################################################################
# Method: Host$infect()
###############################################################################

test_that("Host$infect() errors when a non infectious host tries to infect someone", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)
    infector <- pop$susceptible_hosts()[[2]]
    infectee <- pop$susceptible_hosts()[[1]]
    time <- 3

    expect_error(infector$infect(infectee, time))
})

test_that("Host$infect() errors when invalid infectee object given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)
    infector <- pop$infectious_hosts()[[1]]

    Class <- R6::R6Class("Class")
    class <- Class$new()
    time <- 3

    expect_error(infector$infect(class, time))
})

test_that("Host$infect() errors when infectee hasn't been prepared for infection", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)
    infector <- pop$infectious_hosts()[[1]]
    infectee <- pop$susceptible_hosts()[[1]]
    time <- 3

    expect_error(infector$infect(infectee, time))
})

test_that("Host$infect() errors when infectee hasn't been prepared for infection", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain, init_inf = 2)
    infector <- pop$infectious_hosts()[[1]]
    infectee <- pop$infectious_hosts()[[2]]
    time <- 3

    expect_error(infector$infect(infectee, time))
})

test_that("Host$infect() errors when invalid time given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain, init_inf = 2)
    infector <- pop$infectious_hosts()[[1]]
    infectee <- pop$susceptible_hosts()[[1]]
    time <- 3
    infectee$prepare_for_infection(time)

    expect_error(infector$infect(infectee, "a"))
    expect_error(infector$infect(infectee, 1.3))
    expect_error(infector$infect(infectee, c(time, 2)))
})

test_that("Host$infect() completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain, init_inf = 2)
    infector <- pop$infectious_hosts()[[1]]
    infectee <- pop$susceptible_hosts()[[1]]
    time <- 3
    infectee$prepare_for_infection(time)

    infector$infect(infectee, time)

    expect_equal(infector$realisation_time, time)
    ancestors <- lapply(
        infectee$strains,
        function(strain) strain$ancestor
    )
    expect_setequal(ancestors, infector$strains)
    expect_true(any(
        vapply(
            infector$infectees,
            identical,
            logical(1L),
            infectee
        )
    ))
})

###############################################################################
# Method: Host$is_exposed()
###############################################################################

test_that("Host$is_exposed() errors invalid time given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)
    host <- pop$susceptible_hosts()[[2]]

    expect_error(host$is_exposed("a"))
    expect_error(host$is_exposed(c(1, 2)))
})

test_that("Host$is_exposed() runs successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain, init_inf = 2)
    infector <- pop$infectious_hosts()[[1]]
    infectee <- pop$susceptible_hosts()[[1]]
    
    time <- 15
    expect_false(infectee$is_exposed(time))

    # prepare for infection
    infectee$prepare_for_infection(time)
    expect_false(infectee$is_exposed(time))

    # infect
    infector$infect(infectee, time)
    expect_true(infectee$is_exposed(time))
    expect_true(infectee$is_exposed(infectee$infectious_time - 1))
    expect_false(infectee$is_exposed(infectee$infectious_time))
})

###############################################################################
# Method: Host$is_infectious()
###############################################################################

test_that("Host$is_infectious() errors invalid time given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)
    host <- pop$susceptible_hosts()[[2]]

    expect_error(host$is_infectious("a"))
    expect_error(host$is_infectious(c(1, 2)))
})

test_that("Host$is_infectious() runs successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain, inc_shape = 0)
    infector <- pop$infectious_hosts()[[1]]
    infectee <- pop$susceptible_hosts()[[1]]
    
    time <- 15
    expect_false(infectee$is_infectious(time))

    # prepare for infection
    infectee$prepare_for_infection(time)
    expect_false(infectee$is_infectious(time))

    # infect
    infector$infect(infectee, time)
    expect_true(infectee$is_infectious(time))
    expect_true(infectee$is_infectious(infectee$recovery_time - 1))
    expect_false(infectee$is_infectious(infectee$recovery_time))
})

###############################################################################
# Method: Host$is_recovered()
###############################################################################

test_that("Host$is_recovered() errors invalid time given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)
    host <- pop$susceptible_hosts()[[2]]

    expect_error(host$is_recovered("a"))
    expect_error(host$is_recovered(c(1, 2)))
})

test_that("Host$is_recovered() runs successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain, inc_shape = 0)
    infector <- pop$infectious_hosts()[[1]]
    infectee <- pop$susceptible_hosts()[[1]]
    
    time <- 15
    expect_false(infectee$is_recovered(time))

    # prepare for infection
    infectee$prepare_for_infection(time)
    expect_false(infectee$is_recovered(time))

    # infect
    infector$infect(infectee, time)
    expect_false(infectee$is_recovered(infectee$recovery_time - 1))
    expect_true(infectee$is_recovered(infectee$recovery_time))
})

###############################################################################
# Method: Host$is_sampling_due()
###############################################################################

test_that("Host$is_sampling_due() errors invalid time given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)
    host <- pop$susceptible_hosts()[[2]]

    expect_error(host$is_sampling_due("a"))
    expect_error(host$is_sampling_due(c(1, 2)))
})

test_that("Host$is_sampling_due() runs successfully", {

    withr::with_seed(12434, {

        ref_strain <- ReferenceStrain$new("ref_strain")
        sample_freq <- 5
        pop <- Population$new(
            1, ref_strain, 
            inc_shape = 0,
            rec_shape = 20, rec_rate = 1,
            sample_schedule = "calendar", sample_freq = sample_freq
        )
        infector <- pop$infectious_hosts()[[1]]
        infectee <- pop$susceptible_hosts()[[1]]
        
        time <- 6
        expected_sample_time <- 10
        expect_false(infectee$is_sampling_due(expected_sample_time))

        # # prepare for infection
        infectee$prepare_for_infection(time)
        expect_false(infectee$is_sampling_due(expected_sample_time))

        # infect
        infector$infect(infectee, time)
        expect_false(infectee$is_sampling_due(expected_sample_time - 1))
        expect_true(infectee$is_sampling_due(expected_sample_time))
        expect_false(infectee$is_sampling_due(expected_sample_time + 1))
    })
})

###############################################################################
# Method: Host$is_susceptible()
###############################################################################

test_that("Host$is_susceptible() errors invalid time given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)
    host <- pop$susceptible_hosts()[[2]]

    expect_error(host$is_susceptible("a"))
    expect_error(host$is_susceptible(c(1, 2)))
})

test_that("Host$is_susceptible() runs successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain, inc_shape = 0)
    infector <- pop$infectious_hosts()[[1]]
    infectee <- pop$susceptible_hosts()[[1]]
    
    time <- 15
    expect_true(infectee$is_susceptible(time))

    # prepare for infection
    infectee$prepare_for_infection(time)
    expect_true(infectee$is_susceptible(time))

    # infect
    infector$infect(infectee, time)
    expect_true(infectee$is_susceptible(time - 1))
    expect_false(infectee$is_susceptible(time))
})

###############################################################################
# Method: Host$prepare_for_infection()
###############################################################################

test_that("Host$prepare_for_infection() errors invalid time given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)
    host <- pop$susceptible_hosts()[[2]]

    expect_error(host$prepare_for_infection("a"))
    expect_error(host$prepare_for_infection(1.2))
    expect_error(host$prepare_for_infection(c(1, 2)))
})

test_that("Host$prepare_for_infection() errors invalid initial given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)
    host <- pop$susceptible_hosts()[[2]]

    expect_error(host$prepare_for_infection("a"))
    expect_error(host$prepare_for_infection(c(TRUE, FALSE)))
})

test_that("Host$prepare_for_infection() errors when we attempt to prepare a non-susceptible host for infection", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)

    expect_error(host$prepare_for_infection(1))
})

test_that("Host$prepare_for_infection() runs successfully with initial = FALSE", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    TestP <- R6::R6Class("TestP", inherit = Population,
        active = list(
            interval_stack = function() return(private$interval_stack_)
        )
    )
    pop <- TestP$new(1, ref_strain)
    host <- pop$susceptible_hosts()[[1]]

    time <- 5
    intervals <- pop$interval_stack[1, ]
    infector_interval <- host$prepare_for_infection(time)

    expect_equal(host$exposure_time, time)
    expect_equal(host$realisation_time, time)
    expect_equal(host$infectious_time, time + ceiling(intervals$inc))
    expect_false(is.na(host$recovery_time))
    expect_true(host$recovery_time >= host$infectious_time)
    expect_false(is.na(host$sample_time))
    expect_true(host$sample_time >= host$infectious_time)
    expect_true(host$sample_time < host$recovery_time)
    expect_equal(infector_interval, intervals$infector_interval)
})

test_that("Host$prepare_for_infection() runs successfully with initial = TRUE", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    TestP <- R6::R6Class("TestP", inherit = Population,
        active = list(
            interval_stack = function() return(private$interval_stack_)
        )
    )
    pop <- TestP$new(1, ref_strain)
    host <- pop$susceptible_hosts()[[1]]

    time <- 5
    intervals <- pop$interval_stack[1, ]
    infector_interval <- host$prepare_for_infection(time, TRUE)

    expect_equal(host$exposure_time, time)
    expect_equal(host$realisation_time, time)
    expect_equal(host$infectious_time, time)
    expect_false(is.na(host$recovery_time))
    expect_true(host$recovery_time >= host$infectious_time)
    expect_false(is.na(host$sample_time))
    expect_true(host$sample_time >= host$infectious_time)
    expect_true(host$sample_time < host$recovery_time)
    expect_equal(infector_interval, NA)
})

###############################################################################
# Method: Host$realise()
###############################################################################

test_that("Host$realise() errors invalid time given", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)
    host <- pop$susceptible_hosts()[[2]]

    expect_error(host$realise("a"))
    expect_error(host$realise(1.2))
    expect_error(host$realise(c(1, 2)))
})

test_that("Host$realise() runs successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)
    host <- pop$susceptible_hosts()[[2]]
    
    time <- 15
    host$realise(time)

    expect_equal(host$realisation_time, time)
})

###############################################################################
# Method: Host$realise()
###############################################################################

test_that("Host$realise() runs successfully with random sample_schedule", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)
    host <- pop$infectious_hosts()[[1]]

    sample_time <- host$sample_time
    host$update_sample_time()

    expect_equal(host$sample_time, sample_time)
})

test_that("Host$realise() runs successfully with calendar sample_schedule", {

    withr::with_seed(1234, {
        ref_strain <- ReferenceStrain$new("ref_strain")
        sample_freq <- 5
        pop <- Population$new(
            1, ref_strain,
            sample_schedule = "calendar", sample_freq = sample_freq,
            rec_shape = 11
        )
        host <- pop$infectious_hosts()[[1]]

        sample_time <- host$sample_time
        host$update_sample_time()
        expect_equal(host$sample_time, sample_time + sample_freq)
        host$update_sample_time()
        expect_equal(host$sample_time, Inf)

    })
})

test_that("Host$realise() runs successfully with individual sample_schedule", {

    withr::with_seed(1234, {
        ref_strain <- ReferenceStrain$new("ref_strain")
        sample_freq <- 5
        pop <- Population$new(
            1, ref_strain,
            sample_schedule = "individual", sample_freq = sample_freq,
            rec_shape = 11
        )
        host <- pop$infectious_hosts()[[1]]

        sample_time <- host$sample_time
        host$update_sample_time()
        expect_equal(host$sample_time, sample_time + sample_freq)
        host$update_sample_time()
        expect_equal(host$sample_time, Inf)

    })
})

###############################################################################
# Active Binding: Host$is_index()
###############################################################################

test_that("Host$is_index() completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)
    infector <- pop$infectious_hosts()[[1]]
    infectee <- pop$susceptible_hosts()[[1]]
    host <- pop$susceptible_hosts()[[2]]

    time <- 3
    infectee$prepare_for_infection(time)
    infector$infect(infectee, time)

    expect_true(infector$is_index)
    expect_false(infectee$is_index)
    expect_false(host$is_index)
})  

###############################################################################
# Private member: Host$calendar_sample_time()
###############################################################################

test_that("Host$calendar_sample_time() completes successfully on an index case", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    TestH <- R6::R6Class("TestH", inherit = Host,
        public = list(
            calendar_sample_time_ = function() return(private$calendar_sample_time())
        )
    )
    sample_freq <- 5
    pop <- Population$new(
        1, ref_strain,
        sample_freq = sample_freq,
        host_class = TestH
    )
    host <- pop$infectious_hosts()[[1]]

    expect_equal(host$calendar_sample_time_(), sample_freq)
})  

test_that("Host$calendar_sample_time() completes successfully on a non-index case", {

    withr::with_seed(1234, {
        ref_strain <- ReferenceStrain$new("ref_strain")
        TestH <- R6::R6Class("TestH", inherit = Host,
            public = list(
                calendar_sample_time_ = function() return(private$calendar_sample_time())
            )
        )
        sample_freq <- 5
        pop <- Population$new(
            1, ref_strain,
            inc_shape = 0,
            sample_freq = sample_freq,
            host_class = TestH
        )
        host <- pop$susceptible_hosts()[[1]]
        time <- sample_freq * 2 - 1
        host$prepare_for_infection(time)

        expect_equal(host$calendar_sample_time_(), sample_freq * 2)
    })  
})  

###############################################################################
# Private member: Host$generate_recovery_time()
###############################################################################

test_that("Host$generate_recovery_time() completes successfully on an index case", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    TestH <- R6::R6Class("TestH", inherit = Host,
        public = list(
            generate_recovery_time_ = function() return(private$generate_recovery_time())
        )
    )
    pop <- Population$new(
        1, ref_strain,
        host_class = TestH
    )
    host <- pop$infectious_hosts()[[1]]

    expect_gte(host$generate_recovery_time_(), host$infectious_time)
})  

###############################################################################
# Private member: Host$generate_sample_time()
###############################################################################

test_that("Host$generate_sample_time() completes successfully with individual sample_schedule", {

    withr::with_seed(1234, {
        ref_strain <- ReferenceStrain$new("ref_strain")
        TestH <- R6::R6Class("TestH", inherit = Host,
            public = list(
                generate_sample_time_ = function() return(private$generate_sample_time()),
                individual_sample_time_ = function() return(private$individual_sample_time())
            )
        )
        sample_freq <- 2
        pop <- Population$new(
            1, ref_strain,
            host_class = TestH,
            rec_shape = 5,
            sample_schedule = "individual", sample_freq = sample_freq
        )

        host <- pop$infectious_hosts()[[1]]
        expect_equal(host$generate_sample_time_(), host$individual_sample_time_())

        new_sample_freq <- 200
        pop$sample_freq <- new_sample_freq
        expect_equal(host$generate_sample_time_(), Inf)
    })  
})  

test_that("Host$generate_sample_time() completes successfully with calendar sample_schedule", {

    withr::with_seed(1234, {
        ref_strain <- ReferenceStrain$new("ref_strain")
        TestH <- R6::R6Class("TestH", inherit = Host,
            public = list(
                generate_sample_time_ = function() return(private$generate_sample_time()),
                calendar_sample_time_ = function() return(private$calendar_sample_time())
            )
        )
        sample_freq <- 2
        pop <- Population$new(
            1, ref_strain,
            host_class = TestH,
            rec_shape = 5,
            sample_schedule = "calendar", sample_freq = sample_freq
        )

        host <- pop$infectious_hosts()[[1]]
        expect_equal(host$generate_sample_time_(), host$calendar_sample_time_())
    })  
})  

test_that("Host$generate_sample_time() completes successfully with random sample_schedule", {

    withr::with_seed(1234, {
        ref_strain <- ReferenceStrain$new("ref_strain")
        TestH <- R6::R6Class("TestH", inherit = Host,
            public = list(
                generate_sample_time_ = function() return(private$generate_sample_time()),
                random_sample_time_ = function() return(private$random_sample_time())
            )
        )
        pop <- Population$new(
            1, ref_strain,
            host_class = TestH,
            rec_shape = 5,
            sample_schedule = "random"
        )

        host <- pop$infectious_hosts()[[1]]
        expect_equal(host$generate_sample_time_(), host$random_sample_time_())
    })  
})  

###############################################################################
# Private member: Host$individual_sample_time()
###############################################################################

test_that("Host$individual_sample_time() completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    TestH <- R6::R6Class("TestH", inherit = Host,
        public = list(
            individual_sample_time_ = function() return(private$individual_sample_time())
        )
    )
    sample_freq <- 5
    pop <- Population$new(
        1, ref_strain,
        sample_freq = sample_freq,
        host_class = TestH
    )
    host <- pop$infectious_hosts()[[1]]

    expect_equal(host$individual_sample_time_(), host$infectious_time + sample_freq)
})  

###############################################################################
# Private member: Host$random_sample_time()
###############################################################################

test_that("Host$random_sample_time() completes successfully", {

    ref_strain <- ReferenceStrain$new("ref_strain")
    TestH <- R6::R6Class("TestH", inherit = Host,
        public = list(
            random_sample_time_ = function() return(private$random_sample_time())
        )
    )
    sample_freq <- 5
    pop <- Population$new(
        1, ref_strain,
        sample_freq = sample_freq,
        host_class = TestH
    )
    host <- pop$infectious_hosts()[[1]]
    sample_time <- host$random_sample_time_()

    expect_true(host$infectious_time <= sample_time)
    expect_true(sample_time < host$recovery_time)
})  

