###############################################################################
# Method: Strain$new()
###############################################################################

test_that("Strain$new() errors when invalid ref_strain object given", {
    
    Class <- R6::R6Class("Class")
    class <- Class$new()

    expect_error(Strain$new(class))
})

test_that("Strain$new() errors when invalid time given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Strain$new(ref_strain, time = 'a'))
    expect_error(Strain$new(ref_strain, time = 1.5))
    expect_error(Strain$new(ref_strain, time = c(2, 3)))
    expect_error(Strain$new(ref_strain, time = -2))
})

test_that("Strain$new() errors when invalid distance given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")

    expect_error(Strain$new(ref_strain, distance = 'a'))
    expect_error(Strain$new(ref_strain, distance = 1.5))
    expect_error(Strain$new(ref_strain, distance = c(2, 3)))
    expect_error(Strain$new(ref_strain, distance = -2))
})

test_that("Strain$new() sets up successfully", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    time <- 10
    strain <- Strain$new(ref_strain, time)

    expect_identical(strain$ref_strain, ref_strain)
    expect_identical(strain$spawn_time, time)
    expect_identical(strain$realisation_time, time)
})

###############################################################################
# Method: Strain$genome()
###############################################################################

test_that("Strain$genome() errors when invalid as_character given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    strain <- Strain$new(ref_strain)

    expect_error(strain$genome("a"))    
    expect_error(strain$genome(c(TRUE, FALSE)))    
})

test_that("Strain$genome() returns correct characters", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    distance <- 5
    strain <- Strain$new(ref_strain, distance = distance)

    genome <- strain$genome()

    expect_true(is.character(genome))
    expect_length(genome, ref_strain$g_len)
    expect_equal(genome[strain$loci], c("A", "C", "G", "T")[strain$nucleotides])
    expect_equal(sum(ref_strain$genome() != genome), distance)
    expect_equal(
        which(ref_strain$genome() != genome), 
        strain$loci[order(strain$loci)]
    )
})

test_that("Strain$genome() returns correct integers", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    distance <- 5
    strain <- Strain$new(ref_strain, distance = distance)

    genome <- strain$genome(FALSE)

    expect_true(is.integer(genome))
    # other aspects already tested
})

###############################################################################
# Method: Strain$realise()
###############################################################################

test_that("Strain$realise() errors when invalid time given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    strain <- Strain$new(ref_strain)

    expect_error(strain$realise("a"))    
    expect_error(strain$realise(1.5))    
    expect_error(strain$realise(c(2, 3)))    
})

test_that("Strain$realise() successfully completes with zero time elapsed since last realisation", {
    
    ref_strain <- ReferenceStrain$new("ref_strain", mut_rate = 1)
    strain <- Strain$new(ref_strain)
    strain$realise(0) # there should be no change

    expect_length(strain$loci, 0)
    expect_equal(strain$realisation_time, 0)        
})

test_that("Strain$realise() successfully completes with time elapsed since last realisation", {withr::with_seed(1234, {
    
    ref_strain <- ReferenceStrain$new("ref_strain", mut_rate = 0.001)
    strain <- Strain$new(ref_strain)
    time <- 10
    strain$realise(time)

    expected_loci <- c(64, 499, 726, 962, 223, 891, 626, 244, 866, 536, 57, 799, 71)
    expected_nucleotides <- c(1, 3, 2, 3, 4, 1, 3, 2, 1, 3, 2, 4, 3)

    expect_equal(strain$loci, expected_loci)       
    expect_equal(strain$nucleotides, expected_nucleotides)       
})})

###############################################################################
# Method: Strain$reproduce()
###############################################################################

test_that("Strain$reproduce() errors when invalid time given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    strain <- Strain$new(ref_strain)

    expect_error(strain$reproduce("a"))    
    expect_error(strain$reproduce(1.5))    
    expect_error(strain$reproduce(c(2, 3)))    
})

test_that("Strain$reproduce() successfully completes", {
    #implicitly this tests deep_clone
    
    ref_strain <- ReferenceStrain$new("ref_strain", mut_rate = 0)
    strain <- Strain$new(ref_strain)
    time <- 10
    d_strain <- strain$reproduce(time) 

    expect_equal(strain$realisation_time, time)
    expect_identical(strain$descendents, c(d_strain))

    expect_identical(class(strain), class(d_strain))
    expect_identical(d_strain$ancestor, strain)
    expect_equal(d_strain$spawn_time, time)
})

###############################################################################
# Private Method: Strain$mutate()
###############################################################################

test_that("Strain$mutate() successfull complettion", {
    
    TestS <- R6::R6Class("TestS",
        inherit = Strain,
        public = list(mutate_ = function(...) private$mutate(...))
    )
    ref_strain <- ReferenceStrain$new("ref_strain")
    strain <- TestS$new(ref_strain)
    
    distance <- 5
    strain$mutate_(distance)

    expect_length(strain$loci, distance)
    expect_true(is.integer(strain$loci))
    expect_length(strain$nucleotides, distance)
    expect_true(is.integer(strain$nucleotides))
})

###############################################################################
# Private Method: Strain$mutate_nucleotides()
###############################################################################

test_that("Strain$mutate_nucleotides() successfull complettion", {
    
    TestS <- R6::R6Class("TestS",
        inherit = Strain,
        public = list(mutate_nucleotides_ = function(...) private$mutate_nucleotides(...))
    )
    ref_strain <- ReferenceStrain$new("ref_strain")
    strain <- TestS$new(ref_strain)
    
    old_nucs <- sample(c(1, 2, 3, 4), 20, replace = TRUE)
    new_nucs <- strain$mutate_nucleotides_(old_nucs)

    expect_length(new_nucs, length(old_nucs))
    expect_false(any(new_nucs == old_nucs))
    expect_true(all(new_nucs > 0))
    expect_true(all(new_nucs < 5))
})

###############################################################################
# Private Method: Strain$substitution_number()
###############################################################################

test_that("Strain$substitution_number() successfull complettion", {
    
    TestS <- R6::R6Class("TestS",
        inherit = Strain,
        public = list(substitution_number_ = function(...) private$substitution_number(...))
    )
    ref_strain <- ReferenceStrain$new("ref_strain", mut_rate = 0.001)
    strain <- TestS$new(ref_strain)
    time <- 10

    output <- NA
    expected_output <- NA

    seed <- 1234
    withr::with_seed(seed, {
        mean <- ref_strain$mut_rate * ref_strain$g_len * time
        expected_output <- stats::rpois(1, mean) 
    })
    withr::with_seed(seed, {
        output <- strain$substitution_number_(time)
    })

    expect_equal(output, expected_output)
})