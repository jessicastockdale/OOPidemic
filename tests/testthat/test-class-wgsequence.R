###############################################################################
# Method: WGSequence$new()
###############################################################################

test_that("WGSequence$new() errors when invalid sample_num given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)
    host <- pop$infectious_hosts()[[1]]
    strain <- host$strains[[1]]

    expect_error(WGSequence$new(sample_num = "a", strain, host))
    expect_error(WGSequence$new(sample_num = 1.5, strain, host))
    expect_error(WGSequence$new(sample_num = 0, strain, host))
    expect_error(WGSequence$new(sample_num = c(1, 2, strain, host)))
})

test_that("WGSequence$new() errors when invalid strain given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)
    host <- pop$infectious_hosts()[[1]]

    Class <- R6::R6Class("Class")
    class <- Class$new()

    expect_error(WGSequence$new(1, class, host))
})

test_that("WGSequence$new() errors when invalid strain given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(1, ref_strain)
    host <- pop$infectious_hosts()[[1]]
    strain <- host$strains[[1]]

    Class <- R6::R6Class("Class")
    class <- Class$new()

    expect_error(WGSequence$new(1, strain, class))
})

test_that("WGSequence$new() completes successfully", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(
        1, ref_strain,
        init_inf = 2,
        min_init_dist = 1, max_init_dist = 20
    )
    host <- pop$infectious_hosts()[[1]]
    strain <- host$strains[[1]]
    wgs <- WGSequence$new(1, strain, host)

    expect_setequal(strain$loci, wgs$loci)
    expect_setequal(strain$nucleotides, wgs$nucleotides)
    expect_identical(strain$ref_strain, ref_strain)
    expect_equal(wgs$time, host$realisation_time)
})

###############################################################################
# Method: WGSequence$genome()
###############################################################################

test_that("WGSequence$genome() errors when invalid as_character given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(
        1, ref_strain,
        init_inf = 2,
        min_init_dist = 1, max_init_dist = 20
    )
    host <- pop$infectious_hosts()[[1]]
    strain <- host$strains[[1]]
    wgs <- WGSequence$new(1, strain, host)

    expect_error(wgs$genome("a"))    
    expect_error(wgs$genome(c(TRUE, FALSE)))    
})

test_that("WGSequence$genome() returns correct characters", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(
        1, ref_strain,
        init_inf = 2,
        min_init_dist = 1, max_init_dist = 20
    )
    host <- pop$infectious_hosts()[[1]]
    strain <- host$strains[[1]]
    wgs <- WGSequence$new(1, strain, host)

    genome <- wgs$genome()

    expect_true(is.character(genome))
    expect_length(genome, ref_strain$g_len)
    expect_equal(genome[wgs$loci], c("A", "C", "G", "T")[wgs$nucleotides])
    expect_equal(sum(ref_strain$genome() != genome), length(wgs$loci))
    expect_equal(
        which(ref_strain$genome() != genome), 
        wgs$loci[order(wgs$loci)]
    )
})

test_that("WGSequence$genome() returns correct integers", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(
        1, ref_strain,
        init_inf = 2,
        min_init_dist = 1, max_init_dist = 20
    )
    host <- pop$infectious_hosts()[[1]]
    strain <- host$strains[[1]]
    wgs <- WGSequence$new(1, strain, host)

    genome <- wgs$genome(FALSE)

    expect_true(is.integer(genome))
    # other aspects already tested
})

###############################################################################
# Method: WGSequence$name()
###############################################################################

test_that("WGSequence$name() errors when invalid as_character given", {
    
    ref_strain <- ReferenceStrain$new("ref_strain")
    pop <- Population$new(
        1, ref_strain,
        init_inf = 2,
        min_init_dist = 1, max_init_dist = 20
    )
    host <- pop$infectious_hosts()[[1]]
    strain <- host$strains[[1]]
    sample_num <- 1
    wgs <- WGSequence$new(sample_num, strain, host)

    expect_equal(
        wgs$name, 
        paste0(
            sprintf("%02d", pop$id), "_",
            sprintf("%04d", host$id), "_",
            sprintf("%04d", host$realisation_time), "_",
            sprintf("%02d", sample_num)
        )
    )
})

