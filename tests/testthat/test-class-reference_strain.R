test_that("ReferenceStrain$new() errors when invalid name given", {
    expect_error(ReferenceStrain$new(name = 3))
    expect_error(ReferenceStrain$new(name = c("ref", "strain")))
})

test_that("ReferenceStrain$new() errors when invalid genome given", {
    expect_error(ReferenceStrain$new(name = "ref_strain", genome = 3))
    expect_error(ReferenceStrain$new(name = "ref_strain", genome = character(0L)))
})

test_that("ReferenceStrain$new() errors when invalid g_len given", {
    expect_error(ReferenceStrain$new(name = "ref_strain", g_len = "a"))
    expect_error(ReferenceStrain$new(name = "ref_strain", g_len = 12.3))
    expect_error(ReferenceStrain$new(name = "ref_strain", g_len = 0))
    expect_error(ReferenceStrain$new(name = "ref_strain", g_len = c(3, 4)))
})

test_that("ReferenceStrain$new() errors when invalid mut_rate given", {
    expect_error(ReferenceStrain$new(name = "ref_strain", mut_rate = "a"))
    expect_error(ReferenceStrain$new(name = "ref_strain", mut_rate = c(0.2, 0.3)))
    expect_error(ReferenceStrain$new(name = "ref_strain", mut_rate = -1))
    expect_error(ReferenceStrain$new(name = "ref_strain", mut_rate = 1.5))
})

test_that("ReferenceStrain$new() errors when invalid dna given", {
    expect_error(ReferenceStrain$new(name = "ref_strain", dna = "a"))
    expect_error(ReferenceStrain$new(name = "ref_strain", dna = c(TRUE, FALSE)))
})

test_that("ReferenceStrain$new() errors when invalid dna genome given", {
    
    genome <- c("A", "U")
    expect_error(ReferenceStrain$new(name = "ref_strain", genome = genome))    
})

test_that("ReferenceStrain$new() errors when invalid rna genome given", {
    
    genome <- c("A", "T")
    expect_error(ReferenceStrain$new(
        name = "ref_strain", 
        genome = genome, 
        dna = FALSE
    ))    
})

test_that("ReferenceStrain$new() sets up successfully when no genome given", {
    
    g_len <- 25
    name <- "ref_strain"
    mut_rate <- 0.1
    ref_strain <- ReferenceStrain$new(
        name = name, 
        g_len = g_len,
        mut_rate = mut_rate,
    )
    
    expect_equal(ref_strain$g_len, g_len)
    expect_false(is.null(ref_strain$genome()))
    expect_equal(ref_strain$name, name)
    expect_equal(ref_strain$mut_rate, mut_rate)
    expect_true(ref_strain$is_dna)
})

test_that("ReferenceStrain$new() sets up successfully when valid genome given", {
    
    genome <- c("A", "T")
    ref_strain <- ReferenceStrain$new(name = "ref_strain", genome = genome)
    
    expect_equal(ref_strain$g_len, length(genome))
    expect_equal(ref_strain$genome(FALSE), genome_to_integers(genome))
})

test_that("ReferenceStrain$genome() errors when invalid as_character given", {
    
    ref_strain <- ReferenceStrain$new(name = "ref_strain")
    
    expect_error(ref_strain$genome("a"))
    expect_error(ref_strain$genome(c(TRUE, FALSE)))
})

test_that("ReferenceStrain$genome() returns correct characters", {
    
    g <- c("A", "T")
    ref_strain <- ReferenceStrain$new(name = "ref_strain", genome = g)
    genome <- ref_strain$genome()
    
    expect_true(is.character(genome))
    expect_equal(genome, g)
})

test_that("ReferenceStrain$genome() returns correct integers", {
    
    g <- c("A", "T")
    ref_strain <- ReferenceStrain$new(name = "ref_strain", genome = g)
    genome <- ref_strain$genome(FALSE)
    
    expect_true(is.integer(genome))
    expect_equal(genome, genome_to_integers(g))
})

test_that("ReferenceStrain$nucleotides_at_loci() errors when invalid loci given", {
    
    g_len <- 20
    ref_strain <- ReferenceStrain$new(name = "ref_strain", g_len = g_len)
    
    expect_error(ref_strain$nucleotides_at_loci(c(0)))
    expect_error(ref_strain$nucleotides_at_loci(c(g_len + 1)))
    expect_error(ref_strain$nucleotides_at_loci("a"))
    expect_error(ref_strain$nucleotides_at_loci(integer(0L)))
})

test_that("ReferenceStrain$nucleotides_at_loci() errors when invalid as_character given", {
    
    ref_strain <- ReferenceStrain$new(name = "ref_strain")
    
    expect_error(ref_strain$nucleotides_at_loci(c(1), "a"))
    expect_error(ref_strain$nucleotides_at_loci(c(1), c(TRUE, FALSE)))
})

test_that("ReferenceStrain$nucleotides_at_loci() returns correct integers", {
    
    genome <- c("A", "A", "T", "G")
    g_ints <- genome_to_integers(genome)
    ref_strain <- ReferenceStrain$new(name = "ref_strain", genome)
    
    loci <- sample(length(genome), 2)
    nucs <- ref_strain$nucleotides_at_loci(loci)
    expect_true(is.integer(nucs))
    expect_equal(nucs, g_ints[loci])
})

test_that("ReferenceStrain$nucleotides_at_loci() returns correct characters", {
    
    genome <- c("A", "A", "T", "G")
    ref_strain <- ReferenceStrain$new(name = "ref_strain", genome)
    
    loci <- sample(length(genome), 2)
    nucs <- ref_strain$nucleotides_at_loci(loci, TRUE)
    expect_true(is.character(nucs))
    expect_equal(nucs, genome[loci])
})

test_that("ReferenceStrain$randomise_genome() *private* returns valid genome", {
    
    TestRS <- R6::R6Class(
        "TestRS",
        inherit = ReferenceStrain,
        public = list(
            randomise_genome_ = function(...) private$randomise_genome(...)
        )
    )

    ref_strain <- TestRS$new(name = "ref_strain")
    g_len <- 20

    rand_genome <- ref_strain$randomise_genome_(g_len)

    expect_true(is.integer(rand_genome))
    expect_length(rand_genome, g_len)
    expect_true(all(rand_genome >= 1))
    expect_true(all(rand_genome <= 4))
})