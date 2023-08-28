###############################################################################
# function: genome_to_integers()
###############################################################################

test_that("genome_to_integers() errors when a non character input is given for char_genome", {
	
	length <- 10
	g_ints <- sample(c(1:4), length, replace = TRUE)

	expect_error(genome_to_integers(g_ints))
})

test_that("genome_to_integers() errors when non logical input is given for dna", {
	
	expect_error(genome_to_integers(c("A", "C", "G", "T"), dna = 3))
})

test_that("genome_to_integers() errors when a object of length greater then 1 is given for dna", {
	
	expect_error(genome_to_integers(c("A", "C", "G", "T"), dna = c(TRUE, TRUE)))
})

test_that("genome_to_integers() errors when a DNA genome contains 'U'", {
	
	expect_error(genome_to_integers(c("A", "C", "G", "U")))
})

test_that("genome_to_integers() errors when a DNA genome contains 'U'", {
	
	expect_error(genome_to_integers(c("A", "C", "G", "T"), dna = FALSE))
})

test_that("genome_to_integers() returns the correct output for DNA", {

	length <- 10
	g_ints <- sample(c(1:4), length, replace = TRUE)
	base_chars <- c("A", "C", "G", "T")
	g_chars <- base_chars[g_ints]

	expect_equal(genome_to_integers(g_chars), g_ints)
})

test_that("genome_to_integers() returns the correct output for RNA", {

	length <- 10
	g_ints <- sample(c(1:4), length, replace = TRUE)
	base_chars <- c("A", "C", "G", "U")
	g_chars <- base_chars[g_ints]

	expect_equal(genome_to_integers(g_chars, FALSE), g_ints)
})

###############################################################################
# function: genome_to_characters()
###############################################################################

test_that("genome_to_characters() errors when a non character input is given for int_genome", {
	
	length <- 10
	g_chars <- sample(c("A", "C", "G", "U"), length, replace = TRUE)

	expect_error(genome_to_characters(g_chars))
})

test_that("genome_to_characters() errors when invalid integer input is given for int_genome", {
	
	expect_error(genome_to_characters(c(0, 1)))
	expect_error(genome_to_characters(c(5, 1)))
})


test_that("genome_to_characters() errors when non logical input is given for dna", {
	
	expect_error(genome_to_characters(c(1, 3, 2, 4), dna = 3))
})

test_that("genome_to_characters() errors when a object of length greater then 1 is given for dna", {
	
	expect_error(genome_to_characters(c(1, 3, 2, 4), dna = c(TRUE, TRUE)))
})

test_that("genome_to_characters() returns the correct output for DNA", {

	length <- 10
	g_ints <- sample(c(1:4), length, replace = TRUE)
	base_chars <- c("A", "C", "G", "T")
	g_chars <- base_chars[g_ints]

	expect_equal(genome_to_characters(g_ints), g_chars)
})

test_that("genome_to_characters() returns the correct output for RNA", {

	length <- 10
	g_ints <- sample(c(1:4), length, replace = TRUE)
	base_chars <- c("A", "C", "G", "U")
	g_chars <- base_chars[g_ints]

	expect_equal(genome_to_characters(g_ints, FALSE), g_chars)
})
