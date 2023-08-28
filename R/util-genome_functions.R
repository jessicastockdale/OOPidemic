base_chars <- c("A", "C", "G")
base_ints <- c(1L, 2L, 3L, 4L)

#' Convert a genome represeneted by characters to a genome represented by integers
#'
#' @param char_genome A character vector composed of "A", "C", "G" & "T" or "U"
#' @param dna A bool indicating if the genome is DNA or RNA based
#'
#' @return A integer vector
#'
genome_to_integers <- function(char_genome, dna = TRUE) {

    # catch invalid input
    stopifnot(
        is.character(char_genome),
        is.logical(dna),
        length(dna) == 1,
        any(
            all(dna, !"U" %in% char_genome),
            all(!dna, !"T" %in% char_genome)
        )
    )

    if (dna) {base_chars <- c(base_chars, "T")}
    else {base_chars <- c(base_chars, "U")}

    names(base_ints) <- base_chars
    genome <- unlist(char_genome)
    unname(base_ints[genome])
}

#' Convert a genome represeneted by integers to a genome represented by characters
#'
#' @param int_genome A integer vector compose of 1, 2, 3, & 4
#' @param dna A bool indicating if the genome is DNA or RNA based
#'
#' @return A character vector
#'
genome_to_characters <- function(int_genome, dna = TRUE) {

    stopifnot(
        is.integer(int_genome),
        all(0 < int_genome),
        all(int_genome < 5),
        is.logical(dna),
        length(dna) == 1
    )

    if (dna) {base_chars <- c(base_chars, "T")}
    else {base_chars <- c(base_chars, "U")}

    base_chars[int_genome]
}
