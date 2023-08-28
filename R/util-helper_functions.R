#' Function to recursively test if class_1 is the same class as or a ancestor class of class_2
#'
#' @param class_1,class_2 R6ClassGenerator objects
#'
#' @return A boolean 
#' @export
#'
is_ancestor_r6class <- function(class_1, class_2) {

    stopifnot(
        any(R6::is.R6Class(class_1), is.null(class_1)),
        R6::is.R6Class(class_2)
    )

    if (class_1$classname == class_2$classname) {
        return(TRUE)

    } else if (is.null(class_2$inherit)) {
        # if we've progressed through all ancestors of class_1 without a match class_1$inherit should be null
        return(FALSE)

    } else {
        return(is_ancestor_r6class(class_1, get(class_2$inherit)))

    }
}
