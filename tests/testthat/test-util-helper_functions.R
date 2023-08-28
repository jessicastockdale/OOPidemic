###############################################################################
# function: is_ancestor_r6class()
###############################################################################

test_that("is_ancestor_r6class() errors when inputs are not R6Classes", {
    
    expect_error(is_ancestor_r6class(Reference_strain, 2))
    expect_error(is_ancestor_r6class("a", Reference_strain))
})

test_that("is_ancestor_r6class() returns FALSE if class_1 is null", {
    
    expect_error(is_ancestor_r6class(NULL, ReferenceStrain))
})

test_that("is_ancestor_r6class() returns TRUE if class_1 is the same as class_2", {
    
    expect_true(is_ancestor_r6class(ReferenceStrain, ReferenceStrain))
})

test_that("is_ancestor_r6class() returns TRUE if class_1 is an ancester of  class_2", {
    
    # define test R6Class
    DescendentClass <- R6::R6Class("DescendentClass", inherit = ReferenceStrain)
    expect_true(is_ancestor_r6class(ReferenceStrain, DescendentClass))
})

test_that("is_ancestor_r6class() returns False if class_1 is not an ancester of  class_2", {
    
    # define test R6Class
    DescendentClass <- R6::R6Class("DescendentClass", inherit = ReferenceStrain)
    expect_false(is_ancestor_r6class(DescendentClass, ReferenceStrain))
})
