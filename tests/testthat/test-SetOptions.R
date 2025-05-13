test_that("SetOptions returns all arguments in the list", {
  # Get the formals (arguments) of the SetOptions function
  function_args <- names(formals(SetOptions))

  # Get the names of the elements in the returned list
  returned_list <- SetOptions()
  returned_list_names <- names(returned_list)

  # Check if all function arguments are present in the returned list
  expect_true(
    all(function_args %in% returned_list_names),
    info = "Some arguments in SetOptions are not returned in the list."
  )

  # Check if there are no extra elements in the returned list
  expect_true(
    all(returned_list_names %in% function_args),
    info = "The returned list contains elements not present in the function arguments."
  )
})
