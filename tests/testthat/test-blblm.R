test_that("return object in blblm class", {
  blblm_fit <- blblm(speed ~ dist, data = cars, m = 5, B = 100)
  expect_s3_class(blblm_fit, "blblm")
})
test_that("check dimension of coef returned", {
  blblm_fit <- blblm(speed ~ dist, data = cars, m = 5, B = 100)
  co <- coef(blblm_fit)
  expect_equal(length(co), 2)
})
