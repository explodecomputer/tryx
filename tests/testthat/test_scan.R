context("initialise")
library(tryx)

test_that("initialise", 
{
	a <- extract_instruments("ieu-a-300")
	b <- extract_outcome_data(a$SNP, "ieu-a-7", access_token=NULL)
	dat <- harmonise_data(a,b)
	X <- Tryx$new(dat)
	expect_true(nrow(X$output$dat) == sum(dat$mr_keep))
})


