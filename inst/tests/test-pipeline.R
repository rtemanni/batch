require(pasilla)
datafile=system.file("extdata/pasilla_gene_counts.tsv", package="pasilla")
counts = read.table(datafile, header=TRUE, row.names=1)
design = data.frame(row.names=colnames(counts), 
                    condition=c("untreated","untreated","untreated","untreated","treated","treated","treated"),
                    libType=c("single-end","single-end","paired-end","paired-end","single-end","paired-end","paired-end"))

test_that("pipeline works", {
  mod=model.matrix(~1+condition, data=design)
  res=batchSEQ(counts, mod, design$batch)
  expect_that(res$elist, is_a("EList"))
  expect_false(is.null(res$combatEstimates$gamma.star))
  expect_false(is.null(res$combatEstimates$delta.star))
})
