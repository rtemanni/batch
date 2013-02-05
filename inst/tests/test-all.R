require(pasilla)
datafile=system.file("extdata/pasilla_gene_counts.tsv", package="pasilla")
counts = read.table(datafile, header=TRUE, row.names=1)
design = data.frame(row.names=colnames(counts), 
                    condition=c("untreated","untreated","untreated","untreated","treated","treated","treated"),
                    libType=c("single-end","single-end","paired-end","paired-end","single-end","paired-end","paired-end"))

test_that("qnorm works", {
  qcounts=qNorm(counts)
  expect_that(dim(qcounts), equals(dim(counts)))
})

test_that("svd works", {
  qcounts=qNorm(counts)
  res=makeSVD(qcounts)
  expect_false(is.null(res$v))
  expect_false(is.null(res$d))
  expect_equal(ncol(res$v), length(res$d))
})

test_that("pcres works", {
  qcounts=qNorm(counts)
  res=makeSVD(qcounts)
  res=pcRes(res$v,res$d, design$condition,design$libType)
  expect_that(names(res), equals(c("propVar","cumPropVar","cond.R2","batch.R2")))
})

test_that("plotpc works", {
  qcounts=qNorm(counts)
  res=makeSVD(qcounts)
  plotPC(res$v,res$d, col=design$condition, pch=ifelse(design$libType=="single-end", 21, 22))
})

test_that("logcpm works", {
  qcounts=qNorm(counts)
  res=log2CPM(qcounts)
  expect_that(dim(res$y), equals(dim(qcounts)))
  expect_that(names(res$y), equals(names(qcounts)))
  expect_that(res$lib.size, equals(colSums(qcounts)))
})

test_that("combatmod works", {
  qcounts=qNorm(counts)
  res=log2CPM(qcounts)
  y=res$y
  lib.size=res$lib.size
  
  res=combatMod(y, design$libType, design$condition)
  expect_false(is.null(res$bayesdata))
  expect_that(dim(res$bayesdata), equals(dim(y)))
  expect_false(is.null(res$gamma.star))
  expect_false(is.null(res$delta.star))
})

test_that("voomMod works", {
  qcounts=qNorm(counts)
  res=log2CPM(qcounts)
  y=res$y
  lib.size=res$lib.size
  res=combatMod(y, design$libType, design$condition)
  y=res$bayesdata
  
  mod=model.matrix(~1+condition, data=design)
  voomRes = voomMod(y, mod, lib.size, plot=FALSE)
  expect_that(voomRes, is_a("EList"))
})

