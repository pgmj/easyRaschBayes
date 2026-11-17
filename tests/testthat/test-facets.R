# Tests for many-facet helpers (R/facets.R) and contrast-aware item
# effect extraction (R/item_draws.R).

# --- contr_sum_zero() ----------------------------------------------------

test_that("contr_sum_zero() is orthonormal and sums to zero", {
  for (k in 2:7) {
    h <- contr_sum_zero(k)
    expect_equal(dim(h), c(k, k - 1L))
    expect_equal(unname(crossprod(h)), diag(k - 1L))
    expect_equal(unname(colSums(h)), rep(0, k - 1L))
    # Implied prior covariance of the level effects: I - J/k
    expect_equal(unname(h %*% t(h)), diag(k) - 1 / k)
  }
  expect_equal(rownames(contr_sum_zero(c("a", "b", "c"))), c("a", "b", "c"))
  expect_error(contr_sum_zero(1), "at least two")
  # Callable by name as a contrasts function (fourth review, Finding 2)
  x <- factor(c("a", "b", "c", "a"))
  contrasts(x) <- "contr_sum_zero"
  expect_equal(unname(contrasts(x)), unname(contr_sum_zero(3)))
  old <- options(contrasts = c("contr_sum_zero", "contr.poly"))
  mm <- model.matrix(~ x)
  options(old)
  expect_equal(unname(mm[, -1]), unname(contr_sum_zero(3)[x, ]))
  expect_equal(contr_sum_zero(c("a", "b"), contrasts = FALSE),
               matrix(c(1, 0, 0, 1), 2, dimnames = list(c("a", "b"),
                                                        c("a", "b"))))
  expect_error(contr_sum_zero(3, sparse = TRUE), "sparse")
  expect_error(contr_sum_zero(2.5), "whole number")
})

test_that("contr.sum gives the last level a larger prior variance", {
  # The reason for contr_sum_zero(): iid coefficient priors
  k <- 4
  v_sum <- diag(contr.sum(k) %*% t(contr.sum(k)))
  v_ozs <- diag(contr_sum_zero(k) %*% t(contr_sum_zero(k)))
  expect_equal(unname(v_sum[k] / v_sum[1]), k - 1)
  expect_equal(unname(v_ozs), rep((k - 1) / k, k))
})

# --- contrast lookup ------------------------------------------------------

test_that(".factor_contrast_matrix() reads stored and default contrasts", {
  lev <- c("A", "B", "C")
  x <- factor(c("A", "B", "C"))
  expect_equal(unname(.factor_contrast_matrix(x, lev)),
               unname(contr.treatment(3)))
  contrasts(x) <- contr.sum(3)
  expect_equal(unname(.factor_contrast_matrix(x, lev)), unname(contr.sum(3)))
  y <- factor(c("A", "B", "C"))
  contrasts(y) <- "contr.helmert"
  expect_equal(unname(.factor_contrast_matrix(y, lev)),
               unname(contr.helmert(3)))
  z <- factor(c("A", "B", "C"), ordered = TRUE)
  expect_equal(unname(.factor_contrast_matrix(z, lev)), unname(contr.poly(3)))
  expect_equal(rownames(.factor_contrast_matrix(z, lev)), lev)
  # Character columns are converted with the default (treatment)
  expect_equal(unname(.factor_contrast_matrix(c("A", "B", "C"), lev)),
               unname(contr.treatment(3)))
})

# --- extraction under different codings, checked against brms ------------

# For each coding: build the brms design matrix, check that its item columns
# are the contrast rows (the assumption behind the mapping), then name fake
# draws after the brms coefficients and check the recovered item effects.
check_coding <- function(set_contrasts, labels) {
  d <- expand.grid(id = 1:3, item = labels, stringsAsFactors = FALSE)
  d$item <- factor(d$item, levels = labels)
  d$y <- rep(0:1, length.out = nrow(d))
  if (!is.null(set_contrasts)) contrasts(d$item) <- set_contrasts
  sd <- brms::make_standata(y ~ 1 + item + (1 | id), data = d,
                            family = brms::bernoulli())
  X <- sd$X[, colnames(sd$X) != "Intercept", drop = FALSE]
  lev <- .item_levels(list(data = d), "item")
  cm <- .item_contrast_candidates(lev)[[1L]]
  expect_equal(unname(X), unname(cm[as.character(d$item), , drop = FALSE]))

  k <- length(labels)
  # True item effects in the span of the contrasts (sum to zero, or first
  # level 0 for treatment coding)
  e <- cm %*% seq(0.3, by = -0.4, length.out = k - 1L)
  beta <- solve(crossprod(cm), crossprod(cm, e))
  draws <- as.data.frame(matrix(rep(beta, each = 2), nrow = 2))
  names(draws) <- paste0("b_", colnames(X))
  draws$b_Intercept <- 0.1
  tau <- .item_tau_draws(draws, lev, "item", ordinal = FALSE)
  expect_equal(vapply(tau, function(m) m[1, 1], numeric(1)),
               setNames(-(0.1 + as.numeric(e)), labels))
}

test_that("item effects are recovered under treatment coding", {
  check_coding(NULL, c("I1", "I2", "I3", "I4"))
})

test_that("item effects are recovered under contr.sum and contr_sum_zero", {
  check_coding(contr.sum(4), c("I1", "I2", "I3", "I4"))
  check_coding(contr_sum_zero(4), c("I1", "I2", "I3", "I4"))
  check_coding(contr.helmert(3), c("a", "b", "c"))
})

test_that("non-plain labels are no longer mapped as treatment coding", {
  # Before 0.4.0.9001 these went silently wrong with contr.sum
  check_coding(contr.sum(3), c("item 1", "item-2", "item.3"))
  check_coding(contr_sum_zero(3), c("item 1", "item-2", "item.3"))
})

test_that("treatment coding with labels that brms renames", {
  check_coding(NULL, c("item 1", "item-2", "item (3)"))
})

test_that("fixed effects of another variable with the same prefix are ignored", {
  lev <- c("R1", "R2", "R3")
  attr(lev, "contrast_spec") <- list(stored = NULL, ordered = FALSE)
  d <- data.frame(b_Intercept = 0, b_raterR2 = 0.5, b_raterR3 = -0.5,
                  b_rater_exp = 9)
  expect_equal(unname(.item_effect_draws(d, lev, "rater")[1, ]),
               c(0, 0.5, -0.5))
})

test_that("coefficients from other contrasts than the stored ones are an error", {
  # Fitted under options(contrasts = c("contr.sum", ...)): brms names the
  # coefficients item1, item2. Read back without stored contrasts, the
  # treatment names would be itemitem2, itemitem3.
  lev <- c("item 1", "item 2", "item 3")
  attr(lev, "contrast_spec") <- list(stored = NULL, ordered = FALSE)
  d <- data.frame(b_Intercept = 0, b_item1 = 0.2, b_item2 = -0.1)
  expect_error(.item_effect_draws(d, lev, "item"), "options\\(contrasts")
})

# Evaluate code with options(contrasts) set to sum coding
with_sum_option <- function(code) {
  op <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(op))
  force(code)
}

test_that("contrasts option of the reading session is tried after the default", {
  lev <- c("item 1", "item 2", "item 3")
  attr(lev, "contrast_spec") <- list(stored = NULL, ordered = FALSE)
  # Fitted with the default treatment coding, read under contr.sum
  trt <- data.frame(b_Intercept = 0, b_itemitem2 = 0.5, b_itemitem3 = -0.5)
  with_sum_option({
    expect_equal(unname(.item_effect_draws(trt, lev, "item")[1, ]),
                 c(0, 0.5, -0.5))
  })
  # Fitted under options(contrasts = contr.sum) and read under it
  sum_d <- data.frame(b_Intercept = 0, b_item1 = 0.2, b_item2 = -0.1)
  with_sum_option({
    expect_equal(unname(.item_effect_draws(sum_d, lev, "item")[1, ]),
                 c(0.2, -0.1, -0.1))
  })
})

test_that("names that fit both treatment and the session option are an error", {
  # Labels 0..3: treatment names item1..item3 (levels 1..3), contr.sum
  # names item1..item3 (columns 1..3)
  lev <- as.character(0:3)
  attr(lev, "contrast_spec") <- list(stored = NULL, ordered = FALSE)
  d <- data.frame(b_Intercept = 0, b_item1 = 1, b_item2 = 2, b_item3 = 3)
  expect_equal(unname(.item_effect_draws(d, lev, "item")[1, ]), 0:3)
  with_sum_option({
    expect_error(.item_effect_draws(d, lev, "item"), "match both")
  })
})

test_that("ordinal common thresholds with sum-to-zero item effects", {
  d <- data.frame(item = factor(rep(c("A", "B", "C"), 2)))
  contrasts(d$item) <- contr_sum_zero(3)
  lev <- .item_levels(list(data = d), "item")
  e <- c(0.5, -0.2, -0.3)
  beta <- as.numeric(crossprod(contr_sum_zero(3), e))
  draws <- data.frame("b_Intercept[1]" = -1, "b_Intercept[2]" = 1,
                      b_item1 = beta[1], b_item2 = beta[2],
                      check.names = FALSE)
  tau <- .item_tau_draws(draws, lev, "item", ordinal = TRUE)
  expect_equal(unname(tau$A[1, ]), c(-1, 1) - 0.5)
  expect_equal(unname(tau$C[1, ]), c(-1, 1) + 0.3)
})

test_that("random items never evaluate the item contrasts", {
  # contr.poly() fails above 95 levels, so building contrasts for an
  # ordered item factor in a (1 | item) model would stop
  d <- data.frame(item = factor(sprintf("i%03d", 1:100), ordered = TRUE))
  expect_no_error(lev <- .item_levels(list(data = d), "item"))
  draws <- as.data.frame(matrix(0.1, nrow = 2, ncol = 101))
  names(draws) <- c("b_Intercept", paste0("r_item[", lev, ",Intercept]"))
  expect_no_error(tau <- .item_tau_draws(draws, lev, "item",
                                         ordinal = FALSE))
  expect_length(tau, 100L)
})

# --- .resolve_facet() ----------------------------------------------------

# A brmsfit shell without sampling: enough for formula and data.
facet_model <- function(formula, data, family = brms::bernoulli()) {
  suppressMessages(brms::brm(formula, data = data, family = family,
                             empty = TRUE))
}

mfrm_data <- function(persons_by_rater, items = c("c1", "c2", "c3")) {
  d <- do.call(rbind, lapply(names(persons_by_rater), function(r) {
    expand.grid(id = persons_by_rater[[r]], item = items, rater = r,
                stringsAsFactors = FALSE)
  }))
  d$y <- rep(0:1, length.out = nrow(d))
  d
}

test_that("connected design: cells, roles, no warning", {
  d <- mfrm_data(list(R1 = 1:6, R2 = 4:9, R3 = c(1, 9)))
  m <- facet_model(y ~ 1 + item + (1 | rater) + (1 | id), d)
  expect_no_warning(f <- .resolve_facet(m, "id", "item", "rater", "test"))
  expect_equal(f$rank_deficiency, 0L)
  expect_equal(levels(f$cell)[1:3], c("c1:R1", "c1:R2", "c1:R3"))
  expect_equal(nlevels(f$cell), 9L)
  expect_equal(f$cell_item[4], "c2")
  expect_equal(f$cell_facet[4], "R1")
  expect_equal(as.character(f$cell[1]), paste(d$item[1], d$rater[1], sep = ":"))
  expect_equal(f$roles, list(fixed = FALSE, random = TRUE, thres = FALSE))
  expect_false(f$halo)
})

test_that("fixed facet with contrasts, halo term, thresholds by facet", {
  d <- mfrm_data(list(R1 = 1:6, R2 = 4:9))
  d$rater <- factor(d$rater)
  contrasts(d$rater) <- contr_sum_zero(2)
  m <- facet_model(y ~ 1 + item + rater + (1 | id) + (1 | rater:id), d)
  f <- .resolve_facet(m, "id", "item", "rater", "test")
  expect_true(f$roles$fixed)
  expect_true(f$halo)

  d$yo <- rep(1:3, length.out = nrow(d))
  m2 <- facet_model(yo | thres(gr = rater) ~ 1 + item + (1 | id), d,
                    family = brms::acat())
  f2 <- .resolve_facet(m2, "id", "item", "rater", "test")
  expect_equal(f2$roles, list(fixed = FALSE, random = FALSE, thres = TRUE))
})

test_that("a prefix of another variable name is not a role", {
  d <- mfrm_data(list(R1 = 1:6, R2 = 4:9))
  d$r <- d$rater
  m <- facet_model(y ~ 1 + item + (1 | rater) + (1 | id), d)
  # 'r' is not in model$data, since brms keeps formula variables only
  expect_error(.resolve_facet(m, "id", "item", "r", "test"), "not found")
  expect_false(.formula_roles(m, "r")$random)
})

test_that("nested and gr() grouping terms are read as brms expands them", {
  expect_equal(.group_factors(c("id/rater", "gr(rater, by = g)", "rater:id")),
               c("id", "id:rater", "rater", "rater:id"))
  d <- mfrm_data(list(R1 = 1:6, R2 = 4:9))
  m <- facet_model(y ~ 1 + item + (1 | id/rater), d)
  f <- .resolve_facet(m, "id", "item", "rater", "test")
  expect_true(f$halo)
  expect_equal(f$halo_term, "id:rater")
  # by = g needs g constant within rater
  d$g <- ifelse(d$rater == "R1", "a", "b")
  m2 <- facet_model(y ~ 1 + item + (1 | gr(rater, by = g)) + (1 | id), d)
  expect_true(.formula_roles(m2, "rater")$random)
})

test_that("the halo message names the term and the right infit figure", {
  fr <- list(halo = TRUE, halo_term = "id:rater")
  expect_message(.facet_halo_message(fr, "q3_statistic"),
                 "(1 | id:rater)", fixed = TRUE)
  expect_message(
    .facet_halo_message(fr, "infit_statistic",
                        .infit_halo_note("facet", "item", "rater")),
    "1-2% of fitting levels of rater"
  )
  expect_match(.infit_halo_note("cell", "item", "rater"), "2% of fitting cells")
  expect_silent(.facet_halo_message(list(halo = FALSE), "q3_statistic"))
  expect_message(.facet_halo_message(fr, "q3_statistic", method = "theta"),
                 "expected values from the posterior include its effects")
})

test_that("colliding cell labels stop with a clear error", {
  d <- rbind(
    expand.grid(id = 1:4, item = c("a:b", "a"), rater = c("c", "b:c"),
                stringsAsFactors = FALSE)
  )
  d$y <- rep(0:1, length.out = nrow(d))
  m <- facet_model(y ~ 1 + item + (1 | rater) + (1 | id), d)
  expect_error(.resolve_facet(m, "id", "item", "rater", "test"),
               "same cell label 'a:b:c'")
})

test_that("raters who never rate the same persons are reported", {
  d <- mfrm_data(list(R1 = 1:4, R2 = 1:4, R3 = 5:8))
  m <- facet_model(y ~ 1 + item + (1 | rater) + (1 | id), d)
  expect_warning(f <- .resolve_facet(m, "id", "item", "rater", "test"),
                 "\\{R1, R2\\}, \\{R3\\}")
  expect_equal(f$rank_deficiency, 1L)
})

test_that("items nested in raters are reported", {
  # Each rater uses different items for the same persons: rater and item
  # set effects are confounded although the persons connect the raters
  d <- rbind(
    expand.grid(id = 1:5, item = c("c1", "c2"), rater = "R1",
                stringsAsFactors = FALSE),
    expand.grid(id = 1:5, item = c("c3", "c4"), rater = "R2",
                stringsAsFactors = FALSE)
  )
  d$y <- rep(0:1, length.out = nrow(d))
  m <- facet_model(y ~ 1 + item + (1 | rater) + (1 | id), d)
  expect_warning(f <- .resolve_facet(m, "id", "item", "rater", "test"),
                 "only observed together")
  expect_equal(f$rank_deficiency, 1L)
})

test_that("a facet outside the model is not checked for connectivity", {
  d <- mfrm_data(list(R1 = 1:4, R3 = 5:8))
  m <- facet_model(y ~ 1 + item + (1 | id) + (0 + rater | item), d)
  # rater is in the data (used in a slope term) but not as a facet role
  expect_no_warning(f <- .resolve_facet(m, "id", "item", "rater", "test"))
  expect_true(is.na(f$rank_deficiency))
})

test_that("thres(gr = ) is matched on the whole facet name", {
  d <- mfrm_data(list(R1 = 1:6, R2 = 4:9))
  d$rater.grp <- d$rater
  d$yo <- rep(1:3, length.out = nrow(d))
  m <- facet_model(yo | thres(gr = rater.grp) ~ 1 + item + (1 | id) +
                     (1 | rater), d, family = brms::acat())
  expect_false(.formula_roles(m, "rater")$thres)
  expect_true(.formula_roles(m, "rater.grp")$thres)
  m2 <- facet_model(yo | thres(2, gr = rater) ~ 1 + item + (1 | id), d,
                    family = brms::acat())
  expect_true(.formula_roles(m2, "rater")$thres)
})

test_that("an item or facet variable called cell is rejected", {
  d <- mfrm_data(list(R1 = 1:4, R2 = 1:4))
  names(d)[names(d) == "item"] <- "cell"
  m <- facet_model(y ~ 1 + cell + (1 | rater) + (1 | id), d)
  expect_error(.resolve_facet(m, "id", "cell", "rater", "test"),
               "cannot be called 'cell'")
})

test_that("the three variables must differ", {
  d <- mfrm_data(list(R1 = 1:4, R2 = 1:4))
  m <- facet_model(y ~ 1 + item + (1 | rater) + (1 | id), d)
  expect_error(.resolve_facet(m, "id", "item", "item", "test"),
               "three different")
})

test_that("rank deficiency agrees with the rank of the full design", {
  set.seed(3)
  for (i in 1:5) {
    d <- expand.grid(id = 1:8, item = paste0("c", 1:4),
                     rater = paste0("R", 1:4), stringsAsFactors = FALSE)
    d <- d[runif(nrow(d)) < 0.15, ]
    X <- model.matrix(~ factor(id) + item + rater, data = d)
    expected <- (length(unique(d$id)) + length(unique(d$item)) +
                   length(unique(d$rater)) - 2L) - qr(X)$rank
    expect_equal(.design_rank_deficiency(d$id, d$item, d$rater), expected)
  }
})

# --- facet_var in the fit functions (mocked predictions) -----------------

# A brmsfit shell with posterior_epred() replaced by fixed probabilities
# from an additive dichotomous MFRM, so the facet wrappers run end to end
# without sampling. Every draw has the same parameters.
mock_mfrm <- function(env = parent.frame(),
                      persons_by_rater = list(R1 = 1:30, R2 = 16:45,
                                              R3 = c(1:10, 31:45))) {
  set.seed(21)
  d <- mfrm_data(persons_by_rater, items = c("c1", "c2", "c3"))
  theta <- rnorm(max(unlist(persons_by_rater)))
  b_item <- c(c1 = -0.5, c2 = 0, c3 = 0.5)
  sev <- c(R1 = -0.3, R2 = 0, R3 = 0.3)
  p <- plogis(theta[d$id] - b_item[d$item] - sev[d$rater])
  d$y <- rbinom(nrow(d), 1, p)
  m <- facet_model(y ~ 1 + item + (1 | rater) + (1 | id), d)
  testthat::local_mocked_bindings(
    posterior_epred = function(object, draw_ids = NULL, ...) {
      matrix(p, nrow = length(draw_ids), ncol = length(p), byrow = TRUE)
    },
    posterior_predict = function(object, draw_ids = NULL, ...) {
      matrix(stats::rbinom(length(draw_ids) * length(p), 1,
                           rep(p, each = length(draw_ids))),
             nrow = length(draw_ids))
    },
    ndraws = function(x) 4L,
    variables = function(x, ...) {
      c("b_Intercept", "sd_id__Intercept", "sd_rater__Intercept")
    },
    .package = "brms", .env = env
  )
  list(model = m, data = d, p = p)
}

test_that("infit_statistic() with facet_var: cells and pooled by item/facet", {
  mk <- mock_mfrm()
  set.seed(1)
  cell <- infit_statistic(mk$model, facet_var = rater)
  expect_equal(names(cell)[1:4], c("cell", "item", "rater", "draw"))
  expect_equal(dplyr::n_distinct(cell$cell), 9L)
  expect_equal(nrow(cell), 9L * 4L)
  expect_true(all(paste(cell$item, cell$rater, sep = ":") ==
                    as.character(cell$cell)))

  by_item <- infit_statistic(mk$model, facet_var = rater, by = "item")
  by_rater <- infit_statistic(mk$model, facet_var = rater, by = "facet")
  expect_equal(names(by_item)[1], "item")
  expect_equal(names(by_rater)[1], "rater")
  expect_equal(levels(by_rater$rater), c("R1", "R2", "R3"))
  expect_equal(as.character(unique(by_rater$rater)), c("R1", "R2", "R3"))
  # Both methods give the same type and order
  th_rater <- infit_statistic(mk$model, facet_var = rater, by = "facet",
                              method = "theta", ndraws_use = 2)
  expect_identical(levels(th_rater$rater), levels(by_rater$rater))
  expect_equal(as.character(unique(th_rater$rater)), c("R1", "R2", "R3"))

  # Observed pooled infit equals pooling the engine's cell sums directly
  d <- mk$data
  cellf <- factor(paste(d$item, d$rater, sep = ":"),
                  levels = levels(.resolve_facet(mk$model, "id", "item",
                                                 "rater", "t")$cell))
  pid <- match(d$id, sort(unique(d$id)))
  iid <- as.integer(cellf)
  X <- matrix(NA_real_, max(pid), nlevels(cellf))
  X[cbind(pid, iid)] <- d$y
  grp <- match(sub(".*:", "", levels(cellf)), c("R1", "R2", "R3"))
  ref <- .cond_infit_draw(cbind(1 - mk$p, mk$p), pid, iid, X, group = grp)
  obs <- by_rater$infit[by_rater$draw == 1]
  expect_equal(obs[order(by_rater$rater[by_rater$draw == 1])], ref$infit)

  post <- infit_post(cell)
  expect_equal(names(post$summary)[1:3], c("cell", "item", "rater"))
  expect_equal(names(post$hdi)[1:3], c("cell", "item", "rater"))
  expect_error(infit_statistic(mk$model, by = "item"), "facet_var")
})

test_that("item_restscore_statistic() and q3_statistic() with facet_var", {
  mk <- mock_mfrm()
  set.seed(2)
  rs <- item_restscore_statistic(mk$model, facet_var = rater)
  expect_equal(names(rs$result)[1:3], c("cell", "item", "rater"))
  expect_equal(names(rs$draws)[1:3], c("cell", "item", "rater"))
  expect_equal(nrow(rs$result), 9L)
  rsp <- item_restscore_post(rs)
  expect_equal(names(rsp$summary)[1:3], c("cell", "item", "rater"))
  # Plot rows (bottom to top) in the order of the table
  expect_equal(rev(levels(rsp$plot$data$item)),
               as.character(rsp$summary$cell))

  q3 <- q3_statistic(mk$model, facet_var = rater)
  types <- table(q3$pair_type[q3$draw == 1])
  # 9 cells: 3 raters x choose(3, 2) same-rater pairs, 3 items x 3
  # same-item pairs, the rest neither
  expect_equal(as.numeric(types[c("same rater", "same item", "neither")]),
               c(9, 9, 18))
  q3p <- q3_post(q3)
  expect_true("pair_type" %in% names(q3p$summary))
  expect_equal(names(q3p$hdi)[1:4],
               c("item_pair", "item_1", "item_2", "pair_type"))
})

test_that("item_restscore_post() plots factor cells in table order", {
  lev <- c("c2:R1", "c10:R1", "c1:R1")
  cell <- factor(lev, levels = lev)
  rs <- list(
    result = tibble::tibble(cell = cell, item = c("c2", "c10", "c1"),
                            rater = "R1", gamma_obs = c(.1, .2, .3),
                            gamma_exp = .2, gamma_diff = 0, ppp = .5,
                            p_adj = 1),
    draws = tibble::tibble(cell = rep(cell, each = 4), item = "x",
                           rater = "R1", draw = rep(1:4, 3),
                           gamma_rep = seq(0, 1, length.out = 12))
  )
  rsp <- item_restscore_post(rs)
  expect_equal(as.character(rsp$summary$cell), lev)
  expect_equal(rev(levels(rsp$plot$data$item)), lev)
})

test_that("q3_statistic() with facet_var leaves out unshared cell pairs", {
  # R1 and R3 rate no person in common, so their "neither" and "same item"
  # pairs have no Q3
  mk <- mock_mfrm(persons_by_rater = list(R1 = 1:20, R2 = 11:30, R3 = 21:40))
  q3 <- q3_statistic(mk$model, facet_var = rater, method = "theta",
                     ndraws_use = 2)
  expect_false(anyNA(q3$q3[q3$draw == 1]))
  types <- table(q3$pair_type[q3$draw == 1])
  # Same rater: 3 x 3. Same item: R1-R2 and R2-R3 for 3 items. Neither:
  # 6 per sharing rater pair, 2 sharing pairs.
  expect_equal(as.numeric(types[c("same rater", "same item", "neither")]),
               c(9, 6, 12))
})

test_that("conditional Q3 with facet_var computes shared cell pairs only", {
  mk <- mock_mfrm(persons_by_rater = list(R1 = 1:20, R2 = 11:30, R3 = 21:40))
  q3 <- q3_statistic(mk$model, facet_var = rater, ndraws_use = 2)
  types <- table(q3$pair_type[q3$draw == 1])
  expect_equal(as.numeric(types[c("same rater", "same item", "neither")]),
               c(9, 6, 12))
})

test_that(".q3_shared_pairs() follows the row-major upper triangle", {
  # Person 1 answers items 1, 2; person 2 answers 2, 4. Pairs in order
  # 1:2, 1:3, 1:4, 2:3, 2:4, 3:4
  expect_equal(.q3_shared_pairs(c(1, 1, 2, 2), c(1, 2, 2, 4), 2, 4),
               c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE))
})

test_that("facet functions reject several responses per cell", {
  mk <- mock_mfrm()
  m <- mk$model
  m$data <- rbind(m$data, m$data[1, ])
  expect_error(infit_statistic(m, facet_var = rater), "item x rater")
})

test_that(".label_cells() puts cell, item and facet first", {
  fr <- list(cell = factor(c("a:R1", "b:R1"), levels = c("a:R1", "b:R1")),
             cell_item = c("a", "b"), cell_facet = c("R1", "R1"),
             item_levels = c("a", "b"), facet_levels = "R1")
  df <- tibble::tibble(.cell = c("b:R1", "a:R1"), draw = 1:2, x = 3:4)
  out <- .label_cells(df, fr, "item", "rater")
  expect_equal(names(out), c("cell", "item", "rater", "draw", "x"))
  expect_equal(as.character(out$item), c("b", "a"))
})

# --- facet_parameters() --------------------------------------------------

# Fake draws for a shell model: one row per draw.
fake_draws <- function(cols, n = 6) {
  d <- as.data.frame(lapply(cols, function(v) rep(v, length.out = n)),
                     check.names = FALSE)
  d$.chain <- rep(1:2, each = n / 2)
  d
}

test_that("facet locations: severity sign, centred, coding-invariant", {
  d <- mfrm_data(list(R1 = 1:6, R2 = 4:9, R3 = c(1, 9)))
  r <- c(R1 = 0.4, R2 = 0, R3 = -0.7)          # added to the linear predictor
  base <- list(b_Intercept = 0.2, b_itemc2 = -0.3, b_itemc3 = 0.5,
               "r_id[1,Intercept]" = 0.1, "r_id[2,Intercept]" = -0.1)

  m_rand <- facet_model(y ~ 1 + item + (1 | rater) + (1 | id), d)
  dr <- fake_draws(c(base, setNames(as.list(r),
                                    paste0("r_rater[", names(r), ",Intercept]")),
                     list(sd_rater__Intercept = 0.5)))
  fp_r <- .facet_parameters_core(m_rand, dr, "rater", "item", "id",
                                 draws = TRUE, prob = 0.9)
  expected <- -(r - mean(r))
  expect_equal(fp_r$summary$location, unname(round(expected, 4)))
  expect_equal(unname(fp_r$draws_matrix[, 1]), unname(expected))
  expect_equal(fp_r$facet_sd$mean, 0.5)
  expect_null(fp_r$thresholds)

  # Same effects as a fixed facet with sum-to-zero contrasts
  d2 <- d
  d2$rater <- factor(d2$rater)
  contrasts(d2$rater) <- contr_sum_zero(3)
  m_fix <- facet_model(y ~ 1 + item + rater + (1 | id), d2)
  beta <- as.numeric(crossprod(contr_sum_zero(3), r - mean(r)))
  base_f <- base
  base_f$b_Intercept <- 0.2 + mean(r)
  df <- fake_draws(c(base_f, list(b_rater1 = beta[1], b_rater2 = beta[2])))
  fp_f <- .facet_parameters_core(m_fix, df, "rater", "item", "id",
                                 draws = FALSE, prob = 0.9)
  expect_equal(fp_f$summary$location, fp_r$summary$location)
  expect_equal(fp_f$summary$fair_avg, fp_r$summary$fair_avg)
  expect_null(fp_f$facet_sd)
})

test_that("fair average is the mean expected rating at the mean person", {
  d <- mfrm_data(list(R1 = 1:6, R2 = 4:9))
  m <- facet_model(y ~ 1 + item + (1 | rater) + (1 | id), d)
  dr <- fake_draws(list(b_Intercept = 0.2, b_itemc2 = -0.3, b_itemc3 = 0.5,
                        "r_id[1,Intercept]" = 0.3, "r_id[2,Intercept]" = -0.1,
                        "r_rater[R1,Intercept]" = 0.4,
                        "r_rater[R2,Intercept]" = -0.4))
  fp <- .facet_parameters_core(m, dr, "rater", "item", "id", FALSE, 0.9)
  theta <- 0.1
  eta_items <- 0.2 + c(0, -0.3, 0.5)
  manual <- c(mean(plogis(theta + eta_items + 0.4)),
              mean(plogis(theta + eta_items - 0.4)))
  expect_equal(fp$summary$fair_avg, round(manual, 3))
  expect_equal(fp$summary$n, as.integer(table(d$rater)))
  expect_equal(fp$summary$observed_avg,
               round(as.numeric(tapply(d$y, d$rater, mean)), 3))
})

test_that("rater-specific thresholds give a thresholds table", {
  d <- mfrm_data(list(R1 = 1:6, R2 = 4:9))
  d$yo <- rep(1:3, length.out = nrow(d))
  m <- facet_model(yo | thres(gr = rater) ~ 1 + item + (1 | id), d,
                   family = brms::acat())
  dr <- fake_draws(list("b_Intercept[R1,1]" = -1, "b_Intercept[R1,2]" = 1,
                        "b_Intercept[R2,1]" = -0.2, "b_Intercept[R2,2]" = 0.6,
                        b_itemc2 = 0.3, b_itemc3 = -0.1,
                        "r_id[1,Intercept]" = 0))
  fp <- .facet_parameters_core(m, dr, "rater", "item", "id", FALSE, 0.9)
  # Locations: threshold means 0 and 0.2, centred -> -0.1, 0.1
  expect_equal(fp$summary$location, c(-0.1, 0.1))
  expect_equal(fp$thresholds$location, c(-1, 1, -0.2, 0.6) - 0.1)
  expect_equal(names(fp$thresholds)[1:2], c("rater", "threshold"))
  # Both fair averages lie inside the category range
  expect_true(all(fp$summary$fair_avg > 1 & fp$summary$fair_avg < 3))
})

test_that("different threshold counts per rater stop (review 4, F1)", {
  d <- mfrm_data(list(R1 = 1:6, R2 = 4:9))
  d$yo <- rep(1:3, length.out = nrow(d))
  m <- facet_model(yo | thres(gr = rater) ~ 1 + item + (1 | id), d,
                   family = brms::acat())
  dr <- fake_draws(list("b_Intercept[R1,1]" = -1, "b_Intercept[R1,2]" = 1,
                        "b_Intercept[R2,1]" = -0.2,
                        b_itemc2 = 0.3, b_itemc3 = -0.1,
                        "r_id[1,Intercept]" = 0))
  expect_error(.facet_parameters_core(m, dr, "rater", "item", "id", FALSE,
                                      0.9),
               "R1: 2, R2: 1.*thres\\(2, gr = rater\\)")
})

test_that("a covariate named like the facet is not its effect (review 4, F4)", {
  d <- mfrm_data(list(R1 = 1:6, R2 = 4:9))
  d$yo <- rep(1:3, length.out = nrow(d))
  d$rater_exp <- rep(c(0, 1), length.out = nrow(d))
  m <- facet_model(yo | thres(gr = rater) ~ 1 + item + rater_exp + (1 | id),
                   d, family = brms::acat())
  dr <- fake_draws(list("b_Intercept[R1,1]" = -1, "b_Intercept[R1,2]" = 1,
                        "b_Intercept[R2,1]" = -0.2, "b_Intercept[R2,2]" = 0.6,
                        b_itemc2 = 0.3, b_itemc3 = -0.1, b_rater_exp = 0.5,
                        "r_id[1,Intercept]" = 0))
  fp <- .facet_parameters_core(m, dr, "rater", "item", "id", FALSE, 0.9)
  expect_equal(fp$summary$location, c(-0.1, 0.1))
  # Without the variable names the old check still stops
  expect_error(.item_effect_draws(dr, c("R1", "R2"), "rater"),
               "Could not match")
  expect_null(.item_effect_draws(dr, c("R1", "R2"), "rater",
                                 other_vars = names(m$data)))
})

test_that("thresholds by facet with item effects not found is an error", {
  d <- mfrm_data(list(R1 = 1:6, R2 = 4:9))
  d$yo <- rep(1:3, length.out = nrow(d))
  m <- facet_model(yo | thres(gr = rater) ~ 1 + item + (1 | id), d,
                   family = brms::acat())
  # Item coefficients under names that match no coding
  dr <- fake_draws(list("b_Intercept[R1,1]" = -1, "b_Intercept[R1,2]" = 1,
                        "b_Intercept[R2,1]" = -0.2, "b_Intercept[R2,2]" = 0.6,
                        b_itemX = 0.3, b_itemY = -0.1, b_itemZ = 0,
                        b_itemW = 0, "r_id[1,Intercept]" = 0))
  expect_error(.facet_parameters_core(m, dr, "rater", "item", "id", FALSE,
                                      0.9),
               "could not find the parameters of 'item'")
})

test_that("vectorised cumulative, sratio, cratio expected ratings", {
  set.seed(5)
  theta <- rnorm(6)
  tau <- t(apply(matrix(rnorm(18), 6, 3), 1L, sort))
  for (fam in c("cumulative", "sratio", "cratio")) {
    for (li in list(stats::plogis, stats::pnorm)) {
      got <- .expected_rating(theta, tau, fam, li)
      ref <- vapply(1:6, function(s) {
        p <- .compute_cat_probs(theta[s], tau[s, ], FALSE,
                                fam == "cumulative", fam == "sratio",
                                fam == "cratio", linkinv = li)
        sum(p * seq_along(p))
      }, numeric(1))
      expect_equal(got, ref)
    }
  }
})

test_that("vectorised acat expected rating matches the scalar version", {
  set.seed(4)
  theta <- rnorm(5)
  tau <- matrix(rnorm(15), 5, 3)
  got <- .expected_rating(theta, tau, "acat", stats::plogis)
  ref <- vapply(1:5, function(s) {
    p <- .compute_cat_probs(theta[s], tau[s, ], TRUE, FALSE, FALSE, FALSE)
    sum(p * seq_along(p))
  }, numeric(1))
  expect_equal(got, ref)
  expect_equal(.expected_rating(theta, tau[, 1, drop = FALSE], "acat",
                                stats::plogis),
               1 + plogis(theta - tau[, 1]))
})

test_that("facet_parameters() rejects hurdle families", {
  d <- mfrm_data(list(R1 = 1:4, R2 = 3:6))
  m <- facet_model(y ~ 1 + item + (1 | id) + (1 | rater), d)
  m$formula$family$family <- "hurdle_acat"
  expect_error(.facet_parameters_core(m, fake_draws(list(b_Intercept = 0)),
                                      "rater", "item", "id", FALSE, 0.9),
               "hurdle")
})

test_that("facet_parameters() needs the facet in the model", {
  d <- mfrm_data(list(R1 = 1:4, R2 = 3:6))
  m <- facet_model(y ~ 1 + item + (1 | id) + (0 + rater | item), d)
  expect_error(.facet_parameters_core(m, fake_draws(list(b_Intercept = 0)),
                                      "rater", "item", "id", FALSE, 0.9),
               "not in the model")
})

test_that("plot_residual_pca() with facet_var counts person x facet units", {
  mk <- mock_mfrm()
  m <- mk$model
  # Drop one rating, so one id x rater combination is incomplete
  m$data <- m$data[-1, ]
  p <- mk$p[-1]
  testthat::local_mocked_bindings(
    posterior_epred = function(object, draw_ids = NULL, ...) {
      matrix(p, nrow = length(draw_ids), ncol = length(p), byrow = TRUE)
    },
    posterior_predict = function(object, draw_ids = NULL, ...) {
      matrix(stats::rbinom(length(draw_ids) * length(p), 1,
                           rep(p, each = length(draw_ids))),
             nrow = length(draw_ids))
    },
    .package = "brms"
  )
  n_units <- nrow(unique(mk$data[c("id", "rater")]))
  # The shell has no draws for the later steps, so only the message is
  # checked
  expect_message(
    try(plot_residual_pca(m, facet_var = rater, ndraws_use = 2),
        silent = TRUE),
    paste0("uses the ", n_units - 1L, " id x rater combinations")
  )
})
