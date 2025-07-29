AUCBoot
function (pred, target, nboot = 10000, seed = NA, digits = NULL)
{
  assert_lengths(pred, target)
  assert_noNA(pred)
  assert_01(target)
  y <- as.logical(target)
  ord <- order(pred, y)
  pred <- pred[ord]
  y <- y[ord]
  if (!is.na(seed)) {
    old <- .Random.seed
    on.exit({
      .Random.seed <<- old
    })
    set.seed(seed)
  }
  repl <- boot_auc_sorted_tab(pred, y, nboot)
  if (nbNA <- sum(is.na(repl)))
    warning2("%d/%d bootstrap replicates were mono-class.",
             nbNA, nboot)
  res <- c(Mean = mean(repl, na.rm = TRUE), stats::quantile(repl,
                                                            c(0.025, 0.975), na.rm = TRUE), Sd = stats::sd(repl,
                                                                                                           na.rm = TRUE))
  round2(res, digits)
}
<bytecode: 0x55e390286188>
  <environment: namespace:bigstatsr>
  