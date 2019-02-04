ihh2ihs <- function(res_ihh,
                    freqbin = 0.025,
                    minmaf = 0.05) {
  if (!all(c("CHR", "POSITION", "freq_A", "iHH_A", "iHH_D") %in% colnames(res_ihh))) {
    stop("Data does not contain the columns 'CHR', 'POSITION', 'freq_A', iHH_A' and 'iHH_D'")
  }
  res_ihh = res_ihh[(res_ihh$freq_A > minmaf) &
                      (res_ihh$freq_A < (1 - minmaf)),]

  summary_class_colnames = c("#mrk",
                             "mean(log(iHHA/iHHD))",
                             "sd(log(iHHA/iHHD))")
  if (freqbin > 0) {
    if (freqbin >= 1) {
      freqbin = (1 - 2 * minmaf) / round(freqbin)
    }
    freq_class = seq(minmaf, (1 - minmaf), freqbin)
    summary_class = matrix(0, length(freq_class) - 1, 3)
    colnames(summary_class) = summary_class_colnames
    tmp.classnames = rep(NA, nrow(summary_class))
    ihs = log(res_ihh$iHH_A / res_ihh$iHH_D)
    ihs[ihs == "Inf" | ihs == "-Inf"] = NA
    for (c in 1:(length(freq_class) - 1)) {
      lim_inf = freq_class[c]
      lim_sup = freq_class[c + 1]
      mrk_sel = ((res_ihh$freq_A >= lim_inf) &
                   (res_ihh$freq_A < lim_sup))
      bin_size = sum(mrk_sel)
      if (bin_size < 10) {
        warning(
          paste(
            "The number of SNPs with allele frequencies comprised between ",
            lim_inf,
            " and ",
            lim_sup,
            " is less than 10: you should probably increase freqbin\n",
            sep = ""
          )
        )
      }
      tmp.classnames[c] = paste(lim_inf, "-", lim_sup)
      summary_class[c, 1] = bin_size
      summary_class[c, 2] = mean(ihs[mrk_sel], na.rm = TRUE)
      summary_class[c, 3] = sd(ihs[mrk_sel], na.rm = TRUE)
      ihs[mrk_sel] = (ihs[mrk_sel] - summary_class[c, 2]) / summary_class[c, 3]
    }
  } else {
    freq_class = unique(res_ihh$freq_A)
    freq_class = freq_class[order(freq_class)]
    summary_class = matrix(0, length(freq_class), 3)
    tmp.classnames = rep(NA, nrow(summary_class))
    colnames(summary_class) = summary_class_colnames
    ihs = log(res_ihh$iHH_A / res_ihh$iHH_D)
    for (f in 1:length(freq_class)) {
      mrk_sel = (res_ihh$freq_A == freq_class[f])
      tmp.classnames[f] = freq_class[f]
      summary_class[f, 1] = sum(mrk_sel)
      summary_class[f, 2] = mean(ihs[mrk_sel], na.rm = TRUE)
      summary_class[f, 3] = sd(ihs[mrk_sel], na.rm = TRUE)
      ihs[mrk_sel] = (ihs[mrk_sel] - summary_class[f, 2]) / summary_class[f, 3]
    }
  }
  rownames(summary_class) = tmp.classnames

  tmp_pval = -log10(1 - 2 * abs(pnorm(ihs) - 0.5))
  tmp_pval2 = tmp_pval
  tmp_pval2[tmp_pval2 == "Inf"] = NA
  tmp_pval[tmp_pval == "Inf"] = max(tmp_pval2, na.rm = TRUE) + 1

  res_ihs = data.frame(
    res_ihh$CHR,
    res_ihh$POSITION,
    ihs,
    tmp_pval,
    stringsAsFactors = FALSE,
    row.names = rownames(res_ihh)
  )
  colnames(res_ihs) = c("CHR", "POSITION", "iHS", "-log10(p-value)")
  return(list(iHS = res_ihs, frequency.class = summary_class))
}
