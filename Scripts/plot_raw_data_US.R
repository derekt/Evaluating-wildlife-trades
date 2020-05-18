# Print out the raw and transformed data
pdf(
  'c:/users/derekt/work/research/dropbox/legalillegalanalogues/raw_data_US.pdf', family = "ArialMT"
)

par(mfrow = c(4, 2), mar = c(3, 3, 2, 3))

for (i in 1:length(DAT)) {
  N <- DAT[[i]]$N
  min_year = 2014 - N
  years = (min_year + 1):2014
  legal <- DAT[[i]]$legal
  illegal <- DAT[[i]]$illegal
  illegal_original <- DAT[[i]]$illegal_original
  legal_original = DAT[[i]]$legal_original
  legal_original_non_robust = DAT[[i]]$legal_original_non_robust
  rightplot_ylimit = c(min(illegal_original), max(illegal_original))
  
  # Plot the legal and illegal data
  plot(years, legal_original, type = 'p', pch = 16,  col = "grey", ylim = c(min(c(legal_original_non_robust)), max(c(legal_original_non_robust))))
  lines(years, legal_original, col = "grey")
  points(years, legal_original_non_robust, pch = 16, cex = 0.8)
  #  points(years, legal_original_non_robust, pch = 16, col = "grey")
  
  if (i == 1)
    legend(
      "topright",
      legend = c("Legal volume (robust)", "Legal volume (non-robust)"),
      col = c("grey", "black"),
      lty = c(1, NA),
      pch  = c(16, 16)
    )
  
  mtext(bquote(paste(italic(.(DAT[[i]]$taxa)))), line = 1.03, cex = 0.7)
  if (DAT[[i]]$unit == "KIL")
    mtext(bquote(paste("CITES Appendix: ", .(get_cites_appendix(DAT[[i]])), " (",.(DAT[[i]]$terms), " [kg])", sep="")), cex = 0.7, line = -0.1)
  else
    mtext(bquote(paste("CITES Appendix: ", .(get_cites_appendix(DAT[[i]])), " (",.(DAT[[i]]$terms), .(DAT[[i]]$unit), ")", sep="")), cex = 0.7, line = -0.1)
  
  
  plot(years, 
       illegal,
       col = 'orange',
       type = 'p',
       pch = 16,
       ylab = 'n',
       xlab = 'n', ylim = rightplot_ylimit
  )
  lines(years, illegal, col = 'orange')
  points(years, illegal_original, col ="red", pch = 16, cex = 0.8)
  mtext(bquote(paste(italic(.(DAT[[i]]$taxa)))), line = 1.03, cex = 0.7)
  if (DAT[[i]]$unit == "KIL")
    mtext(bquote(paste("CITES Appendix: ", .(get_cites_appendix(DAT[[i]])), " (",.(DAT[[i]]$terms), " [kg])", sep="")), cex = 0.7, line = -0.1)
  else
    mtext(bquote(paste("CITES Appendix: ", .(get_cites_appendix(DAT[[i]])), " (",.(DAT[[i]]$terms), .(DAT[[i]]$unit), ")", sep="")), cex = 0.7, line = -0.1)
  
  
  if (i == 1)
    legend(
      "topright",
      legend = c("Seizures (robust)", "Seizures (non-robust)"),
      col = c("orange", "red"),
      lty = c(1, NA),
      pch  = c(16, 16)
    )
  
}
dev.off()