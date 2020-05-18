# A function to extract appropriate fits to a data-frame and save results
extract_save_fits <- function(FITS, DAT, model_selected, filename)
{
  TaxaName = numeric(0)
  Class = numeric(0)
  Product = numeric(0)
  Unit = numeric(0)
  Poisson25 = numeric(0)
  PoissonMed = numeric(0)
  Poisson975 = numeric(0)
  PoissonSignif = numeric(0)
  MeanLegalLevel = numeric(0)
  MeanIllegalLevel = numeric(0)
  MaxLegalLevel = numeric(0)
  MaxIllegalLevel =numeric(0)
  SumLegal = numeric(0)
  SumIllegal = numeric(0)
  CitesListing = numeric(0)
  
  for (ii in 1:length(DAT))
  {
      TaxaName = c(TaxaName, DAT[[ii]]$taxa)
      Class = c(Class, DAT[[ii]]$class)
      Product = c(Product, DAT[[ii]]$terms)
      Unit = c(Unit, DAT[[ii]]$unit)
      MeanIllegalLevel = c(MeanIllegalLevel, mean(DAT[[ii]]$illegal_original))
      MeanLegalLevel = c(MeanLegalLevel, mean(DAT[[ii]]$legal_original))
      MaxIllegalLevel = c(MaxIllegalLevel, max(DAT[[ii]]$illegal_original))
      MaxLegalLevel = c(MaxLegalLevel, max(DAT[[ii]]$legal_original))
      SumLegal = c(SumLegal, sum(DAT[[ii]]$legal_original))
      SumIllegal = c(SumIllegal, sum(DAT[[ii]]$illegal_original))
      
      Poisson25 = c(Poisson25, quantile(extract(FITS$nb)$slope[,ii], 0.025))
      PoissonMed = c(PoissonMed, mean(extract(FITS$nb)$slope[,ii]))
      Poisson975 = c(Poisson975, quantile(extract(FITS$nb)$slope[,ii], 0.975))
      CitesListing = c(CitesListing, get_cites_appendix(DAT[[ii]]))
      
      if (Poisson25[length(Poisson25)] < 0)
      {
        if (Poisson975[length(Poisson975)] < 0)
        {
          PoissonSignif = c(PoissonSignif, TRUE)
        } else
        {
          PoissonSignif = c(PoissonSignif, FALSE)
        }
        
      } else
      {
        PoissonSignif = c(PoissonSignif, TRUE)
      }
  }

  new_Data_frame = data.frame(
    TaxaName = TaxaName,
    Class = Class,
    Product = Product,
    Unit = Unit,
    MeanLegalLevel = MeanLegalLevel,
    MeanIllegalLevel = MeanIllegalLevel,
    MaxLegalLevel = MaxLegalLevel,
    MaxIllegalLevel = MaxIllegalLevel,
    SumLegal = SumLegal,
    SumIllegal = SumIllegal,
    CitesListing = CitesListing,
    Poisson25 = Poisson25,
    PoissonMed = PoissonMed,
    Poisson975 = Poisson975,
    PoissonSignif = PoissonSignif
  )
  
  write.csv(new_Data_frame, file = filename)
  
  return(new_Data_frame)
}

# Concatenate text in multiple colours
concat.text<-function(x,y,txt,col) {
  thisx<-x
  for(txtstr in 1:length(txt)) {
    text(thisx,y,txt[txtstr],col=col[txtstr], cex = 0.4, pos = 4)
    thisx<-thisx+strwidth(txt[txtstr], cex = 0.4)
  }
}


# Make caterpillar plot of results
plot_caterpillar <- function(to_plot, models_to_include, overall_slope, filename, colour_bars = FALSE, show_appendix = FALSE, show_trade_levels = FALSE)
{
  
  dev.new()
  pdf(paste('c:/users/derekt/work/research/dropbox/legalillegalanalogues/', filename, sep=""), width = 8.5, height = 11, family = "ArialMT", useDingbats = FALSE)
  par(mar = c(4,4,2,4))
  names = as.character(to_plot$TaxaName)
  products = as.character(to_plot$Product)
  classes = as.character(to_plot$Class)
  units = as.character(to_plot$Unit)
  pt25 = to_plot$Poisson25
  pt975 = to_plot$Poisson975
  meds = to_plot$PoissonMed
  MeanIllegalLevel = to_plot$MeanIllegalLevel
  MeanLegalLevel = to_plot$MeanLegalLevel
  MaxIllegalLevel = to_plot$MaxIllegalLevel
  MaxLegalLevel = to_plot$MaxLegalLevel
  SumLegal = to_plot$SumLegal
  SumIllegal = to_plot$SumIllegal
  CitesListing = as.character(to_plot$CitesListing)
  minp = min(pt25)
  maxp = max(pt975)
  color_list = list()
  
  if (!colour_bars)
  {
    for (ii in 1:length(meds))
    {
      alpha_level = 1
      color_list[ii] = rgb(0,0,0,alpha_level)     
    }

  } else
  {
  for (ii in 1:length(meds))
  {
    alpha_level = 1
    if (classes[ii] == "AVES")
    {
      color_list[ii] = rgb(0.8,0,0,alpha_level)
    } else if (classes[ii] == "MAMMALIA")
    {
      color_list[ii] = rgb(0,0,0.8,alpha_level)
    } else if (classes[ii] == "REPTILIA")
    {
      color_list[ii] = rgb(0,0.8,0,alpha_level)
    } else
    {
      color_list[ii] = rgb(0,0,0,alpha_level)
    }
  }
  }

  plot(-500,-500, xlim = c(-4, 1.5), ylim = c(-1, length(meds)),axes = F, xlab = "Slope value", ylab = "", cex = 0.8)
  axis(side = 1, at = c(-1.5,-1,-0.5,0,0.5,1,1.5),labels = c(-1.5,-1,-0.5,0,0.5,1,1.5), xlab = "Slope value")
  #box()
  for (ii in 1:length(models_to_include))
  {
    print(names[models_to_include[ii]])
    print(meds[models_to_include[ii]])
    if (names[models_to_include[ii]] != "Acipenseriformes")
    {
      points(meds[models_to_include[ii]], ii, pch = ifelse(to_plot$PoissonSignif[models_to_include[ii]] == TRUE,ifelse(CitesListing[models_to_include[ii]] == "I",17,16), ifelse(CitesListing[models_to_include[ii]] == "I",2,1)), 
           cex = ifelse(to_plot$PoissonSignif[models_to_include[ii]] == TRUE,0.8, 0.5), col = ifelse(CitesListing[models_to_include[ii]] == "I", "red", "black"))#unlist(color_list[[ii]]))
      #axis(side = 1)
      lines(c(pt25[models_to_include[ii]],pt975[models_to_include[ii]]),c(ii,ii), col = ifelse(CitesListing[models_to_include[ii]] == "I", "red", "black"))#unlist(color_list[[ii]]))
    } else
    {
      points(meds[models_to_include[ii]], ii, pch = ifelse(to_plot$PoissonSignif[models_to_include[ii]] == TRUE,ifelse(CitesListing[models_to_include[ii]] == "I",17,17), ifelse(CitesListing[models_to_include[ii]] == "I",2,2)), 
             cex = ifelse(to_plot$PoissonSignif[models_to_include[ii]] == TRUE,0.8, 0.5), col = ifelse(CitesListing[models_to_include[ii]] == "I", "red", "black"))#unlist(color_list[[ii]]))
      #axis(side = 1)
      lines(c(pt25[models_to_include[ii]],pt975[models_to_include[ii]]),c(ii,ii), col = ifelse(CitesListing[models_to_include[ii]] == "I", "red", "black"))#unlist(color_list[[ii]]))
      
    }
    if (show_trade_levels)
    {
      ctext <- c(names[models_to_include[ii]], " ", products[models_to_include[ii]], " ", units[models_to_include[ii]], "   ", "Means: ", round(MeanLegalLevel[models_to_include[ii]], 1), 
                  " ", round(MeanIllegalLevel[models_to_include[ii]], 1), "   Maximums: ", round(MaxLegalLevel[models_to_include[ii]],0), " ", 
                 round(MaxIllegalLevel[models_to_include[ii]],0), " L/I ratio: ", round(SumLegal[models_to_include[ii]]/SumIllegal[models_to_include[ii]], 1))
       concat.text(-4, ii, ctext, col = c("black", "black", "black", "black", "black", "black", "black", "blue", "black", "red", "black", "blue", "black", "red", "black", "black"))
      # 
    } else
    {
      
      
      if (show_appendix)
      {
        if (units[models_to_include[ii]] == "")
        {
          text(-4,ii,bquote(paste(italic(.(names[models_to_include[ii]])), " ", " (", .(CitesListing[models_to_include[ii]]), ")", " (",.(products[models_to_include[ii]]), ") ", sep="")), cex = 0.8, pos = 4)
        } else
        {
          if (units[models_to_include[ii]] == "KIL")
          {
            text(-4,ii,bquote(paste(italic(.(names[models_to_include[ii]])), " ", " (", .(CitesListing[models_to_include[ii]]), ")", " (",.(products[models_to_include[ii]]), " kg", ") ", sep="")), cex = 0.8, pos = 4)
          } else
          {
            
            text(-4,ii,bquote(paste(italic(.(names[models_to_include[ii]])), " ", " (", .(CitesListing[models_to_include[ii]]), ")", " (",.(products[models_to_include[ii]]), ") ", sep="")), cex = 0.8, pos = 4)
          }
        }
    #    text(-4,ii,paste(names[models_to_include[ii]], products[models_to_include[ii]], units[models_to_include[ii]], "CITES listing: ", CitesListing[models_to_include[ii]], sep =" "), cex = 0.4, pos =4)
      } else
      {
        if (units[models_to_include[ii]] == "")
        {
      #    text(-4,ii,bquote(paste(italic(.(names[models_to_include[ii]])), " (",.(products[models_to_include[ii]]), ") ", sep="")), cex = 0.8, pos = 4)
          text(-4,ii + 0.2,bquote(paste(italic(.(names[models_to_include[ii]])), sep="")), cex = 0.8, pos = 4)
          text(-4,ii-0.2,bquote(paste( "   (",.(products[models_to_include[ii]]), ") ", sep="")), cex = 0.8, pos = 4)
        } else
        {
          if (units[models_to_include[ii]] == "KIL")
          {
         #   text(-6,ii,bquote(paste(italic(.(names[models_to_include[ii]])), " (",.(products[models_to_include[ii]]), " kg", ") ", sep="")), cex = 0.8, pos = 4)
            text(-4,ii + 0.2,bquote(paste(italic(.(names[models_to_include[ii]])),  sep="")), cex = 0.8, pos = 4)
            text(-4,ii-0.2,bquote(paste("   (",.(products[models_to_include[ii]]), " kg", ") ", sep="")), cex = 0.8, pos = 4)
          } else
          {
            #text(-6,ii,bquote(paste(italic(.(names[models_to_include[ii]])), " (",.(products[models_to_include[ii]]), ") ", sep="")), cex = 0.8, pos = 4)
            text(-4,ii + 0.2,bquote(paste(italic(.(names[models_to_include[ii]])),  sep="")), cex = 0.8, pos = 4)
            text(-4,ii-0.2,bquote(paste("   (",.(products[models_to_include[ii]]), ") ", sep="")), cex = 0.8, pos = 4)
          }
        }
      }
    }
  }
  
  # Plot meta-analytic mean
  meta_analytic_list = list(meta_analytic_min = quantile(overall_slope, 0.025), 
                            meta_analytic_med = quantile(overall_slope, 0.5), 
                            meta_analytic_max = quantile(overall_slope, 0.975))
  if (meta_analytic_list$meta_analytic_min < 0)
  {
    if (meta_analytic_list$meta_analytic_max < 0)
    {
      TempSig = TRUE
    } else
    {
      TempSig = FALSE
    }
    
  } else
  {
    TempSig = TRUE
  }
  
  points(meta_analytic_list$meta_analytic_med, -1, pch = ifelse(TempSig == TRUE,16, 0.8),
         cex = ifelse(TempSig == TRUE,1, 0.5))
  lines(c(meta_analytic_list$meta_analytic_min, meta_analytic_list$meta_analytic_max), c(-1,-1))
  text(-4,-1,"Meta-analytic median", cex = 0.8, pos =4)
  
  abline(v = 0, lty = 2)
  
  

  
  dev.off()
}
  #-------------------------------------------------------------------------------------
  
  
  
  
  
  
  # PLOT THE RATIOS OF LEGAL TO ILLEGAL DATA
  
  # Make caterpillar plot of results
  plot_caterpillar_ratio <- function(DAT, filename)
  {
    
    dev.new()
    pdf(paste('c:/users/derekt/work/research/dropbox/legalillegalanalogues/', filename, sep=""), width = 8.5, height = 11, family="ArialMT")
    par(mar = c(4,4,2,4))
    
    TaxaName = numeric(0)
    Class = numeric(0)
    Product = numeric(0)
    Unit = numeric(0)
    SumLegal = numeric(0)
    SumIllegal = numeric(0)
    Appendix = numeric(0)
    
    # Extract the important data
    for (ii in 1:length(DAT))
    {
      TaxaName = c(TaxaName, DAT[[ii]]$taxa)
      Class = c(Class, DAT[[ii]]$class)
      Product = c(Product, DAT[[ii]]$terms)
      Unit = c(Unit, DAT[[ii]]$unit)
      SumLegal = c(SumLegal, sum(DAT[[ii]]$legal_original_non_robust))
      SumIllegal = c(SumIllegal, sum(DAT[[ii]]$illegal_original))
      Appendix = c(Appendix, get_cites_appendix(DAT[[ii]]))
    }
    
    Ratio = SumIllegal / SumLegal * 100
    
    first_sort_order = sort(Ratio, index.return = T)
    Ratio = Ratio[first_sort_order$ix]
    TaxaName = TaxaName[first_sort_order$ix]
    Class = Class[first_sort_order$ix]
    Product = Product[first_sort_order$ix]
    Unit = Unit[first_sort_order$ix]
    SumLegal = SumLegal[first_sort_order$ix]
    SumIllegal = SumIllegal[first_sort_order$ix]
    Appendix = Appendix[first_sort_order$ix]
    
    color_list = list()
    pch_list = list()
    
      for (ii in 1:length(Ratio))
      {
        alpha_level = 1
        color_list[ii] = ifelse(Appendix[ii] == "I", "red", "black")
        pch_list[ii] = ifelse(Appendix[ii] == "I", 17, 19)
      }

    
    plot(-500,-500,  xlim = c(0.000001, 1000), log = "x", ylim = c(-1, length(Ratio)),
              axes = F, xlab = "Ratio of seizure volume to legal trade volume", ylab = "", cex = 0.8)
    axis(side = 1, at = c(0.01, 0.1, 1, 10, 100, 1000, 10000, 100000), labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100", "1000"), xlab = "Ratio of seizures to legal trade")
    for (ii in 1:length(Ratio))
    {
      if (TaxaName[ii] != "Acipenseriformes")
      points(Ratio[ii], ii, pch = unlist(pch_list[ii]), 
             cex = 1, col = unlist(color_list[[ii]]))
      else
        points(Ratio[ii], ii, pch = 17, 
               cex = 1, col = 'black')

        if (Unit[ii] == "")
        {
         text(0.000001,ii + 0.2,bquote(paste(italic(.(TaxaName[ii])), sep="")), cex = 0.7, pos = 4)
          text(0.000001,ii - 0.2,bquote(paste("   (",.(format(round(SumLegal[ii] + SumIllegal[ii]), big.mark = ",",scientific = FALSE)), " ", .(Product[ii]),") ", sep="")), cex = 0.7, pos = 4)
        
          } else
        {
          if (Unit[ii] == "KIL")
          {
            # text(0.000001,ii,bquote(paste(italic(.(TaxaName[ii])), " (", .(format((SumLegal[ii] + SumIllegal[ii]), big.mark = ",",scientific = FALSE)), " kg ",.(Product[ii]),")",  sep="")), cex = 0.7, pos = 4)
            text(0.000001,ii + 0.2,bquote(paste(italic(.(TaxaName[ii])), sep="")), cex = 0.7, pos = 4)
            text(0.000001,ii - 0.2,bquote(paste("   (", .(format(round(SumLegal[ii] + SumIllegal[ii]), big.mark = ",",scientific = FALSE)), " kg ",.(Product[ii]),")",  sep="")), cex = 0.7, pos = 4)
          } else
          {
            
          #  text(0.000001,ii,bquote(paste(italic(.(TaxaName[ii])), " (", .(format((SumLegal[ii] + SumIllegal[ii]), big.mark = ",",scientific = FALSE)), " ",.(Product[ii]),  " (", .(Unit[ii]),")", sep="")), cex = 0.7, pos = 4)
            text(0.000001,ii + 0.2,bquote(paste(italic(.(TaxaName[ii])),  sep="")), cex = 0.7, pos = 4)
            text(0.000001,ii - 0.2,bquote(paste("   (", .(format(round(SumLegal[ii] + SumIllegal[ii]), big.mark = ",",scientific = FALSE)), " ",.(Product[ii]),  " (", .(Unit[ii]),")", sep="")), cex = 0.7, pos = 4)
          }
        }

    }
    
    #abline(v = exp(1), lty = 2, untf =FALSE)
    #abline(v =exp(5), lty = 3, untf = FALSE)
    lines(c(100,100),c(-10,50), lty = 2)
    print(mean(Ratio))
    
    
    points(mean(Ratio), 0, pch = 16)
    text(0.000001,0,"Mean", cex = 0.7, pos =4)
  
    #box()
    
    dev.off()
  }