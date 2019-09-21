  setwd(getwd())
  
  if(!require("openxlsx")){install.packages("openxlsx");library(openssl)}
  if (!require("ggplot2")) {install.packages("ggplot2"); library(ggplot2)}
  if (!require("reshape2")) {install.packages("reshape2"); library(reshape2)}
  if (!require("plotly")) {install.packages("plotly"); library(plotly)}
  if (!require("autoplotly")) {install.packages("autoplotly"); library(autoplotly)}
  if (!require("gridExtra")) {install.packages("gridExtra"); library(gridExtra)}
  
  #This file is the bin file for each sample that after the fitting
  muscle_data <- read.csv("data_nmr_leucine_snr.csv", header = TRUE)[c(1:7, 14:20),]
  muscle_data$Treatment <- c(rep("Ctrl", 7),
                             rep("LT", 7))
  tx <- factor(muscle_data$Treatment)
  
  
  
  norm_factor <- read.xlsx("normalization_factor.xlsx", sheet = 1, startRow = 1, colNames = TRUE, rowNames = FALSE)
    norm_factor <- norm_factor[c(1:7, 14:20),]
  
  
  muscle <- data.frame(muscle_data[,1:3], muscle_data[,4:ncol(muscle_data)]*norm_factor[,4])
  
  muscle$Biopsy <-factor(muscle$Biopsy)
  bp <- muscle$Biopsy
  
  
  #This is the file that contains after-fitting bin information (center, min, max and width)
  bin_range <- read.csv("bin_range.csv", header = TRUE) 
    bin_range$ppm <- as.numeric(gsub("_", ".", gsub("B", "", bin_range$name)))
  
  #This is the bin file for each sample using the small bucket to make high-resolution graph
  graph_data <- read.csv("uniforme_snr.csv", header = TRUE)
  graph_data <- graph_data[c(1:7, 14:20),]
  
  
  
  #Line Plot for whole spectra UNIFORME bin_tx#####
  graph_data$Treatment <- c(rep("control", 7),
                            rep("leucine+threonine", 7))
  graph_data$Treatment <- factor(graph_data$Treatment)
  
  
  graph_long <- melt(graph_data)
  graph_long$variable <- as.numeric(gsub("_", ".", 
                                         gsub("B", "", 
                                              gsub("[.]", "", 
                                                   graph_long$variable
                                              )
                                         )
  )
  )
  
  Intensity <- graph_long$variable
  ppm <- graph_long$value
  Treatment <- graph_long$Treatment
  
  graph <- ggplot() + 
    geom_line(aes(x = Intensity, 
                  y = ppm, 
                  color = Treatment)) + 
    scale_x_reverse(breaks = round(seq(min(graph_long$variable), 
                                       max(graph_long$variable), by = 0.05),2))+
    labs(title = "Spectra_muscle", x = "ppm", y = "Intensity")
  
  
  
  #output spectra of all treatments
  # htmlwidgets::saveWidget(ggplotly(graph, dynamicTicks = TRUE), "spectra_tx.html")  
  
  #Stat_tx####
  sig <- list()
  i <- NULL
  for (i in 4:ncol(muscle)) {
    ppm <- as.numeric(gsub("_", ".", gsub("B", "", colnames(muscle)[i])))
    Ctrl <- muscle[grep("Ctrl", muscle[, 2]),i]
    LT <- muscle[grep("LT", muscle[,2]),i]
    pval <- data.frame(ppm,
                       ppm_chenomx = ppm - 0.017,
                       p.wilcox = wilcox.test(muscle[[i]] ~ tx)[["p.value"]],
                       mean.all = round(mean(muscle[,i]),2),
                       
                       mean.Ctrl = paste(round(mean(Ctrl),2),
                                         "±", round(sd(Ctrl),2)),
                       mean.LT = paste(round(mean(LT),2),
                                           "±", round(sd(LT),2))
    )
    if (pval$p.wilcox >= 0.1) {
      pval$i <- i
      pval$sig.level = ""
      sig[[i]] <- pval
    }
    if (pval$p.wilcox < 0.1 & pval$p.wilcox >= 0.05) {
      pval$i <- i
      pval$sig.level = "Tendency"
      sig[[i]] <- pval
    }
    
    if (pval$p.wilcox < 0.05 & pval$p.wilcox >= 0.01) {
      pval$i <- i
      pval$sig.level = "Signif"
      sig[[i]] <- pval
    }
    if (pval$p.wilcox < 0.01) {
      pval$i <- i
      pval$sig.level = "Very Signif"
      sig[[i]] <- pval
    }
  }
  
  signif <-  do.call(rbind.data.frame, sig)
  signif <- signif[ , !(names(signif) %in% "i")]
  rownames(signif) <- NULL
  rm(sig)
  
  #Match significant range_tx######
  sig <- list()
  for (i in 1:nrow(signif)){
    for (j in 1:nrow(bin_range)){
      if (signif$ppm[i] == bin_range$ppm[j]){
        pval <- cbind.data.frame(signif[i,],bin_range[j,3:5])
      }
    }
    pval$i <- i
    sig[[i]] <- pval
  }
  significant <-  do.call(rbind.data.frame, sig)
  significant <- significant[ , !(names(significant) %in% "i")]
  significant <- significant[,c(1, 8:10, 2:7)]
  
  rm(sig, pval)
  #Identify the chemical shifts
  
  rm(signif)
  # signif <- read.csv("Assigned_signif.csv", header = TRUE, fileEncoding= "UTF-8-BOM")
  
  
  #Match significant with uniforme data_tx#####
  graph_mean <- cbind.data.frame(mean = colMeans(graph_data[,4:ncol(graph_data)]),
                                 ppm = as.numeric(gsub("_", ".", 
                                                       gsub("B", "", 
                                                            gsub("[.]", "", 
                                                                 colnames(graph_data)[4:ncol(graph_data)])))))
  
  rownames(graph_mean) <- NULL
  i <- NULL
  j <- NULL
  sig <- list()
  for (j in 1:nrow(significant )) {{
    for (i in 1:nrow(graph_mean))
      if (significant$min[j] <= graph_mean$ppm[i] & significant$max[j] >= graph_mean$ppm[i]){
        signif_graph <- graph_mean[i,]
        signif_graph$pval <- significant$p.wilcox[j]
        signif_graph$sig.lev <- significant$sig.lev[j]
        signif_graph$i <- i
        sig[[i]] <- signif_graph
      }     
  }
  }
  signiff <-  do.call(rbind.data.frame, sig)
  missing.value <- c(1:nrow(graph_mean))[!c(1:nrow(graph_mean)) %in% signiff$i]
  i <- NULL
  for (i in missing.value[1: length(missing.value)]) {
    signif_graph <- graph_mean[i,]
    signif_graph$pval <- NA
    signif_graph$sig.lev <- NA
    signif_graph$i <- i
    sig[[i]] <- signif_graph
  }
  
  
  signiff <-  do.call(rbind.data.frame, sig)
  rm(sig)
  
  write.xlsx(significant, "significant_tx.xlsx")


#Make significant_graph_tx####
tendency  <- signiff
tendency[which(is.na(tendency$sig.lev)),]$mean <- NA
tendency[which(tendency$sig.lev != "Tendency"),]$mean <- NA

sigg <- signiff
sigg[which(is.na(sigg$sig.lev)),]$mean <- NA
sigg[which(sigg$sig.lev != "Signif"),]$mean <- NA

very <- signiff
very[which(is.na(very$sig.lev)),]$mean <- NA
very[which(very$sig.lev != "Very Signif"),]$mean <- NA


p <- ggplot()+geom_line(aes(y = signiff$mean, x = signiff$ppm))
pp <- p + 
  geom_line(aes(y = sigg$mean, x = sigg$ppm, color = "Significant")) +
  geom_line(aes(y = very$mean, x = very$ppm, color = "Very Significant")) +
  scale_x_reverse()+
  scale_color_manual(values = c("yellow",  "red"))+
  labs(title="1H NMR Mean Muscle Spectrum With wilcox Significance between Ctrl and LT Highlight",
       x="Chemical Shift", y = "SNR")

htmlwidgets::saveWidget(ggplotly(pp, dynamicTicks = TRUE), "spectra_sig_tx.html")  

