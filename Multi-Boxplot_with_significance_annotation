#Original scientific paper: https://www.frontiersin.org/articles/10.3389/fvets.2019.00123/full

setwd("/home/kuai/Dropbox/2018-19/Newborn parameter/Boxplot")
nb<-read.csv("All_nb.csv", header = TRUE)

animal<-nb$Animal
age<-nb$Age
exp<-nb$Experiment
tx<-nb$Treatment

for (i in 5:ncol(nb)) {
  nb[[i]] <- suppressWarnings(as.numeric(as.character(nb[[i]])))
}
i <-NULL


#########################################
########Boxplot with Tukey HSD###########
#########################################

if (!require("dplyr")) {install.packages("dplyr"); library(dplyr)}
if (!require("multcompView")) {install.packages("multcompView"); library(multcompView)}
if (!require("ggplot2")) {install.packages("ggplot2"); library(ggplot2)}
if (!require("ggpubr")) {install.packages("ggpubr"); library(ggpubr)}
if (!require("gridExtra")) {install.packages("gridExtra"); library(gridExtra)}

names(nb)[5:ncol(nb)] <- c("ALT (U/L)", "AST (U/L)",
               "Cholesterol (mg/dL)", "Creatinine (mg/dL)", "GGT (U/L)", "Glucose (mg/dL)",
               "GPx (U/L)", "Haptoglobin (mg/mL)", "NEFAs (mmol/L)", "TP (g/dL)", "SOD (U/mL)",
               "TGs (mg/dL)", "Urea (mg/dL)")

  ###ggplot
  p <- list()
  posthoc <- list()
  i <- NULL
  for (i in colnames(nb)[5:ncol(nb)]){
    #Create posthoc list (if is necessary)
    tukey <- TukeyHSD(aov((nb[[i]])~age), conf.level = 0.95)$age  #Show TukeyHSD result
     Bonferroni <- p.adjust(tukey[,4], "bonferroni")  #Bonferroni correction
     names(Bonferroni) <- c("24h-2w", "24h-5w", "24h-7w", "2w-5w",  "2w-7w",  "5w-7w")
      bflist <- cbind.data.frame(tukey, Bonferroni)  #Make list combining TukeyHSD and Bonferroni
        names(bflist) <- c("Mean Difference", "Lower Bound", "Upper Bound", 
                           "Tukey p-Value", "Bonferroni p-value Adjusted") #Name list title
    bflist$i <- i  #Give list's corrsponding analyte's name
     posthoc[[i]] <- format(round(bflist[,1:5], digits = 2), nsmall = 2)  # Take analyte's name and drop the names in the list
    #Reorder Bonferroni for boxplot (if is necessary)
    #Create letter annotation
     annotation <- data.frame(multcompLetters(Bonferroni)['Letters'])
      bxpt <- boxplot(nb[[i]] ~ age, range=3, plot = FALSE)  #Create bxpt
       top <- bxpt$stats[nrow(bxpt$stats),]  #find the top places in boxplot
        names(top) <- bxpt$names  #Add age names
         highest <- data.frame(top)
    #Create data.frame that countains age, highest point and annotation
    df <- data.frame(x = rownames(highest), y = highest[,1], label = annotation[,1])
    #Create plot
    ggbox <- ggplot(data=nb) + 
     geom_boxplot(aes(x = age, y = nb[[i]], fill = NULL), coef = 3, size =1, notch = TRUE)+
      geom_text(data = df, aes(x = x, y = y, label = label), 
                size=6, vjust=0.4, hjust=-0.5)+
       ggtitle(i) +
        theme_bw()+
         theme(panel.border = element_rect(colour = "black", fill=NA, size=2),
                plot.title = element_text(hjust = 0.5, family = "Times New Roman", size = 15, face = "bold"),
                 axis.text.x = element_text(family = "Times New Roman", size=14),
                  axis.title.y = element_blank(),
                   axis.text.y = element_text(family = "Times New Roman", size=14),
                    axis.title.x = element_blank())
    p[[i]] = ggplotGrob(ggbox)
  }

####Haptoglobin_rm####
  #give which case (its numeric order in the vector) is inside the fence
  outlier.fence <- function(x, fence){ 
    if(quantile(x, na.rm = TRUE)[4] + #75%quartile
       fence*IQR(x, na.rm = TRUE) >  #3 times IQR
       max(x, na.rm = TRUE)
    ){ #up fence higher than max
      uplim <- max(x, na.rm = TRUE)
    } else{
      uplim <- quantile(x, na.rm = TRUE)[4] + 
        fence*IQR(x, na.rm = TRUE)
    }
    
    if(quantile(x, na.rm = TRUE)[2] - #25%quartile
       fence*IQR(x, na.rm = TRUE) <  #3 times IQR
       min(x, na.rm = TRUE)
    ){ #up fence higher than max
      lwlim <- min(x, na.rm = TRUE)
    } else{
      lwlim <- quantile(x, na.rm = TRUE)[2] - 
        fence*IQR(x, na.rm = TRUE)
    }
    print(c(lwlim, uplim))
  }  

  
  ##plot
  tukey <- TukeyHSD(aov((nb$"Haptoglobin (mg/mL)")~age), conf.level = 0.95)$age  #Show TukeyHSD result
  Bonferroni <- p.adjust(tukey[,4], "bonferroni")  #Bonferroni correction
  names(Bonferroni) <- c("24h-2w", "24h-5w", "24h-7w", "2w-5w",  "2w-7w",  "5w-7w")
 
  #Create letter annotation
  annotation <- data.frame(multcompLetters(Bonferroni)['Letters'])
  bxpt <- boxplot((nb$"Haptoglobin (mg/mL)") ~ age, range=3, plot = FALSE)  #Create bxpt
  top <- bxpt$stats[nrow(bxpt$stats),]  #find the top places in boxplot
  names(top) <- bxpt$names  #Add age names
  highest <- data.frame(top)
  #Create data.frame that countains age, highest point and annotation
  df <- data.frame(x = rownames(highest), y = highest[,1], label = annotation[,1])
  #Create plot
  ggbox <- ggplot(data=nb) + 
    geom_boxplot(aes(x = age, y = (nb$"Haptoglobin (mg/mL)"), fill = NULL), coef = 3, size =1, notch = TRUE, outlier.shape = NA)+
    geom_text(data = df, aes(x = x, y = y, label = label), 
              size=6, vjust=0.4, hjust=-0.5)+
    ylim(outlier.fence(nb$`Haptoglobin (mg/mL)`, 6))+
    ggtitle("Haptoglobin (mg/mL) Without Outliers") +
    theme_bw()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=2),
          plot.title = element_text(hjust = 0.5, family = "Times New Roman", size = 15, face = "bold"),
          axis.text.x = element_text(family = "Times New Roman", size=14),
          axis.title.y = element_blank(),
          axis.text.y = element_text(family = "Times New Roman", size=14),
          axis.title.x = element_blank())
  p[["Haptoglobin (mg/mL) Without Outliers"]] = ggplotGrob(ggbox)
  
  ####ggt####
  i = colnames(nb)[9]
    #Create posthoc list (if is necessary)
    tukey <- TukeyHSD(aov((nb[[i]])~age), conf.level = 0.95)$age  #Show TukeyHSD result
    Bonferroni <- p.adjust(tukey[,4], "bonferroni")  #Bonferroni correction
    names(Bonferroni) <- c("24h-2w", "24h-5w", "24h-7w", "2w-5w",  "2w-7w",  "5w-7w")
    bflist <- cbind.data.frame(tukey, Bonferroni)  #Make list combining TukeyHSD and Bonferroni
    names(bflist) <- c("Mean Difference", "Lower Bound", "Upper Bound", 
                       "Tukey p-Value", "Bonferroni p-value Adjusted") #Name list title
    bflist$i <- i  #Give list's corrsponding analyte's name
    posthoc[[i]] <- format(round(bflist[,1:5], digits = 2), nsmall = 2)  # Take analyte's name and drop the names in the list
    #Reorder Bonferroni for boxplot (if is necessary)
    #Create letter annotation
    annotation <- data.frame(multcompLetters(Bonferroni)['Letters'])
    bxpt <- boxplot(nb[[i]] ~ age, range=3, plot = FALSE)  #Create bxpt
    top <- bxpt$stats[nrow(bxpt$stats),]  #find the top places in boxplot
    names(top) <- bxpt$names  #Add age names
    highest <- data.frame(top)
    #Create data.frame that countains age, highest point and annotation
    df <- data.frame(x = rownames(highest), y = highest[,1], label = annotation[,1])
    df$y <- c(3770, 348, 250, 250)
    #Create plot
    ggbox <- ggplot(data=nb) + 
      geom_boxplot(aes(x = age, y = nb[[i]], fill = NULL), coef = 3, size =1, notch = TRUE)+
      geom_text(data = df, aes(x = x, y = y, label = label), 
                size=6, vjust=0.4, hjust=-0.5)+
      ggtitle(i) +
      theme_bw()+
      theme(panel.border = element_rect(colour = "black", fill=NA, size=2),
            plot.title = element_text(hjust = 0.5, family = "Times New Roman", size = 15, face = "bold"),
            axis.text.x = element_text(family = "Times New Roman", size=14),
            axis.title.y = element_blank(),
            axis.text.y = element_text(family = "Times New Roman", size=14),
            axis.title.x = element_blank())
    p[[i]] = ggplotGrob(ggbox)
  
  
  
  ####ggt_facet######
  
  bxpt_ggt_facet <- ggplot(data=nb, aes(x = age, y = nb$`GGT (U/L)`)) + 
    geom_boxplot(aes(x = age, y = nb$`GGT (U/L)`, fill = NULL), coef = 3, size =1, notch = TRUE)+
    facet_wrap(age, scales = "free", ncol = 4)+
    ggtitle("GGT (U/L) In Different Scales") +
    theme_bw()+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=2),
          plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size=14),
          axis.title.x = element_blank())
  
  #plot
  p[["GGT (U/L) facet"]] = ggplotGrob(bxpt_ggt_facet)
  
 # names(p)
  
  #[1] "ALT (U/L)"                        2    "AST (U/L)"                   3         "Cholesterol (mg/dL)"                 
  #[4] "Creatinine (mg/dL)"               5    "GGT (U/L)"                   6         "Glucose (mg/dL)"                     
  #[7] "GPx (U/L)"                       8     "Haptoglobin (mg/mL)"            9      "NEFAs (mmol/L)"                      
  #[10] "TP (g/dL)"                       11     "SOD (U/mL)"                   12        "TGs (mg/dL)"                         
  #[13] "Urea (mg/dL)"                 14        "Haptoglobin (mg/mL) Without Outliers" 15  "GGT (U/L) facet"        
  
  
#Save multiple plots
main <- grid.arrange(grobs= p[c(  1, 2, 7,
                                  5, 15, 11,
                                  9, 12, 3,
                                  6, 4, 10,
                                  13, 8, 14)] , ncol=3)
                               
                     #           p$`ALT (U/L)`,       p$`AST (U/L)`,           p$`GPx (U/L)`,
                      #          p$`GGT (U/L)`,       p$`GGT (U/L) facet`,     p$`SOD (U/mL)`,
                       #         p$`NEFAs (mmol/L)`,  p$`TGs (mg/dL)`,         p$`Cholesterol (mg/dL)`,
                        #        p$`Glucose (mg/dL)`, p$`Creatinine (mg/dL)`,  p$`TP (g/dL)`,
                         #       p$`Urea (mg/dL)`,    p$`Haptoglobin (mg/mL)`, p$`Haptoglobin (mg/mL) Without Outliers`
  

 

dev.off()

 ggsave("bxpt_calf.jpg", plot = main,
       height = (length(p)/3)*9,
        width = 15*3,
         dpi = 300, units = "cm", limitsize = FALSE)



 
                   ####Table####
if (!require("tableHTML")) {install.packages("tableHTML"); library(tableHTML)}
###Create posthoc longlist
DF <- do.call(rbind.data.frame, posthoc)
 Comparison <- rep(c("2w-24h","5w-24h","7w-24h","5w-2w","7w-2w","7w-5w"),length(posthoc))
  tukey_tab <- cbind.data.frame(Comparison,DF)
   rownames(tukey_tab) <- NULL
hsd_calf_tab <- tableHTML(tukey_tab, 
                          widths = rep(140, 7), theme = "default", rownames = FALSE,
                          row_groups = list(rep(6,13), 
                                            c("ALT (U/L)", "AST (U/L)",
                                              "Cholesterol (mg/dL)", "Creatinine (mg/dL)", "GGT (U/L)", "Glucose (mg/dL)",
                                              "GPx (U/L)", "Haptoglobin (mg/mL)", "NEFAs (mmol/L)", "TP (g/dL)", "SOD (U/mL)",
                                              "TGs (mg/dL)", "Urea (mg/dL)"
                                            )),
                          caption = 'tukeyHSD posthoc test with Bonferroni correction-calf') %>%
  add_css_table(css = list('text-align', 'center')) 
  
write_tableHTML(hsd_calf_tab, "C:/Users/kuaiy/Dropbox/2018-19/Newborn parameter/manuscript_nb/hsd_calf_tab.HTML")
