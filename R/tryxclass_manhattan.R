Tryx$set("public", "manhattan_plot", function(dat = self$output$candidate_outcome_mr, id_remove=NULL, y_scale=NULL, lable = TRUE){
  
  cpg <- require(ggplot2)
  if(!cpg)
  {
    stop("Please install the ggplot2 package")
  }
  cpg <- require(ggrepel)
  if(!cpg)
  {
    stop("Please install the ggrepel package")
  }
  
  #Open & clean data
  #mr outcome: candidate traits-outcome / candidate traits-exposure / exposure-candidate traits
  temp <- subset(dat, !id.exposure %in% id_remove)
  
  #X-axis preparing: 
  #Each trait needs to be shown in different phenotype groups, along the X axis 
  #Sort the traits based on their p-value for MR results
  temp <- temp[order(temp$pval_adj),]
  
  #Numbering row according to phenotype group where numeric variable is required.
  temp <- temp %>% mutate(id = row_number())
  
  #To compute the (mimic) cumulative position for the traits
  #pos = temp %>% 
  #      #group_by(cat) %>% 
  #      summarize(center=( max(id) + min(id) ) / 2 )
  
  #
  name <- tidyr::separate(temp, exposure, sep="\\|", c("exposure", "blank", "id"))
  temp <- subset(name, select= -c(blank))
  
  #manhattan plot
  p <- ggplot(temp, aes(x=id, y=as.numeric(-log10(pval_adj)*sign(b)))) +
    geom_point(alpha=0.8, size=3) +
    geom_point(colour = "snow", size = 1.5) +
    
    # custom X axis
    # scale_x_discrete("Phenotype", breaks=gd$center, labels = c(gd$cat))+
    #scale_x_continuous(breaks=pos$center) +
    scale_y_continuous(limits = y_scale) +
    # scale_y_continuous(expand = c(0, 0) ) +   
    
    # Add highlighted points
    geom_point(data=subset(temp, sig), size=3) +
    
    # Add line
    geom_hline(yintercept=0,colour="black", alpha=I(2/3)) +

    # Custom the titles
    #ggtitle("Effect of the candidate traits on the outcome") +
    # xlab("cat") + 
    ylab("-Log10(Adjusted P-value) x sign(beta)") + 
    xlab(paste("Association of the candidate traits and ", temp$outcome[1], sep=""))+
    
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x=element_text(hjust = 0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
      
    )     
  
  # Add label using ggrepel to avoid overlapping
  if(label)
  {
    p + geom_label_repel(data=subset(temp, sig), aes(label=exposure), size=2)
  }
  
  self$output$plots$manhattan_plot <- p
  invisible(self$output)
}
)


