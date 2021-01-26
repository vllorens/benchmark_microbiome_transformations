
source("R/corrmod.R")



plot_comparisons <- function(dunn, title=NULL){
    
    dunn_effsize <- dunn %>% mutate(n=n1+n2) %>% mutate(effsize=statistic/sqrt(n)) %>% dplyr::select(group1, group2, effsize) %>% pivot_wider(names_from = group2, values_from = effsize)
    dunn_adj <- dunn %>% mutate(n=n1+n2) %>% mutate(effsize=statistic/sqrt(n)) %>% dplyr::select(group1, group2, p.adj) %>% pivot_wider(names_from = group2, values_from = p.adj)
    
    dunn_effsize <- as.data.frame(dunn_effsize)
    rownames(dunn_effsize) <- dunn_effsize[,1]
    dunn_effsize <- dunn_effsize[,-1]
    dunn_effsize <- as.matrix(dunn_effsize)
    
    dunn_adj <- as.data.frame(dunn_adj)
    rownames(dunn_adj) <- dunn_adj[,1]
    dunn_adj <- dunn_adj[,-1]
    dunn_adj <- as.matrix(dunn_adj)
    
    dunn_effsize <- cbind(Seq=NA, dunn_effsize)
    dunn_effsize <- rbind(dunn_effsize, ACS=NA)
    dunn_adj <- cbind(Seq=NA, dunn_adj)
    dunn_adj <- rbind(dunn_adj, ACS=NA)

    
    ggcorrmod(dunn_effsize, method = "square", type = "upper", 
              colors = c("#354B99", "white", "#A50026"),
            #  colors = c("#6D9EC1", "white", "#E46726"), 
              p.mat=dunn_adj, 
                           legend.title = "Effect size", outline.color = "black",
                           insig="blank",  ggtheme = theme_void(), title=title) + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
                                                           vjust = 0, hjust = 0))+ 
        scale_x_discrete(position="top") + 
        theme(axis.title = element_text(face = "bold"), 
              plot.title = element_text(size = 10, 
                                        face = "bold"), legend.title = element_text(face = "bold"))

}
