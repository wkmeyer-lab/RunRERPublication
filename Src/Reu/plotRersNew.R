plotRersNew = function (rermat = NULL, index = NULL, phenv = NULL, rers = NULL, method = "k", xlims = NULL, plot = 1, xextend = 0.2, sortrers = F, bgcols = "black", fgcols = "blue", hadjust = 1, vadjust = 0.5, sizeVal = 4.5) {
  {
    if (!is.null(phenv) && length(unique(phenv[!is.na(phenv)])) > 
        2) {
      categorical = TRUE
      if (method != "aov") {
        method = "kw"
      }
    }
    else {
      categorical = FALSE
    }
    if (is.null(rers)) {
      e1 = rermat[index, ][!is.na(rermat[index, ])]
      colids = !is.na(rermat[index, ])
      e1plot <- e1
      if (exists("speciesNames")) {
        names(e1plot) <- speciesNames[names(e1), ]
      }
      if (is.numeric(index)) {
        gen = rownames(rermat)[index]
      }
      else {
        gen = index
      }
    }
    else {
      e1plot = rers
      gen = "rates"
    }
    names(e1plot)[is.na(names(e1plot))] = ""
    if (!is.null(phenv)) {
      phenvid = phenv[colids]
      if (categorical) {
        fgdcor = getAllCor(rermat[index, , drop = F], phenv, 
                           method = method)[[1]]
      }
      else {
        fgdcor = getAllCor(rermat[index, , drop = F], phenv, 
                           method = method)
      }
      plottitle = paste0(gen, ": rho = ", round(fgdcor$Rho, 
                                                4), ", p = ", round(fgdcor$P, 4))
      if (categorical) {
        n = length(unique(phenvid))
        if (n > length(palette())) {
          pal = colorRampPalette(palette())(n)
        }
        else {
          pal = palette()[1:n]
        }
      }
      if (categorical) {
        df <- data.frame(species = names(e1plot), rer = e1plot, 
                         stringsAsFactors = FALSE) %>% mutate(mole = as.factor(phenvid))
      }
      else {
        df <- data.frame(species = names(e1plot), rer = e1plot, 
                         stringsAsFactors = FALSE) %>% mutate(mole = as.factor(ifelse(phenvid > 
                                                                                        0, 2, 1)))
      }
    }
    else {
      plottitle = gen
      df <- data.frame(species = names(e1plot), rer = e1plot, 
                       stringsAsFactors = FALSE) %>% mutate(mole = as.factor(ifelse(0, 
                                                                                    2, 1)))
    }
    if (sortrers) {
      df = filter(df, species != "") %>% arrange(desc(rer))
    }
    if (is.null(xlims)) {
      ll = c(min(df$rer) * 1.1, max(df$rer) + xextend)
    }
    else {
      ll = xlims
    }
  }
  if (categorical) {
    g <- ggplot(df, aes(x = rer, y = factor(species, levels = unique(ifelse(rep(sortrers, 
                                                                                nrow(df)), species[order(rer)], sort(unique(species))))), 
                        col = mole, label = species)) + scale_size_manual(values = c(1, 
                                                                                     1, 1, 1)) + geom_point(aes(size = mole)) + scale_color_manual(values = pal) + 
      scale_x_continuous(limits = ll) + geom_text(hjust = 1, 
                                                  size = 2) + ylab("Branches") + xlab("relative rate") + 
      ggtitle(plottitle) + geom_vline(xintercept = 0, linetype = "dotted") + 
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
            legend.position = "none", panel.background = element_blank(), 
            axis.text = element_text(size = 18, face = "bold", 
                                     colour = "black"), axis.title = element_text(size = 24, 
                                                                                  face = "bold"), plot.title = element_text(size = 24, 
                                                                                                                            face = "bold")) + theme(axis.line = element_line(colour = "black", 
                                                                                                                                                                             size = 1)) + theme(axis.line.y = element_blank())
  }
  else {
    g <- ggplot(df, 
                aes(x = rer, 
                    y = factor(species, levels = unique(ifelse(rep(sortrers, nrow(df)), species[order(rer)], sort(unique(species))))), 
                    col = mole, 
                    label = species)
    ) + 
      scale_size_manual(values = c(1, 1, 1, 1)) + 
      geom_point(aes(size = mole)) + 
      scale_color_manual(values = c(bgcols, fgcols)) + 
      scale_x_continuous(limits = ll) + 
      geom_text(hjust = hadjust,  vjust = vadjust, size = sizeVal) + 
      ylab("Branches") + 
      xlab("relative rate") + 
      ggtitle(plottitle) + 
      geom_vline(xintercept = 0, linetype = "dotted") + 
      theme(
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position = "none", 
        panel.background = element_blank(), 
        axis.text = element_text(size = 18, face = "bold", colour = "black"), 
        axis.title = element_text(size = 24, face = "bold"), 
        plot.title = element_text(size = 24, face = "bold")) + 
      theme(axis.line = element_line(colour = "black", size = 1)) + 
      theme(axis.line.y = element_blank())
  }
  if (plot) {
    print(g)
  }
  else {
    g
  }
}
