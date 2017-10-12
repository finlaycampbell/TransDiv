##====== PLOTTING FUNCTIONS =====##

## Plot a discrete violin plot of transmission divergence distributions
vis.gensig <- function(o.store, p.store, alpha = 0.5) {
  
  mk.df <- function(store) {
    
    store$gensig$disease <- factor(store$gensig$disease,
                                   levels = sort.gensig(store))
    
    ## Calculate proportions manually
    df <- by(store$gensig$gensig,
             store$gensig$disease,
             function(sig) prop.table(table(sig)))
    df <- plyr::ldply(df, data.frame)
                                        #df$sig <- as.numeric(df$sig)
    df$sig <- factor(df$sig, levels = as.character(0:max(as.numeric(df$sig))))

    df$.id <- factor(df$.id, levels = rev(sort.gensig(store)))
    df$Freq <- df$Freq/2
    
    return(df)

  }

  o.df <- mk.df(o.store)
  o.df$mod <- "ob"

  p.df <- mk.df(p.store)
  p.df$mod <- 'phyb'

  df <- rbind(o.df, p.df)
  df$fac <- factor(paste0(df$.id, ".", df$mod),
                   levels = paste0(rep(rev(sort.gensig(o.store)), each = 2),
                                   rep(c(".ob", ".phyb"), length(sort.gensig(o.store)))))

  df2 <- df
  df2$Freq <- -df2$Freq

  ggplot(df, aes(x = sig, y = Freq, group = .id, fill = mod)) +
    geom_bar(stat = "identity", position = 'dodge', alpha = alpha) +
    geom_bar(data = df2,
             mapping = aes(x = sig, y = Freq, group = .id, fill = mod),
             stat = 'identity', position = 'dodge', alpha = alpha) +
    facet_wrap( ~ .id + mod,
               ncol = 2*length(unique(df$fac)),
               strip.position = 'bottom',
               labeller = label_wrap_gen(multi_line = FALSE)) +
    xlim(as.character(0:15)) +
    coord_flip() +
    scale_fill_manual(name = 'Model',
                      labels = c('outbreaker2', 'phybreak'),
                      values = c('#7fc97f', '#beaed4')) +
  labs(x = "Transmission divergence", subtitle = 'A') +
    theme_minimal(base_size = 24) +
    theme(#legend.position = 'none',
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_text(colour = 'white'),
      panel.spacing = unit(-0.45, "lines"),
      legend.direction = 'horizontal',
      plot.subtitle = element_text(size = 35))

}

## Violin plot of prop for each pathogen for o and p
vis.uniq <- function(o.store, p.store, alpha = 0.2) {

  o.store$acc$uniq <- o.store$acc$uniq
  o.store$acc$disease <- factor(o.store$acc$disease, levels = rev(sort.gensig(o.store)))
  o.store$acc$mod <- 'ob'
  
  p.store$acc$uniq <- p.store$acc$uniq
  p.store$acc$disease <- factor(p.store$acc$disease, levels = rev(sort.gensig(p.store)))
  p.store$acc$mod <- 'phyb'
  
  df <- rbind(o.store$acc, p.store$acc)
  
  ggplot(df, aes(x = 1, y = uniq)) +
    geom_violin(aes(fill = mod, colour = mod), alpha = alpha, adjust = 1.5) +
    #geom_violin(data = p.store$acc, aes(1, uniq, fill = mod), alpha = alpha, adjust = 1.5) +
    #geom_point(data = df2, aes(1, uniq), colour = 'black', size = 4) +
    facet_wrap( ~ disease,
               labeller = as_labeller(create.lab()),
               ncol = length(unique(df$disease)),
               strip.position = 'bottom') +
  create.factheme() +
  labs(y = create.axlab('uniq'), subtitle = 'B') +
    scale_fill_manual(values = c(c('#7fc97f', '#beaed4'))) +
    scale_colour_manual(values = c(c('#7fc97f', '#beaed4'))) +
    scale_x_discrete(labels = as_labeller(create.lab())) +
    scale_y_continuous(breaks = seq(0, 1, 0.25),
                       minor_breaks = seq(0, 1, 0.25/2))

} 

## vis.lm but facet_wrapped over ob and phyb
vis.lm <- function(o.store, p.store, a, b, mod = 'exp', pol = 1,
                        y.lim = NULL, mod2 = mod, lab = NULL) {
  
  if(is.null(y.lim)) {
    if(b %in% c('gensig', 'ent', 'dna.ent', 'nodna.ent')) {
      y.lim <- ylim(0, NA)
    } else if(b %in% 'improv') {
      y.lim <- ylim(NA, 1)
    } else if(b == 'improv.ent') {
      y.lim <- NULL
    }else {
      y.lim <- ylim(0, 1)
    }
  }

  o.df <- create.lm(o.store, a, b, mod, pol, mut = 1e-5)
  o.df$mod <- 'ob'
  p.df <- create.lm(p.store, a, b, mod2, pol, mut = 1e-5)
  p.df$mod <- 'phyb'

  df <- rbind(o.df, p.df)

  facet_labeller <- as_labeller(c('ob' = 'outbreaker2', 'phyb' = 'phybreak'))
  
  o.store$acc$mod <- 'ob'
  p.store$acc$mod <- 'phyb'

  o.store$acc <- rbind(o.store$acc, p.store$acc)
  
  p <- ggplot(df, aes(x, value)) +
    geom_point(data = o.store$acc, aes_string(a, b, colour = 'disease')) +
    facet_grid(. ~ mod, labeller = facet_labeller) +
    geom_line(aes(group = variable, size = variable, linetype = variable)) +
    scale_size_manual(values = c(1.5, 1, 1)) +
    scale_linetype_manual(values = c('solid', 'dashed', 'dashed')) +
    labs(x = create.axlab(a), y = create.axlab(b), subtitle = lab) +
    vis.lab(o.store, a, b, double = TRUE) +
    create.lmtheme() +
    y.lim
  
  return(p)
  
}

## Violin plot of prop for each pathogen for o and p
vis.raw.acc <- function(o.store, p.store, alpha = 0.2) {

  o.store$acc$disease <- factor(o.store$acc$disease, levels = rev(sort.gensig(o.store)))
  o.store$acc$mod <- 'ob'
  
  p.store$acc$disease <- factor(p.store$acc$disease, levels = rev(sort.gensig(p.store)))
  p.store$acc$mod <- 'phyb'
  
  df <- rbind(o.store$acc, p.store$acc)

  facet_labeller <- as_labeller(c('ob' = 'outbreaker2', 'phyb' = 'phybreak'))
  
  ggplot(df) +
    geom_violin(aes(x = disease, y = nodna, fill = 'No WGS'), alpha = alpha, adjust = 1.5) +
    geom_violin(aes(x = disease, y = dna, fill = 'WGS'), alpha = alpha, adjust = 1.5) +
    #geom_point(data = df2, aes(1, uniq), colour = 'black', size = 4) +
    facet_wrap( ~ mod, nrow = 2, labeller = facet_labeller) +
    theme_minimal(base_size = 24) +
    scale_x_discrete(name = NULL, labels = create.lab()) +
    theme(legend.position = 'bottom',
          strip.text.x = element_text(colour = "black"),
          strip.text.y = element_text(colour = "black"),
          strip.background = element_rect(colour = "darkgrey", fill = "grey90")) +
                                        #labs(y = create.axlab('uniq')) +
    labs(y = create.axlab('dna')) +
  scale_fill_manual(name = NULL, values = c('#ffff99', '#fdc086'))

}

## Violin plot of prop for each pathogen for o and p
vis.raw.ent <- function(o.store, p.store, alpha = 0.2) {

  o.store$acc$disease <- factor(o.store$acc$disease, levels = rev(sort.gensig(o.store)))
  o.store$acc$mod <- 'ob'
  
  p.store$acc$disease <- factor(p.store$acc$disease, levels = rev(sort.gensig(p.store)))
  p.store$acc$mod <- 'phyb'
  
  df <- rbind(o.store$acc, p.store$acc)

  facet_labeller <- as_labeller(c('ob' = 'outbreaker2', 'phyb' = 'phybreak'))
  
  ggplot(df) +
    geom_violin(aes(x = disease, y = nodna.ent, fill = 'No WGS'), alpha = alpha, adjust = 1.5) +
    geom_violin(aes(x = disease, y = dna.ent, fill = 'WGS'), alpha = alpha, adjust = 1.5) +
    #geom_point(data = df2, aes(1, uniq), colour = 'black', size = 4) +
    facet_wrap( ~ mod, nrow = 2, labeller = facet_labeller) +
    theme_minimal(base_size = 24) +
    scale_x_discrete(name = NULL, labels = create.lab()) +
    theme(legend.position = 'bottom',
          strip.text.x = element_text(colour = "black"),
          strip.text.y = element_text(colour = "black"),
          strip.background = element_rect(colour = "darkgrey", fill = "grey90")) +
                                        #labs(y = create.axlab('uniq')) +
    labs(y = create.axlab('dna.ent')) +
  scale_fill_manual(name = NULL, values = c('#ffff99', '#fdc086'))

} 



##===== SAVING =====##

## Wrapper for ggsave to put everything in the right directory
vis.save <- function(p, name, ext = 'svg', dpi = 250,
                     width = 15, height = 10, ...) {

  ggsave(paste0("../figs/", name, ".", ext),
         p, dpi = dpi, width = width, height = height,
         ...)
}

## Update plots in the final manuscript folder
create.figs <- function() {

  p <- vis.gensig(o.store, p.store, 1)
  q <- vis.uniq(o.store, p.store, 1)
  
  r <- vis.lm(o.store, p.store, 'gensig', 'improv', mod = 'exp', lab = "A")
  s <- vis.lm(o.store, p.store, 'gensig', 'improv.ent', mod = 'neg.exp', lab = "B")
  
  t <- vis.lm(o.store, p.store, 'uniq', 'improv', mod = 'pol', lab = "A")
  u <- vis.lm(o.store, p.store, 'uniq', 'improv.ent', mod = 'pol', lab = "B")
  
  v <- vis.raw.acc(o.store, p.store, 0.5)
  w <- vis.raw.ent(o.store, p.store, 0.5)
  
  vis.save(joint.legend(p, q), 'Fig1', ext = 'png', width = 20, height = 20)
  vis.save(grid.arrange(r, s), 'Fig2', ext = 'png', width = 20, height = 21)
  vis.save(grid.arrange(t, u), 'Fig3', ext = 'png', width = 20, height = 21)
  vis.save(v, 'S1_Fig', ext = 'png', width = 18, height = 20)
  vis.save(w, 'S2_Fig', ext = 'png', width = 18, height = 20)

}



##===== HELPER FUNCTIONS =====##

## Returns the dataframe output of a linear model
create.lm <- function(store, a, b, mod, pol, mut) {

  ## Initialise the linear model
  dat <<- data.frame(x = store$acc[['gensig']],
                     y = store$acc[['prop']])

  lin <<- nls(y ~ 1 - (1 - 0.001)^(fac*x/0.001), data = dat,
              start = list(fac = 1),
              trace = TRUE)

  cof <<- coef(lin)

  lin <<- as.lm(lin)

  if(mod == 'exp') {
    dat <<- data.frame(x = store$acc[[a]],
                       y = store$acc[[b]])

    lin <<- nls(y ~ i - i*exp(-a*x), data = dat,
                start = list(a = 1, i = 1),
                trace = TRUE, control = list(maxiter = 200))

    lin <<- as.lm(lin)
    
  } else if(mod == 'neg.exp') {
    dat <<- data.frame(x = store$acc[[a]],
                       y = store$acc[[b]])

    lin <<- nls(y ~ i*exp(-a*x) - i, data = dat,
                start = list(a = 1, i = 3),
                trace = TRUE, control = list(maxiter = 200))

    lin <<- as.lm(lin)
    
  } else if(mod == 'grow') {
    dat <<- data.frame(x = store$acc[[a]],
                       y = store$acc[[b]])

    lin <<- nls(y ~ exp(a*x) - i, data = dat,
                start = list(a = 1, i = 1),
                trace = TRUE, control = list(maxiter = 200))

    lin <<- as.lm(lin)

  } else if(mod == 'asymp') {

    dat <<- data.frame(x = store$acc[[a]],
                      y = store$acc[[b]])

    lin <<- nls(y ~ 1 - (1 - mut)^(fac*x/mut), data = dat,
               start = list(fac = 1),
               trace = TRUE, control = list(maxiter = 200))

    cof <<- coef(lin)

    lin <<- as.lm(lin)

  } else if(mod == 'neg') {

    dat <<- data.frame(x = store$acc[[a]],
                      y = store$acc[[b]])

    lin <<- nls(y ~ i*exp(-a*x), data = dat,
                start = list(a = 1, i = 3),
                trace = TRUE, control = list(maxiter = 200))

    cof <<- coef(lin)

    lin <<- as.lm(lin)

  } else if(mod == 'pol') {

    var <- store$acc[[b]]
    q <- store$acc[[a]]

    lin <<- lm(var ~ poly(q, pol))

  }

  df <- as.data.frame(predict(lin, data.frame(x = store$acc[[a]]),
                              interval = 'prediction',
                              level = 0.95)) %>%
    mutate(x = store$acc[[a]]) %>%
    melt(id.vars = 'x') 
  
  return(df)
  
}

## A labeller for pathogen names
create.lab <- function(input) {

  c(ebola = "EBOV",
    sars = "SARS-CoV",
    mers = "MERS-CoV",
    ifz = "Influenza A",
    mrsa = "MRSA",
    klebs = "K. pneumoniae",
    strep = "S. pneumoniae",
    shig = "S. sonnei",
    tb = "M. tuberculosis",
    cdif = "C. difficile")

}

## Create a labeller for x and y axes from variable names
create.axlab <- function(input) {

  reference <- c(gensig = 'Average transmission divergence of outbreak',
                 improv = 'Change in accuracy of outbreak reconstruction',
                 dna = 'Accuracy of outbreak reconstruction',
                 nodna = 'Accuracy of outbreak reconstruction',
                 prop = 'Proportion of genetically distinct transmission pairs',
                 prop.id = 'Proportion of genetically identical transmission pairs',
                 uniq = 'Number of unique sequences / outbreak size',
                 dna.ent = 'Posterior entropy',
                 nodna.ent = 'Posterior entropy',
                 improv.ent = 'Change in posterior entropy')

  if(!input %in% names(reference)) {
    return(input)
  } else {
    return(reference[[input]])
  }
  
}

## Create theme for uniq/gensig
create.factheme <- function() {

    theme_minimal(base_size = 24) +
    theme(#legend.position = 'none',
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_blank(),
          #axis.text.x = element_text(size = 40),
      panel.spacing = unit(-0.45, "lines"),
      legend.direction = 'horizontal',
      plot.subtitle=element_text(size = 35))

}

## Create theme for vis.lm
create.lmtheme <- function() {

  theme_light(base_size = 24) +
    theme(legend.position = 'none',
          strip.text.x = element_text(colour = "black"),
          strip.text.y = element_text(colour = "black"),
          strip.background = element_rect(colour = "darkgrey", fill = "grey90"),
          plot.subtitle=element_text(size = 35))

}

## Add label
vis.lab <- function(store, a, b, size = 5, double = FALSE) {

  df <- mk.lab(store, a, b, double = double)

  geom_label_repel(data = df,
                   aes(a, b, colour = disease,
                       label = lab),
                   segment.alpha = 0,
                   size = size,
                   show.legend = FALSE)
  
}

## Create df for vis.lab
mk.lab <- function(store, a, b, double = FALSE) {

  if(double) {
  
    df <- store$acc %>%
      group_by(disease, mod) %>%
      summarise_(a = interp(~mean(var), var = as.name(a)),
                 b = interp(~mean(var), var = as.name(b)))

  } else {
    
    df <- store$acc %>%
      group_by(disease) %>%
      summarise_(a = interp(~mean(var), var = as.name(a)),
                 b = interp(~mean(var), var = as.name(b)))

  }
    
  if(double) {
    df$lab <- rep(create.lab()[sort.gensig(store)], each = 2)
  } else {
    df$lab <- create.lab()[sort.gensig(store)]
  }

  return(df)

}

## Extracts a ggplot legend
extract.legend <- function(gplot){

    tmp <- ggplot_gtable(ggplot_build(gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)

}

## Arranges two plots horizontally with a joint legend at the bottom (from p1)
joint.legend <- function(p1, p2, pos = 'bottom') {

    leg <- extract.legend(p1)

    if(pos == 'bottom') {
    p <- grid.arrange(arrangeGrob(p1 + theme(legend.position = 'none'),
                                   p2 + theme(legend.position = 'none'),
                                   nrow = 2),
                      leg, nrow = 2, ncol = 1, heights = c(10, 0.5))
    }

    if(pos == 'right') {
        p <- grid.arrange(arrangeGrob(p1 + theme(legend.position = 'none'),
                                      p2 + theme(legend.position = 'none'),
                                      nrow = 2),
                          leg, nrow = 1, ncol = 2, widths = c(10, 0.5))
    }

    return(p)

}

## Random package not downloaded properly
as.proto.list <- function(x, envir, parent, all.names = FALSE, ...,
                          funEnvir = envir, SELECT = function(x) TRUE) {
  if (missing(envir)) {
    if (missing(parent))
      parent <- parent.frame()
    envir <- if (is.proto(parent))
               parent$proto(...)
             else
               proto(parent, ...)
  }
  for (s in names(x))
    if (SELECT(x[[s]])) {
      assign(s, x[[s]], envir = envir)
      if (is.function(x[[s]]) && !identical(funEnvir, FALSE))
        environment(envir[[s]]) <- funEnvir
    }
  if (!missing(parent))
    parent.env(envir) <- parent
  as.proto.environment(envir)  # force refresh of .that and .super
}
as.proto.environment <- function(x, ...) {
  assign(".that", x, envir = x)
  assign(".super", parent.env(x), envir = x)
  structure(x, class = c("proto", "environment"))
}
as.lm <- function(object, ...) UseMethod("as.lm")
as.lm.nls <- function(object, ...) {
  if (!inherits(object, "nls")) {
    w <- paste("expected object of class nls but got object of class:",
               paste(class(object), collapse = " "))
    warning(w)
  }

  gradient <- object$m$gradient()
  if (is.null(colnames(gradient))) {
    colnames(gradient) <- names(object$m$getPars())
  }

  response.name <- if (length(formula(object)) == 2) "0" else
                                                           as.character(formula(object)[[2]])

  lhs <- object$m$lhs()
  L <- data.frame(lhs, gradient)
  names(L)[1] <- response.name

  fo <- sprintf("%s ~ %s - 1", response.name,
		paste(colnames(gradient), collapse = "+"))
  fo <- as.formula(fo, env = as.proto.list(L))

  do.call("lm", list(fo, offset = substitute(fitted(object))))

}
