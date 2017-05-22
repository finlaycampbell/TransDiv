#====== Plots =====#

## Plot a discrete violin plot of gensig distributions
vis.gensig <- function(store) {

    ## Calculate proportions manually
    df <- by(store$gensig$gensig,
             store$gensig$disease,
             function(sig) prop.table(table(sig)))
    df <- plyr::ldply(df, data.frame)
    #df$sig <- as.numeric(df$sig)
    df$sig <- factor(df$sig, levels = as.character(0:max(as.numeric(df$sig))))

    df$.id <- factor(df$.id, levels = rev(sort.gensig(store)))
    df$Freq <- df$Freq/2
    df2 <- df
    df2$Freq <- -df$Freq

    p <- ggplot(df, aes(x = sig, y = Freq, fill = .id, colour = .id)) +
        geom_bar(stat = "identity") +
        geom_bar(data = df2,
                 mapping = aes(x = sig, y = Freq, fill = .id, colour = .id),
                 stat = 'identity') +
        geom_vline(xintercept = 1.5, linetype = 'dashed') +
        theme_minimal(base_size = 24) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              strip.background = element_blank(),
              panel.spacing = unit(-0.45, "lines"),
              legend.position = "none") +
        facet_wrap( ~ .id,
                   labeller = as_labeller(create.lab()),
                   ncol = length(unique(df$.id)),
                   strip.position = 'bottom') +
        labs(x = "Transmission divergence") +
        geom_label(data = data.frame(sig = 13.5,
                                    Freq = 0,
                                    .id = factor('shig', levels = sort.gensig(store)),
                                    lab = 'A'),
                   label = 'A', size = 15, colour = 'black',
                   fill = 'white', label.r = unit(0, 'cm'), label.size = unit(0, 'cm')) +
        coord_flip()

    p

}

## Plot the analytical solution to the genetic signature
vis.anylt <- function(store) {

    param <- create.param()

    df <- sapply(names(param), function(i) with(param[[i]], mut*seql*w.mean)) %>%
        plyr::ldply(data.frame)



    p <- store$gensig %>%
        group_by(disease) %>%
        summarise(mean = mean(gensig)) %>%
        mutate(anlyt = df$X..i..[match(disease, df$.id)]) %>%
        gather(model, value, -disease) %>%
        ggplot(aes(disease, value, fill = disease, alpha = model)) +
        geom_bar(stat = 'identity', position = 'dodge') +
        labs(x = 'Pathogen', y = 'Genetic signature') +
        scale_x_discrete(labels = as_labeller(create.lab())) +
        scale_alpha_manual(values = c(0.5, 1),
                           name = NULL,
                           labels = c('Analytical', 'Simulated')) +
        theme(legend.position = 'bottom') +
        guides(fill = FALSE)
    p

}

## prop ~ improv
vis.atleast <- function(store) {

    #pred <- store$acc$prop
    #pred[length(pred)] <- 1

    vis.lm(store, 'prop', 'improv', mod = 'pol', pol = 1, add.line = TRUE) +
        vis.lab(store, 'prop', 'improv') +
        scale_y_continuous(breaks = seq(0, 1, 0.25))

}

## Plot the distribution of inference accuracy (dna vs nodna)
vis.acc <- function(store) {

    store$acc$disease <- factor(store$acc$disease, levels = sort.gensig(store))
    df <- melt(store$acc, measure.vars = c('dna', 'nodna'))

    p <- ggplot(df, aes(variable, value, fill = disease, alpha = variable)) +
        geom_violin(scale = 'width') +
        facet_wrap( ~ disease, labeller = as_labeller(create.lab())) +
        scale_alpha_manual(name = 'Data', values = c(1, 0.2)) +
        theme(legend.position = 'none') +
        labs(y = 'Accuracy of outbreak reconstruction', x = NULL)
    p

}

## gensig ~ improv
vis.rel <- function(store) {

    vis.lm(store, 'gensig', 'dna', mod = 'asymp') +
        scale_y_continuous(limits = c(0, 1), breaks = seq(-0.2, 1, 0.2))

}

## Violin plot of prop for each pathogen
vis.prop <- function(store) {

    store$acc$prop <- 1 - store$acc$prop

    ggplot(store$acc, aes(x = 1, y = prop, fill = disease)) +
        geom_violin() +
        facet_wrap( ~ disease,
                   labeller = as_labeller(create.lab()),
                   ncol = length(unique(store$acc$disease)),
                   strip.position = 'bottom') +
        theme_minimal(base_size = 24) +
        theme(legend.position = 'none',
              axis.title.x = element_blank(),
              #axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              strip.background = element_blank(),
              axis.text.x = element_text(size = 40),
              panel.spacing = unit(-0.45, "lines")) +
        labs(x = 'Pathogen', y = 'Proportion of genetically identical transmission pairs') +
        scale_x_discrete(labels = as_labeller(create.lab())) +
        scale_y_continuous(breaks = seq(0, 1, 0.25),
                           minor_breaks = seq(0, 1, 0.25/2)) +
        geom_label(data = data.frame(x = 1,
                                    prop = 1.1,#0.9375,
                                    disease = factor('shig', levels = sort.gensig(store)),
                                    lab = 'B'),
                  label = 'B', size = 15, colour = 'black',
                  fill = 'white', label.r = unit(0, 'cm'), label.size = unit(0, 'cm'))

}

## Visualise size of outbreaks
vis.n <- function(store) {

    ggplot(store$sim, aes(n, fill = disease)) +
        geom_histogram() +
        facet_wrap( ~ disease, labeller = as_labeller(create.lab())) +
        theme(legend.position = 'none')

}

## Visualise generation time distributions
vis.w <- function() {

    param <- create.param()
    df <- data.frame(x = seq_len(length(param[[1]]$w)))
    for(i in names(param)) df[[i]] <- param[[i]]$w
    df <- gather(df, disease, value, -x)

    p <- ggplot(df, aes(x, value, colour = disease)) + geom_line()
    p

}

## Fit a polynomial model
vis.lm <- function(store, a, b, pol = 3, mod = 'asymp',
                   mut = 0.00001, add.line = FALSE) {

    ## Initialise the linear model
    dat <<- data.frame(x = store$acc[['gensig']],
                      y = store$acc[['prop']])

    lin <<- nls(y ~ 1 - (1 - 0.001)^(fac*x/0.001), data = dat,
               start = list(fac = 1),
               trace = TRUE)

    cof <<- coef(lin)

    lin <<- as.lm(lin)

    if(mod == 'exp') {
        dat <- data.frame(x = store$acc[[a]],
                          y = store$acc[[b]])

        lin <- nls(y ~ 1 - exp(-a*x), data = dat,
                   start = list(a = 1), trace = TRUE)

        lin <- as.lm(lin)

    }

    if(mod == 'asymp') {

        dat <- data.frame(x = store$acc[[a]],
                          y = store$acc[[b]])

        lin <- nls(y ~ 1 - (1 - mut)^(fac*x/mut), data = dat,
                   start = list(fac = 1),
                   trace = TRUE)

        cof <- coef(lin)

        lin <- as.lm(lin)

    }

    if(mod == 'pol') {

        var <- store$acc[[b]]
        q <- store$acc[[a]]

        lin <- lm(var ~ poly(q, pol))

    }

    if(add.line) unity.line <- geom_segment(x = 0, y = 0, xend = 1, yend = 1, size = 1.5,
                                            linetype = 'solid', colour = 'gray65')
    else unity.line <- NULL

    as.data.frame(predict(lin, data.frame(x = store$acc[[a]]),
                          interval = 'prediction',
                          level = 0.95)) %>%
        mutate(x = store$acc[[a]]) %>%
        melt(id.vars = 'x') %>%
        ggplot(aes(x, value)) +
        geom_point(data = store$acc, aes_string(a, b, colour = 'disease')) +
        unity.line +
        geom_line(aes(group = variable, size = variable, linetype = variable)) +
        scale_size_manual(values = c(1.5, 1, 1)) +
        scale_linetype_manual(values = c('solid', 'dashed', 'dashed')) +
        labs(x = create.axlab(a), y = create.axlab(b)) +
        vis.lab(store, a, b) +
        theme_minimal(base_size = 18) +
        theme(legend.position = 'none')

}

## Add label
vis.lab <- function(store, a, b, size = 5) {

    df <- store$acc %>%
        group_by(disease) %>%
        summarise_(a = interp(~mean(var), var = as.name(a)),
                   b = interp(~mean(var), var = as.name(b)))

    geom_label_repel(data = df,
                     aes(a, b, colour = disease,
                         label = create.lab()[rev(sort.gensig(store))]),
                     segment.alpha = 0,
                     size = size,
                     show.legend = FALSE)
}

## Generalised function to plot correlate variables a and b
vis.corr <- function(store, a, b, xlab = NULL, ylab = NULL) {

    if(is.null(xlab)) xlab <- a
    if(is.null(ylab)) ylab <- b

    df <- store$acc %>%
        group_by(disease) %>%
        summarise_(a = interp(~mean(var), var = as.name(a)),
                   b = interp(~mean(var), var = as.name(b)))

    ggplot(store$acc, aes_string(a, b)) +
        geom_point(aes(colour = disease)) +
        geom_smooth() +
        geom_label_repel(data = df,
                         aes(a, b, colour = disease,
                             label = create.lab()[rev(sort.gensig(store))]),
                         segment.alpha = 0,
                         show.legend = FALSE) +
        theme(legend.position = 'none') +
        labs(x = xlab, y = ylab)


}
