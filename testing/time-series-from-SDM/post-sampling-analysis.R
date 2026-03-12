# if you are running this after the fact: be sure to load the RData file

catout <- function(printtext) {
    cat(paste0("\n", printtext, "\n"))
}

# check run time and some MH acceptance rate summaries --------------------

print(endtime - starttime)

# check the acceptance rates
# gibbsPostBurn <- round(gibbsIts * burnin):gibbsIts
gibbsPostBurn <- (gibbsIts*burnin / t_thin + 1):num_samp
acc_rate_PostBurn <- (gibbsIts * burnin):gibbsIts

catout("sample-wise summary of Ukl acc. rates, all and post burn-in")
apply(accCount_s[, , 2:gibbsIts], c(1, 2), mean) |> t() |> summary()
apply(accCount_s[, , acc_rate_PostBurn], c(1, 2), mean) |> t() |> summary()

catout("overall summary of Ukl acc. rates, post burn-in")
apply(accCount_s[, , acc_rate_PostBurn], c(1, 2), mean) |>
    as.vector() |> 
    quantile(probs = c(0, 0.025, .25, .5, .75, .975, 1))

# check Sigmal acceptance rates
catout("Sigmal acc. rates, all and post burn-in")
rowMeans(accCount_Sigma_s) |> 
    quantile(probs = c(0, 0.025, .25, .5, .75, .975, 1))
rowMeans(accCount_Sigma_s[, acc_rate_PostBurn]) |> 
    quantile(probs = c(0, 0.025, .25, .5, .75, .975, 1))

# evaluate Lambdas --------------------------------------------------------

upper_q <- apply(Lambdak_l_s[, , , gibbsPostBurn], 
                 c(1, 2, 3), quantile, probs = 0.975)
lower_q <- apply(Lambdak_l_s[, , , gibbsPostBurn], 
                 c(1, 2, 3), quantile, probs = 0.025)
Lambda_means <- apply(Lambdak_l_s[, , , gibbsPostBurn], c(1, 2, 3), mean)

# compare quantiles and means for all frequencies of each k
for (k in 1:K) {
    plot(upper_q[1, k, ], type = "l", lty = 2, 
         ylim = c(0, max(Lambdakl0[1, k, 1:num_freqs])),
         main = paste0("post. mean and true Lambda, k = ", k),
         ylab = "Lambda")
    lines(lower_q[1, k, ], type = "l", lty = 2)
    lines(Lambda_means[1, k, ])
    lines(Lambdakl0[1, k, 1:num_freqs], lty = 3)

    lines(upper_q[2, k, ], type = "l", lty = 2, col = "red")
    lines(lower_q[2, k, ], type = "l", lty = 2, , col = "red")
    lines(Lambda_means[2, k, ], , col = "red")
    lines(Lambdakl0[2, k, 1:num_freqs], lty = 3, , col = "red")

    save_plot_pdf(file.path(result_dir, "post-Lambda-and-true",
                paste0("post-Lambda-and-true-", k, ".pdf")))
}

# distances to true Ukl0 for all Ukl --------------------------------------

ds_to_true <- array(NA, c(K, num_freqs))

for (k in 1:K) {
    for (l in 1:num_freqs) {
        avgUkl <- apply(U_kls_all[, , k, l, gibbsPostBurn],
                        c(1, 2), mean)
        avgUkl <- avgUkl %*% 
            solve( EigenR::Eigen_sqrt( t(Conj(avgUkl)) %*% avgUkl ) )
        
        ds_to_true[k, l] <- fast_evec_Frob_stat(U_kl0[, , k, l], avgUkl)
    }
}

# summary(as.vector(ds_to_true))
catout("summary of axis Frobenius distances from sample mean Ukl to true Ukl0")
quantile(as.vector(ds_to_true), probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1))
catout("how many distances are >= sqrt(2): ")
sum(ds_to_true >= sqrt(2))
# which(ds_to_true >= sqrt(2), arr.ind = TRUE)

# get the k and l indices for the Ukl with max and min distances to truth
catout("k and l indices of Ukl with max and min dist. to true Ukl0")
which.max(ds_to_true) |> arrayInd(.dim = dim(ds_to_true))
which.min(ds_to_true) |> arrayInd(.dim = dim(ds_to_true))

# plot the max and min distances per frequency
# create a temporary data frame first
rbind(data.frame(distance = apply(ds_to_true, 2, min), datalabel = "min"),
      data.frame(distance = apply(ds_to_true, 2, max), datalabel = "max")) |> 
    # ggplot
    ggplot(aes(x = rep(1:num_freqs, times = 2), 
               y = distance, color = datalabel)) +
    geom_line() +
    labs(x = "freq. index", y = "axis Frob. dist.", 
         title = "max. and min. dist. of post. mean Ukl to true Ukl0") +
    theme(legend.position = c(0, 1), legend.justification = c(0, 1),
          legend.background = element_rect(fill = alpha("white", 0.6)))

save_plot_pdf(file.path(result_dir, 
                        paste0("Ukl-dist-to-true-max-min.pdf")))

# calculate distances and plot ESS for all k and l ------------------------

all_Ukl_dists <- array(NA, c(K, num_freqs, num_samp))
ESS_Ukl_dists <- array(NA, c(K, num_freqs))

for (k in 1:K) {
    for (l in 1:num_freqs) {
        avgUkl <- apply(U_kls_all[, , k, l, gibbsPostBurn],
                        c(1, 2), mean)
        avgUkl <- avgUkl %*% solve(EigenR::Eigen_sqrt(t(Conj(avgUkl)) %*% 
                                                          avgUkl))
        for (s in 1:num_samp) {
            all_Ukl_dists[k, l, s] <- 
                fast_evec_Frob_stat(U_kls_all[, , k, l, s], avgUkl)
        }
        
        ESS_Ukl_dists[k, l] <- 
            posterior::ess_basic(all_Ukl_dists[k, l, gibbsPostBurn])
    }
}

# plot of min and max ESS at each frequency
# create a temporary data frame first
rbind(data.frame(ESS = apply(ESS_Ukl_dists, 2, min), datalabel = "min"),
      data.frame(ESS = apply(ESS_Ukl_dists, 2, max), datalabel = "max")
      ) |> 
    # ggplot
    ggplot(aes(x = rep(1:num_freqs, times = 2), 
               y = ESS, color = datalabel)) +
    geom_line() +
    labs(x = "freq. index", title = "max. and min. ESS of Ukl dist.") +
    theme(legend.position = c(0, 1), legend.justification = c(0, 1),
          legend.background = element_rect(fill = alpha("white", 0.6)))

save_plot_pdf(file.path(result_dir, 
                        paste0("Ukl-dist-summary-ESS-max-min.pdf")))

# create distance summary trace plots for closest and furthest Ukl --------

# get the k and l indices for the Ukl with max and min distances to truth
ind_close_Ukl <- which.min(ds_to_true) |> arrayInd(.dim = dim(ds_to_true))
ind_far_Ukl <- which.max(ds_to_true) |> arrayInd(.dim = dim(ds_to_true))

# distance trace plot of closest Ukl
plot(all_Ukl_dists[ind_close_Ukl[1], ind_close_Ukl[2], ], 
     type = "l",
     main = paste0("Trace plot, (closest) Ukl axis dist. to mean, k = ", 
                   ind_close_Ukl[1], ", l = ", ind_close_Ukl[2]),
     ylab = "axis Frobenius dist.")
save_plot_pdf(file.path(result_dir, 
                        paste0("post-Ukl-trace-closest-axis-dist-k-",
                               ind_close_Ukl[1], "-l-", ind_close_Ukl[2], 
                               ".pdf")))

# distance trace plot of furthest Ukl
plot(all_Ukl_dists[ind_far_Ukl[1], ind_far_Ukl[2], ], 
     type = "l",
     main = paste0("Trace plot, (furthest) Ukl axis dist. to mean, k = ", 
                   ind_far_Ukl[1], ", l = ", ind_far_Ukl[2]),
     ylab = "axis Frobenius dist.")
save_plot_pdf(file.path(result_dir, 
                        paste0("post-Ukl-trace-furthest-axis-dist-k-", 
                               ind_far_Ukl[1], "-l-", ind_far_Ukl[2], ".pdf")))

# evaluate Sigmal ---------------------------------------------------------

quantile(n_Sig, c(0, .025, .25, .5, .75, .975, 1))

l <- ind_far_Ukl[2]

# check Sigmal acceptance rates
catout("total and post burn-in acc. rates of Sigmal at l for furthest Ukl")
mean(accCount_Sigma_s[l, ])
mean(accCount_Sigma_s[l, gibbsPostBurn])

avgSigmal <- apply(Sigmal_s[, , l, gibbsPostBurn], c(1, 2), mean)

d_to_avgSigmal <- rep(NA, num_samp)
for (s in 1:num_samp) {
    d_to_avgSigmal[s] <- frob_dist(Sigmal_s[, , l, s], avgSigmal)
}

plot(d_to_avgSigmal, type = "l",
     main = paste0("Trace plot, Sigmal dist. to mean, l = ", l),
     ylab = "Frobenius dist.")
save_plot_pdf(file.path(result_dir, 
    paste0("post-Sigmal-dist-trace-l-", l, ".pdf")))

par(mfrow = c(2, 2), mar = c(2, 2, 2, 1))
for(j in 1:4) {
    plot(Re(Sigmal_s[j, j, l, ]), type = "l", 
         main = paste0("Trace plot, Sigmal[j,j], l = ", l, ", j = ", j),
         xlab = "", ylab = "")
}
save_plot_pdf(file.path(result_dir, 
    paste0("post-Sigmal-diag-trace-l-", l, ".pdf")))

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

catout("how many and what fraction of Sigmal have 0 acceptance rate?")
sum(rowMeans(accCount_Sigma_s[, gibbsPostBurn]) == 0)
# what fraction have 0 acceptance rate?
sum(rowMeans(accCount_Sigma_s[, gibbsPostBurn]) == 0) / num_freqs
# show the indices of which ones have 0 acceptance rate
# which(rowMeans(accCount_Sigma_s[, gibbsPostBurn]) == 0)

# smoothing parameter summaries -------------------------------------------

plot(taujkl2_s[1, 1, 1, gibbsPostBurn], type = "l")

catout("median of taujkl2 samples at some diff j, k, and l")
median(taujkl2_s[1, 1, 1, gibbsPostBurn])

median(taujkl2_s[1, 1, 150, gibbsPostBurn])
median(taujkl2_s[2, 1, 150, gibbsPostBurn])

median(taujkl2_s[1, 1, 256, gibbsPostBurn])
median(taujkl2_s[1, 1, 350, gibbsPostBurn])

# evaluate sigmak2 -------------------------------------------------------

catout("true sigmak02 and posterior quantiles for some k")

k <- 1
sigmak02[k]
quantile(sigmak2_s[k, gibbsPostBurn], c(.025, .5, .975))

plot(sigmak2_s[k, ], type = "l",
     main = paste0("Trace plot, sigmak2, k = ", k),
     ylab = "sigmak2")
abline(h = sigmak02[k])
save_plot_pdf(file.path(result_dir, paste0("post-sigmak2-k-", k, ".pdf")))

k <- 2
sigmak02[k]
quantile(sigmak2_s[k, gibbsPostBurn], c(.025, .5, .975))

plot(sigmak2_s[k, ], type = "l",
     main = paste0("Trace plot, sigmak2, k = ", k),
     ylab = "sigmak2")
abline(h = sigmak02[k])
save_plot_pdf(file.path(result_dir, paste0("post-sigmak2-k-", k, ".pdf")))

# evaluate full SDM estimates ---------------------------------------------

# sigmakl2 * (Ukl %*% Lambdakl %*% t(Conj(Ukl)) + I_P)
posterior_dists <- array(NA, c(K, num_freqs))
multitaper_dists <- array(NA, c(K, num_freqs))
post_mean_SDMs <- array(NA, c(P, P, K, num_freqs))

post_SDMs <- array(NA, c(P, P, K, num_freqs, length(gibbsPostBurn)))
post_sq_cohe_phase <- data.frame()

for (k in 1:K) {
    cat(paste0(k, " of ", K, ": post. SDM, sq. coherence, and phase\n"))
    
    # calculate individual and posterior mean SDMs
    for (l in 1:num_freqs) {
        meanSDM <- matrix(0 + 0i, P, P)
        for (s in gibbsPostBurn) {
            thisSDM <- sigmak2_s[k, s] * 
                (U_kls_all[, , k, l, s] %*% diag(Lambdak_l_s[, k, l, s])
                 %*% t(Conj(U_kls_all[, , k, l, s])) + diag(P))
            
            t_idx <- s - (gibbsIts*burnin)/t_thin
            
            post_SDMs[, , k, l, t_idx] <- thisSDM
            meanSDM <- meanSDM + thisSDM
        }
        
        meanSDM <- meanSDM / length(gibbsPostBurn)
        post_mean_SDMs[, , k, l] <- meanSDM
        
        # compare distances
        posterior_dists[k, l] <- frob_dist(meanSDM, fkTR[, , k, l])
        multitaper_dists[k, l] <- frob_dist(data_list_w[[l]][[k]] / LL, 
                                            fkTR[, , k, l])
        
        # calculate posterior mean squared coherence and phase
        theseCohe <- data.frame()
        for (ir in 1:(P-1)) {
            for (ic in (ir+1):P) {
                thisCohe <- mean(Mod(post_SDMs[ir, ic, k, l, ])^2 / 
                                     (post_SDMs[ir, ir, k, l, ] * 
                                          post_SDMs[ic, ic, k, l, ]))
                
                
                thisPhase <- mean(Arg(post_SDMs[ir, ic, k, l, ]))
                
                theseCohe <- 
                    rbind(theseCohe,
                          data.frame(ir = ir, ic = ic, k = k,
                                     sq_cohe = Re(thisCohe), phase = thisPhase,
                                     freq_index = l, datalabel = "posterior"))
            }
        }
        
        post_sq_cohe_phase <- rbind(post_sq_cohe_phase, theseCohe)
    }
}

# Mean Squared Error
catout("true sigmak02")
sigmak02

catout("individual MSE for estimated SDMs")
cbind(
    posterior = rowMeans((posterior_dists/P)^2),
    multitaper = rowMeans((multitaper_dists/P)^2),
    scaled_post = rowMeans((posterior_dists/P)^2 / sigmak02),
    scaled_multi = rowMeans((multitaper_dists/P)^2 / sigmak02)
)

# AMSE, unscaled and scaled matrices
catout("AMSE, posterior and multitaper")
c(rowMeans((posterior_dists/P)^2) |> mean(),
  rowMeans((multitaper_dists/P)^2) |> mean())

catout("AMSE, scaled posterior vs. scaled multitaper")
c(rowMeans((posterior_dists/P)^2 / sigmak02) |> mean(),
  rowMeans((multitaper_dists/P)^2 / sigmak02) |> mean())

for (k in 1:K) {
    cat(paste0("\nk = ", k, "\n"))
    quantile(as.vector(posterior_dists[k, ]), 
             probs = c(0, .025, .5, .975, 1)) |> print()
    quantile(as.vector(multitaper_dists[k, ]), 
             probs = c(0, .025, .5, .975, 1)) |> print()
    
    plot(density(as.vector(posterior_dists[k, ]), from = 0), col = 1,
         main = paste0("densities of SDM estimate distances, k = ", k),
         xlab = "Frob. distance", ylab = "density")
    lines(density(as.vector(multitaper_dists[k, ]), from = 0), col = 2)
    legend(x = "topright", legend = c("posterior", "multitaper"),
           col = c(1, 2), lwd = 2)
    
    save_plot_pdf(file.path(result_dir, "post-SDM-est-dist-density",
                            paste0("post-SDM-est-dist-density-k-", k, ".pdf")))
}

# compare true, multitaper, and posterior mean SDMs -----------------------

SDMests_df <- data.frame()
trueSDM_df <- data.frame()
post_ests_df <- data.frame()

for (k in 1:K) {
    for (j in 1:P) {
        SDMests_df <- rbind(SDMests_df, data.frame(
            power = Re(SDMests[[k]][j, j, 1:num_freqs]), 
            k = k, j = j, l = 1:num_freqs, datalabel = "multitaper"))
        
        trueSDM_df <- rbind(trueSDM_df, data.frame(
            power = Re(fkTR[j, j, k, 1:num_freqs]), 
            k = k, j = j, l = 1:num_freqs, datalabel = "true"))
        
        post_ests_df <- rbind(post_ests_df, data.frame(
            power = Re(post_mean_SDMs[j, j, k, ]),
            k = k, j = j, l = 1:num_freqs, datalabel = "post. mean"))
    }
}

SDMests_df <- dplyr::mutate(SDMests_df, across(c(k, j, datalabel), as.factor))
trueSDM_df <- dplyr::mutate(trueSDM_df, across(c(k, j, datalabel), as.factor))
post_ests_df <- dplyr::mutate(post_ests_df, 
                              across(c(k, j, datalabel), as.factor))

all_df <- rbind(SDMests_df, trueSDM_df, post_ests_df)
all_df$datalabel <- factor(all_df$datalabel, 
                           levels = c("multitaper", "post. mean", "true"))

for (kk in 1:K) {
    plotp <- dplyr::filter(all_df, k == kk) |> 
        ggplot(aes(x = l, y = power)) +
        geom_line(aes(color = j, linetype = datalabel), linewidth = .3) +
        ggtitle(paste0("True and est. SDMs, k = ", kk))
    
    print(plotp)
    save_plot_pdf(file.path(result_dir, "compare-SDM-powers",
                            paste0("compare-SDM-powers-k-", kk, ".pdf")))
}

# squared coherence and phase ---------------------------------------------

# fkTR[, , k, l] is a true SDM
# post_mean_SDMs[, , k, l] is a posterior mean SDM

multi_SDMs <- array(NA, c(P, P, K, num_freqs))
for (k in 1:K) {
    for (l in 1:num_freqs) {
        multi_SDMs[, , k, l] <- data_list_w[[l]][[k]] / LL
    }
}

for (k in 1:K) {
    all_cohe_phase <- data.frame()
    
    for (ir in 1:(P-1)) {
        for(ic in (ir+1):P) {
            true_sq_cohe <- 
                Re(Mod(fkTR[ir, ic, k, ])^2 /
                       (fkTR[ir, ir, k, ] * fkTR[ic, ic, k, ]))
            # post_sq_cohe <- 
            #     Re(Mod(post_mean_SDMs[ir, ic, k, ])^2 / 
            #            (post_mean_SDMs[ir, ir, k, ] * post_mean_SDMs[ic, ic, k, ]))
            multi_sq_cohe <- 
                Re(Mod(multi_SDMs[ir, ic, k, ])^2 / 
                       (multi_SDMs[ir, ir, k, ] * multi_SDMs[ic, ic, k, ]))
            
            true_phase <- Arg(fkTR[ir, ic, k, 1:num_freqs])
            # post_phase <- Arg(post_mean_SDMs[ir, ic, k, ])
            multi_phase <- Arg(multi_SDMs[ir, ic, k, 1:num_freqs])
            
            all_cohe_phase <- 
                rbind(all_cohe_phase,
                      data.frame(
                          ir = ir, ic = ic, k = k,
                          sq_cohe = c(true_sq_cohe[1:num_freqs], 
                                      multi_sq_cohe[1:num_freqs]),
                          phase = c(true_phase, multi_phase),
                          freq_index = 1:num_freqs,
                          datalabel = rep(c("true", "multitaper"), 
                                          each = num_freqs))
                )
        }
    }
    
    # combine with the posterior mean coherence and phase calculated earlier
    all_cohe_phase <- 
        post_sq_cohe_phase[post_sq_cohe_phase$k == k, ] |> 
        rbind(all_cohe_phase)
    
    # plot squared coherence
    plot_sq_cohe <- all_cohe_phase |> 
        # filter and just plot the first few, for now. need to expand this for larger P
        dplyr::filter(ir <= 3 & ic <= 4) |> 
        ggplot(aes(x = freq_index, y = sq_cohe)) +
        ggplot2::facet_wrap(vars(ir, ic), scales = "free_y", 
                            strip.position = "left", 
                            labeller = \(x) label_both(x, multi_line = FALSE)) +
        geom_line(aes(linetype = datalabel, color = datalabel), 
                  linewidth = .2) +
        scale_color_manual(values = c("true" = "red", "multitaper" = "black", 
                                      "posterior" = "black")) +
        labs(title = paste0("Estimated squared coherence for k = ", k)) +
        theme(legend.position = c(0, 1), legend.justification = c(0, 1),
              legend.background = element_rect(fill = alpha("white", 0.6)),
              legend.title = element_blank())
    
    print(plot_sq_cohe)
    save_plot_pdf(file.path(result_dir, "SDM-est-coherence-and-phase",
                            paste0("SDM-est-squared-coherences-k-", k, ".pdf")))
    
    # plot phase
    plot_phase <- all_cohe_phase |> 
        # filter and just plot the first few, for now. need to expand this for larger P
        dplyr::filter(ir <= 3 & ic <= 4) |>
        ggplot(aes(x = freq_index, y = phase)) +
        ggplot2::facet_wrap(vars(ir, ic), scales = "free_y", 
                            strip.position = "left", 
                            labeller = \(x) label_both(x, multi_line = FALSE)) +
        geom_line(aes(linetype = datalabel, color = datalabel), 
                  linewidth = .2) +
        scale_color_manual(values = c("true" = "red", "multitaper" = "black", 
                                      "posterior" = "black")) +
        labs(title = paste0("Estimated phase for k = ", k)) +
        theme(legend.position = c(0, 1), legend.justification = c(0, 1),
              legend.background = element_rect(fill = alpha("white", 0.6)),
              legend.title = element_blank())
    
    print(plot_phase)
    save_plot_pdf(file.path(result_dir, "SDM-est-coherence-and-phase",
                            paste0("SDM-est-phases-k-", k, ".pdf")))
} # end of for loop

# plot likelihood of data -------------------------------------------------

calc_logLikeliS <- function (k) {
    logLikeliS <- rep(NA, num_samp)
    for (s in 1:num_samp) {
        lsum <- 0
        for (l in 1:num_freqs) {
            Omegakl <- diag(1 / ( 1/Lambdak_l_s[, k, l, s] + 1 ))
            invGamma <- 1/sigmak2_s[k, s] * 
                (diag(P) - U_kls_all[, , k, l, s] %*% 
                     Omegakl %*% t(Conj(U_kls_all[, , k, l, s])))
            tracekl <- Re(sum(t(invGamma) * data_list_w[[l]][[k]]))
            
            logdetData <- log(EigenR::Eigen_det(data_list_w[[l]][[k]]))
            logdetGamma <- 
                log(prod(sigmak2_s[k, s] * 
                             c(Lambdak_l_s[, k, l, s] + 1, rep(1, P-d))))
            
            lsum <- lsum + (LL-P)*logdetData - LL*logdetGamma - tracekl
        }
        logLikeliS[s] <- lsum
    }
    return(logLikeliS)
}

# this_logLike <- calc_logLikeliS(1)

if (useMclapply) {
    stopCluster(cluster)
}

cluster <- makeCluster(n_cores)
registerDoParallel(cluster)

# took a couple minutes for P = 4, d = 2, K = 10, 5000 Gibbs iterations
all_logLike <- mclapply(1:K, calc_logLikeliS, mc.cores = n_cores)
# all_logLike <- lapply(1:5, calc_logLikeliS)

total_logLike <- rep(0, num_samp)

for (s in 1:num_samp) {
    for (k in 1:K) {
        total_logLike[s] <- total_logLike[s] + Re(all_logLike[[k]][s])
    }
}

data.frame(index = 1:num_samp, logLike = total_logLike) |> 
    ggplot(aes(x = index, y = logLike)) + 
    geom_line() +
    labs(y = "log likelihood", title = "log likelihood trace plot",
         subtitle = paste0("P = ", P, ", d = ", d, ", K = ", K, ", T = ", Tt))

save_plot_pdf(file.path(result_dir,
                        paste0("log-likelihood-trace-plot.pdf")))
