# if you are running this after the fact: be sure to load the RData file

catout <- function(printtext) {
    cat(paste0("\n", printtext, "\n"))
}

ess_prop_rows <- list()
ess_s_rows <- list()

# check run time and some MH acceptance rate summaries --------------------

print(endtime - starttime)
sampling_time_sec <- as.numeric(endtime - starttime, units = "secs")

# check the acceptance rates
# gibbsPostBurn <- round(gibbsIts * burnin):gibbsIts
gibbsPostBurn <- (gibbsIts*burnin / t_thin + 1):num_samp
num_post_samp <- length(gibbsPostBurn)
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

# evaluate Lambda curve estimates -----------------------------------------

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

# evaluate Lambda ESS -----------------------------------------------------

Lambda_ESS <- array(NA, c(d, K, num_freqs))

for (j in 1:d) {
    for (k in 1:K) {
        for (l in 1:num_freqs) {
            Lambda_ESS[j, k, l] <- posterior::ess_basic(Lambdak_l_s[j, k, l, gibbsPostBurn])
        }
    }
}

# these are d x num_freqs
max_Lambda_ESS <- apply(Lambda_ESS, c(1, 3), max)
med_Lambda_ESS <- apply(Lambda_ESS, c(1, 3), median)
min_Lambda_ESS <- apply(Lambda_ESS, c(1, 3), min)

for (j in 1:d) {
    plot(max_Lambda_ESS[j, ] / num_post_samp, type = "l", 
         main = paste0("min, median, max ESS prop. for Lambda, j = ", j),
         ylim = c(0, max(max_Lambda_ESS[j, ] / num_post_samp)),
         ylab = "ESS proportion")
    lines(med_Lambda_ESS[j, ] / num_post_samp, col = "green")
    lines(min_Lambda_ESS[j, ] / num_post_samp, col = "red")
    
    save_plot_pdf(file.path(result_dir,
                            paste0("Lambda-min-med-max-ESS-j-", j, ".pdf")))
}


# which.max(ds_to_true) |> arrayInd(.dim = dim(ds_to_true))
which.min(min_Lambda_ESS) |> arrayInd(.dim = dim(min_Lambda_ESS))
min(min_Lambda_ESS)

plot(Lambdak_l_s[2, 1, 28, ], type = "l")
plot(Lambdak_l_s[2, 2, 28, ], type = "l")

# calculate ESS proportion and per second summaries
ess_prop_rows[["Lambda"]] <- quantile(as.vector(Lambda_ESS) / num_post_samp,
         probs = c(0, .025, .25, .5, .75, .975, 1))
ess_s_rows[["Lambda"]] <- quantile(as.vector(Lambda_ESS) / sampling_time_sec, 
         probs = c(0, .025, .25, .5, .75, .975, 1))

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
      data.frame(distance = apply(ds_to_true, 2, median), datalabel = "median"),
      data.frame(distance = apply(ds_to_true, 2, max), datalabel = "max")) |> 
    # ggplot
    ggplot(aes(x = rep(1:num_freqs, times = 3), 
               y = distance, color = datalabel)) +
    geom_line() +
    labs(x = "freq. index", y = "axis Frob. dist.", 
         title = "max, median, and min dist. of post. mean Ukl to true Ukl0") +
    theme(legend.position = c(0, 1), legend.justification = c(0, 1),
          legend.background = element_rect(fill = alpha("white", 0.6)))

save_plot_pdf(file.path(result_dir, 
                        paste0("Ukl-dist-to-true-max-median-min.pdf")))

# calculate Ukl distances and plot ESS for all k and l --------------------

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

# plot of min, median, and max ESS at each frequency
# create a temporary data frame first
rbind(data.frame(ESS = apply(ESS_Ukl_dists, 2, min) / num_post_samp, 
                 datalabel = "min"),
      data.frame(ESS = apply(ESS_Ukl_dists, 2, median) / num_post_samp, 
                 datalabel = "median"),
      data.frame(ESS = apply(ESS_Ukl_dists, 2, max) / num_post_samp, 
                 datalabel = "max")
      ) |> 
    # ggplot
    ggplot(aes(x = rep(1:num_freqs, times = 3), 
               y = ESS, color = datalabel)) +
    geom_line() +
    labs(x = "freq. index", y = "ESS proportion", 
         title = "max, median, min ESS prop. of Ukl dist.") +
    theme(legend.position = c(0, 1), legend.justification = c(0, 1),
          legend.background = element_rect(fill = alpha("white", 0.6)))

save_plot_pdf(file.path(result_dir, 
                        paste0("Ukl-dist-summary-ESS-max-min.pdf")))

# calculate ESS proportion and per second summaries
ess_prop_rows[["Ukl_dist"]] <- 
    quantile(as.vector(ESS_Ukl_dists) / num_post_samp, 
             probs = c(0, .025, .25, .5, .75, .975, 1))
ess_s_rows[["Ukl_dist"]] <- 
        quantile(as.vector(ESS_Ukl_dists) / sampling_time_sec,
                 probs = c(0, .025, .25, .5, .75, .975, 1))

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

# evaluate sampling of individual elements of Ukl -------------------------

Id_Pd <- diag(1, P, d)

std_U_kls_all <- array(NA, c(P, d, K, num_freqs, num_samp))

# calculate a standardized Ukl for each sample
# the new Ukl is the closest to the d columns of PxP Identity
for (k in 1:K) {
    for (l in 1:num_freqs) {
        for (s in 1:num_samp) {
            newU <- evec_Frob_stat(
                Id_Pd, U_kls_all[, , k, l, s], returnMats = TRUE)
            std_U_kls_all[, , k, l, s] <- newU$Yopt
        }
    }
}

# the complex version ends up having a lot of values very close to 0
ESS_Ukl_real <- apply(Re(std_U_kls_all[, , , , gibbsPostBurn]), 
                      c(1, 2, 3, 4), posterior::ess_basic)
max_ESS_Ukl_real <- apply(ESS_Ukl_real, c(4), max)
min_ESS_Ukl_real <- apply(ESS_Ukl_real, c(4), min)
med_ESS_Ukl_real <- apply(ESS_Ukl_real, c(4), median)

# aggregate the ESS at each frequency (i.e. over P, d, K)
summ_ESS_Ukl_real <- data.frame(
    values = c(max_ESS_Ukl_real, min_ESS_Ukl_real, med_ESS_Ukl_real),
    stat = rep(c("max", "min", "median"), each = num_freqs),
    freq_index = rep(1:num_freqs, 3)
)

# plot ESS per frequency
summ_ESS_Ukl_real |> 
    ggplot(aes(x = freq_index, y = values/num_post_samp, color = stat)) +
    geom_line() +
    labs(x = "freq. index", y = "ESS proportion", 
         title = "Ukl (real values) ESS prop. aggregated over P, d, K")
save_plot_pdf(file.path(result_dir, 
                        paste0("post-Ukl-element-ESS-by-freq.pdf")))

catout("indices of Re(Ukl) with the worse ESS:")
ind_worst_ReUkl <- which.min(ESS_Ukl_real) |> arrayInd(.dim = dim(ESS_Ukl_real))
print(ind_worst_ReUkl)

worst_ReUkl <- Re(std_U_kls_all[ind_worst_ReUkl[1], ind_worst_ReUkl[2],
                            ind_worst_ReUkl[3], ind_worst_ReUkl[4], ])

# trace plot of the entry with worst ESS
plot(worst_ReUkl, type = "l",
     main = paste0("Trace plot of Re(Ukl) with (p, j, k, l) = (",
                   paste(ind_worst_ReUkl, collapse = ", " ), ")"),
     xlab = "freq. index", ylab = "ESS")
save_plot_pdf(file.path(result_dir, 
                        paste0("post-Ukl-element-worstESS-traceplot.pdf")))

# calculate ESS proportion and per second summaries
ess_prop_rows[["Ukl_real"]] <- 
    quantile(as.vector(ESS_Ukl_real) / num_post_samp, 
             probs = c(0, .025, .25, .5, .75, .975, 1))
ess_s_rows[["Ukl_real"]] <- 
        quantile(as.vector(ESS_Ukl_real) / sampling_time_sec,
                 probs = c(0, .025, .25, .5, .75, .975, 1))

# evaluating ESS of projection matrices -----------------------------------

Proj_kl_s <- array(NA, c(P, P, K, num_freqs, num_samp))

# calculate each sampled projection matrix
for (k in 1:K) {
    for (l in 1:num_freqs) {
        for (s in 1:num_samp) {
            Proj_kl_s[, , k, l, s] <- 
                U_kls_all[, , k, l, s] %*% t(Conj(U_kls_all[, , k, l, s]))
        }
    }
}

# ESS: for each k and l:
# - one Real for each diagonal entry
# - one Real for each off-diagonal entry
# - one Complex for each off-diagonal entry
Proj_ESS <- matrix(NA, K*num_freqs*(P + 2*choose(P, 2)), 5)
Proj_Re_or_Im <- rep(NA, K*num_freqs*(P + 2*choose(P, 2)))

plot(Re(Proj_kl_s[1, 1, 1, 1, ]), type = "l")

# calculate ESS for real and complex, diagonal and upper-triangular entries
rowi <- 1
for (k in 1:K) {
    for (l in 1:num_freqs) {
        for (ii in 1:P) {
            for (j in (ii):P) {
                # Complex
                if (j != ii) {
                    Proj_ESS[rowi, ] <- c(posterior::ess_basic(
                        Im(Proj_kl_s[ii, j, k, l, gibbsPostBurn])),
                        k, l, ii, j)
                    Proj_Re_or_Im[rowi] <- "Im"
                    rowi <- rowi + 1
                }
                
                # Real
                Proj_ESS[rowi, ] <- c(posterior::ess_basic(
                    Re(Proj_kl_s[ii, j, k, l, gibbsPostBurn])),
                    k, l, ii, j)
                Proj_Re_or_Im[rowi] <- "Re"
                rowi <- rowi + 1
            }
        }
    }
}

# put together into a data frame
Proj_ESS <- data.frame(Proj_ESS)
Proj_ESS$type <- Proj_Re_or_Im
names(Proj_ESS) <- c("ESS", "k", "l", "p", "j", "type")

# what are the indices of the entry with the worst ESS?
Proj_ESS |> dplyr::group_by(type) |> dplyr::slice_min(ESS, n = 1)

# save the aggregated min, med, max in long form for easier plotting
summ_ESS <- Proj_ESS |> 
    dplyr::summarise(min = min(ESS), median = median(ESS), max = max(ESS),
                     .by = c(l, type)) |> 
    tidyr::pivot_longer(cols = !c(l, type),
                        names_to = "stat",
                        values_to = "ESS")

# plot ESS by frequency for Real values
summ_ESS |> 
    dplyr::filter(type == "Re") |> 
    ggplot(aes(x = l, y = ESS / num_post_samp, color = stat)) +
    geom_line() +
    labs(title = "min, median, max ESS prop. for Real proj. matrix entries",
         x = "freq. index", y = "ESS proportion")

# plot ESS by frequency for Complex values
summ_ESS |> 
    dplyr::filter(type == "Im") |> 
    ggplot(aes(x = l, y = ESS / num_post_samp, color = stat)) +
    geom_line() +
    labs(title = "min, median, max ESS prop. for Complex proj. matrix entries",
         x = "freq. index", y = "ESS proportion")

# calculate ESS proportion and per second summaries
ess_prop_rows[["Proj"]] <- quantile(
    Proj_ESS$ESS / num_post_samp, 
    probs = c(0, .025, .25, .5, .75, .975, 1))
ess_s_rows[["Proj"]] <- quantile(
    Proj_ESS$ESS / sampling_time_sec, 
    probs = c(0, .025, .25, .5, .75, .975, 1))

# Finding worst combinations of ESS for Projection matrices ---------------

# what are the 25% worst combinations of entries?
# Proj_ESS_Re |>
#     dplyr::filter(ESS <= quantile(ESS, 0.25)) |>
#     dplyr::count(p, j) |> dplyr::arrange(desc(n))
# Proj_ESS_Im |>
#     dplyr::filter(ESS <= quantile(ESS, 0.25)) |>
#     dplyr::count(p, j) |> dplyr::arrange(desc(n))
# 
# # what are the 25% worst observations?
# # I think that this is maybe not as meaningful
# Proj_ESS_Re |>
#     dplyr::filter(ESS <= quantile(ESS, 0.25)) |>
#     dplyr::count(k) |> dplyr::arrange(desc(n))
# 
# # a trace plot of one projection matrix entry
# plot(Re(Proj_kl_s[3, 4, 10, 224, ]), type = "l")
# plot(Im(Proj_kl_s[3, 4, 8, 259, ]), type = "l")

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

# Evaluate Sigmal ESS -----------------------------------------------------

Sigmal_ESS <- matrix(NA, num_freqs*(P + 2*choose(P, 2)), 4)
Sigmal_Re_or_Im <- rep(NA, num_freqs*(P + 2*choose(P, 2)))

# calculate ESS for real and complex, diagonal and upper-triangular entries
rowi <- 1
for (l in 1:num_freqs) {
    for (ii in 1:P) {
        for (j in (ii):P) {
            if (j != ii) {
                Sigmal_ESS[rowi, ] <- c(
                    posterior::ess_basic(Im(Sigmal_s[ii, j, l, gibbsPostBurn])),
                    l, ii, j)
                Sigmal_Re_or_Im[rowi] <- "Im"
                rowi <- rowi + 1
            }
            
            Sigmal_ESS[rowi, ] <- c( 
                posterior::ess_basic(Re(Sigmal_s[ii, j, l, gibbsPostBurn])),
                l, ii, j)
            Sigmal_Re_or_Im[rowi] <- "Re"
            rowi <- rowi + 1
        }
    }
}

# put into data frame
Sigmal_ESS <- cbind(data.frame(Sigmal_ESS), Sigmal_Re_or_Im)
names(Sigmal_ESS) <- c("ESS", "l", "p", "j", "type")

Sigmal_ESS |> dplyr::group_by(type) |> dplyr::slice_min(ESS, n = 1)

# summarize min, max, median ESS in a long form for easier plotting
summ_Sigmal_ESS <- Sigmal_ESS |>  
    dplyr::summarise(min = min(ESS), median = median(ESS), max = max(ESS),
                     .by = c(l, type)) |>
    tidyr::pivot_longer(cols = !c(l, type),
                        names_to = "stat",
                        values_to = "ESS")

# plot ESS for Real values
summ_Sigmal_ESS |> 
    dplyr::filter(type == "Re") |> 
    ggplot(aes(x = l, y = ESS / num_post_samp, color = stat)) +
    geom_line() +
    labs(title = "ESS prop. for Re Sigmal aggregated over P, d",
         x = "freq. index", y = "ESS proportion")

# plot ESS for Complex values
summ_Sigmal_ESS |> 
    dplyr::filter(type == "Im") |> 
    ggplot(aes(x = l, y = ESS / num_post_samp, color = stat)) +
    geom_line() +
    labs(title = "ESS prop. for Im Sigmal aggregated over P, d",
         x = "freq. index", y = "ESS proportion")

# calculate ESS proportion and per second summaries
ess_prop_rows[["Sigmal"]] <- quantile(
    Sigmal_ESS$ESS / num_post_samp, 
    probs = c(0, .025, .25, .5, .75, .975, 1))
ess_s_rows[["Sigmal"]] <- quantile(
    Sigmal_ESS$ESS / sampling_time_sec,
    probs = c(0, .025, .25, .5, .75, .975, 1))

# smoothing parameter summaries -------------------------------------------

if (Lambda_method == "bspline") {
    catout("median of tau2 samples at some diff j and k")
    catout("(j, k) = (1, 1) and (2, 1)")
    median(tau2_bs_s[1, 1, gibbsPostBurn])
    median(tau2_bs_s[2, 1, gibbsPostBurn])
    catout("(j, k) = (1, 2) and (2, 2)")
    median(tau2_bs_s[1, 2, gibbsPostBurn])
    median(tau2_bs_s[2, 2, gibbsPostBurn])
} else {
    plot(taujkl2_s[1, 1, 1, gibbsPostBurn], type = "l")
    
    catout("median of taujkl2 samples at some diff j, k, and l")
    median(taujkl2_s[1, 1, 1, gibbsPostBurn])
    
    median(taujkl2_s[1, 1, 150, gibbsPostBurn])
    median(taujkl2_s[2, 1, 150, gibbsPostBurn])
    
    median(taujkl2_s[1, 1, 256, gibbsPostBurn])
    median(taujkl2_s[1, 1, 350, gibbsPostBurn])    
}

# evaluate sigmak2 trace and CIs ------------------------------------------

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

# evaluate sigmak2 ESS ----------------------------------------------------

sigmak2_ESS <- apply(sigmak2_s[, gibbsPostBurn], c(1), 
                     posterior::ess_basic)

# calculate ESS proportion and per second summaries
ess_prop_rows[["sigmak2"]] <- quantile(sigmak2_ESS / num_post_samp, 
         probs = c(0, .025, .25, .5, .75, .975, 1))
ess_s_rows[["sigmak2"]] <- quantile(sigmak2_ESS / sampling_time_sec, 
         probs = c(0, .025, .25, .5, .75, .975, 1))

# create table of ESS stats -----------------------------------------------

ess_prop_df <- do.call(rbind, ess_prop_rows)
ess_s_df <- do.call(rbind, ess_s_rows)

catout("ESS proportion summary for each parameter")
print(ess_prop_df)
catout("ESS/s summary for each parameter")
print(ess_s_df)

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
