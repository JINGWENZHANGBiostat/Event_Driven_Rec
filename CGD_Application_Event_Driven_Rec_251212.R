###############################################
# Application to the CGD dataset
# - Proposed event-driven design
# - Fixed design, target number of events as 39
###############################################

library(survival)
library(tidyr)
library(dplyr)
library(purrr)

################################################
# 1. Data preparation for CGD dataset
################################################

# Load CGD data
data(cgd, package = "survival")

# 1.1 Randomization dates and study time scale
cgd0$rand_str <- as.character(cgd0$random)
cgd0$rand_str <- ifelse(
        nchar(cgd0$rand_str) == 5,
        paste0("0", cgd0$rand_str),
        cgd0$rand_str
)

cgd0$random_date <- as.Date(cgd0$rand_str, format = "%m%d%y")
study_start      <- min(cgd0$random_date)

# Calendar time of randomization (days since study start)
cgd0$s_i <- as.numeric(cgd0$random_date - study_start)

# 1.2 End of follow-up and total follow-up time
cgd0$end_followup_date <- cgd0$random_date + cgd0$futime
cgd0$T_i <- as.numeric(cgd0$end_followup_date - study_start)

# 1.3 Convert wide event times to long format
long_events <- cgd0 %>%
        pivot_longer(
                cols      = starts_with("etime"),
                names_to  = "event",
                values_to = "time"
        ) %>%
        mutate(time = as.numeric(time)) %>%
        filter(!is.na(time)) %>%
        arrange(id, time)

# 1.4 List of relative event times per patient
event_times_rel <- long_events %>%
        group_by(id) %>%
        summarize(events_rel = list(sort(time)), .groups = "drop") %>%
        complete(id = cgd$id, fill = list(events_rel = list(numeric(0)))) %>%
        arrange(id) %>%
        select(id, events_rel)

cgd0 <- cgd0 %>%
        left_join(event_times_rel, by = "id")

# 1.5 Calendar-time event times
cgd0 <- cgd0 %>%
        mutate(events_cal = map2(s_i, events_rel, ~ .x + .y))

# 1.6 Final patient-level dataset for analysis
patients_df <- cgd0 %>%
        select(
                id, treat, futime,
                random_date, s_i,
                end_followup_date, T_i,
                events_rel, events_cal
        )

################################################
# Helper: build Andersen–Gill dataset
################################################

build_ag_data <- function(eligible_idx, t_monitor, patients_df) {
        
        analysis_data_list <- lapply(seq_along(eligible_idx), function(idx) {
                
                i <- eligible_idx[idx]
                
                # relative follow-up time at monitoring time t
                t_rel_end <- min(
                        t_monitor - patients_df$s_i[i],
                        patients_df$T_i[i] - patients_df$s_i[i]
                )
                
                ev_times <- patients_df$events_rel[[i]]
                ev_times <- ev_times[ev_times <= t_rel_end]
                
                # no events
                if (length(ev_times) == 0) {
                        return(data.frame(
                                id    = i,
                                start = 0,
                                stop  = t_rel_end,
                                event = 0,
                                x     = patients_df$treat[i]
                        ))
                }
                
                # events exactly at t_rel_end
                if (abs(t_rel_end - tail(ev_times, 1)) <= 1e-6) {
                        intervals   <- c(0, ev_times)
                        n_intervals <- length(intervals) - 1
                        
                        return(data.frame(
                                id    = rep(i, n_intervals),
                                start = head(intervals, -1),
                                stop  = tail(intervals, -1),
                                event = rep(1, length(ev_times)),
                                x     = rep(patients_df$treat[i], n_intervals)
                        ))
                }
                
                # events occur before t_rel_end
                intervals   <- c(0, ev_times, t_rel_end)
                n_intervals <- length(intervals) - 1
                
                data.frame(
                        id    = rep(i, n_intervals),
                        start = head(intervals, -1),
                        stop  = tail(intervals, -1),
                        event = c(rep(1, length(ev_times)), 0),
                        x     = rep(patients_df$treat[i], n_intervals)
                )
        })
        
        analysis_data <- do.call(rbind, analysis_data_list)
        analysis_data <- na.omit(analysis_data)
        analysis_data <- analysis_data[(analysis_data$stop - analysis_data$start) > 1e-6, ]
        
        analysis_data
}

################################################
# monitoring grid: common for the fixed design and the blind continuous monitoring procedure
################################################

study_period <- as.numeric(
        max(patients_df$end_followup_date) -
                min(patients_df$random_date)
)

monitor_start <- 7
monitor_times <- seq(monitor_start, study_period, by = 1)

################################################
# 2. Proposed procedure: blind continuous monitoring procedure
################################################

# 2.1 Target v^2 from desired power and effect size
target.power <- 0.8
z_alpha_2    <- qnorm(1 - 0.05 / 2)
beta_use     <- log(0.3)

find_v <- function(v) {
        target.power - (
                pnorm(-z_alpha_2 - beta_use / v) +
                        1 - pnorm(z_alpha_2 - beta_use / v)
        )
}

solution <- uniroot(find_v, c(0.01, 10))
v        <- solution$root
v2       <- v^2

# 2.2 Monitoring loop for proposed procedure

beta_est_prop   <- NA
beta_robse_prop <- NA
p_value_prop    <- NA
L_prop          <- NA
t_final_prop    <- NA
statistic_prop  <- NA
pred_power_prop <- NA

for (t in monitor_times) {
        
        eligible_idx <- which(patients_df$s_i <= t)
        if (length(eligible_idx) == 0) next
        
        # relative follow-up at time t
        t_rel_end_vec <- pmin(
                t - patients_df$s_i[eligible_idx],
                patients_df$T_i[eligible_idx] - patients_df$s_i[eligible_idx]
        )
        
        # pooled event times
        cum_events <- unlist(lapply(eligible_idx, function(i) {
                ev <- patients_df$events_rel[[i]]
                ev[ev <= min(
                        t - patients_df$s_i[i],
                        patients_df$T_i[i] - patients_df$s_i[i]
                )]
        }))
        
        L <- length(cum_events)
        if (L == 0) next
        
        # Nelson–Aalen estimator based on pooled data
        unique_times_rel <- sort(unique(cum_events))
        d_counts         <- as.numeric(table(cum_events))[order(unique_times_rel)]
        
        cum_hazards <- numeric(length(unique_times_rel))
        cum_hazard  <- 0
        
        for (j in seq_along(unique_times_rel)) {
                t_rel <- unique_times_rel[j]
                Y     <- sum(t_rel <= t_rel_end_vec)  # number at risk
                
                if (Y == 0) next
                
                cum_hazard  <- cum_hazard + d_counts[j] / Y
                cum_hazards[j] <- cum_hazard
        }
        
        haz_func_rel <- stepfun(unique_times_rel, c(0, cum_hazards))
        
        # Sum of squared residuals (observed - expected)^2
        sum_sq <- 0
        for (i in eligible_idx) {
                t_rel_end <- min(
                        t - patients_df$s_i[i],
                        patients_df$T_i[i] - patients_df$s_i[i]
                )
                observed <- sum(patients_df$events_rel[[i]] <= t_rel_end)
                expected <- haz_func_rel(t_rel_end)
                sum_sq   <- sum_sq + (observed - expected)^2
        }
        
        # Estimated v^2 from pooled data
        statistic <- (sum_sq * 4) / (L^2)
        
        # Decision rule for final analysis under proposed method
        if ((L > 0 && !is.na(statistic) && statistic > 0 && statistic <= v2) ||
            (L > 0 && t == study_period)) {
                
                t_final_prop   <- t
                L_prop         <- L
                statistic_prop <- statistic
                
                pred_power_prop <- pnorm(-z_alpha_2 - beta_use / sqrt(statistic)) +
                        1 - pnorm(z_alpha_2 - beta_use / sqrt(statistic))
                
                analysis_data <- build_ag_data(eligible_idx, t, patients_df)
                
                if (nrow(analysis_data) > 0) {
                        tryCatch({
                                fit <- coxph(
                                        Surv(start, stop, event) ~ x,
                                        data = analysis_data,
                                        cluster = id,
                                        robust  = TRUE
                                )
                                
                                fit_summary      <- summary(fit)
                                beta_est_prop    <- fit_summary$coefficients["x", "coef"]
                                beta_robse_prop  <- fit_summary$coefficients["x", "robust se"]
                                p_value_prop     <- fit_summary$waldtest["pvalue"]
                                
                        }, error = function(e) {
                                cat("Proposed method: Cox model failed:", conditionMessage(e), "\n")
                                p_value_prop <- NA
                        })
                } else {
                        cat("Proposed method: empty analysis_data, skipping Cox model\n")
                        p_value_prop <- NA
                }
                
                break
        }
}

################################################
# 3. Fixed design analysis
################################################

beta_est_fix   <- NA
beta_robse_fix <- NA
p_value_fix    <- NA
L_fix          <- NA
t_final_fix    <- NA

for (t in monitor_times) {
        
        eligible_idx <- which(patients_df$s_i <= t)
        if (length(eligible_idx) == 0) next
        
        # pooled event times
        cum_events <- unlist(lapply(eligible_idx, function(i) {
                ev <- patients_df$events_rel[[i]]
                ev[ev <= min(
                        t - patients_df$s_i[i],
                        patients_df$T_i[i] - patients_df$s_i[i]
                )]
        }))
        
        L <- length(cum_events)
        if (L == 0) next
        
        # decision: stop when total events reach the fixed target 39
        if (L >= 39 || (L > 0 && t == study_period)) {
                
                t_final_fix <- t
                L_fix       <- L
                
                analysis_data <- build_ag_data(eligible_idx, t, patients_df)
                
                if (nrow(analysis_data) > 0) {
                        tryCatch({
                                fit <- coxph(
                                        Surv(start, stop, event) ~ x,
                                        data = analysis_data,
                                        cluster = id,
                                        robust  = TRUE
                                )
                                
                                fit_summary     <- summary(fit)
                                beta_est_fix    <- fit_summary$coefficients["x", "coef"]
                                beta_robse_fix  <- fit_summary$coefficients["x", "robust se"]
                                p_value_fix     <- fit_summary$robscore["pvalue"]
                                
                        }, error = function(e) {
                                cat("Fixed design: Cox model failed:", conditionMessage(e), "\n")
                                p_value_fix <- NA
                        })
                } else {
                        cat("Fixed design: empty analysis_data, skipping Cox model\n")
                        p_value_fix <- NA
                }
                
                break
        }
}

################################################
# 4. Print results
################################################

cat("===== Proposed event-driven design =====\n")
cat("Monitoring time t       =", t_final_prop, "\n")
cat("Total events L          =", L_prop, "\n")
cat("Estimated v^2 statistic =", statistic_prop, "\n")
cat("Predicted power         =", pred_power_prop, "\n")
cat("beta                    =", beta_est_prop, "\n")
cat("robust SE               =", beta_robse_prop, "\n")
cat("p-value                 =", p_value_prop, "\n\n")

cat("===== Conventional fixed design =====\n")
cat("Monitoring time t       =", t_final_fix, "\n")
cat("Total events L          =", L_fix, "\n")
cat("beta                    =", beta_est_fix, "\n")
cat("robust SE               =", beta_robse_fix, "\n")
cat("p-value                 =", p_value_fix, "\n")



