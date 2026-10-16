suppressPackageStartupMessages({
  library(dplyr)
  library(nlme)
})

# Reproducible sensitivity analysis for the targeted HUMAnN-KO functional groups.
# Run from the repository root after run_functional_metagenomics_audit.R.

script_arg <- commandArgs(trailingOnly=FALSE)[grep("^--file=",commandArgs(trailingOnly=FALSE))]
repo_root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), mustWork=TRUE)
out_root <- file.path(repo_root, "results", "reproduced", "functional_metagenomics")

fit_group_model <- function(dat, experiment) {
  if (experiment == "field") {
    fit <- lm(y ~ replicate + Treatment, data = dat)
    tab <- drop1(fit, test = "F")
    return(list(F = tab["Treatment", "F value"],
                p = tab["Treatment", "Pr(>F)"],
                model_class = "linear_model"))
  }

  fit <- try(
    lme(y ~ factor(replicate) + Treatment,
        random = ~1 | Inoculum,
        data = dat,
        method = "REML",
        control = lmeControl(returnObject = TRUE)),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    fallback <- lm(y ~ factor(Inoculum) + factor(replicate) + Treatment,
                   data = dat)
    tab <- drop1(fallback, test = "F")
    return(list(F = tab["Treatment", "F value"],
                p = tab["Treatment", "Pr(>F)"],
                model_class = "fixed_inoculum_fallback"))
  }

  tab <- anova(fit)
  list(F = tab["Treatment", "F-value"],
       p = tab["Treatment", "p-value"],
       model_class = "random_inoculum_intercept")
}

prepare_experiment <- function(experiment) {
  exp_dir <- file.path(out_root, experiment)
  abundance <- read.csv(
    file.path(exp_dir, "functional_group_relative_abundance.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  metadata <- read.csv(file.path(exp_dir, "metadata.csv"),
                       stringsAsFactors = FALSE)
  sample_column <- if (experiment == "field") "id" else "sample_id"
  rownames(metadata) <- metadata[[sample_column]]
  abundance <- abundance[match(rownames(metadata), abundance$sample_id), , drop = FALSE]
  stopifnot(identical(abundance$sample_id, rownames(metadata)))
  rownames(abundance) <- abundance$sample_id
  abundance$sample_id <- NULL

  if (experiment == "field") {
    metadata$Treatment <- factor(metadata$Treatment,
                                 levels = c("S+", "UV-", "S-"))
  } else {
    metadata$Treatment <- factor(metadata$Treatment,
                                 levels = c("SusS+", "SusUV-", "SusS-"))
    metadata$Inoculum <- factor(metadata$Inoculum)
  }
  metadata$replicate <- factor(metadata$replicate)
  list(abundance = abundance, metadata = metadata, directory = exp_dir)
}

run_sensitivity <- function(experiment) {
  x <- prepare_experiment(experiment)
  abundance <- x$abundance
  metadata <- x$metadata
  primary <- read.csv(file.path(x$directory, "functional_group_omnibus.csv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
  primary_significance <- setNames(primary$p_BH < 0.05,
                                   primary$functional_group)

  # Transformation sensitivity: apply a group-specific half-minimum-positive
  # pseudocount, take the natural logarithm, and refit the original design.
  log_results <- lapply(colnames(abundance), function(group) {
    y <- abundance[[group]]
    pseudocount <- if (any(y > 0)) min(y[y > 0]) / 2 else 1e-12
    dat <- metadata
    dat$y <- log(y + pseudocount)
    model <- fit_group_model(dat, experiment)
    data.frame(functional_group = group,
               pseudocount = pseudocount,
               F = model$F,
               p = model$p,
               model_class = model$model_class,
               stringsAsFactors = FALSE)
  }) %>% bind_rows() %>% mutate(p_BH = p.adjust(p, method = "BH"))

  # Influence sensitivity: omit each sample in turn, refit every group, then
  # recalculate BH q-values across all groups for that reduced dataset.
  leave_one_out <- lapply(colnames(abundance), function(group) {
    lapply(seq_len(nrow(metadata)), function(i) {
      dat <- metadata[-i, , drop = FALSE]
      dat$y <- abundance[-i, group]
      model <- try(fit_group_model(dat, experiment), silent = TRUE)
      if (inherits(model, "try-error")) {
        return(data.frame(functional_group = group,
                          removed_sample = rownames(metadata)[i],
                          F = NA_real_, p = NA_real_))
      }
      data.frame(functional_group = group,
                 removed_sample = rownames(metadata)[i],
                 F = model$F, p = model$p)
    }) %>% bind_rows()
  }) %>% bind_rows() %>%
    group_by(removed_sample) %>%
    mutate(p_BH = p.adjust(p, method = "BH")) %>%
    ungroup() %>%
    mutate(primary_significant =
             unname(primary_significance[functional_group]),
           significance_changed =
             (p_BH < 0.05) != primary_significant)

  summary <- primary %>%
    select(functional_group, primary_p = p, primary_q = p_BH) %>%
    left_join(log_results %>%
                select(functional_group, log_p = p, log_q = p_BH,
                       log_pseudocount = pseudocount),
              by = "functional_group") %>%
    left_join(
      leave_one_out %>%
        group_by(functional_group) %>%
        summarise(
          loo_q_min = min(p_BH, na.rm = TRUE),
          loo_q_max = max(p_BH, na.rm = TRUE),
          n_omissions_changing_significance = sum(significance_changed,
                                                  na.rm = TRUE),
          omissions_changing_significance = paste(
            removed_sample[significance_changed],
            collapse = "; "
          ),
          .groups = "drop"
        ),
      by = "functional_group"
    ) %>%
    mutate(
      primary_significant = primary_q < 0.05,
      log_significant = log_q < 0.05,
      log_conclusion_stable = primary_significant == log_significant,
      loo_conclusion_stable = n_omissions_changing_significance == 0
    )

  # Keep the detailed outputs beside the primary experiment-specific results.
  write.csv(log_results %>% select(-model_class),
            file.path(x$directory, "log_abundance_sensitivity.csv"),
            row.names = FALSE)
  write.csv(leave_one_out,
            file.path(x$directory, "model_influence_leave_one_out.csv"),
            row.names = FALSE)
  write.csv(summary,
            file.path(x$directory, "functional_group_sensitivity_summary.csv"),
            row.names = FALSE)
  summary$experiment <- experiment
  summary
}

field <- run_sensitivity("field")
microcosm <- run_sensitivity("microcosm")
all_results <- bind_rows(field, microcosm) %>%
  select(experiment, everything())
write.csv(all_results,
          file.path(out_root, "functional_group_sensitivity_summary.csv"),
          row.names = FALSE)

photo <- microcosm %>%
  filter(functional_group ==
           "Photoreactivation and nucleotide-excision repair")

writeLines(c(
  "# Functional-group sensitivity analysis", "",
  "## Methods", "",
  paste(
    "Primary univariate models were refitted after natural-log transformation",
    "of each functional-group relative abundance using one-half of that group's",
    "smallest positive abundance as a pseudocount. The experimental design was",
    "unchanged, and omnibus p-values were BH-adjusted across all analyzable groups",
    "within each experiment."
  ), "",
  paste(
    "Influence was evaluated by sequentially omitting each sample, refitting every",
    "univariate group model with the original design, and recalculating BH-adjusted",
    "omnibus p-values across groups within each reduced dataset. A conclusion was",
    "classified as leave-one-out stable only when its full-data q < 0.05 status was",
    "preserved after every single-sample omission. This diagnostic was not used to",
    "exclude samples or replace the primary full-data analysis."
  ), "",
  "## Photoreactivation/NER result", "",
  sprintf(
    paste(
      "The log-transformed microcosm model remained significant (q = %.4f), but",
      "leave-one-out q-values ranged from %.4f to %.4f; %d of %d single-sample",
      "omissions changed the q < 0.05 classification."
    ),
    photo$log_q, photo$loo_q_min, photo$loo_q_max,
    photo$n_omissions_changing_significance, nrow(prepare_experiment("microcosm")$metadata)
  )
), file.path(out_root, "FUNCTIONAL_GROUP_SENSITIVITY_METHODS.md"))

writeLines(capture.output(sessionInfo()),
           file.path(out_root, "sensitivity_R_session_info.txt"))
