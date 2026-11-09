# Set up #######################################################################

library(tidyverse)
library(ggplot2)
library(patchwork)

make_circadian_data <- function(data) {
  circadian_data <- data
  
  CT_ZT0 <- data[data$Diet == "CT" & data$ZT == 0, ]
  if (nrow(CT_ZT0) > 0) {
    CT_ZT24 <- CT_ZT0
    CT_ZT24$ZT <- 24
    circadian_data <- rbind(circadian_data, CT_ZT24)
  }
  
  DF_ZT0 <- data[data$Diet == "DF" & data$ZT == 0, ]
  if (nrow(DF_ZT0) > 0) {
    DF_ZT24 <- DF_ZT0
    DF_ZT24$ZT <- 24
    circadian_data <- rbind(circadian_data, DF_ZT24)
  }
  
  return(circadian_data)
}

DosageTable <- read.table(
  "~/circanutr/data/Dosage.csv",
  header = TRUE,
  sep = ";",
  dec = ","
) |>
  separate(
    Condition,
    into = c("Diet", "ZT"),
    sep = "\\.ZT",
    remove = TRUE
  ) |>
  transform(ZT = as.numeric(ZT), Diet = as.factor(Diet))

DosageTable.mean <- data.frame(
  DosageTable |>
    group_by(Diet, ZT) |>
    summarise(
      AGL.mean = mean(AGL, na.rm = TRUE),
      Albumine.mean = mean(Albumine, na.rm = TRUE),
      T_bilirubine.mean = mean(T.bilirubine, na.rm = TRUE),
      Gluc_PAP.mean = mean(Gluc.PAP, na.rm = TRUE),
      Lact.mean = mean(Lact, na.rm = TRUE),
      Trigly.mean = mean(Trigly, na.rm = TRUE),
      Insu.mean = mean(Insuline, na.rm = TRUE)
    ) |>
    make_circadian_data()
)

DosageTable.mean <- DosageTable.mean[order(DosageTable.mean$ZT), ]
bioch.var <- names(DosageTable.mean)[-(1:2)]



# Venn AL rhythmic genes #######################################################
tpm <- read.table(file = "~/circanutr/data/salmon.merged.gene_tpm.tsv", dec = ".", header = T)
rownames(tpm) <- tpm[, 1]   # Row names are gene_id
tpm <- tpm[, -1]

gn_total <- as.vector(rownames(tpm))
# A vector of all genes names from the ENSEMBL GRCm39 Mus musculus genome
# This is the "Mouse" set

OneCount <- read.table(file = "~/circanutr/data/OneCountTable.csv", sep = ";")
OneCount <- OneCount[-1, ]
# Number of mouse genes detected in "liver"Mouse" with >= 1 count in AdLib

FilterByCount <- data.frame(t(rbind(gn_total, OneCount$V2)))
FilterByCount$V2 <- as.logical(FilterByCount$V2)

Liver <- subset(FilterByCount, V2)
# This is the "Liver" set

liver_detected <- as.matrix(
  read.csv2(
    "~/circanutr/data/FilterCountTable.csv",
    stringsAsFactors = FALSE,
    row.names = 1
  )
)
# Number of mouse genes detected in the liver with >= 10 Counts in AdLib
# This is the "Detected" set

plot(
  eulerr::euler(
    list(
      `Mouse` = gn_total,
      `Liver` = Liver$gn_total,
      `Detected` = rownames(liver_detected)
    ),
    shape = "ellipse",
  ),
  quantities = T,
  fills = c("#828282", "#D8D8D8", "#FFFFFF")
)
# Mouse = 78298 = 45295+12657+11481+8865  (GRCm39)
# Liver = 33003 = 12657+11481+8865        (>=1 Count in AdLib)
# Detected = 20346 = 11481+8865           (>=10 Counts in AdLib)



# Modeling #####################################################################

## CT ##########################################################################

CT.list.models <- lapply(bioch.var, function(y_var) {
  formulas <- list(
    LM1 = as.formula(paste(
      y_var, "~ cos(2*pi*ZT/24) + sin(2*pi*ZT/24)"
    )),
    LM2 = as.formula(paste(y_var, "~ 1")),
    LM3 = as.formula(paste(y_var, "~ ZT")),
    LM4 = as.formula(paste(y_var, "~ exp(ZT)"))
  )
  
  names(formulas) <- c("LM1", "LM2", "LM3", "LM4")
  
  # Fit the models
  lapply(formulas, function(form) {
    lm(formula = form, data = DosageTable.mean[DosageTable.mean$Diet == "CT", ])
  })
})

names(CT.list.models) <- bioch.var

CT.list.BIC <- list()

for (y_var in bioch.var) {
  # Set the three models for this response variable
  model1 <- CT.list.models[[y_var]][["LM1"]]
  model2 <- CT.list.models[[y_var]][["LM2"]]
  model3 <- CT.list.models[[y_var]][["LM3"]]
  model4 <- CT.list.models[[y_var]][["LM4"]]
  
  # Get BIC for each model
  BIC_values <- BIC(model1, model2, model3, model4)
  
  # Get the results with the variable name
  CT.list.BIC[[y_var]] <- list(
    LM1 = BIC_values[1, 2],
    LM2 = BIC_values[2, 2],
    LM3 = BIC_values[3, 2],
    LM4 = BIC_values[4, 2]
  )
}
rm(BIC_values, model1, model2, model3, model4, y_var)

CT.lowest.BIC.summary <- data.frame(
  Diet = "CT",
  Variable = names(CT.list.BIC),
  Best.model = sapply(CT.list.BIC, function(x)
    names(x)[which.min(unlist(x))]),
  formula = sapply(CT.list.BIC, function(x) {
    best_model <- names(x)[which.min(unlist(x))]
    switch(
      best_model,
      "LM1" = "y ~cos(2*pi*x/24) + sin(2*pi*x/24)",
      "LM2" = "y ~ 1",
      "LM3" = "y ~ x",
      "LM4" = "y ~ exp(x)",
      NA_character_
    )
  }),
  Min.BIC = sapply(CT.list.BIC, function(x)
    min(unlist(x))),
  row.names = NULL
)


## DF ##########################################################################

DF.list.models <- lapply(bioch.var, function(y_var) {
  formulas <- list(
    LM1 = as.formula(paste(
      y_var, "~ cos(2*pi*ZT/24) + sin(2*pi*ZT/24)"
    )),
    LM2 = as.formula(paste(y_var, "~ 1")),
    LM3 = as.formula(paste(y_var, "~ ZT")),
    LM4 = as.formula(paste(y_var, "~ exp(ZT)"))
  )
  
  names(formulas) <- c("LM1", "LM2", "LM3", "LM4")
  
  # Fit the models
  lapply(formulas, function(form) {
    lm(formula = form, data = DosageTable.mean[DosageTable.mean$Diet == "DF", ])
  })
})

names(DF.list.models) <- bioch.var

DF.list.BIC <- list()

for (y_var in bioch.var) {
  # Set the three models for this response variable
  model1 <- DF.list.models[[y_var]][["LM1"]]
  model2 <- DF.list.models[[y_var]][["LM2"]]
  model3 <- DF.list.models[[y_var]][["LM3"]]
  model4 <- DF.list.models[[y_var]][["LM4"]]
  
  # Get BIC for each model
  BIC_values <- BIC(model1, model2, model3, model4)
  
  # Get the results with the variable name
  DF.list.BIC[[y_var]] <- list(
    LM1 = BIC_values[1, 2],
    LM2 = BIC_values[2, 2],
    LM3 = BIC_values[3, 2],
    LM4 = BIC_values[4, 2]
  )
}
rm(BIC_values, model1, model2, model3, model4, y_var)

DF.lowest.BIC.summary <- data.frame(
  Diet = "DF",
  Variable = names(DF.list.BIC),
  Best.model = sapply(DF.list.BIC, function(x)
    names(x)[which.min(unlist(x))]),
  formula = sapply(DF.list.BIC, function(x) {
    best_model <- names(x)[which.min(unlist(x))]
    switch(
      best_model,
      "LM1" = "y ~cos(2*pi*x/24) + sin(2*pi*x/24)",
      "LM2" = "y ~ 1",
      "LM3" = "y ~ x",
      "LM4" = "y ~ exp(x)",
      NA_character_
    )
  }),
  Min.BIC = sapply(DF.list.BIC, function(x)
    min(unlist(x))),
  row.names = NULL
)

## FA ##########################################################################

FA.DosageTable.mean <- DosageTable.mean[DosageTable.mean$Diet == "FA", ] |>
  rbind(DosageTable.mean[DosageTable.mean$Diet == "CT" &
                           DosageTable.mean$ZT == 0, ]) |>
  mutate((across('Diet', str_replace, 'CT', 'FA')))

FA.DosageTable.mean <- FA.DosageTable.mean[order(FA.DosageTable.mean$ZT), ]

FA.list.models <- lapply(bioch.var, function(y_var) {
  formulas <- list(
    LM1 = as.formula(paste(
      y_var, "~ cos(2*pi*ZT/24) + sin(2*pi*ZT/24)"
    )),
    LM2 = as.formula(paste(y_var, "~ 1")),
    LM3 = as.formula(paste(y_var, "~ ZT")),
    LM4 = as.formula(paste(y_var, "~ exp(ZT)"))
  )
  
  names(formulas) <- c("LM1", "LM2", "LM3", "LM4")
  
  # Fit the models
  lapply(formulas, function(form) {
    lm(formula = form, data = FA.DosageTable.mean)
  })
})

names(FA.list.models) <- bioch.var

FA.list.BIC <- list()

for (y_var in bioch.var) {
  # Set the three models for this response variable
  model1 <- FA.list.models[[y_var]][["LM1"]]
  model2 <- FA.list.models[[y_var]][["LM2"]]
  model3 <- FA.list.models[[y_var]][["LM3"]]
  model4 <- FA.list.models[[y_var]][["LM4"]]
  
  
  # Get BIC for each model
  BIC_values <- BIC(model1, model2, model3, model4)
  
  # Get the results with the variable name
  FA.list.BIC[[y_var]] <- list(
    LM1 = BIC_values[1, 2],
    LM2 = BIC_values[2, 2],
    LM3 = BIC_values[3, 2],
    LM4 = BIC_values[4, 2]
  )
}
rm(BIC_values, model1, model2, model3, model4, y_var)

FA.lowest.BIC.summary <- data.frame(
  Diet = "FA",
  Variable = names(FA.list.BIC),
  Best.model = sapply(FA.list.BIC, function(x)
    names(x)[which.min(unlist(x))]),
  formula = sapply(FA.list.BIC, function(x) {
    best_model <- names(x)[which.min(unlist(x))]
    switch(
      best_model,
      "LM1" = "y ~cos(2*pi*x/24) + sin(2*pi*x/24)",
      "LM2" = "y ~ 1",
      "LM3" = "y ~ x",
      "LM4" = "y ~ exp(x)",
      NA_character_
    )
  }),
  Min.BIC = sapply(FA.list.BIC, function(x)
    min(unlist(x))),
  row.names = NULL
)

## FR ##########################################################################

FR.DosageTable.mean <- DosageTable.mean[DosageTable.mean$Diet == "FR", ] |>
  rbind(DosageTable.mean[DosageTable.mean$Diet == "CT" &
                           DosageTable.mean$ZT == 0, ]) |>
  mutate((across('Diet', str_replace, 'CT', 'FR')))

FR.DosageTable.mean <- FR.DosageTable.mean[order(FR.DosageTable.mean$ZT), ]

FR.list.models <- lapply(bioch.var, function(y_var) {
  formulas <- list(
    LM1 = as.formula(paste(
      y_var, "~ cos(2*pi*ZT/24) + sin(2*pi*ZT/24)"
    )),
    LM2 = as.formula(paste(y_var, "~ 1")),
    LM3 = as.formula(paste(y_var, "~ ZT")),
    LM4 = as.formula(paste(y_var, "~ exp(ZT)"))
  )
  
  names(formulas) <- c("LM1", "LM2", "LM3", "LM4")
  
  # Fit the models
  lapply(formulas, function(form) {
    lm(formula = form, data = FR.DosageTable.mean)
  })
})

names(FR.list.models) <- bioch.var

FR.list.BIC <- list()

for (y_var in bioch.var) {
  # Set the three models for this response variable
  model1 <- FR.list.models[[y_var]][["LM1"]]
  model2 <- FR.list.models[[y_var]][["LM2"]]
  model3 <- FR.list.models[[y_var]][["LM3"]]
  model4 <- FR.list.models[[y_var]][["LM4"]]
  
  # Get BIC for each model
  BIC_values <- BIC(model1, model2, model3, model4)
  
  # Get the results with the variable name
  FR.list.BIC[[y_var]] <- list(
    LM1 = BIC_values[1, 2],
    LM2 = BIC_values[2, 2],
    LM3 = BIC_values[3, 2],
    LM4 = BIC_values[4, 2]
  )
}
rm(BIC_values, model1, model2, model3, model4, y_var)

FR.lowest.BIC.summary <- data.frame(
  Diet = "FR",
  Variable = names(FR.list.BIC),
  Best.model = sapply(FR.list.BIC, function(x)
    names(x)[which.min(unlist(x))]),
  formula = sapply(FR.list.BIC, function(x) {
    best_model <- names(x)[which.min(unlist(x))]
    switch(
      best_model,
      "LM1" = "y ~cos(2*pi*x/24) + sin(2*pi*x/24)",
      "LM2" = "y ~ 1",
      "LM3" = "y ~ x",
      "LM4" = "y ~ exp(x)",
      NA_character_
    )
  }),
  Min.BIC = sapply(FR.list.BIC, function(x)
    min(unlist(x))),
  row.names = NULL
)


# Graph data ###################################################################

# Combine mean values after accounting for non-circadian diets
Adjusted.DosageTable.mean <- rbind(DosageTable.mean[DosageTable.mean$Diet == "CT", ],
                                   FR.DosageTable.mean,
                                   FA.DosageTable.mean,
                                   DosageTable.mean[DosageTable.mean$Diet == "DF", ]) |>
  pivot_longer(!c(Diet, ZT), names_to = "var_name", values_to = "Concentration") |>
  rename(Variable = var_name) |>
  left_join(CT.lowest.BIC.summary, by = join_by(Diet, Variable)) |>
  rows_update(FR.lowest.BIC.summary, by = c("Diet", "Variable")) |>
  rows_update(FA.lowest.BIC.summary, by = c("Diet", "Variable")) |>
  rows_update(DF.lowest.BIC.summary, by = c("Diet", "Variable"))


var_results <- list()

for (var in bioch.var) {
  # Extract rows for this variable from each diet
  var_results[[var]] <- rbind(
    CT.lowest.BIC.summary[CT.lowest.BIC.summary$Variable == var, ],
    FR.lowest.BIC.summary[FR.lowest.BIC.summary$Variable == var, ],
    FA.lowest.BIC.summary[FA.lowest.BIC.summary$Variable == var, ],
    DF.lowest.BIC.summary[DF.lowest.BIC.summary$Variable == var, ]
  )
  
  # Add a column to indicate which summary came from which diet
  var_results[[var]]$Diet <- c("CT", "FR", "FA", "DF")
}

# unused
SEM.table <- cbind(
  DosageTable |>
    group_by(Diet, ZT) |>
    summarise(
      AGL.SEM = sd(AGL, na.rm = T) / sqrt(n() - 1),
      Albumine.SEM = sd(Albumine, na.rm = T) / sqrt(n() - 1),
      T_bilirubine.SEM = sd(T.bilirubine, na.rm = T) / sqrt(n() - 1),
      Gluc_PAP.SEM = sd(Gluc.PAP, na.rm = T) / sqrt(n() - 1),
      Lact.SEM = sd(Lact, na.rm = T) / sqrt(n() - 1),
      Trigly.SEM = sd(Trigly, na.rm = T) / sqrt(n() - 1),
      Insu.SEM = sd(Insuline, na.rm = T) / sqrt(n() - 1)
    )
)


## Graph per variables #########################################################

### Glycémie ###################################################################
Adjusted.DosageTable.mean |>
  filter(Variable == "Gluc_PAP.mean") |>
  arrange(ZT) |>
  mutate(Diet = factor(Diet, levels = c("CT", "FR", "FA", "DF"))) |>
  ggplot(aes(
    x = ZT,
    y = Concentration,
    color = Diet,
    group = Diet
  )) +
  stat_smooth(
    data = . %>%
      filter(Diet == "CT"),
    formula = as.formula(subset(var_results[["Gluc_PAP.mean"]], Diet == "CT")$formula),
    method = "lm",
    se = FALSE,
    aes(fill = Diet, group = Diet)
  ) +
  stat_smooth(
    data = . %>%
      filter(Diet == "FR"),
    formula = as.formula(subset(var_results[["Gluc_PAP.mean"]], Diet == "FR")$formula),
    method = "lm",
    se = FALSE,
    aes(fill = Diet, group = Diet)
  )  +
  stat_smooth(
    data = . %>%
      filter(Diet == "FA"),
    formula = as.formula(subset(var_results[["Gluc_PAP.mean"]], Diet == "FA")$formula),
    method = "lm",
    se = FALSE,
    aes(fill = Diet, group = Diet)
  ) +
  stat_smooth(
    data = . %>%
      filter(Diet == "DF"),
    formula = as.formula(subset(var_results[["Gluc_PAP.mean"]], Diet == "DF")$formula),
    method = "lm",
    se = FALSE,
    aes(fill = Diet, group = Diet)
  ) +
  geom_point(aes(shape = Diet, fill = Diet)) +
  scale_color_manual(values = c(
    "CT" = "#67987E",
    "FR" = "#007FFF",
    "FA" = "#E0115F",
    "DF" = "#FF7F00"
  )) +
  scale_shape_manual(values = c(
    "CT" = 19,
    "FR" = 19,
    "FA" = 19,
    "DF" = 19
  )) +
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, by = 4)) +
  theme_classic() +
  theme(
    text = element_text(size = 8),
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    )
  ) +
  facet_wrap( ~ Diet,
              ncol = 4,
              scales = "fixed",
              axes = "all_y")

ggsave(
  filename = "Glycemie.svg",
  path = "~/circanutr/output/",
  width = 15,
  height = 3,
  dpi = 600,
  units = "cm"
)

### Insuline####################################################################
Adjusted.DosageTable.mean |>
  filter(Variable == "Insu.mean") |>
  arrange(ZT) |>
  mutate(Diet = factor(Diet, levels = c("CT", "FR", "FA", "DF"))) |>
  ggplot(aes(
    x = ZT,
    y = Concentration,
    color = Diet,
    group = Diet
  )) +
  stat_smooth(
    data = . %>%
      filter(Diet == "CT"),
    formula = as.formula(subset(var_results[["Insu.mean"]], Diet == "CT")$formula),
    method = "lm",
    se = FALSE,
    aes(fill = Diet, group = Diet)
  ) +
  stat_smooth(
    data = . %>%
      filter(Diet == "FR"),
    formula = as.formula(subset(var_results[["Insu.mean"]], Diet == "FR")$formula),
    method = "lm",
    se = FALSE,
    aes(fill = Diet, group = Diet)
  )  +
  stat_smooth(
    data = . %>%
      filter(Diet == "FA"),
    formula = as.formula(subset(var_results[["Insu.mean"]], Diet == "FA")$formula),
    method = "lm",
    se = FALSE,
    aes(fill = Diet, group = Diet)
  ) +
  stat_smooth(
    data = . %>%
      filter(Diet == "DF"),
    formula = as.formula(subset(var_results[["Insu.mean"]], Diet == "DF")$formula),
    method = "lm",
    se = FALSE,
    aes(fill = Diet, group = Diet)
  ) +
  geom_point(aes(shape = Diet, fill = Diet)) +
  scale_color_manual(values = c(
    "CT" = "#67987E",
    "FR" = "#007FFF",
    "FA" = "#E0115F",
    "DF" = "#FF7F00"
  )) +
  scale_shape_manual(values = c(
    "CT" = 19,
    "FR" = 19,
    "FA" = 19,
    "DF" = 19
  )) +
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, by = 4)) +
  theme_classic() +
  theme(
    text = element_text(size = 8),
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    )
  ) +
  facet_wrap( ~ Diet,
              ncol = 4,
              scales = "fixed",
              axes = "all_y")

ggsave(
  filename = "Insuline.svg",
  path = "~/circanutr/output/",
  width = 15,
  height = 3,
  dpi = 600,
  units = "cm"
)


### FFA ########################################################################
Adjusted.DosageTable.mean |>
  filter(Variable == "AGL.mean") |>
  arrange(ZT) |>
  mutate(Diet = factor(Diet, levels = c("CT", "FR", "FA", "DF"))) |>
  ggplot(aes(
    x = ZT,
    y = Concentration,
    color = Diet,
    group = Diet
  )) +
  stat_smooth(
    data = . %>%
      filter(Diet == "CT"),
    formula = as.formula(subset(var_results[["AGL.mean"]], Diet == "CT")$formula),
    method = "lm",
    se = FALSE,
    aes(fill = Diet, group = Diet)
  ) +
  stat_smooth(
    data = . %>%
      filter(Diet == "FR"),
    formula = as.formula(subset(var_results[["AGL.mean"]], Diet == "FR")$formula),
    method = "lm",
    se = FALSE,
    aes(fill = Diet, group = Diet)
  )  +
  stat_smooth(
    data = . %>%
      filter(Diet == "FA"),
    formula = as.formula(subset(var_results[["AGL.mean"]], Diet == "FA")$formula),
    method = "lm",
    se = FALSE,
    aes(fill = Diet, group = Diet)
  ) +
  stat_smooth(
    data = . %>%
      filter(Diet == "DF"),
    formula = as.formula(subset(var_results[["AGL.mean"]], Diet == "DF")$formula),
    method = "lm",
    se = FALSE,
    aes(fill = Diet, group = Diet)
  ) +
  geom_point(aes(shape = Diet, fill = Diet)) +
  scale_color_manual(values = c(
    "CT" = "#67987E",
    "FR" = "#007FFF",
    "FA" = "#E0115F",
    "DF" = "#FF7F00"
  )) +
  scale_shape_manual(values = c(
    "CT" = 19,
    "FR" = 19,
    "FA" = 19,
    "DF" = 19
  )) +
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, by = 4)) +
  ylab(label = "mmol/L") +
  theme_classic() +
  theme(
    text = element_text(size = 8),
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    )
  ) +
  facet_wrap( ~ Diet,
              ncol = 4,
              scales = "fixed",
              axes = "all_y")

ggsave(
  filename = "FFA.svg",
  path = "~/circanutr/output/",
  width = 15,
  height = 3,
  dpi = 600,
  units = "cm"
)


### TG #########################################################################
Adjusted.DosageTable.mean |>
  filter(Variable == "Trigly.mean") |>
  arrange(ZT) |>
  mutate(Diet = factor(Diet, levels = c("CT", "FR", "FA", "DF"))) |>
  ggplot(aes(
    x = ZT,
    y = Concentration,
    color = Diet,
    group = Diet
  )) +
  stat_smooth(
    data = . %>%
      filter(Diet == "CT"),
    formula = as.formula(subset(var_results[["Trigly.mean"]], Diet == "CT")$formula),
    method = "lm",
    se = FALSE,
    aes(fill = Diet, group = Diet)
  ) +
  stat_smooth(
    data = . %>%
      filter(Diet == "FR"),
    formula = as.formula(subset(var_results[["Trigly.mean"]], Diet == "FR")$formula),
    method = "lm",
    se = FALSE,
    aes(fill = Diet, group = Diet)
  )  +
  stat_smooth(
    data = . %>%
      filter(Diet == "FA"),
    formula = as.formula(subset(var_results[["Trigly.mean"]], Diet == "FA")$formula),
    method = "lm",
    se = FALSE,
    aes(fill = Diet, group = Diet)
  ) +
  stat_smooth(
    data = . %>%
      filter(Diet == "DF"),
    formula = as.formula(subset(var_results[["Trigly.mean"]], Diet == "DF")$formula),
    method = "lm",
    se = FALSE,
    aes(fill = Diet, group = Diet)
  ) +
  geom_point(aes(shape = Diet, fill = Diet)) +
  scale_color_manual(values = c(
    "CT" = "#67987E",
    "FR" = "#007FFF",
    "FA" = "#E0115F",
    "DF" = "#FF7F00"
  )) +
  scale_shape_manual(values = c(
    "CT" = 19,
    "FR" = 19,
    "FA" = 19,
    "DF" = 19
  )) +
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, by = 4)) +
  theme_classic() +
  theme(
    text = element_text(size = 8),
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    )
  ) +
  facet_wrap( ~ Diet,
              ncol = 4,
              scales = "fixed",
              axes = "all_y")

ggsave(
  filename = "Trigly.svg",
  path = "~/circanutr/output/",
  width = 15,
  height = 3,
  dpi = 600,
  units = "cm"
)

### Bilirubine #################################################################
Adjusted.DosageTable.mean |>
  filter(Variable == "T_bilirubine.mean") |>
  arrange(ZT) |>
  mutate(Diet = factor(Diet, levels = c("CT", "FR", "FA", "DF"))) |>
  ggplot(aes(
    x = ZT,
    y = Concentration,
    color = Diet,
    group = Diet
  )) +
  stat_smooth(
    data = . %>%
      filter(Diet == "CT"),
    formula = as.formula(subset(var_results[["T_bilirubine.mean"]], Diet == "CT")$formula),
    method = "lm",
    se = FALSE,
    aes(fill = Diet, group = Diet)
  ) +
  stat_smooth(
    data = . %>%
      filter(Diet == "FR"),
    formula = as.formula(subset(var_results[["T_bilirubine.mean"]], Diet == "FR")$formula),
    method = "lm",
    se = FALSE,
    aes(fill = Diet, group = Diet)
  )  +
  stat_smooth(
    data = . %>%
      filter(Diet == "FA"),
    formula = as.formula(subset(var_results[["T_bilirubine.mean"]], Diet == "FA")$formula),
    method = "lm",
    se = FALSE,
    aes(fill = Diet, group = Diet)
  ) +
  stat_smooth(
    data = . %>%
      filter(Diet == "DF"),
    formula = as.formula(subset(var_results[["T_bilirubine.mean"]], Diet == "DF")$formula),
    method = "lm",
    se = FALSE,
    aes(fill = Diet, group = Diet)
  ) +
  geom_point(aes(shape = Diet, fill = Diet)) +
  scale_color_manual(values = c(
    "CT" = "#67987E",
    "FR" = "#007FFF",
    "FA" = "#E0115F",
    "DF" = "#FF7F00"
  )) +
  scale_shape_manual(values = c(
    "CT" = 19,
    "FR" = 19,
    "FA" = 19,
    "DF" = 19
  )) +
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, by = 4)) +
  theme_classic() +
  theme(
    text = element_text(size = 8),
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    )
  ) +
  facet_wrap( ~ Diet,
              ncol = 4,
              scales = "fixed",
              axes = "all_y")

ggsave(
  filename = "Bilirubine.svg",
  path = "~/circanutr/output/",
  width = 15,
  height = 3,
  dpi = 600,
  units = "cm"
)



# qPCR Clock genes ############################################################
ClockTable <- read.table(
  "~/circanutr/data/qPCR_clock_genes.csv",
  header = TRUE,
  sep = ";",
  dec = ","
) |>
  transform(ZT = as.numeric(ZT), Diet = as.factor(Diet))
# Ready to use file, no need to rearragne stuff inside nor to compute mean, sd, sem


## Clock graph data ###########################################################
clock_theme <- list(
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, by = 4)),
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 5)),
  theme_classic(),
  theme(
    text = element_text(size = 10),
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    )
  )
)


### CT ########################################################################
g.CT.clock <- ClockTable |>
  filter(ClockTable$Diet == "CT") |>
  ggplot(aes(x = ZT)) +
  geom_point(aes(y = Bmal1.MEAN),
             shape = 16,
             color = "#67987E") +
  geom_point(aes(y = Nr1d1.MEAN),
             shape = 1,
             color = "#67987E") +
  stat_smooth(
    aes(y = Bmal1.MEAN),
    formula = as.formula("y ~cos(2*pi*x/24) + sin(2*pi*x/24)"),
    method = "lm",
    se = FALSE,
    lty = 1,
    color = "#67987E"
  ) +
  stat_smooth(
    aes(y = Nr1d1.MEAN),
    formula = as.formula("y ~cos(2*pi*x/24) + sin(2*pi*x/24)"),
    method = "lm",
    se = FALSE,
    lty = 2,
    color = "#67987E"
  ) +
  clock_theme

### NF ########################################################################
g.NF.clock <- ClockTable |>
  filter(ClockTable$Diet == "NF") |>
  ggplot(aes(x = ZT)) +
  geom_point(aes(y = Bmal1.MEAN),
             shape = 16,
             color = "#007FFF") +
  geom_point(aes(y = Nr1d1.MEAN),
             shape = 1,
             color = "#007FFF") +
  stat_smooth(
    aes(y = Bmal1.MEAN),
    formula = as.formula("y ~cos(2*pi*x/24) + sin(2*pi*x/24)"),
    method = "lm",
    se = FALSE,
    lty = 1,
    color = "#007FFF"
  ) +
  stat_smooth(
    aes(y = Nr1d1.MEAN),
    formula = as.formula("y ~cos(2*pi*x/24) + sin(2*pi*x/24)"),
    method = "lm",
    se = FALSE,
    lty = 2,
    color = "#007FFF"
  ) +
  clock_theme

### FA ########################################################################
g.FA.clock <- ClockTable |>
  filter(ClockTable$Diet == "FA") |>
  ggplot(aes(x = ZT)) +
  geom_point(aes(y = Bmal1.MEAN),
             shape = 16,
             color = "#E0115F") +
  geom_point(aes(y = Nr1d1.MEAN),
             shape = 1,
             color = "#E0115F") +
  stat_smooth(
    aes(y = Bmal1.MEAN),
    formula = as.formula("y ~cos(2*pi*x/24) + sin(2*pi*x/24)"),
    method = "lm",
    se = FALSE,
    lty = 1,
    color = "#E0115F"
  ) +
  stat_smooth(
    aes(y = Nr1d1.MEAN),
    formula = as.formula("y ~cos(2*pi*x/24) + sin(2*pi*x/24)"),
    method = "lm",
    se = FALSE,
    lty = 2,
    color = "#E0115F"
  ) +
  clock_theme

### DF ########################################################################
g.DF.clock <- ClockTable |>
  filter(ClockTable$Diet == "DF") |>
  ggplot(aes(x = ZT)) +
  geom_point(aes(y = Bmal1.MEAN),
             shape = 16,
             color = "#FF7F00") +
  geom_point(aes(y = Nr1d1.MEAN),
             shape = 1,
             color = "#FF7F00") +
  stat_smooth(
    aes(y = Bmal1.MEAN),
    formula = as.formula("y ~cos(2*pi*x/24) + sin(2*pi*x/24)"),
    method = "lm",
    se = FALSE,
    lty = 1,
    color = "#FF7F00"
  ) +
  stat_smooth(
    aes(y = Nr1d1.MEAN),
    formula = as.formula("y ~cos(2*pi*x/24) + sin(2*pi*x/24)"),
    method = "lm",
    se = FALSE,
    lty = 2,
    color = "#FF7F00"
  ) +
  clock_theme


## Arrange plots #############################################################
ggarrange(g.CT.clock,
          g.NF.clock,
          g.FA.clock,
          g.DF.clock,
          ncol = 4,
          nrow = 1)

ggsave(
  filename = "clock_genes.svg",
  path = "~/circanutr/output/",
  width = 16,
  height = 2.2,
  dpi = 600,
  units = "cm"
)
