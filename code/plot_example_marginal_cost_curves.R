library(here)
options(warn=2, scipen=10)
`%>%` <- magrittr::`%>%`

x = seq(0.75, 1, length.out=4000)
methane_kg_per_mcf <- 18.8916
CO2_SCC <- 100 # very rough
methane_scm_per_ton <- CO2_SCC * 29.8 # very rough
gas_frac_methane <- 0.95
ton_per_kg <- 1 / 1000
SCM_PER_MCF <- methane_scm_per_ton * ton_per_kg * methane_kg_per_mcf * gas_frac_methane
PRICE_PER_MCF <- 2.5

df <- purrr::map_dfr(c(0.005, 0.01, 0.05, 0.1, 0.2, 0.5),
  ~tibble::tibble(x = x, alpha = ., y= alpha / (1 - x))) %>%
  dplyr::mutate(alpha = factor(alpha))

plt <- ggplot2::ggplot(df, ggplot2::aes(x, y, color=alpha)) +
  ggplot2::geom_hline(yintercept=SCM_PER_MCF, linetype="dashed", alpha=0.7) +
  ggplot2::geom_hline(yintercept=PRICE_PER_MCF, linetype="dashed", alpha=0.7) +
  ggplot2::geom_line() +
  ggplot2::theme_bw() +
  ggplot2::scale_color_brewer(palette = "Blues") +
  ggplot2::labs(
    x="Prob no leak",
    y="Marginal cost ($/mcf)",
    color="α",
    title="Marginal costs for different α",
    subtitle = glue::glue("Approx price (${PRICE_PER_MCF}/mcf) and SMC (${CO2_SCC}/ton CO2e) in dashed lines")
  ) +
  ggplot2::ylim(NA, SCM_PER_MCF * 1.05)

filename <- here::here("graphics/example_marginal_cost_curves.pdf")
save_plot(plt, filename, reproducible=TRUE)
