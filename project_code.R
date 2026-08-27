# install packages
packages <- c("dplyr", "ggplot2", "lubridate", "readr", "tidyr", "purrr",
              "scales", "maps", "viridis","SAPP","mclust","splines")
for (p in packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)}}
# read data
data <- readLines("SearchResults.txt", warn = FALSE)
data <- data[grepl("^\\d{4}/\\d{2}/\\d{2}", data)]
catalog <- read.table(text = data,header = FALSE,stringsAsFactors = FALSE)
names(catalog) <- c("date", "time", "event_type", "geo_type","mag", "mag_type", 
                    "lat", "lon", "depth","quality", "evid", "nph", "ngrm")
catalog <- catalog %>%
  mutate(datetime = ymd_hms(paste(date, time)),mag = as.numeric(mag),
         lat = as.numeric(lat),lon = as.numeric(lon),depth = as.numeric(depth)) %>%
  filter(!is.na(mag), !is.na(datetime), !is.na(lat), !is.na(lon), !is.na(depth)) %>%
  arrange(datetime)
# check data and some exploratory analysis
cat("Number of events:", nrow(catalog), "\n")
cat("Time range:", as.character(min(catalog$datetime)), "to",
    as.character(max(catalog$datetime)), "\n")
cat("Magnitude range:", min(catalog$mag), "to", max(catalog$mag), "\n")

head(catalog)
summary(catalog$mag)
# ============================================================
# STEP 1: EDA and model comparison
# ============================================================
catalog_main <- catalog %>%
  mutate(date_only = as.Date(datetime),year = year(datetime),
         month = floor_date(datetime, "month"),year_month = format(datetime, "%Y-%m")) %>%
  arrange(datetime)
# Annual event counts
annual_counts <- catalog_main %>%
  group_by(year) %>%
  summarise(n_events = n(),max_mag = max(mag, na.rm = TRUE),
            mean_mag = mean(mag, na.rm = TRUE),median_mag = median(mag, na.rm = TRUE),
            .groups = "drop")

print(annual_counts)
write.csv(annual_counts,"annual_event_counts.csv",row.names = FALSE)

p_annual <- ggplot(annual_counts, aes(x = year, y = n_events)) +
  geom_col() +
  labs(x = "Year",y = "Number of events") +
  theme_minimal()

print(p_annual)
ggsave("annual_event_counts.png",p_annual,width = 8,height = 5,dpi = 300)

# Monthly event counts
monthly_counts <- catalog_main %>%
  group_by(month) %>%
  summarise(n_events = n(),max_mag = max(mag, na.rm = TRUE),
            mean_mag = mean(mag, na.rm = TRUE),
            .groups = "drop")

print(head(monthly_counts))
write.csv(monthly_counts,"monthly_event_counts.csv",row.names = FALSE)

p_monthly <- ggplot(monthly_counts, aes(x = month, y = n_events)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1) +
  labs(title = "Monthly Number of Earthquakes",x = "Time",y = "Number of events") +
  theme_minimal()

print(p_monthly)
ggsave("monthly_event_counts.png",p_monthly,width = 9,height = 5,dpi = 300)

# Cumulative event count through time
catalog_cum <- catalog_main %>%
  arrange(datetime) %>%
  mutate(cumulative_events = row_number())

p_cumulative <- ggplot(catalog_cum, aes(x = datetime, y = cumulative_events)) +
  geom_line(linewidth = 0.7) +
  labs(title = "Cumulative Number of Earthquakes over Time",
       subtitle = "Steeper segments indicate periods of increased seismic activity",
       x = "Time",y = "Cumulative number of events") +
  theme_minimal()

print(p_cumulative)
ggsave("cumulative_event_count.png",p_cumulative,width = 9,height = 5,dpi = 300)

# Daily Event Counts around the Largest Event
largest_event <- catalog_main %>%
  arrange(desc(mag)) %>%
  slice(1)
print(largest_event)
mainshock_time <- largest_event$datetime[1]
window_days <- 30
catalog_around_mainshock <- catalog_main %>%
  filter(datetime >= mainshock_time - days(window_days),
         datetime <= mainshock_time + days(window_days)) %>%
  mutate(days_from_mainshock = as.numeric(difftime(datetime, mainshock_time, units = "days")))
daily_around_mainshock <- catalog_around_mainshock %>%
  mutate(relative_day = floor(days_from_mainshock)) %>%
  group_by(relative_day) %>%
  summarise(n_events = n(), max_mag = max(mag, na.rm = TRUE),
            .groups = "drop")
p_around_daily <- ggplot(daily_around_mainshock, aes(x = relative_day, y = n_events)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Daily Event Counts around the Largest Event",
       subtitle = paste0("Largest event: M = ",largest_event$mag),
       x = "Days from largest event",
       y = "Number of events") +
  theme_minimal()
print(p_around_daily)
ggsave("daily_counts_around_largest_event_Mc_2.png",p_around_daily,width = 8,
       height = 5,dpi = 300)

# Spatial distribution map
catalog_main_spatial <- catalog_main %>%
  mutate(period_group = ifelse(year(datetime) == 2019,"2019","non-2019"))

usa_map <- map_data("state")
california_map <- usa_map %>%
  filter(region == "california")

p_spatial_map <- ggplot() +
  geom_polygon(data = california_map,
               aes(x = long, y = lat,group = group),
               fill = "gray95",color = "gray60") +
  geom_point(data = catalog_main_spatial,
             aes(x = lon, y = lat,size = mag,color = period_group),
             alpha = 0.45) +
  coord_fixed(xlim = c(min(catalog_main_spatial$lon, na.rm = TRUE) - 0.2,
                       max(catalog_main_spatial$lon, na.rm = TRUE) + 0.2),
              ylim = c(min(catalog_main_spatial$lat, na.rm = TRUE) - 0.2,
                       max(catalog_main_spatial$lat, na.rm = TRUE) + 0.2)) +
  scale_size_continuous(range = c(0.5, 4)) +
  labs(title = "Spatial Distribution of Earthquakes",
       subtitle = "Point size represents magnitude; color distinguishes 2019 from other years",
    x = "Longitude", y = "Latitude",size = "Magnitude",color = "Period") +
  theme_minimal()

print(p_spatial_map)
ggsave("spatial_distribution_2019_vs_non2019.png",
       p_spatial_map,width = 8,height = 6,dpi = 300)
catalog_non2019 <- catalog_main_spatial %>%
  filter(period_group == "non-2019")

catalog_2019 <- catalog_main_spatial %>%
  filter(period_group == "2019")

p_spatial_map <- ggplot() +
  geom_polygon(data = california_map,
               aes(x = long,y = lat,group = group),
               fill = "gray95",
               color = "gray60") +
  geom_point(data = catalog_non2019,
             aes(x = lon,y = lat,size = mag,color = period_group),alpha = 0.25) +
  geom_point(data = catalog_2019,
             aes(x = lon,y = lat,size = mag,color = period_group),alpha = 0.65) +
  coord_fixed(xlim = c(min(catalog_main_spatial$lon, na.rm = TRUE) - 0.2,
                       max(catalog_main_spatial$lon, na.rm = TRUE) + 0.2),
              ylim = c(min(catalog_main_spatial$lat, na.rm = TRUE) - 0.2,
                       max(catalog_main_spatial$lat, na.rm = TRUE) + 0.2)) +
  scale_size_continuous(range = c(0.5, 4)) +
  labs(title = "Spatial Distribution of Earthquakes",
       subtitle = "Magnitude shown by point size; 2019 events highlighted",
       x = "Longitude",y = "Latitude",size = "Magnitude",color = "Period") +
  theme_minimal()

print(p_spatial_map)
ggsave("spatial_distribution_2019_vs_non2019.png",p_spatial_map,width = 8,height = 6,dpi = 300)
# ============================================================
# step 2: Estimate Mc
# ============================================================
# Frequency-Magnitude Distribution
bin_width <- 0.1
make_fmd <- function(mags, bin_width = 0.1) {
  mags <- mags[!is.na(mags)]
  min_m <- floor(min(mags) / bin_width) * bin_width
  max_m <- ceiling(max(mags) / bin_width) * bin_width
  bins <- seq(min_m, max_m, by = bin_width)
  fmd <- data.frame(mag_bin = bins) %>%
    mutate(count = sapply(mag_bin, function(m) {
      sum(mags >= m & mags < m + bin_width)}),
      cumulative_count = sapply(mag_bin, function(m) {
        sum(mags >= m)})) %>%
    filter(cumulative_count > 0)
  return(fmd)}
fmd <- make_fmd(catalog$mag, bin_width)

# Incremental FMD
p1 <- ggplot(fmd, aes(x = mag_bin, y = count)) +
  geom_point() +
  geom_line() +
  scale_y_log10() +
  labs(title = "Incremental Frequency-Magnitude Distribution",
       x = "Magnitude",y = "Number of events") +
  theme_minimal()
# Cumulative FMD
p2 <- ggplot(fmd, aes(x = mag_bin, y = cumulative_count)) +
  geom_point() +
  geom_line() +
  scale_y_log10() +
  labs(title = "Cumulative Frequency-Magnitude Distribution",
       x = "Magnitude",y = "Cumulative number of events") +
  theme_minimal()

print(p1)
print(p2)

#  b-value MLE
estimate_b_value <- function(mags, mc, bin_width = 0.1) {
  x <- mags[mags >= mc]
  n <- length(x)
  if (n < 30) {
    return(data.frame(mc = mc,n = n,b = NA,delta_b = NA))}
  mean_mag <- mean(x)
  b <- log10(exp(1)) / (mean_mag - (mc - bin_width / 2))
  delta_b <- 2.3 * b^2 * sqrt(sum((x - mean_mag)^2) / (n * (n - 1)))
  return(data.frame(mc = mc,n = n,b = b,delta_b = delta_b))}

# MBS-WW
estimate_mbs_ww <- function(mags, bin_width = 0.1, delta_m = 0.5) {
  mags <- mags[!is.na(mags)]
  candidates <- seq(floor(min(mags) / bin_width) * bin_width,
                    floor(max(mags) / bin_width) * bin_width - delta_m,
                    by = bin_width)
  b_table <- purrr::map_dfr(candidates,
                            ~ estimate_b_value(mags, .x, bin_width))
  b_table <- b_table %>%
    filter(!is.na(b))
  mbs_table <- b_table %>%
    rowwise() %>%
    mutate(b_avg = mean(b_table$b[b_table$mc >= mc & b_table$mc <= mc + delta_m],
                        na.rm = TRUE),
           abs_diff = abs(b - b_avg),
           stable = abs_diff <= delta_b) %>%
    ungroup()

  mc_est <- mbs_table %>%
    filter(stable == TRUE) %>%
    slice(1) %>%
    pull(mc)
  
  if (length(mc_est) == 0) mc_est <- NA
  
  return(list(mc = mc_est,table = mbs_table))}

mbs <- estimate_mbs_ww(catalog$mag, bin_width, delta_m = 0.5)
mc_mbs <- mbs$mc
cat("MBS-WW Mc =", mc_mbs, "\n")

# Plot b-value stability
ggplot(mbs$table, aes(x = mc, y = b)) +
  geom_line() +
  geom_point() +
  geom_line(aes(y = b_avg), linetype = "dashed") +
  geom_errorbar(aes(ymin = b - delta_b, ymax = b + delta_b),width = 0.03,alpha = 0.4) +
  geom_vline(xintercept = mc_mbs, linetype = "dashed") +
  labs(title = "MBS-WW b-value Stability Method",x = "Candidate Mc",y = "b-value") +
  theme_minimal()

# GFT check
gft_at_mc <- function(mags, mc, bin_width = 0.1) {
  mags <- mags[!is.na(mags)]
  x <- mags[mags >= mc]
  if (length(x) < 50) {
    warning("Too few events above this Mc.")
    return(data.frame(mc = mc,n = length(x),b = NA,gft = NA))}
  b_info <- estimate_b_value(mags, mc, bin_width)
  b <- b_info$b
  
  if (is.na(b)) {
    warning("b-value could not be estimated.")
    return(data.frame(mc = mc,n = length(x),b = NA,gft = NA))}
  
  mag_bins <- seq(mc, max(mags), by = bin_width)
  obs <- sapply(mag_bins, function(m) {sum(mags >= m)})
  pred <- obs[1] * 10^(-b * (mag_bins - mc))
  gft <- 100 - (sum(abs(obs - pred)) / sum(obs)) * 100
  result <- data.frame(mc = mc,n = length(x),b = b,gft = gft)
  
  return(result)
}

# Use MBS-WW result
mc_mbs <- 1.7
gft_check_mbs <- gft_at_mc(mags = catalog$mag,mc = mc_mbs,bin_width = bin_width)
print(gft_check_mbs)
cat("GFT value at MBS-WW Mc =", gft_check_mbs$gft, "\n")

if (gft_check_mbs$gft >= 95) {
  cat("Result: Mc = 1.7 passes the GFT-95% criterion.\n")
} else if (gft_check_mbs$gft >= 90) {
  cat("Result: Mc = 1.7 passes the GFT-90% criterion, but not GFT-95%.\n")
} else {
  cat("Result: Mc = 1.7 does not pass GFT-90%; consider using a higher Mc.\n")
}

# Plot observed vs predicted cumulative FMD at Mc = 1.7
plot_gft_at_mc <- function(mags, mc, bin_width = 0.1) {
  mags <- mags[!is.na(mags)]
  b_info <- estimate_b_value(mags, mc, bin_width)
  b <- b_info$b
  mag_bins <- seq(mc, max(mags), by = bin_width)
  obs <- sapply(mag_bins, function(m) {sum(mags >= m)})
  pred <- obs[1] * 10^(-b * (mag_bins - mc))
  plot_data <- data.frame(mag_bin = mag_bins,observed = obs,predicted = pred)
  
  ggplot(plot_data, aes(x = mag_bin)) +
    geom_point(aes(y = observed), size = 2) +
    geom_line(aes(y = observed), linewidth = 0.8) +
    geom_line(aes(y = predicted), linetype = "dashed", linewidth = 0.9) +
    scale_y_log10() +
    geom_vline(xintercept = mc, linetype = "dotted") +
    labs(title = paste0("GFT Check at MBS-WW Mc = ", mc),
         subtitle = paste0("Observed cumulative FMD vs Gutenberg-Richter prediction; b = ",
        round(b, 3)),
        x = "Magnitude",y = "Cumulative number of events") +
    theme_minimal()
  }

plot_gft_at_mc(mags = catalog$mag,mc = 1.7,bin_width = bin_width)

gft_monte_carlo_envelope <- function(mags,mc = 1.7,bin_width = 0.1,n_sim = 1000,
                                     conf = 0.95,seed = 123) {
  set.seed(seed)
  mags <- mags[!is.na(mags)]
  mags_mc <- mags[mags >= mc]
  n <- length(mags_mc)
  b_info <- estimate_b_value(mags, mc, bin_width)
  b <- b_info$b
  mag_bins <- seq(mc, max(mags_mc), by = bin_width)
  observed <- sapply(mag_bins, function(m) {sum(mags_mc >= m)})
  predicted <- n * 10^(-b * (mag_bins - mc))

  sim_cum <- matrix(NA, nrow = length(mag_bins), ncol = n_sim)
  
  for (s in 1:n_sim) {
    u <- runif(n)
    sim_mags <- mc - log10(u) / b
    
    sim_cum[, s] <- sapply(mag_bins, function(m) {
      sum(sim_mags >= m)
    })
  }
  
  lower_prob <- (1 - conf) / 2
  upper_prob <- 1 - lower_prob
  envelope <- data.frame(mag_bin = mag_bins,observed = observed,
                         predicted = predicted,
                         lower = apply(sim_cum, 1, quantile, probs = lower_prob, na.rm = TRUE),
                         upper = apply(sim_cum, 1, quantile, probs = upper_prob, na.rm = TRUE))

  envelope <- envelope %>%
    mutate(inside_envelope = observed >= lower & observed <= upper,
           expected_count = predicted)
  
  summary <- envelope %>%
    summarise(mc = mc,n = n,b = b,n_bins = n(),bins_inside = sum(inside_envelope),
              proportion_inside = mean(inside_envelope),bins_outside = sum(!inside_envelope))
  return(list(envelope = envelope,summary = summary,b_info = b_info))
}


mc_env_1_7 <- gft_monte_carlo_envelope(mags = catalog$mag,mc = 1.7,bin_width = bin_width,
                                       n_sim = 1000,conf = 0.95,seed = 123)
print(mc_env_1_7$summary)
# ensure if the large number of earthquakes in 2019 will affect the choice of Mc
# Split catalogue into 2019 and non-2019 groups
catalog_fmd_grouped <- catalog %>%
  mutate(year = lubridate::year(datetime),
         period_group = ifelse(year == 2019, "2019", "non-2019"))
group_count_summary <- catalog_fmd_grouped %>%
  count(period_group) %>%
  mutate(percentage = 100 * n / sum(n))
print(group_count_summary)
write.csv(group_count_summary,"FMD_group_event_counts_2019_vs_non2019.csv",
          row.names = FALSE)

# Incremental FMD comparison: 2019 vs non-2019
make_fmd_by_group <- function(df, bin_width = 0.1) {
  df <- df %>%
    filter(!is.na(mag), !is.na(period_group))
  
  mag_min <- floor(min(df$mag) / bin_width) * bin_width
  mag_max <- ceiling(max(df$mag) / bin_width) * bin_width
  
  mag_bins <- seq(mag_min, mag_max, by = bin_width)
  
  fmd_grouped <- df %>%
    mutate(mag_bin = floor(mag / bin_width) * bin_width) %>%
    count(period_group, mag_bin, name = "incremental_count") %>%
    tidyr::complete(period_group,mag_bin = mag_bins,fill = list(incremental_count = 0)) %>%
    group_by(period_group) %>%
    arrange(mag_bin, .by_group = TRUE) %>%
    mutate(cumulative_count = purrr::map_dbl( mag_bin,
                                              ~ sum(incremental_count[mag_bin >= .x]))) %>%
    ungroup()
  return(fmd_grouped)
}

fmd_2019_compare <- make_fmd_by_group(df = catalog_fmd_grouped,bin_width = bin_width)
p_incremental_2019 <- ggplot(fmd_2019_compare,
                             aes(x = mag_bin, y = incremental_count, linetype = period_group)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 1.7, linetype = "dotted") +
  scale_y_log10() +
  labs(title = "Incremental FMD: 2019 vs non-2019",
       subtitle = "Dotted line indicates Mc = 1.7",x = "Magnitude",
       y = "Incremental number of events",linetype = "Period") +
  theme_minimal()
print(p_incremental_2019)

ggsave("Incremental_FMD_2019_vs_non2019.png",p_incremental_2019,width = 8,height = 5,
       dpi = 300)

# Cumulative FMD comparison: 2019 vs non-2019
p_cumulative_2019 <- ggplot(fmd_2019_compare,
                            aes(x = mag_bin, y = cumulative_count, linetype = period_group)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 1.7, linetype = "dotted") +
  scale_y_log10() +
  labs(title = "Cumulative FMD: 2019 vs non-2019",
       subtitle = "Dotted line indicates Mc = 1.7",x = "Magnitude",
       y = "Cumulative number of events",linetype = "Period") +
  theme_minimal()
print(p_cumulative_2019)
ggsave("Cumulative_FMD_2019_vs_non2019.png",p_cumulative_2019,width = 8,height = 5,dpi = 300)

# Normalized cumulative FMD comparison
fmd_2019_compare_norm <- fmd_2019_compare %>%
  group_by(period_group) %>%
  mutate(cumulative_norm = cumulative_count / max(cumulative_count, na.rm = TRUE)) %>%
  ungroup()

p_cumulative_norm_2019 <- ggplot(fmd_2019_compare_norm,
                                 aes(x = mag_bin, y = cumulative_norm, linetype = period_group)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 1.7, linetype = "dotted") +
  scale_y_log10() +
  labs(title = "Normalized Cumulative FMD: 2019 vs non-2019",
       subtitle = "Curves are normalized by each group's total event count",
       x = "Magnitude",y = "Normalized cumulative proportion",linetype = "Period") +
  theme_minimal()
print(p_cumulative_norm_2019)
ggsave("Normalized_Cumulative_FMD_2019_vs_non2019.png",p_cumulative_norm_2019,width = 8,
       height = 5,dpi = 300)

# b-value and GFT comparison at Mc = 1.7
mc_main <- 1.7
fmd_group_quality <- catalog_fmd_grouped %>%
  group_by(period_group) %>%
  group_modify(~ {
    mags_group <- .x$mag
    b_info <- estimate_b_value(mags = mags_group, mc = mc_main,bin_width = bin_width)
    gft_info <- gft_at_mc(mags = mags_group, mc = mc_main,bin_width = bin_width)
    data.frame(n_total = length(mags_group),
               n_above_mc = sum(mags_group >= mc_main, na.rm = TRUE),
               mc = mc_main,b = b_info$b,delta_b = b_info$delta_b,gft = gft_info$gft)
  }) %>%
  ungroup()
print(fmd_group_quality)
write.csv(fmd_group_quality,"FMD_quality_2019_vs_non2019_at_Mc_1_7.csv",row.names = FALSE)

# Estimate MBS-WW Mc separately for 2019 and non-2019
mc_mbs_by_group <- catalog_fmd_grouped %>%
  group_by(period_group) %>%
  group_modify(~ {
    mags_group <- .x$mag
    mbs_group <- estimate_mbs_ww(mags = mags_group,bin_width = bin_width,delta_m = 0.5)
    data.frame(n_total = length(mags_group), mc_mbs = mbs_group$mc)
  }) %>%
  ungroup()
print(mc_mbs_by_group)
write.csv( mc_mbs_by_group,"MBS_WW_Mc_2019_vs_non2019.csv",row.names = FALSE)

# Temporal Mc analysis using rolling event windows
# check if mc changed a lot during 2019
window_size <- 1000 
step_size   <- 250     
bin_width   <- 0.1
mc_check_values <- c(1.7, 2.0)
catalog_temporal <- catalog %>%
  filter(!is.na(datetime), !is.na(mag)) %>%
  arrange(datetime)
n_events <- nrow(catalog_temporal)
window_starts <- seq(from = 1,to = n_events - window_size + 1,by = step_size)
temporal_mc_results <- purrr::map_dfr(window_starts,function(start_idx) {
    end_idx <- start_idx + window_size - 1
    temp <- catalog_temporal[start_idx:end_idx, ]
    mags <- temp$mag
    start_time  <- min(temp$datetime)
    end_time    <- max(temp$datetime)
    center_time <- start_time + (end_time - start_time) / 2
    mc_maxc_temp <- estimate_maxc(mags = mags,bin_width = bin_width)
    mbs_temp <- estimate_mbs_ww(mags = mags,bin_width = bin_width,delta_m = 0.5)
    mc_mbs_temp <- mbs_temp$mc
    
    b_17 <- estimate_b_value(mags = mags,mc = 1.7, bin_width = bin_width)

    b_20 <- estimate_b_value(mags = mags,mc = 2.0,bin_width = bin_width)
    
    data.frame(start_time = start_time,center_time = center_time,end_time = end_time,
               n_total = length(mags), mc_maxc = mc_maxc_temp, mc_mbs = mc_mbs_temp,
               n_above_1_7 = b_17$n, b_at_1_7 = b_17$b,delta_b_1_7 = b_17$delta_b,
               n_above_2_0 = b_20$n,b_at_2_0 = b_20$b, delta_b_2_0 = b_20$delta_b)
  }
)
print(head(temporal_mc_results))

# Plot temporal Mc
temporal_mc_long <- temporal_mc_results %>%
  select( center_time, mc_mbs) %>%
  pivot_longer( cols = c(mc_mbs),names_to = "method",values_to = "mc")

p_temporal_mc <- ggplot(temporal_mc_long,
                        aes(x = center_time, y = mc, linetype = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 1.7,linetype = "dotted" ) +
  geom_hline( yintercept = 2.0, linetype = "dashed") +
  labs(title = "Temporal Variation of Catalogue Completeness Magnitude",
       subtitle = paste0("Rolling window: ", window_size," events; step: ", step_size,
                         " events" ),
       x = "Time",
       y = "Estimated Mc",
       linetype = "Method") +
  theme_minimal()

print(p_temporal_mc)
p_temporal_mc_2019 <- ggplot(temporal_mc_long,
                             aes(x = center_time,y = mc,linetype = method )) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.4) +
  annotate( "rect",xmin = as.POSIXct("2019-01-01", tz = "UTC"),
            xmax = as.POSIXct("2019-12-31 23:59:59", tz = "UTC"), ymin = -Inf, ymax = Inf,
            alpha = 0.08) +
  geom_hline(yintercept = 1.7,linetype = "dotted") +
  geom_hline( yintercept = 2.0, linetype = "dashed" ) +
  labs(x = "Time",y = "Estimated Mc",linetype = "Method" ) +
  theme_minimal()

print(p_temporal_mc_2019)

# Temporal b-value plot
p_b_temporal <- ggplot(temporal_mc_results,
                       aes( x = center_time,y = b_at_2_0)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = b_at_2_0 - delta_b_2_0, ymax = b_at_2_0 + delta_b_2_0),
                width = 10, alpha = 0.4) +
  annotate( "rect",xmin = as.POSIXct("2019-01-01", tz = "UTC"),
            xmax = as.POSIXct("2019-12-31 23:59:59", tz = "UTC"),ymin = -Inf,
            ymax = Inf,alpha = 0.08) +
  labs( title = "Temporal Variation of b-value",
        subtitle = "b-value estimated using events with M >= 2.0",x = "Time",y = "b-value" ) +
  theme_minimal()

print(p_b_temporal)
write.csv( temporal_mc_results,"Temporal_Mc_Rolling_Window.csv",row.names = FALSE)
ggsave("Temporal_Mc_Rolling_Window.png", p_temporal_mc,width = 10,height = 5,dpi = 300)
ggsave("Temporal_Mc_2019_Highlight.png", p_temporal_mc_2019, width = 10,height = 5, dpi = 300)
ggsave("Temporal_b_value_Mc_2_0.png",p_b_temporal,width = 10, height = 5,dpi = 300)
temporal_mc_results %>%
  filter( year(center_time) == 2019,mc_maxc >= 1.7 | mc_mbs >= 1.7) %>%
  select(start_time,center_time, end_time,mc_maxc, mc_mbs,n_above_1_7, b_at_1_7, n_above_2_0,
         b_at_2_0) %>%
  arrange(center_time)


# ============================================================
# STEP 3: Declustering and ETAS initialization
# ============================================================
# Methods
# A. Gardner-Knopoff window declustering
# B. Zaliapin-Ben-Zion-style nearest-neighbour cluster declustering
# C. Reasenberg-style interaction-link declustering

`%||%` <- function(x, y) if (is.null(x) || length(x)==0 || all(is.na(x))) y else x

mc_main <- 2.0
bin_width <- 0.1
mc_tag <- gsub("\\.", "_", sprintf("%.1f", mc_main))

catalog_main <- catalog %>%
  filter( mag >= mc_main) %>%
  arrange(datetime) %>%
  mutate(event_index = row_number(), original_order = row_number(), year = year(datetime))

stopifnot(all(catalog_main$mag >= mc_main))
cat("Main declustering catalogue: Mc =", mc_main, "; N =", nrow(catalog_main), "\n")

haversine_km <- function(lon1, lat1, lon2, lat2) {
  R <- 6371.227
  lon1 <- lon1*pi/180; lat1 <- lat1*pi/180
  lon2 <- lon2*pi/180; lat2 <- lat2*pi/180
  dlon <- lon2-lon1; dlat <- lat2-lat1
  a <- sin(dlat/2)^2 + cos(lat1)*cos(lat2)*sin(dlon/2)^2
  R * 2 * atan2(sqrt(a), sqrt(1-a))
}

decimal_year <- function(datetime) {
  tz <- attr(datetime, "tzone") %||% "UTC"
  y <- year(datetime)
  y0 <- as.POSIXct(paste0(y,"-01-01 00:00:00"), tz=tz)
  y1 <- as.POSIXct(paste0(y+1,"-01-01 00:00:00"), tz=tz)
  y + as.numeric(difftime(datetime,y0,units="secs")) /
    as.numeric(difftime(y1,y0,units="secs"))
}

estimate_b_value_local <- function(mags, mc, bin_width=0.1) {
  x <- mags[is.finite(mags) & mags >= mc]
  n <- length(x)
  if (n < 30) return(data.frame(mc=mc,n=n,b=NA,delta_b=NA))
  mean_mag <- mean(x)
  b <- log10(exp(1))/(mean_mag-(mc-bin_width/2))
  delta_b <- 2.3*b^2*sqrt(sum((x-mean_mag)^2)/(n*(n-1)))
  data.frame(mc=mc,n=n,b=b,delta_b=delta_b)
}

b_info_declustering <- estimate_b_value_local(catalog_main$mag, mc_main, bin_width)
b_catalogue <- b_info_declustering$b
cat("b-value estimated at Mc=2.0:", b_catalogue, "\n")

save_csv_mc <- function(x, stem) {
  f <- paste0(stem,"_Mc_",mc_tag,".csv")
  write.csv(x,f,row.names=FALSE); invisible(f)
}
save_plot_mc <- function(p, stem, width=9, height=5, dpi=300) {
  f <- paste0(stem,"_Mc_",mc_tag,".png")
  ggsave(f,p,width=width,height=height,dpi=dpi); invisible(f)
}

#  A. GARDNER-KNOPOFF
gk_space_window_km <- function(M) 10^(0.1238*M + 0.983)
gk_time_window_days <- function(M) ifelse(M < 6.5,
                                          10^(0.5409*M - 0.547), 10^(0.032*M + 2.7389))

gardner_knopoff_declustering <- function(df, fs_time_prop=1.0, time_cutoff_days=NULL) {
  stopifnot(fs_time_prop>=0, fs_time_prop<=1)
  d0 <- df %>% arrange(datetime) %>%
    mutate(original_order=row_number(),
           gk_space_window_km=gk_space_window_km(mag),
           gk_time_window_days=gk_time_window_days(mag))
  if (!is.null(time_cutoff_days))
    d0$gk_time_window_days <- pmin(d0$gk_time_window_days,time_cutoff_days)
  
  d <- d0 %>% arrange(desc(mag),datetime)
  n <- nrow(d)
  d$gk_cluster_id <- 0L; d$gk_flag <- 0L
  d$gk_parent_evid <- NA; d$gk_parent_mag <- NA_real_
  cid <- 0L
  
  if (n>=2) for (i in seq_len(n-1)) {
    if (i%%1000==0) cat("GK",i,"/",n,"\n")
    if (d$gk_cluster_id[i]!=0L) next
    dt <- as.numeric(difftime(d$datetime,d$datetime[i],units="days"))
    keep_t <- d$gk_cluster_id==0L &
      dt >= -d$gk_time_window_days[i]*fs_time_prop &
      dt <= d$gk_time_window_days[i]
    idx <- which(keep_t); if (length(idx)<=1) next
    dist <- haversine_km(d$lon[idx],d$lat[idx],d$lon[i],d$lat[i])
    inside <- idx[dist <= d$gk_space_window_km[i]]
    if (length(setdiff(inside,i))==0) next
    cid <- cid+1L
    d$gk_cluster_id[inside] <- cid
    d$gk_flag[inside] <- 1L
    d$gk_flag[inside[dt[inside]<0]] <- -1L
    d$gk_flag[i] <- 0L
    nonmain <- setdiff(inside,i)
    d$gk_parent_evid[nonmain] <- d$evid[i]
    d$gk_parent_mag[nonmain] <- d$mag[i]
  }
  
  d %>% arrange(original_order) %>% mutate(
    gk_label=case_when(
      gk_cluster_id==0L ~ "background",
      gk_cluster_id!=0L & gk_flag==0L ~ "mainshock",
      gk_flag==1L ~ "aftershock",
      gk_flag==-1L ~ "foreshock",
      TRUE ~ "unknown"),
    gk_is_background_like=gk_label %in% c("background","mainshock"),
    gk_background_prob=as.numeric(gk_is_background_like))
}

gk_fs_baseline <- 1.0
declust_gk <- gardner_knopoff_declustering(catalog_main,gk_fs_baseline)
gk_summary <- declust_gk %>% count(gk_label,name="n") %>%
  mutate(percentage=100*n/sum(n),method="Gardner-Knopoff")
print(gk_summary)
save_csv_mc(declust_gk,"declustering_gardner_knopoff")
save_csv_mc(gk_summary,"declustering_gardner_knopoff_summary")

etas_init_gk <- declust_gk %>% select(
  evid,datetime,lat,lon,depth,mag,gk_cluster_id,gk_flag,gk_label,
  gk_is_background_like,gk_background_prob,gk_parent_evid,gk_parent_mag,
  gk_space_window_km,gk_time_window_days)
save_csv_mc(etas_init_gk,"etas_initialization_gardner_knopoff")

# B. ZALIAPIN-BEN-ZION-STYLE NEAREST NEIGHBOUR
compute_nearest_neighbor_eta <- function(df,mc,b_value,d_f=1.6,use_depth=FALSE,
                                         max_lookback_years=Inf,max_previous_events=Inf) {
  d <- df %>% arrange(datetime) %>% mutate(
    event_index=row_number(),decimal_year=decimal_year(datetime),
    nn_parent_index=NA_integer_,nn_parent_evid=NA,
    nn_time_years=NA_real_,nn_distance_km=NA_real_,
    nn_parent_mag=NA_real_,nn_log_eta=NA_real_)
  n <- nrow(d)
  if (n>=2) for (j in 2:n) {
    if (j%%1000==0) cat("NN",j,"/",n,"\n")
    prev <- seq_len(j-1)
    dtall <- d$decimal_year[j]-d$decimal_year[prev]
    keep <- dtall>0
    if (is.finite(max_lookback_years)) keep <- keep & dtall<=max_lookback_years
    prev <- prev[keep]; if (!length(prev)) next
    if (is.finite(max_previous_events) && length(prev)>max_previous_events)
      prev <- tail(prev,as.integer(max_previous_events))
    dt <- d$decimal_year[j]-d$decimal_year[prev]
    dt[dt<=0] <- 1/(365.25*24*3600)
    epi <- haversine_km(d$lon[prev],d$lat[prev],d$lon[j],d$lat[j])
    epi[epi<=0] <- 0.001
    r <- if (use_depth) sqrt(epi^2+(d$depth[j]-d$depth[prev])^2) else epi
    r[r<=0] <- 0.001
    pm <- d$mag[prev]
    logeta <- log10(dt)+d_f*log10(r)-b_value*(pm-mc)
    q <- which.min(logeta); parent <- prev[q]
    d$nn_parent_index[j] <- parent; d$nn_parent_evid[j] <- d$evid[parent]
    d$nn_time_years[j] <- dt[q]; d$nn_distance_km[j] <- r[q]
    d$nn_parent_mag[j] <- pm[q]; d$nn_log_eta[j] <- logeta[q]
  }
  d
}

fit_nn_mixture_threshold <- function(log_eta,seed=123) {
  x <- log_eta[is.finite(log_eta)]
  if (length(x)<100) stop("Too few finite eta values")
  set.seed(seed)
  fit <- mclust::Mclust(x,G=2,modelNames=c("E","V"),verbose=FALSE)
  grid <- seq(quantile(x,.001),quantile(x,.999),length.out=4000)
  pr <- predict(fit,newdata=grid)
  thr <- if (!is.null(pr$z) && ncol(pr$z)==2)
    grid[which.min(abs(pr$z[,1]-pr$z[,2]))] else mean(fit$parameters$mean)
  list(fit=fit,threshold=thr,means=sort(as.numeric(fit$parameters$mean)))
}

classify_nn_clusters <- function(df,threshold) {
  d <- df %>% arrange(datetime) %>% mutate(
    nn_is_clustered_edge=ifelse(is.na(nn_log_eta),FALSE,nn_log_eta<threshold),
    nn_cluster_id=0L,nn_mainshock_evid=NA,nn_label="background")
  n <- nrow(d); adj <- vector("list",n); for(i in seq_len(n)) adj[[i]] <- integer(0)
  for(j in seq_len(n)) {
    p <- d$nn_parent_index[j]
    if (!is.na(p) && d$nn_is_clustered_edge[j]) {
      adj[[j]] <- unique(c(adj[[j]],p)); adj[[p]] <- unique(c(adj[[p]],j))
    }
  }
  visited <- rep(FALSE,n); cid <- 0L
  for(i in seq_len(n)) {
    if (visited[i]) next
    stack <- i; comp <- integer(0)
    while(length(stack)) {
      v <- stack[1]; stack <- stack[-1]
      if (visited[v]) next
      visited[v] <- TRUE; comp <- c(comp,v); stack <- unique(c(stack,adj[[v]]))
    }
    if (length(comp)==1) next
    cid <- cid+1L
    m <- comp[which.max(d$mag[comp])]; mt <- d$datetime[m]; me <- d$evid[m]
    d$nn_cluster_id[comp] <- cid; d$nn_mainshock_evid[comp] <- me
    for(k in comp) d$nn_label[k] <- if(k==m) "mainshock" else if(d$datetime[k]<mt) "foreshock" else "aftershock"
  }
  d %>% mutate(nn_is_background_like=nn_label %in% c("background","mainshock"),
               nn_background_prob=as.numeric(nn_is_background_like))
}

nn_df_baseline <- 1.6
nn_raw <- compute_nearest_neighbor_eta(catalog_main,mc_main,b_catalogue,nn_df_baseline,
                                       use_depth=FALSE,max_lookback_years=Inf,max_previous_events=Inf)
nn_mix <- fit_nn_mixture_threshold(nn_raw$nn_log_eta)
nn_eta_threshold <- nn_mix$threshold
cat("NN: d_f=",nn_df_baseline," b=",b_catalogue," threshold=",nn_eta_threshold,"\n")
nn_component_order <- order(as.numeric(nn_mix$fit$parameters$mean))
nn_mixture_diagnostic <- data.frame(model = nn_mix$fit$modelName,
                                    threshold_log10_eta = nn_eta_threshold,
                                    lower_component_mean =as.numeric(nn_mix$fit$parameters$mean)[nn_component_order[1]],
                                    upper_component_mean =as.numeric(nn_mix$fit$parameters$mean)[nn_component_order[2]],
                                    lower_component_proportion =nn_mix$fit$parameters$pro[nn_component_order[1]],
                                    upper_component_proportion =nn_mix$fit$parameters$pro[nn_component_order[2]])
print(nn_mixture_diagnostic)
save_csv_mc(nn_mixture_diagnostic,"nearest_neighbor_mixture_diagnostic")
declust_nn <- classify_nn_clusters(nn_raw,nn_eta_threshold)
nn_summary <- declust_nn %>% count(nn_label,name="n") %>%
  mutate(percentage=100*n/sum(n),method="Nearest-neighbour")
print(nn_summary)
save_csv_mc(declust_nn,"declustering_nearest_neighbor")
save_csv_mc(nn_summary,"declustering_nearest_neighbor_summary")

etas_init_nn <- declust_nn %>% select(
  evid,datetime,lat,lon,depth,mag,nn_parent_evid,nn_parent_index,
  nn_time_years,nn_distance_km,nn_parent_mag,nn_log_eta,nn_cluster_id,
  nn_mainshock_evid,nn_label,nn_is_background_like,nn_background_prob)
save_csv_mc(etas_init_nn,"etas_initialization_nearest_neighbor")

p_nn_eta <- ggplot(declust_nn %>% filter(is.finite(nn_log_eta)),aes(nn_log_eta))+
  geom_histogram(bins=80)+geom_vline(xintercept=nn_eta_threshold,linetype="dashed")+
  labs(title="Nearest-neighbour proximity distribution",
       subtitle="Two-component Gaussian-mixture threshold",
       x=expression(log[10](eta)),y="Number of events")+theme_minimal()
print(p_nn_eta); save_plot_mc(p_nn_eta,"nearest_neighbor_log_eta_distribution",8,5)

# REASENBERG-STYLE INTERACTION-LINK DECLUSTERING
reas_baseline <- list(tau_min = 1,tau_max = 10,p1 = 0.95,xk = 0.5,xmeff = 4.0,rfact= 10)

haversine_km <- function(lon1, lat1, lon2, lat2) {
  earth_radius_km <- 6371.227
  lon1_rad <- lon1 * pi / 180
  lat1_rad <- lat1 * pi / 180
  lon2_rad <- lon2 * pi / 180
  lat2_rad <- lat2 * pi / 180
  
  dlon <- lon2_rad - lon1_rad
  dlat <- lat2_rad - lat1_rad
  
  a <- sin(dlat / 2)^2 + cos(lat1_rad) *cos(lat2_rad) *sin(dlon / 2)^2
  central_angle <- 2 * atan2(sqrt(a),sqrt(1 - a))
  earth_radius_km * central_angle
}

reasenberg_crack_radius_km <- function(magnitude) {
  10^(0.4 * magnitude - 1.2)
}

reasenberg_interaction_radius_km <- function(magnitude,rfact = 10) {
  rfact *reasenberg_crack_radius_km(magnitude)
}

reasenberg_tau <- function(cluster_max_mag,elapsed_days,tau_min = 1,tau_max = 10,
                           p1 = 0.95,xmeff = 4.0,xk = 0.5) {
  effective_cutoff <- xmeff +xk *pmax(cluster_max_mag - xmeff,0)
  delta_m <- pmax(cluster_max_mag - effective_cutoff,0)
  base_time <- pmax(elapsed_days,tau_min)
  tau_raw <- -log(1 - p1) *base_time /10^(2 * (delta_m - 1) / 3)
  tau_bounded <- pmin(tau_max,pmax(tau_min,tau_raw))
  tau_bounded
}

haversine_km <- function(lon1, lat1, lon2, lat2) {
  R <- 6371.227
  lon1 <- lon1 * pi / 180
  lat1 <- lat1 * pi / 180
  lon2 <- lon2 * pi / 180
  lat2 <- lat2 * pi / 180
  
  dlon <- lon2 - lon1
  dlat <- lat2 - lat1
  
  a <- sin(dlat / 2)^2 +cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  2 * R * atan2(sqrt(a), sqrt(1 - a))
}

reasenberg_crack_radius_km <- function(magnitude) {
  10^(0.4 * magnitude - 1.2)
}

reasenberg_interaction_radius_km <- function(magnitude,rfact = 10) {
  rfact *reasenberg_crack_radius_km(magnitude)
}

reasenberg_tau <- function(cluster_max_mag,elapsed_days,tau_min = 1,tau_max = 10,
                           p1 = 0.95,xmeff = 4.0,xk = 0.5) {
  effective_cutoff <- xmeff +xk * pmax(cluster_max_mag - xmeff,0)
  delta_m <- pmax(cluster_max_mag - effective_cutoff,0)
  base_time <- pmax(elapsed_days,tau_min)
  tau_raw <- -log(1 - p1) *base_time /10^(2 * (delta_m - 1) / 3)
  pmin(tau_max,pmax(tau_min, tau_raw))
}

assign_reasenberg_labels <- function(d,roots) {
  root_counts <- table(roots)
  clustered_roots <- as.integer(names(root_counts[root_counts > 1]))
  cluster_id_map <- setNames(seq_along(clustered_roots),clustered_roots)
  
  d$reasenberg_cluster_id <- ifelse(roots %in% clustered_roots,
                                    unname( cluster_id_map[as.character(roots)]),
                                    0L)
  
  d$reasenberg_label <- "background"
  d$reasenberg_mainshock_evid <- NA
  
  cluster_ids <- sort(unique(d$reasenberg_cluster_id[d$reasenberg_cluster_id > 0]))
  
  for (cluster_id in cluster_ids) {
    cluster_idx <- which(d$reasenberg_cluster_id == cluster_id)
    
    mainshock_idx <- cluster_idx[which.max(d$mag[cluster_idx])]
    mainshock_time <-d$datetime[mainshock_idx]
    mainshock_evid <-d$evid[mainshock_idx]
    d$reasenberg_mainshock_evid[cluster_idx] <- mainshock_evid
    labels <- ifelse(d$datetime[cluster_idx] <mainshock_time,"foreshock",
                     "aftershock")
    
    labels[cluster_idx == mainshock_idx] <- "mainshock"
    
    d$reasenberg_label[cluster_idx] <- labels
  }
  
  d$reasenberg_is_background_like <-
    d$reasenberg_label %in%
    c("background", "mainshock")
  
  d$reasenberg_background_prob <-as.numeric(d$reasenberg_is_background_like)
  
  d
}

reasenberg_style_declustering_fast <- function(df,tau_min = 1,tau_max = 10,p1 = 0.95,
                                               xk = 0.5,xmeff = 4.0, rfact = 10,progress_every = 1000) {
  d <- df %>%
    arrange(datetime)
  n_events <- nrow(d)
  
  if (n_events == 0) {
    return(d)
  }
  
  d$reasenberg_parent_index <- NA_integer_
  d$reasenberg_parent_evid <- NA
  d$reasenberg_distance_km <- NA_real_
  d$reasenberg_time_days <- NA_real_
  
  catalogue_start <- min(d$datetime)
  
  time_days <- as.numeric(difftime(d$datetime,
                                   catalogue_start,
                                   units = "days"))
  uf_parent <- seq_len(n_events)
  uf_size <- rep(1L, n_events)
  
  cluster_max_mag <- d$mag
  cluster_start_day <- time_days
  
  find_root <- function(x) {
    while (uf_parent[x] != x) {
      uf_parent[x] <<-uf_parent[ uf_parent[x] ]
      x <- uf_parent[x]
    }
    x
  }
  union_clusters <- function(a, b) {
    root_a <- find_root(a)
    root_b <- find_root(b)
    if (root_a == root_b) {
      return(root_a)
    }
    if (uf_size[root_a] < uf_size[root_b]) {
      temp <- root_a
      root_a <- root_b
      root_b <- temp
    }
    uf_parent[root_b] <<- root_a
    uf_size[root_a] <<- uf_size[root_a] +uf_size[root_b]
    cluster_max_mag[root_a] <<- max( cluster_max_mag[root_a],cluster_max_mag[root_b])
    cluster_start_day[root_a] <<- min(cluster_start_day[root_a],cluster_start_day[root_b])
    root_a
  }
  
  left_idx <- 1L
  for (j in 2:n_events) {
    
    if (
      progress_every > 0 &&
      j %% progress_every == 0
    ) {
      cat(sprintf("Reasenberg-fast: %d / %d\n",j, n_events))
    }
    while (
      left_idx < j &&
      time_days[j] -
      time_days[left_idx] >
      tau_max
    ) {
      left_idx <- left_idx + 1L
    }
    if (left_idx >= j) {
      next
    }
    candidate_idx <- left_idx:(j - 1L)
    linked_idx <- integer(length(candidate_idx))
    linked_dist <- numeric( length(candidate_idx) )
    n_linked <- 0L
    for (i in candidate_idx) {
      cluster_root <- find_root(i)
      elapsed_days <- time_days[i] -cluster_start_day[ cluster_root]
      tau_i <- reasenberg_tau(cluster_max_mag =cluster_max_mag[ cluster_root],
                              elapsed_days =elapsed_days,tau_min =tau_min,
                              tau_max =tau_max,p1 = p1, xmeff =xmeff,xk =xk)
      dt_ij <- time_days[j] -time_days[i]
      if (dt_ij > tau_i) {
        next
      }
      
      distance_ij <- haversine_km(d$lon[i],d$lat[i], d$lon[j], d$lat[j])
      interaction_radius <-reasenberg_interaction_radius_km(d$mag[i], rfact )
      
      if (
        distance_ij <=interaction_radius
      ) {
        n_linked <- n_linked + 1L
        linked_idx[n_linked] <- i
        linked_dist[n_linked] <- distance_ij
      }
    }
    
    if (n_linked == 0L) {
      next
    }
    
    linked_idx <-linked_idx[seq_len(n_linked)]
    linked_dist <-linked_dist[seq_len(n_linked)]
    for (parent_idx in linked_idx) {
      union_clusters( parent_idx,j)
    }
    
    nearest_pos <- which.min(linked_dist)
    parent_idx <- linked_idx[ nearest_pos]
    
    d$reasenberg_parent_index[j] <-parent_idx
    d$reasenberg_parent_evid[j] <- d$evid[parent_idx]
    d$reasenberg_distance_km[j] <-linked_dist[nearest_pos]
    d$reasenberg_time_days[j] <-time_days[j] - time_days[parent_idx]}
  
  roots <- integer(n_events)
  for (i in seq_len(n_events)) {
    roots[i] <- find_root(i)
  }
  
  d <- assign_reasenberg_labels( d,roots)
  return(d)
}

reas_start_time <- Sys.time()

declust_reasenberg <- reasenberg_style_declustering_fast(df =catalog_main,
                                                         tau_min =reas_baseline$tau_min,
                                                         tau_max =reas_baseline$tau_max,
                                                         p1 =reas_baseline$p1,
                                                         xk =reas_baseline$xk,
                                                         xmeff =reas_baseline$xmeff,
                                                         rfact =reas_baseline$rfact,
                                                         progress_every =1000)

reas_runtime_minutes <- as.numeric(difftime(Sys.time(),reas_start_time,units = "mins"))
cat(sprintf("\nReasenberg runtime: %.2f minutes\n",reas_runtime_minutes))

reasenberg_summary <-declust_reasenberg %>%
  count(reasenberg_label,name = "n") %>%
  mutate(percentage =100 *n /sum(n), method = "Reasenberg-style")
print(reasenberg_summary)


# ETAS initialisation output
etas_init_reasenberg <- declust_reasenberg %>%
  select(evid,datetime,lat,lon,depth, mag,reasenberg_cluster_id,
         reasenberg_parent_index,reasenberg_parent_evid,reasenberg_distance_km,
         reasenberg_time_days,reasenberg_mainshock_evid,reasenberg_label,
         reasenberg_is_background_like,reasenberg_background_prob)

save_csv_mc(declust_reasenberg,"declustering_reasenberg_style")
save_csv_mc(reasenberg_summary, "declustering_reasenberg_style_summary")
save_csv_mc(etas_init_reasenberg, "etas_initialization_reasenberg_style")

# Sanity checks
stopifnot(nrow(declust_reasenberg) ==nrow(catalog_main))
stopifnot(all(declust_reasenberg$mag >=2.0))
cat(sprintf("N events: %d\n",nrow(declust_reasenberg)))
cat(sprintf("Background-like fraction: %.4f\n",
            mean(declust_reasenberg$reasenberg_is_background_like)))

#HARMONISE EVENT-LEVEL OUTPUTS
event_level_compare <- catalog_main %>% select(evid,datetime,mag,lat,lon,year) %>%
  left_join(declust_gk %>% select(evid,gk_label,gk_is_background_like),by="evid") %>%
  left_join(declust_nn %>% select(evid,nn_label,nn_is_background_like),by="evid") %>%
  left_join(declust_reasenberg %>% select(evid,reasenberg_label,reasenberg_is_background_like),by="evid")
stopifnot(nrow(event_level_compare)==nrow(catalog_main))
stopifnot(!anyNA(event_level_compare[c("gk_is_background_like","nn_is_background_like","reasenberg_is_background_like")]))
save_csv_mc(event_level_compare,"declustering_event_level_comparison")

# DIAGNOSTIC 1: OVERALL BACKGROUND FRACTION
overall_background_summary <- bind_rows(
  data.frame(method="Gardner-Knopoff",n_total=nrow(event_level_compare),
             n_background_like=sum(event_level_compare$gk_is_background_like)),
  data.frame(method="Nearest-neighbour",n_total=nrow(event_level_compare),
             n_background_like=sum(event_level_compare$nn_is_background_like)),
  data.frame(method="Reasenberg-style",n_total=nrow(event_level_compare),
             n_background_like=sum(event_level_compare$reasenberg_is_background_like))) %>%
  mutate(n_clustered=n_total-n_background_like,
         background_fraction=n_background_like/n_total,
         background_percentage=100*background_fraction)
print(overall_background_summary)
save_csv_mc(overall_background_summary,"diagnostic_overall_background_fraction")

# DIAGNOSTIC 2: EVENT-LEVEL AGREEMENT
binary_jaccard <- function(x,y) { u<-sum(x|y); if(u==0) NA_real_ else sum(x&y)/u }
cohen_kappa_binary <- function(x,y) {
  x<-as.integer(x); y<-as.integer(y); po<-mean(x==y)
  pe<-mean(x==1)*mean(y==1)+mean(x==0)*mean(y==0)
  if(abs(1-pe)<1e-12) NA_real_ else (po-pe)/(1-pe)
}
agreement_pair <- function(x,y,m1,m2) data.frame(
  method_1=m1,method_2=m2,raw_agreement=mean(x==y),
  jaccard_background=binary_jaccard(x,y),kappa=cohen_kappa_binary(x,y),
  n_disagree=sum(x!=y),disagreement_percentage=100*mean(x!=y))

pairwise_agreement <- bind_rows(
  agreement_pair(event_level_compare$gk_is_background_like,event_level_compare$nn_is_background_like,
                 "Gardner-Knopoff","Nearest-neighbour"),
  agreement_pair(event_level_compare$gk_is_background_like,event_level_compare$reasenberg_is_background_like,
                 "Gardner-Knopoff","Reasenberg-style"),
  agreement_pair(event_level_compare$nn_is_background_like,event_level_compare$reasenberg_is_background_like,
                 "Nearest-neighbour","Reasenberg-style"))
print(pairwise_agreement)
save_csv_mc(pairwise_agreement,"diagnostic_pairwise_event_agreement")

event_level_disagreement <- event_level_compare %>% mutate(
  n_background_votes=as.integer(gk_is_background_like)+as.integer(nn_is_background_like)+as.integer(reasenberg_is_background_like),
  unanimous=n_background_votes %in% c(0,3),any_disagreement=!unanimous)
disagreement_summary <- event_level_disagreement %>% summarise(
  n_total=n(),n_unanimous=sum(unanimous),n_disagreement=sum(any_disagreement),
  disagreement_percentage=100*mean(any_disagreement))
print(disagreement_summary)
save_csv_mc(event_level_disagreement,"diagnostic_event_level_disagreement")
save_csv_mc(disagreement_summary,"diagnostic_event_level_disagreement_summary")

# DIAGNOSTIC 3: 2019 VS NON-2019
period_summary <- event_level_compare %>% mutate(period_group=ifelse(year==2019,"2019","non-2019")) %>%
  group_by(period_group) %>% summarise(
    n_total=n(),
    gk_background=sum(gk_is_background_like),gk_background_fraction=mean(gk_is_background_like),
    nn_background=sum(nn_is_background_like),nn_background_fraction=mean(nn_is_background_like),
    reasenberg_background=sum(reasenberg_is_background_like),reasenberg_background_fraction=mean(reasenberg_is_background_like),
    all_three_background=sum(gk_is_background_like&nn_is_background_like&reasenberg_is_background_like),
    all_three_clustered=sum(!gk_is_background_like&!nn_is_background_like&!reasenberg_is_background_like),
    any_method_disagreement=sum((as.integer(gk_is_background_like)+as.integer(nn_is_background_like)+as.integer(reasenberg_is_background_like)) %in% c(1,2)),
    disagreement_fraction=any_method_disagreement/n_total,.groups="drop")
print(period_summary)
save_csv_mc(period_summary,"diagnostic_2019_vs_non2019_declustering")

agreement_by_period <- event_level_compare %>% mutate(period_group=ifelse(year==2019,"2019","non-2019")) %>%
  group_split(period_group) %>% map_dfr(function(z) bind_rows(
    agreement_pair(z$gk_is_background_like,z$nn_is_background_like,"Gardner-Knopoff","Nearest-neighbour"),
    agreement_pair(z$gk_is_background_like,z$reasenberg_is_background_like,"Gardner-Knopoff","Reasenberg-style"),
    agreement_pair(z$nn_is_background_like,z$reasenberg_is_background_like,"Nearest-neighbour","Reasenberg-style")) %>%
      mutate(period_group=unique(z$period_group)))
print(agreement_by_period)
save_csv_mc(agreement_by_period,"diagnostic_pairwise_agreement_2019_vs_non2019")

#  MONTHLY ETAS INITIAL BACKGROUND RATE (events/day)
monthly_background_rate <- function(df,background_col,method_name) {
  tmp <- df %>% mutate(month=floor_date(datetime,"month"),bg=as.numeric(.data[[background_col]])) %>%
    group_by(month) %>% summarise(total_events=n(),background_events=sum(bg),.groups="drop")
  full <- tibble(month=seq(floor_date(min(catalog_main$datetime),"month"),
                           floor_date(max(catalog_main$datetime),"month"),by="month"))
  full %>% left_join(tmp,by="month") %>% mutate(
    total_events=replace_na(total_events,0),background_events=replace_na(background_events,0),
    days_in_month=days_in_month(month),
    background_rate_per_day=background_events/days_in_month,
    background_count_per_month=background_events,method=method_name)
}

bg_rate_gk_monthly <- monthly_background_rate(declust_gk,"gk_is_background_like","Gardner-Knopoff")
bg_rate_nn_monthly <- monthly_background_rate(declust_nn,"nn_is_background_like","Nearest-neighbour")
bg_rate_reasenberg_monthly <- monthly_background_rate(declust_reasenberg,"reasenberg_is_background_like","Reasenberg-style")
bg_rate_three_methods <- bind_rows(bg_rate_gk_monthly,bg_rate_nn_monthly,bg_rate_reasenberg_monthly)
save_csv_mc(bg_rate_three_methods,"initial_background_rate_three_methods")

p_bg_three <- ggplot(bg_rate_three_methods,aes(month,background_rate_per_day,linetype=method))+
  geom_line(linewidth=.75)+labs(title="Comparison of Declustering-based Initial Background Rates",
                                subtitle=paste0("Common input catalogue: M >= ",mc_main),x="Time",y="Background-like events per day",
                                linetype="Declustering method")+theme_minimal()
print(p_bg_three); save_plot_mc(p_bg_three,"initial_background_rate_three_methods",10,5)


# DIAGNOSTIC 4: QUANTIFY INITIALIZATION DIFFERENCES
bgwide <- bg_rate_three_methods %>% select(month,method,background_rate_per_day) %>%
  pivot_wider(names_from=method,values_from=background_rate_per_day)
init_metrics <- function(x,y,n1,n2) data.frame(method_1=n1,method_2=n2,
                                               MAE=mean(abs(x-y),na.rm=TRUE),RMSE=sqrt(mean((x-y)^2,na.rm=TRUE)),
                                               correlation=cor(x,y,use="complete.obs"))
init_difference_summary <- bind_rows(
  init_metrics(bgwide[["Gardner-Knopoff"]],bgwide[["Nearest-neighbour"]],"Gardner-Knopoff","Nearest-neighbour"),
  init_metrics(bgwide[["Gardner-Knopoff"]],bgwide[["Reasenberg-style"]],"Gardner-Knopoff","Reasenberg-style"),
  init_metrics(bgwide[["Nearest-neighbour"]],bgwide[["Reasenberg-style"]],"Nearest-neighbour","Reasenberg-style"))
print(init_difference_summary)
save_csv_mc(init_difference_summary,"diagnostic_initialisation_difference_metrics")

bg2019 <- bg_rate_three_methods %>% filter(year(month)==2019) %>%
  select(month,method,background_rate_per_day) %>% pivot_wider(names_from=method,values_from=background_rate_per_day)
init_difference_2019 <- bind_rows(
  init_metrics(bg2019[["Gardner-Knopoff"]],bg2019[["Nearest-neighbour"]],"Gardner-Knopoff","Nearest-neighbour"),
  init_metrics(bg2019[["Gardner-Knopoff"]],bg2019[["Reasenberg-style"]],"Gardner-Knopoff","Reasenberg-style"),
  init_metrics(bg2019[["Nearest-neighbour"]],bg2019[["Reasenberg-style"]],"Nearest-neighbour","Reasenberg-style")) %>%
  mutate(period="2019")
print(init_difference_2019)
save_csv_mc(init_difference_2019,"diagnostic_initialisation_difference_metrics_2019")

# ============================================================
# STEP 4: NON-STATIONARY TEMPORAL ETAS
# ============================================================
# model test
# MODEL COMPARISON 1 (Homogeneous Poisson baseline vs ETAS)
catalog_model <- catalog_main 
t0 <- min(catalog_model$datetime)
t_end <- max(catalog_model$datetime)
catalog_model <- catalog_model %>%
  mutate(time_days = as.numeric(difftime(datetime, t0, units = "days")))
N <- nrow(catalog_model)
T_days <- as.numeric(difftime(t_end, t0, units = "days"))
lambda_hat <- N / T_days
logLik_poisson <- N * log(lambda_hat) - lambda_hat * T_days
k_poisson <- 1
AIC_poisson <- -2 * logLik_poisson + 2 * k_poisson
BIC_poisson <- -2 * logLik_poisson + log(N) * k_poisson
poisson_results <- data.frame(model = "Homogeneous Poisson",n_events = N,
                              duration_days = T_days,lambda_per_day = lambda_hat,
                              logLik = logLik_poisson,k = k_poisson,AIC = AIC_poisson,
                              BIC = BIC_poisson)
print(poisson_results)
write.csv(poisson_results,"model_comparison_poisson.csv",row.names = FALSE)
# MODEL COMPARISON 2
# Stationary temporal ETAS 
etas_time <- catalog_model$time_days
etas_mag  <- catalog_model$mag

# Mc
mc_etas <- 2.0

threshold_mag <- mc_etas
reference_mag <- mc_etas


param_initial <- c( 0.5, 0.05,0.01,1.1,1.1)
t_start <- min(etas_time)
t_end   <- max(etas_time)

fit_etas_stationary <- etasap(time = etas_time, mag = etas_mag,
                              threshold = threshold_mag,reference = reference_mag,
                              parami = param_initial,tstart = t_start,zte = t_end,
                              approx = 2, plot = TRUE)
print(fit_etas_stationary)

logLik_stationary <- -fit_etas_stationary$ngmle
# Five free parameters:
# mu, K, c, alpha, p
k_stationary <- 5
AIC_stationary <- -2 * logLik_stationary +2 * k_stationary
BIC_stationary <- -2 * logLik_stationary +log(N) * k_stationary
stationary_results <- data.frame( model = "Stationary ETAS",
                                  logLik = logLik_stationary, k = k_stationary,
                                  AIC = AIC_stationary,BIC = BIC_stationary)
print(stationary_results)
print(fit_etas_stationary$param)
# ============================================================
# AIC / BIC model comparison
# ============================================================
model_comparison_basic <- data.frame(model = c("Homogeneous Poisson", "Stationary ETAS"),
  logLik = c( logLik_poisson, logLik_stationary),
  k = c(k_poisson, k_stationary)) %>%
  mutate(AIC = -2 * logLik + 2 * k, BIC = -2 * logLik + log(N) * k,
         delta_AIC = AIC - min(AIC),delta_BIC = BIC - min(BIC))
print(model_comparison_basic)
catalog_model <- catalog_main %>%
  filter(!is.na(datetime),!is.na(mag)) %>%
  arrange(datetime)

mc_etas <- 2.0
t0 <- min(catalog_model$datetime)
catalog_model <- catalog_model %>%
  mutate(
    time_days = as.numeric(
      difftime(datetime, t0, units = "days")
    )
  )

times <- catalog_model$time_days
mags  <- catalog_model$mag

N <- length(times)
T_end <- max(times)

cat("N =", N, "\n")
cat("Duration =", T_end, "days\n")

trigger_integral <- function(times,mags,T_end, K,c,alpha, p,Mc) {
  dt <- T_end - times
  productivity <- K * exp(alpha * (mags - Mc))
  if (abs(p - 1) < 1e-6) {
    integral <-
      productivity *
      log((dt + c) / c)
    
  } else {
    integral <-productivity * ((dt + c)^(1 - p) - c^(1 - p)) /(1 - p)
  }
  sum(integral)
}
etas_nll_stationary <- function(par, times, mags,T_end, Mc) {
  mu    <- exp(par[1])
  K     <- exp(par[2])
  c     <- exp(par[3])
  alpha <- exp(par[4])
  p     <- exp(par[5])
  N <- length(times)
  lambda <- numeric(N)
  for (j in seq_len(N)) {
    if (j == 1) {
      lambda[j] <- mu
    } else {
      dt <- times[j] - times[1:(j - 1)]
      triggering <- K * exp( alpha *(mags[1:(j - 1)] - Mc) ) / (dt + c)^p
      lambda[j] <- mu +sum(triggering)
    }
  }
  if (any(!is.finite(lambda)) ||any(lambda <= 0)) {
    return(1e100)
  }
  integrated_background <-mu * T_end
  integrated_triggering <-  trigger_integral(times = times, mags = mags,
                                             T_end = T_end,K = K,c = c,
                                             alpha = alpha,p = p,Mc = Mc)
  logLik <- sum(log(lambda)) - integrated_background -integrated_triggering
  return(-logLik)
}

initial_stationary <- log(c(mu = N / T_end * 0.2, K = 0.05,c = 0.01,alpha = 1.0,p = 1.1))
fit_stationary <- optim(par = initial_stationary,fn = etas_nll_stationary,times = times,
                        mags = mags,T_end = T_end,Mc = mc_etas,method = "BFGS",
                        control = list(maxit = 1000,reltol = 1e-8),hessian = TRUE)
fit_stationary$convergence
fit_stationary$value
stationary_parameters <-exp(fit_stationary$par)
names(stationary_parameters) <-c("mu", "K", "c", "alpha","p")
print(stationary_parameters)
logLik_stationary_custom <-  -fit_stationary$value
k_stationary_custom <- 5
segment_years <- 5
segment_days <-segment_years * 365.25
segment_breaks <- seq( from = 0,to = ceiling(T_end / segment_days) *segment_days,
                       by = segment_days)
if (
  tail(segment_breaks, 1) < T_end
) {
  segment_breaks <-c(segment_breaks, T_end)
}
J <- length(segment_breaks) - 1

cat("Number of background-rate segments =",J,"\n")

print(segment_breaks)
segment_id <- findInterval(times, vec = segment_breaks,rightmost.closed = TRUE)
segment_id[segment_id > J] <- J
etas_nll_nonstationary <- function(par,times, mags, T_end,Mc,segment_breaks,segment_id) {
  J <- length(segment_breaks) - 1
  mu_segments <-exp(par[1:J])
  K <- exp(par[J + 1])
  c <- exp(par[J + 2])
  alpha <- exp(par[J + 3])
  p <-exp(par[J + 4])
  N <- length(times)
  lambda <- numeric(N)
  for (j in seq_len(N)) {
    mu_t <- mu_segments[ segment_id[j]]
    if (j == 1) {lambda[j] <- mu_t }
    else {dt <- times[j] -times[1:(j - 1)]
    triggering <- K *exp(alpha * ( mags[1:(j - 1)] - Mc)) /(dt + c)^p
    lambda[j] <-mu_t + sum(triggering)}}
  if (
    any(!is.finite(lambda)) ||
    any(lambda <= 0)
  ) {
    return(1e100)
  }

  segment_lengths <-diff(segment_breaks)
  segment_lengths[J] <-T_end -segment_breaks[J]
  integrated_background <-sum(mu_segments * segment_lengths)

  integrated_triggering <- trigger_integral(times = times,mags = mags,
                                            T_end = T_end,K = K,c = c,alpha = alpha,
                                            p = p,Mc = Mc)
  logLik <-sum(log(lambda)) - integrated_background -integrated_triggering
  return(-logLik)
}

segment_counts <- table( factor(segment_id,levels = 1:J))
segment_lengths <-diff(segment_breaks)
segment_lengths[J] <- T_end -segment_breaks[J]
raw_segment_rates <-as.numeric(segment_counts) /segment_lengths
initial_mu_segments <- pmax(raw_segment_rates * 0.2,1e-6 )
initial_nonstationary <- c(log(initial_mu_segments),log(0.05),log(0.01),log(1.0),log(1.1))
fit_nonstationary <- optim( par = initial_nonstationary,fn = etas_nll_nonstationary,
                            times = times,mags = mags,T_end = T_end, Mc = mc_etas,
                            segment_breaks = segment_breaks,segment_id = segment_id,
                            method = "BFGS",control = list(maxit = 1500, reltol = 1e-8),
                            hessian = TRUE)
fit_nonstationary$convergence
fit_nonstationary$value
J <- length(segment_breaks) - 1
mu_nonstationary <-exp(fit_nonstationary$par[1:J ])
trigger_par_nonstationary <-exp( fit_nonstationary$par[(J + 1):(J + 4) ])
names(trigger_par_nonstationary) <- c("K","c","alpha","p")
print(mu_nonstationary)
print(trigger_par_nonstationary)
logLik_nonstationary <- -fit_nonstationary$value
k_nonstationary <- J + 4
cat( "Stationary ETAS logLik =",logLik_stationary_custom,"\n")
cat("Non-stationary ETAS logLik =",logLik_nonstationary,"\n")
cat("Stationary k =",k_stationary_custom, "\n")
cat( "Non-stationary k =",k_nonstationary,"\n")
etas_model_comparison <- data.frame(model = c("Stationary ETAS","Non-stationary ETAS"),
                                    logLik = c(logLik_stationary_custom,logLik_nonstationary),
                                    k = c(k_stationary_custom,k_nonstationary))
etas_model_comparison <- etas_model_comparison %>%
  mutate(AIC = -2 * logLik + 2 * k, BIC = -2 * logLik +log(N) * k,
         delta_AIC = AIC - min(AIC), delta_BIC =BIC - min(BIC))
print(etas_model_comparison)
write.csv(etas_model_comparison, "stationary_vs_nonstationary_ETAS_Mc_2.csv",row.names = FALSE)
background_rate_results <- data.frame(segment = 1:J, start_day =segment_breaks[1:J],
    
    end_day = pmin( segment_breaks[ 2:(J + 1)],T_end),
    mu = mu_nonstationary)

background_rate_results <-background_rate_results %>%
  mutate(start_date =as.POSIXct( t0 + start_day *24 * 3600,origin = "1970-01-01"),
    end_date = as.POSIXct( t0 +end_day * 24 * 3600,origin = "1970-01-01"))
print(background_rate_results)
fit_stationary$convergence
fit_nonstationary$convergence
mc_main <- 2.0

mu_df_main <- 8
penalty_lambda_main <- 1

max_trigger_days_main <- 3650
maxit_main <- 300

etas_cat <- catalog %>%
  filter(!is.na(datetime),
         !is.na(mag),
         mag >= mc_main) %>%
  arrange(datetime)

t_origin <- min(etas_cat$datetime)
t_end_datetime <- max(etas_cat$datetime)

etas_cat <- etas_cat %>%
  mutate(event_id = row_number(),
         t_days = as.numeric(difftime(datetime, t_origin, units = "days")),
         mag_excess = mag - mc_main)

t_start <- min(etas_cat$t_days)
t_end <- max(etas_cat$t_days)
T_days <- t_end - t_start
N <- nrow(etas_cat)

cat("ETAS catalogue:", N, "events\n")
cat("Mc =", mc_main, "\n")
cat("Duration =", round(T_days, 2), "days\n")

stopifnot(all(etas_cat$mag >= mc_main))
make_spline_basis <- function(etas_cat, mu_df) {
  t_start <- min(etas_cat$t_days)
  t_end <- max(etas_cat$t_days)
  B_event <- splines::ns(etas_cat$t_days,df = mu_df,intercept = TRUE,Boundary.knots = c(t_start, t_end))
  t_grid <- seq(t_start, t_end, by = 1)
  if (tail(t_grid, 1) < t_end) {t_grid <- c(t_grid, t_end)}
  B_grid <- predict( B_event,newx = t_grid)
  list( B_event = B_event,B_grid = B_grid,t_grid = t_grid)
}
basis_main <- make_spline_basis( etas_cat = etas_cat, mu_df = mu_df_main)
B_event <- basis_main$B_event
B_grid <- basis_main$B_grid
t_grid <- basis_main$t_grid
make_initial_beta_from_bg_rate <- function( bg_rate_monthly, B_event, t_origin,
                                            t_end_datetime, pseudo_count = 0.1,ridge = 1e-4) {
  bg <- bg_rate_monthly %>%
    mutate( month_start = as.POSIXct(month), month_end = month_start %m+% months(1),
            month_start_clip = if_else( month_start < t_origin, t_origin, month_start),
            month_end_clip = if_else(month_end > t_end_datetime,t_end_datetime, month_end),
            exposure_days = as.numeric( difftime( month_end_clip,month_start_clip, units = "days")),
      month_mid =
        month_start_clip +
        (month_end_clip - month_start_clip) / 2,
      t_mid_days = as.numeric( difftime(month_mid, t_origin,units = "days")),
      initial_rate_per_day =(background_events + pseudo_count) / pmax(exposure_days, 1),
      log_initial_rate = log(initial_rate_per_day) ) %>%
    filter(exposure_days > 0,
           is.finite(t_mid_days),
           t_mid_days >= 0,
           is.finite(log_initial_rate))
  B_mid <- predict(B_event,newx = bg$t_mid_days)
  
  XtX <- crossprod(B_mid)
  Xty <- crossprod(B_mid,bg$log_initial_rate)
  beta_init <- solve(XtX +ridge * diag(ncol(B_mid)),Xty)
  as.numeric(beta_init)
}
beta_init_gk <- make_initial_beta_from_bg_rate(bg_rate_gk_monthly,B_event,t_origin,t_end_datetime)
beta_init_nn <- make_initial_beta_from_bg_rate(bg_rate_nn_monthly,B_event,t_origin,t_end_datetime)
beta_init_reasenberg <- make_initial_beta_from_bg_rate(bg_rate_reasenberg_monthly,B_event,
                                                       t_origin,t_end_datetime)

trapz_integral <- function(x, y) {
  if (length(x) < 2) {
    return(0)
  }
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}
# ETAS LIKELIHOOD COMPONENTS
etas_components <- function(par,etas_cat,B_event,B_grid,t_grid,max_trigger_days = 3650,penalty_lambda = 1) {
  logK <- par[1]
  log_alpha <- par[2]
  log_c <- par[3]
  log_p_minus_1 <- par[4]
  beta <- par[-(1:4)]
  K <- exp(logK)
  alpha <- exp(log_alpha)
  c_par <- exp(log_c)
  p_par <- 1 + exp(log_p_minus_1)
  t <- etas_cat$t_days
  m_excess <- etas_cat$mag_excess
  n <- length(t)
  mu_event <- as.numeric(exp(B_event %*% beta))
  mu_grid <- as.numeric(exp(B_grid %*% beta))
  bg_integral <- trapz_integral(t_grid,mu_grid)
  trig_event <- numeric(n)
  if (n >= 2) {
    for (i in 2:n) {
      dt <- t[i] - t[1:(i - 1)]
      keep <- ( dt > 0 & dt <= max_trigger_days)
      if (any(keep)) {
        prev_idx <- which(keep)
        trig_event[i] <- sum(K * exp(alpha * m_excess[prev_idx]) *
                               (dt[prev_idx] +c_par)^(-p_par))
      }
    }
  }
  lambda_event <-mu_event +trig_event
  if (
    any(!is.finite(lambda_event)) ||
    any(lambda_event <= 0)
  ) {
    
    return(
      list(
        valid = FALSE,
        penalized_objective = 1e100
      )
    )
  }
  
  event_loglik <-sum(log(lambda_event))
  T_end <- max(t)
  upper_dt <- pmin(max_trigger_days,pmax(T_end - t, 0))
  productivity <-K * exp(alpha * m_excess)
  trigger_integral_each <-productivity * 
    (c_par^(1 - p_par) -(upper_dt + c_par)^(1 - p_par)) /(p_par - 1)
  
  trigger_integral <- sum(trigger_integral_each)
  loglik <-event_loglik -bg_integral -trigger_integral
  penalty <- penalty_lambda * sum(diff(beta,differences = 2)^2)
  
  penalized_objective <- -loglik +penalty
  
  list(valid = is.finite(penalized_objective),loglik = loglik,negloglik = -loglik,
       penalty = penalty,penalized_objective = penalized_objective,K = K,alpha = alpha,
       c = c_par,p = p_par, beta = beta,mu_event = mu_event, mu_grid = mu_grid,
       bg_integral = bg_integral,trigger_integral = trigger_integral)
}


etas_objective <- function(par, etas_cat,B_event, B_grid,t_grid,max_trigger_days,
                           penalty_lambda) {
  comp <- etas_components(par = par,etas_cat = etas_cat,B_event = B_event,
                          B_grid = B_grid,t_grid = t_grid,
                          max_trigger_days = max_trigger_days,
                          penalty_lambda = penalty_lambda)
  if (!isTRUE(comp$valid)) {
    return(1e100)
  }
  comp$penalized_objective
}

fit_etas_from_start <- function(beta_init,method_name,etas_cat,B_event,B_grid,t_grid,
                                max_trigger_days = 3650,penalty_lambda = 1,maxit = 300,
                                trace = 1) {
  K0 <- 0.02
  alpha0 <- 1.0
  c0 <- 0.01
  p0 <- 1.1
  
  par0 <- c(log(K0),log(alpha0),log(c0),log(p0 - 1), beta_init)
  lower <- c( log(1e-6),log(0.05), log(1e-4),log(0.001),rep(-20, length(beta_init)))
  upper <- c(log(10), log(5),log(10),log(5),rep(5, length(beta_init)))
  cat( "\n----------------------------------------\n",
       "Fitting ETAS from: ", method_name, "\n",
       "----------------------------------------\n",
       sep = "")
  start_time <- Sys.time()
  fit <- optim(par = par0,fn = etas_objective,method = "L-BFGS-B",lower = lower,
               upper = upper,control = list(maxit = maxit,trace = trace,REPORT = 5),
               etas_cat = etas_cat,B_event = B_event,B_grid = B_grid,t_grid = t_grid,
               max_trigger_days = max_trigger_days, penalty_lambda = penalty_lambda)
  runtime_min <- as.numeric(difftime(Sys.time(),start_time,units = "mins"))
  final_comp <- etas_components(par = fit$par,etas_cat = etas_cat,B_event = B_event,
                                B_grid = B_grid,t_grid = t_grid,
                                max_trigger_days = max_trigger_days,
                                penalty_lambda = penalty_lambda)
  list(method = method_name,fit = fit,par_hat = fit$par,K = final_comp$K,
       alpha = final_comp$alpha,c = final_comp$c,p = final_comp$p,beta = final_comp$beta,
       logLik = final_comp$loglik, negloglik = final_comp$negloglik,
       penalty = final_comp$penalty, penalized_objective =final_comp$penalized_objective,
       mu_grid = final_comp$mu_grid,bg_integral =final_comp$bg_integral,
       trigger_integral =final_comp$trigger_integral,convergence = fit$convergence,
       message = fit$message,runtime_minutes = runtime_min)
}


# MAIN EXPERIMENT:
fit_gk <- fit_etas_from_start(beta_init = beta_init_gk,method_name = "Gardner-Knopoff",
                              etas_cat = etas_cat,B_event = B_event,B_grid = B_grid,
                              t_grid = t_grid, max_trigger_days = max_trigger_days_main,
                              penalty_lambda = penalty_lambda_main,maxit = maxit_main)
fit_nn <- fit_etas_from_start(beta_init = beta_init_nn,method_name = "Nearest-neighbour",
                              etas_cat = etas_cat,B_event = B_event,B_grid = B_grid,
                              t_grid = t_grid,max_trigger_days = max_trigger_days_main,
                              penalty_lambda = penalty_lambda_main,maxit = maxit_main)
fit_reasenberg <- fit_etas_from_start(beta_init = beta_init_reasenberg,
                                      method_name = "Reasenberg-style",
                                      etas_cat = etas_cat, B_event = B_event,
                                      B_grid = B_grid,t_grid = t_grid,
                                      max_trigger_days = max_trigger_days_main,
                                      penalty_lambda = penalty_lambda_main,
                                      maxit = maxit_main)

#  COMPARE FINAL TRIGGERING PARAMETERS AND LIKELIHOOD
etas_fit_summary <- bind_rows(data.frame( method = fit_gk$method,
                                          convergence = fit_gk$convergence,
                                          logLik = fit_gk$logLik,
                                          negloglik = fit_gk$negloglik,
                                          penalty = fit_gk$penalty,
                                          penalized_objective =fit_gk$penalized_objective,
                                          K = fit_gk$K,alpha = fit_gk$alpha,c = fit_gk$c,
                                          p = fit_gk$p,runtime_minutes =fit_gk$runtime_minutes),
  data.frame(method = fit_nn$method,convergence = fit_nn$convergence,logLik = fit_nn$logLik,
             negloglik = fit_nn$negloglik,penalty = fit_nn$penalty,
             penalized_objective =fit_nn$penalized_objective,K = fit_nn$K,
             alpha = fit_nn$alpha,c = fit_nn$c,p = fit_nn$p,runtime_minutes =fit_nn$runtime_minutes),
  data.frame(method = fit_reasenberg$method,convergence =fit_reasenberg$convergence,
             logLik =fit_reasenberg$logLik,negloglik =fit_reasenberg$negloglik,
             penalty =fit_reasenberg$penalty,penalized_objective =fit_reasenberg$penalized_objective,
             K =fit_reasenberg$K, alpha =fit_reasenberg$alpha,c =fit_reasenberg$c,
             p = fit_reasenberg$p,runtime_minutes =fit_reasenberg$runtime_minutes))
print(etas_fit_summary)
write.csv(etas_fit_summary,"nonstationary_ETAS_fit_summary_three_initializations_Mc_2_0.csv",row.names = FALSE)

mu_final <- bind_rows(data.frame(t_days = t_grid,mu = fit_gk$mu_grid,method = "Gardner-Knopoff"),
                      data.frame(t_days = t_grid, mu = fit_nn$mu_grid,method = "Nearest-neighbour"),
  data.frame(t_days = t_grid,mu = fit_reasenberg$mu_grid, method = "Reasenberg-style")) %>%
  mutate(datetime =t_origin + t_days * 24 * 3600)
write.csv(mu_final,"nonstationary_ETAS_final_background_rates_Mc_2_0.csv",row.names = FALSE)
p_mu_final <- ggplot(mu_final,aes( x = datetime, y = mu,linetype = method)) +
  geom_line(linewidth = 0.75) +
  labs( title ="Final Fitted Non-stationary ETAS Background Rates",
        subtitle ="Three declustering-based starting values; common ETAS model",
        x = "Time",y = "Fitted background rate (events/day)",linetype = "Initialisation") +
  theme_minimal()
print(p_mu_final)
ggsave("nonstationary_ETAS_final_background_rates_Mc_2_0.png", p_mu_final, width = 10,
       height = 5,dpi = 300)

# QUANTIFY FINAL BACKGROUND-RATE DIFFERENCES
trajectory_metrics <- function(x,y,method_1,method_2) {
  data.frame(method_1 = method_1,method_2 = method_2,
             MAE =mean(abs(x - y),na.rm = TRUE),
             RMSE =sqrt(mean((x - y)^2,na.rm = TRUE)),
             correlation =cor(x,y,use = "complete.obs"))
}
mu_final_wide <- mu_final %>%
  select(t_days,method,mu) %>%
  pivot_wider(names_from = method, values_from = mu)
final_mu_difference_summary <- bind_rows(
  trajectory_metrics(mu_final_wide[["Gardner-Knopoff"]],
                     mu_final_wide[["Nearest-neighbour"]],
                     "Gardner-Knopoff","Nearest-neighbour" ),
 trajectory_metrics(mu_final_wide[["Gardner-Knopoff"]],
                    mu_final_wide[["Reasenberg-style"]],
                    "Gardner-Knopoff","Reasenberg-style"),
  trajectory_metrics(mu_final_wide[["Nearest-neighbour"]],
                     mu_final_wide[["Reasenberg-style"]],
                     "Nearest-neighbour","Reasenberg-style"))
print(final_mu_difference_summary)
write.csv(final_mu_difference_summary,
          "nonstationary_ETAS_final_background_difference_metrics_Mc_2_0.csv",
          row.names = FALSE)
# COMPARE SMOOTH STARTING AND FINAL BACKGROUND FUNCTIONS
# Both sets of functions are evaluated at the same daily time points.

mu_initial <- bind_rows(data.frame(t_days = t_grid,
                                   mu = as.numeric(exp(B_grid %*% beta_init_gk)),
                                   method = "Gardner-Knopoff"),
                        data.frame(t_days = t_grid,
                                   mu = as.numeric(exp(B_grid %*% beta_init_nn)),
                                   method = "Nearest-neighbour"),
                        data.frame(t_days = t_grid,
                                   mu = as.numeric(exp(B_grid %*% beta_init_reasenberg)),
                                   method = "Reasenberg-style")) %>%
  mutate(datetime = t_origin + t_days * 24 * 3600)

write.csv(mu_initial,"nonstationary_ETAS_initial_background_rates_Mc_2_0.csv",row.names = FALSE)
mu_initial_wide <- mu_initial %>%
  select(t_days, method, mu) %>%
  pivot_wider(names_from = method, values_from = mu)

initial_mu_difference_summary <- bind_rows(
  trajectory_metrics(
    mu_initial_wide[["Gardner-Knopoff"]],
    mu_initial_wide[["Nearest-neighbour"]],
    "Gardner-Knopoff",
    "Nearest-neighbour"
  ),
  trajectory_metrics(
    mu_initial_wide[["Gardner-Knopoff"]],
    mu_initial_wide[["Reasenberg-style"]],
    "Gardner-Knopoff",
    "Reasenberg-style"
  ),
  trajectory_metrics(
    mu_initial_wide[["Nearest-neighbour"]],
    mu_initial_wide[["Reasenberg-style"]],
    "Nearest-neighbour",
    "Reasenberg-style"
  )
)

print(initial_mu_difference_summary)

write.csv(initial_mu_difference_summary,
          "nonstationary_ETAS_initial_background_difference_metrics_Mc_2_0.csv",row.names = FALSE)

propagation_summary <- initial_mu_difference_summary %>%
  rename(initial_MAE = MAE,initial_RMSE = RMSE,initial_correlation = correlation) %>%
  left_join(final_mu_difference_summary %>%
              rename(final_MAE = MAE,final_RMSE = RMSE,final_correlation = correlation),
            by = c("method_1", "method_2")) %>%
  mutate(MAE_ratio_final_to_initial = final_MAE / initial_MAE,
         RMSE_ratio_final_to_initial = final_RMSE / initial_RMSE,
         MAE_remaining_percentage =100 * MAE_ratio_final_to_initial,
         RMSE_remaining_percentage =100 * RMSE_ratio_final_to_initial,
         MAE_attenuation_percentage =100 * (1 - MAE_ratio_final_to_initial),
         RMSE_attenuation_percentage =100 * (1 - RMSE_ratio_final_to_initial))
print(propagation_summary)
write.csv(propagation_summary,"ETAS_initial_to_final_background_propagation_Mc_2_0.csv",row.names = FALSE)

cat("\nConvergence codes:\n")
print(etas_fit_summary %>%
        select(method,convergence,logLik,penalized_objective))
if (
  any(
    etas_fit_summary$convergence != 0
  )
) {
  warning(
    "At least one ETAS fit did not return convergence code 0."
  )
}

# SMALL ETAS MODEL-SPECIFICATION ROBUSTNESS CHECK
baseline_df <- 8
baseline_penalty <- 1
baseline_cutoff <- 3650

cat("SMALL ETAS SPECIFICATION ROBUSTNESS\n","Reference initialization: Nearest-neighbour\n",
    "Baseline: df=8, penalty=1, cutoff=3650 days\n","========================================\n",
    sep = "")

cat("\n[1/3] Fitting robustness model: df = 6\n")

basis_df6 <- make_spline_basis(etas_cat = etas_cat,mu_df = 6)
beta_init_df6 <- make_initial_beta_from_bg_rate(
  bg_rate_monthly = bg_rate_nn_monthly,
  B_event = basis_df6$B_event,
  t_origin = t_origin,
  t_end_datetime = t_end_datetime)

fit_df6 <- fit_etas_from_start( beta_init = beta_init_df6,
                                method_name = "NN_start_df6",etas_cat = etas_cat,
                                B_event = basis_df6$B_event,B_grid = basis_df6$B_grid,
                                t_grid = basis_df6$t_grid,max_trigger_days = baseline_cutoff,
                                penalty_lambda = baseline_penalty,maxit = 300,trace = 1)
cat("\n[2/3] Fitting robustness model: df = 10\n")
basis_df10 <- make_spline_basis(etas_cat = etas_cat, mu_df = 10)
beta_init_df10 <- make_initial_beta_from_bg_rate(
  bg_rate_monthly = bg_rate_nn_monthly,
  B_event = basis_df10$B_event,
  t_origin = t_origin,
  t_end_datetime = t_end_datetime)

fit_df10 <- fit_etas_from_start(beta_init = beta_init_df10, method_name = "NN_start_df10",
                                etas_cat = etas_cat,B_event = basis_df10$B_event,
                                B_grid = basis_df10$B_grid,t_grid = basis_df10$t_grid,
                                max_trigger_days = baseline_cutoff,penalty_lambda = baseline_penalty,
                                maxit = 300,trace = 1)
cat("\n[3/3] Fitting robustness model: cutoff = 7300 days\n")

basis_cut7300 <- make_spline_basis(etas_cat = etas_cat,mu_df = baseline_df)
beta_init_cut7300 <- make_initial_beta_from_bg_rate(
  bg_rate_monthly = bg_rate_nn_monthly,
  B_event = basis_cut7300$B_event,
  t_origin = t_origin,
  t_end_datetime = t_end_datetime)
fit_cut7300 <- fit_etas_from_start(beta_init = beta_init_cut7300,
                                   method_name = "NN_start_cutoff7300",
                                   etas_cat = etas_cat,B_event = basis_cut7300$B_event,
                                   B_grid = basis_cut7300$B_grid,t_grid = basis_cut7300$t_grid,
                                   max_trigger_days = 7300,penalty_lambda = baseline_penalty,
                                   maxit = 300,trace = 1)

robustness_summary <- bind_rows(
  data.frame(specification = "Baseline",mu_df = 8,penalty_lambda = 1,max_trigger_days = 3650,
             convergence = fit_nn$convergence,logLik = fit_nn$logLik,
             penalized_objective =fit_nn$penalized_objective,K = fit_nn$K,alpha = fit_nn$alpha,c = fit_nn$c,
             p = fit_nn$p,runtime_minutes =fit_nn$runtime_minutes),
  data.frame(specification = "df_6",mu_df = 6,penalty_lambda = 1,max_trigger_days = 3650,
             convergence = fit_df6$convergence,logLik = fit_df6$logLik,
             penalized_objective =fit_df6$penalized_objective,K = fit_df6$K,alpha = fit_df6$alpha,
             c = fit_df6$c,p = fit_df6$p,runtime_minutes =fit_df6$runtime_minutes),
  data.frame(specification = "df_10",mu_df = 10,penalty_lambda = 1,max_trigger_days = 3650,
             convergence = fit_df10$convergence,logLik = fit_df10$logLik,
             penalized_objective =fit_df10$penalized_objective,K = fit_df10$K,alpha = fit_df10$alpha,
             c = fit_df10$c, p = fit_df10$p,runtime_minutes =fit_df10$runtime_minutes),
  data.frame(specification = "cutoff_7300",mu_df = 8,penalty_lambda = 1,max_trigger_days = 7300,
             convergence =fit_cut7300$convergence,logLik =fit_cut7300$logLik,
             penalized_objective =fit_cut7300$penalized_objective,
             K =fit_cut7300$K,alpha =fit_cut7300$alpha,c =fit_cut7300$c,
             p =fit_cut7300$p,runtime_minutes =fit_cut7300$runtime_minutes)
)
cat("ROBUSTNESS SUMMARY\n")
print(robustness_summary)
write.csv(robustness_summary,"ETAS_small_specification_robustness_Mc_2_0.csv",row.names = FALSE)
baseline_row <-robustness_summary %>%
  filter(specification == "Baseline")
robustness_parameter_changes <-robustness_summary %>%
  mutate(K_pct_change = 100 * (K - baseline_row$K) /baseline_row$K,
         alpha_pct_change = 100 * (alpha - baseline_row$alpha) /baseline_row$alpha,
         c_pct_change =100 *(c - baseline_row$c) /baseline_row$c,
         p_pct_change =100 *(p - baseline_row$p) /baseline_row$p,
         logLik_difference = logLik -baseline_row$logLik)
cat("CHANGE RELATIVE TO BASELINE\n")
print(robustness_parameter_changes %>%
        select(specification, convergence,K_pct_change,alpha_pct_change,c_pct_change,p_pct_change,
               logLik_difference))
write.csv(robustness_parameter_changes,"ETAS_small_specification_parameter_changes_Mc_2_0.csv",
          row.names = FALSE)
mu_robustness <- bind_rows(data.frame(t_days = t_grid,mu = fit_nn$mu_grid,specification = "Baseline"),
                           data.frame(t_days = basis_df6$t_grid,mu = fit_df6$mu_grid,specification = "df = 6"),
                           data.frame(t_days = basis_df10$t_grid,mu = fit_df10$mu_grid,specification = "df = 10"),
                           data.frame(t_days = basis_cut7300$t_grid,mu = fit_cut7300$mu_grid,specification = "cutoff = 7300")) %>%
  mutate(datetime =t_origin + t_days * 24 * 3600)
write.csv(mu_robustness,"ETAS_small_specification_background_rates_Mc_2_0.csv",row.names = FALSE)

p_robustness_mu <- ggplot( mu_robustness,aes(x = datetime, y = mu,linetype = specification)) +
  geom_line(linewidth = 0.75) +
  labs(title = "ETAS Background-Rate Specification Robustness",
       subtitle ="Fixed nearest-neighbour initialization", x = "Time",
       y ="Fitted background rate (events/day)",linetype = "Specification") +
  theme_minimal()
print(p_robustness_mu)
ggsave("ETAS_small_specification_background_rates_Mc_2_0.png",p_robustness_mu,width = 10,
       height = 5,dpi = 300)

trajectory_metrics <- function(x,y,reference,comparison) {
  data.frame(reference = reference, comparison = comparison,
             MAE =mean(abs(x - y),na.rm = TRUE),
             RMSE =sqrt(mean((x - y)^2,na.rm = TRUE)),
             correlation =cor(x, y,use = "complete.obs"))
}
mu_wide <-mu_robustness %>%
  select(t_days,specification,mu) %>%
  pivot_wider(names_from = specification,values_from = mu)
mu_robustness_metrics <- bind_rows(trajectory_metrics(mu_wide[["Baseline"]],
                                                      mu_wide[["df = 6"]],
                                                      "Baseline","df = 6"),
                                   trajectory_metrics(mu_wide[["Baseline"]],mu_wide[["df = 10"]],
                                                      "Baseline","df = 10"),
                                   trajectory_metrics(mu_wide[["Baseline"]], mu_wide[["cutoff = 7300"]],
                                                      "Baseline","cutoff = 7300"))

cat("FINAL mu(t) ROBUSTNESS METRICS\n")
print(mu_robustness_metrics)
write.csv(mu_robustness_metrics,"ETAS_small_specification_mu_difference_metrics_Mc_2_0.csv",
          row.names = FALSE)

if (any(robustness_summary$convergence != 0)) {
  warning("At least one robustness model did not converge.")
} else {
  cat("\nAll robustness models returned convergence code 0.\n")
}
cat("\nSmall ETAS specification robustness analysis completed.\n")
