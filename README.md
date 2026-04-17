# MIS410 Culottes Forecasting Dashboard

High-end Shiny dashboard for forecasting and inventory analysis using Google Trends data for **culottes**.

## Project coverage
- Pattern analysis
- Best forecast model using RMSFE
- 6–12 week forecast
- Safety stock and reorder point
- Business insights
- MySQL write-back support

## Files
- `app.R` — main Shiny application
- `gtrends.csv` — input dataset
- `.gitignore` — recommended Git ignore rules

## Required packages
```r
install.packages(c(
  "shiny", "bslib", "forecast", "tseries", "ggplot2", "dplyr",
  "RMariaDB", "DBI", "scales", "readr", "lubridate", "DT", "gridExtra"
))
```

## Run locally
```r
shiny::runApp()
```

## GitHub upload
```bash
git init
git add .
git commit -m "Initial commit - MIS410 culottes forecasting dashboard"
git branch -M main
git remote add origin https://github.com/YOUR-USERNAME/mis410-culottes-dashboard.git
git push -u origin main
```
