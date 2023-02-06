require(downloader)
library(curl)

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
out_dir <- paste0(args[1], "download")

curl_download(
  "https://orcestradata.blob.core.windows.net/prism/profile.sensitivity.PRISM.rds",
  destfile = file.path(out_dir, "profile.sensitivity.PRISM.rds")
)
curl_download(
  "https://orcestradata.blob.core.windows.net/prism/raw.sensitivity.prismii_v3.rds",
  destfile = file.path(out_dir, "raw.sensitivity.prismii_v3.rds")
)
# Download secondary-screen-replicate-collapsed-treatment-info.csv
curl_download(
  "https://ndownloader.figshare.com/files/20237763",
  destfile = file.path(out_dir, "secondary-screen-replicate-collapsed-treatment-info.csv")
)

# Download secondary-screen-cell-line-info.csv
curl_download(
  "https://ndownloader.figshare.com/files/20237769",
  destfile = file.path(out_dir, "secondary-screen-cell-line-info.csv")
)

# Download secondary-screen-replicate-collapsed-logfold-change.csv
curl_download(
  "https://ndownloader.figshare.com/files/20237757",
  destfile = file.path(out_dir, "secondary-screen-replicate-collapsed-logfold-change.csv")
)

# Download secondary-screen-dose-response-curve-parameters.csv
curl_download(
  "https://ndownloader.figshare.com/files/20237739",
  destfile = file.path(out_dir, "secondary-screen-dose-response-curve-parameters.csv")
)
