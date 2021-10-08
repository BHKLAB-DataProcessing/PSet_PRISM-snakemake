# Secondary screen
work_dir = "~/nikta/out/"
setwd(work_dir)
library(reshape2)
library(PharmacoGx)
library(data.table)
options(encoding = "UTF-8")

############################################# Molecular profile #############################################
# There will be no molecular profile as all the cells are from ccle
emptyE <- ExpressionSet()
pData(emptyE)$cellid <- character()
pData(emptyE)$batchid <- character()
fData(emptyE)$BEST <- vector()
fData(emptyE)$Symbol <- character()
annotation(emptyE) <- "PRISM contains no molecular profiles of cell lines. This SE is empty placeholder."

emptySE <- SummarizedExperiment::SummarizedExperiment(
  ## TODO:: Do we want to pass an environment for better memory efficiency?
  assays=S4Vectors::SimpleList(as.list(Biobase::assayData(emptyE))
  ),
  # Switch rearrange columns so that IDs are first, probes second
  rowData=S4Vectors::DataFrame(Biobase::fData(emptyE),
                               rownames=rownames(Biobase::fData(emptyE)) 
  ),
  colData=S4Vectors::DataFrame(Biobase::pData(emptyE),
                               rownames=rownames(Biobase::pData(emptyE))
  ),
  metadata=list("experimentData" = emptyE@experimentData, 
                "annotation" = Biobase::annotation(emptyE), 
                "protocolData" = Biobase::protocolData(emptyE)
  )
)

############################################# Drug object-part1 #############################################
# Drug data (secondary-screen-replicate-collapsed-treatment-info.csv)
drugs <- read.csv("https://ndownloader.figshare.com/files/20237763" , stringsAsFactors = FALSE)
drug.obj <- unique(drugs[!duplicated(drugs$broad_id),- which(colnames(drugs) %in% c("column_name", "dose","screen_id"))])
colnames(drug.obj)[colnames(drug.obj) == "name"] <- "PRISM.drugid"

# lab durgnames
lab.drug.names <- read.csv(paste0(work_dir,"drug_with_ids.csv"), stringsAsFactors = F , na.strings = "") # this file is "drug_with_ids.csv" from pachy-annotations
lab.drug.names <- lab.drug.names[!is.na(lab.drug.names$PRISM.drugid), c("unique.drugid" , "PRISM.drugid","cid","inchikey")]

for (i in seq(nrow(lab.drug.names))){
  name = lab.drug.names$PRISM.drugid[i]
  if(grepl("///" , name)){
    
    all.name = unlist(strsplit(name, "///"))
    print(all.name)
    for(j in 1:length(all.name)){
      lab.drug.names[1+nrow(lab.drug.names) , ] = lab.drug.names[i , ]
      lab.drug.names$PRISM.drugid[nrow(lab.drug.names)] = all.name[j]
    }
    lab.drug.names <- lab.drug.names[-i,]
  }
}

drug.obj <- merge(drug.obj, lab.drug.names, by="PRISM.drugid") 
colnames(drug.obj)[colnames(drug.obj) == "unique.drugid"] <- "drugid"


############################################# Cell object #############################################
# Lab cell names
lab.cell.names <- read.csv(paste0(work_dir,"cell_annotation_all.csv"), stringsAsFactors = F , na.strings = "") # this file is "cell_annotation_all.csv" from pachy annotation 

# Cell-line data (secondary-screen-cell-line-info.csv)
cells <- read.csv("https://ndownloader.figshare.com/files/20237769" , stringsAsFactors = FALSE)
colnames(cells)[colnames(cells) == "ccle_name"] <- "PRISM.cellid"
cells <- cells[!is.na(cells$depmap_id) , ] # Removing NA depmap ids

cell.obj <- merge(lab.cell.names[!is.na(lab.cell.names$PRISM.cellid), c("unique.cellid", "unique.tissueid","PRISM.cellid","PRISM.tissueid")] , cells, by="PRISM.cellid", all.x=T)
any(is.na(cell.obj$unique.cellid)) # FALSE
any(is.na(cell.obj$unique.tissueid)) # FALSE

colnames(cell.obj)[colnames(cell.obj) == "unique.cellid"] <- "cellid"
colnames(cell.obj)[colnames(cell.obj) == "unique.tissueid"] <- "tissueid"
rownames(cell.obj) <- cell.obj$cellid

cell.obj <- cell.obj[,-which(colnames(cell.obj)=="row_name")]
cell.obj [cell.obj == ""] <- NA

############################################# sensitivity data #############################################
# Collapsed log_fold_change (secondary-screen-replicate-collapsed-logfold-change.csv)
# Columns: IDs identifying experimental conditions. See column name from primary_replicate_collapsed_treatment_info
# Rows: IDs identifying cell lines. See row_name from cell_line_info
dose.resp.file <- tempfile()
download.file("https://ndownloader.figshare.com/files/20237757", dose.resp.file )
dose.resp.data <- read.csv(dose.resp.file , stringsAsFactors = F)

dose.resp <- data.frame(t(dose.resp.data))
colnames(dose.resp) <- dose.resp[1,]
dose.resp <- dose.resp[-1,]
dose.resp$column_name <- gsub(".." , "::" , rownames(dose.resp) , fixed = T) # For merging with drugs
dose.resp$column_name <- gsub("." , "-" , dose.resp$column_name , fixed = T) # For merging with drugs
dose.resp$column_name <- gsub("::0-" , "::0." , dose.resp$column_name , fixed = T) # For merging with drugs
dose.resp$column_name <- gsub("::1-" , "::1." , dose.resp$column_name , fixed = T) # For merging with drugs
dose.resp$column_name <- gsub("::2-" , "::2." , dose.resp$column_name , fixed = T) # For merging with drugs
dose.resp$column_name <- gsub("::3-" , "::3." , dose.resp$column_name , fixed = T) # For merging with drugs
dose.resp$column_name <- gsub("::4-" , "::4." , dose.resp$column_name , fixed = T) # For merging with drugs
dose.resp$column_name <- gsub("::5-" , "::5." , dose.resp$column_name , fixed = T) # For merging with drugs
dose.resp$column_name <- gsub("::8-" , "::8." , dose.resp$column_name , fixed = T) # For merging with drugs
dose.resp$column_name <- gsub("::9-" , "::9." , dose.resp$column_name , fixed = T) # For merging with drugs
dose.resp$column_name <- gsub("::10-" , "::10." , dose.resp$column_name , fixed = T) # For merging with drugs
dose.resp$column_name <- gsub("::11-" , "::11." , dose.resp$column_name , fixed = T) # For merging with drugs
dose.resp$column_name <- gsub("::12-" , "::12." , dose.resp$column_name , fixed = T) # For merging with drugs


# To get the cell_line names for creating expid column
dose.resp<-reshape2::melt(dose.resp , id.vars="column_name")
colnames(dose.resp)[colnames(dose.resp) == "variable"] <- "depmap_id"
colnames(dose.resp)[colnames(dose.resp) == "value"] <- "viab" # represent viability
range(as.numeric(dose.resp$viab) ,na.rm =T) # (-12.065364 , 5.007653)


# Merging dose-resp with drugs to get the doses for each experiment
temp <- merge(drugs[, c("column_name", "broad_id", "dose", "screen_id", "name", "compound_plate")] , dose.resp, by="column_name") 
temp$expid <- paste(temp$broad_id , temp$screen_id , sep="::") 
temp$expid <- paste(temp$expid , temp$depmap_id, sep="::") #depmap_id is cell line name
temp$expid <- paste(temp$expid , temp$name, sep="::") # name is drug

######### Sensitivity-raw #########
sub.dose.resp <- setDT(temp[!is.na(temp$viab), ])
sub.dose.resp <- merge(sub.dose.resp[,.N, expid][N>3], sub.dose.resp) # Removing experiments with less than 4 points
sub.dose.resp <- sub.dose.resp[!grepl("FAILED_STR", depmap_id, fixed = T), ] # Removing failed experiments:

str(sub.dose.resp[, .(dose, viab)]) # They both need to be numeric
sub.dose.resp[,viab := as.numeric(viab)]

# test <- sub.dose.resp[N>8,]
# ftable(table(sub.dose.resp$expid))
# 
# raw.sensitivity <- array(data = NA_real_,
#                          dim= c(length(unique(sub.dose.resp$expid)),max(sub.dose.resp$N), 2))
# 
# colnames(raw.sensitivity) <- paste0('dose:', seq_len(NCOL(raw.sensitivity)))
# rownames(raw.sensitivity) <- unique(sub.dose.resp$expid)
# dimnames(raw.sensitivity)[[3]] <- c("Dose", "Viability") 
# setkey(sub.dose.resp, expid) # defining a"key" column for sorting
# setorder(sub.dose.resp, expid, dose)
# 
# for(exp in unique(sub.dose.resp[,expid])){
#   raw.sensitivity[exp, seq_len(sub.dose.resp[expid==exp, unique(N)]),] <- data.matrix(sub.dose.resp[expid==exp, .(dose, viab)])
# }
# 
# saveRDS(raw.sensitivity , "raw.sensitivity.prismii_v3.rds") # v3 is after removing N<4 and str_failed
raw.sensitivity <- readRDS(paste0(work_dir,"raw.sensitivity.prismii_v3.rds")) 

raw.sensitivity[,,"Viability"] <- 2^(raw.sensitivity[ , ,"Viability"])  # Converting log values to decimals (NO need for logs)
#range(raw.sensitivity[,, "Viability"], na.rm = T) #2.333262e-04 3.217019e+01
#summary(raw.sensitivity[,, "Viability"])

raw.sensitivity[,,"Viability"] <- 100 * (raw.sensitivity[, ,"Viability"]) # Multiplying by 100 --> convert to percentage
#range(raw.sensitivity[,, "Viability"], na.rm = T) #2.333262e-02 3.217019e+03
#summary(raw.sensitivity[,, "Viability"])

#Submitting a job on h4h for aac value alculations
# saveRDS(raw.sensitivity , "raw.sensitivity.prismiiH4H.rds")
# # data-slicing - conducted on h4h4 #
# raw.sensitivity <- readRDS("/cluster/projects/bhklab/Data/PRISM/raw.sensitivity.prismiiH4H.rds")
# dir.create("/cluster/projects/bhklab/Data/PRISM/out/")
# dir.create("/cluster/projects/bhklab/Data/PRISM/out/slices/")
# 
# setwd("/cluster/projects/bhklab/Data/PRISM/out/")
# raw.sensitivity.x <- parallel::splitIndices(nrow(raw.sensitivity), floor(nrow(raw.sensitivity)/10000))
# for(i in seq_along(raw.sensitivity.x)){
#   slce <- raw.sensitivity[raw.sensitivity.x[[i]],,]
#   saveRDS(slce, file=paste0("slices/PRISM_raw_sens_", i, ".rds"))
# }
# 
# dim(raw.sensitivity) # 739776,40,2 

# #### h4h job ####
# library(PharmacoGx)
# args <- commandArgs(trailingOnly = TRUE)
# jobIndex <- as.numeric(args[[1]])
# myfn <- list.files("/cluster/projects/bhklab/Data/PRISM/out/slices/", full.names=TRUE)[[jobIndex]]
# print(myfn)
# mybasenm <- basename(myfn)
# slice <- readRDS(myfn)
# res <- PharmacoGx:::.calculateFromRaw(slice)
# saveRDS(res, file=paste0("/cluster/projects/bhklab/Data/PRISM/out/", gsub(mybasenm, pattern = ".rds", replacement="_recomp.rds", fixed=TRUE)))
# 
### output from h4h job ####
# myfn <- list.files("~/nikta/out/out",pattern = ".rds", full.names=TRUE) 
# slices <- list()
# 
# for(fn in myfn){
#   temp <- readRDS(fn)
#   parTable <- do.call(rbind,temp[[3]])
#   #print(head(rownames(parTable)))
#   # print(str(temp[[3]]))
#   n <- cbind("aac_recomputed" = as.numeric(unlist(temp[[1]]))/100, 
#              "ic50_recomputed" = as.numeric(unlist(temp[[2]])), 
#              "HS" = as.numeric(unlist(parTable[,1])),
#              "E_inf" = as.numeric(unlist(parTable[,2])),
#              "EC50" = as.numeric(unlist(parTable[,3]))) 
#   print(head(rownames(n)))
#   rownames(n) <- names(temp[[3]])
#   slices[[fn]] <- n
# }
# 
# profile.sensitivity <- as.data.frame(do.call(rbind, slices))
# saveRDS(profile.sensitivity , "~/nikta/out/profile.sensitivity.PRISM.rds")

profile.sensitivity<- readRDS(paste0(work_dir,"profile.sensitivity.PRISM.rds"))
length(which(is.na(profile.sensitivity$ic50_recomputed))) # 254748
length(which(is.na(profile.sensitivity$aac_recomputed))) #0

# Published sensitivity data:
# secondary-screen-dose-response-curve-parameters.csv
auc.file <- tempfile()
download.file("https://ndownloader.figshare.com/files/20237739", auc.file)
auc.data <- read.csv(auc.file , stringsAsFactors = F)

auc.data$expid <- paste(auc.data$broad_id , auc.data$screen_id , sep="::")
auc.data$expid <- paste(auc.data$expid , auc.data$depmap_id, sep="::")
auc.data$expid <- paste(auc.data$expid , auc.data$name, sep="::")

# Removing the failed experimets:
auc.data <- auc.data[!is.na(auc.data$depmap_id),]
rownames(auc.data) <- auc.data$expid

# Adding published info to sensitivity profile
profile.sensitivity <- merge(profile.sensitivity, auc.data[, c("auc","ec50","ic50")] , by="row.names", all.x=T) #This might take a while to run
rownames(profile.sensitivity) <- profile.sensitivity$Row.names
profile.sensitivity <- profile.sensitivity[,-1]
colnames(profile.sensitivity) <- c("aac_recomputed","ic50_recomputed","HS", "E_inf", "EC50_recomputed","auc_published","ec50_published","ic50_published")

######### Sensitivity-Info #########
info.sensitivity <- data.frame(sub.dose.resp[!duplicated(sub.dose.resp$expid) , c("expid", "broad_id", "screen_id", "compound_plate", "name","depmap_id")]) 

# add "cellid" and tissueid:
#setdiff(info.sensitivity$depmap_id , cell.obj$depmap_id) #character(0)
#setdiff(cell.obj$depmap_id , info.sensitivity$depmap_id) #character(0)
info.sensitivity <- merge(info.sensitivity, cell.obj[, c("cellid" ,"tissueid","depmap_id")] , by="depmap_id", all.x=T)

# add "drugid"
#setdiff(info.sensitivity$broad_id , drug.obj$broad_id)#character(0)
#setdiff(drug.obj$broad_id , info.sensitivity$broad_id)#character(0) 
info.sensitivity <- merge(info.sensitivity, unique(drug.obj[, c("drugid" ,"PRISM.drugid","broad_id")]) , by="broad_id", all.x=T)

info.sensitivity $ Duration <- "5 Days"
rownames(info.sensitivity) <- info.sensitivity$expid

############################################# Drug object-part2 #############################################
#now that drug metadata is added to the drug.obj => its time to concatenate them back so we have unique row.names
# drugid is infact unique.drugid
for( d in unique(drug.obj$drugid)){
  ind = which(drug.obj$drugid == d)
  
  if (length(ind)>1){
    print(ind)
    for(c in colnames(drug.obj)){
      print(c)
      if(c == "drugid"){next}
      id = paste(na.omit(unique(drug.obj[ind, c])), collapse="///")
      print(id)
      drug.obj[ind, c] = id
    }
  }
}

drug.obj <- unique(drug.obj)
rownames(drug.obj) <- drug.obj$drugid
drug.obj[which(duplicated(drug.obj$drugid)), "drugid"]#character(0)
drug.obj[which(duplicated(drug.obj$PRISM.drugid)), "drugid"]#character(0)
drug.obj[which(duplicated(drug.obj$broad_id)), "drugid"]#character(0)
drug.obj[drug.obj == ""] <- NA

########################################### Curation ##########################################################
######### Drug #########
cur.drug <- data.frame(unique.drugid = drug.obj$drugid,
                       PRISM.drugid = drug.obj$PRISM.drugid,
                       row.names = drug.obj$drugid)

######### Cell #########
cur.cell <- data.frame(unique.cellid = cell.obj$cellid,
                       PRISM.cellid = cell.obj$PRISM.cellid,
                       row.names= cell.obj$cellid)

######### Tissue #########
cur.tissue <- data.frame(unique.tissueid = cell.obj$tissueid, 
                         PRISM.tissueid = cell.obj$PRISM.tissueid , 
                         row.names = cell.obj$cellid)

############################################# PSet-curation #############################################
rows<- rownames(raw.sensitivity)

PRISM_PSet <- PharmacoGx::PharmacoSet("PRISM",
                                      
                                      molecularProfiles = list("rna" = emptySE),
                                      cell = cell.obj,
                                      drug = drug.obj,
                                      sensitivityInfo = info.sensitivity[rows , ] ,
                                      sensitivityRaw = raw.sensitivity ,
                                      sensitivityProfiles <- profile.sensitivity [rows , ],
                                      curationDrug = cur.drug,
                                      curationCell = cur.cell,
                                      curationTissue = cur.tissue,
                                      datasetType = "sensitivity",
                                      verify = TRUE)

# Add annotation
PRISM_PSet@annotation$notes <- "This PSet includes drug-dose information from screen II of 'Discovering the anti-cancer potential of non-oncology drugs by systematic viability profiling' paper. Drugs are distinguished by their broad-ids in sensitivity objects. Dose values in the sensitivity objects are reported in micromolar."  
saveRDS(PRISM_PSet, paste0(work_dir, "PRISM.rds"))

