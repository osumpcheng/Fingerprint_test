library(xcms)
library(CAMERA)
library('magrittr')
library(scales)
library("dplyr")
library('stats')
register(SerialParam()) 

file_path <- c('C:/Users/shiche/Box/EcoChem/Student Projects/Lya/forensics_nov_xcms')
file_path1 <- c('C:/Users/shiche/Box/OSU Forensics Data Repository/USGS Data/Year1_mzXML')

list_files <- list.files(file_path, pattern = "*.mzXML", recursive = TRUE, full.names = TRUE)
list_files <- c(list_files, list.files(file_path1, pattern = "*.mzXML", recursive = TRUE, full.names = TRUE))

files_lbl <- list.files(file_path, pattern = "*.mzXML", recursive = TRUE, full.names = FALSE)
files_lbl1 <- list.files(file_path1, pattern = "*.mzXML", recursive = TRUE, full.names = FALSE)

group_lbl <- c(substr(files_lbl, 1,4), substr(files_lbl1, 5,7))

#for (file in files_lbl){
#  group_lbl <- c(group_lbl, substr(file, 1,4))
#}
pd_for <- data.frame(sample_name = c(files_lbl,files_lbl1),
                 sample_group = group_lbl)

raw_for <- readMSData(list_files, pdata = AnnotatedDataFrame(pd_for), mode = "onDisk")

for_colors <- RColorBrewer::brewer.pal(length(unique(group_lbl)), "Set1")

names(for_colors) <- unique(group_lbl)
sampl_colors <- for_colors[raw_for$sample_group]

#plot total ion chromatogram 
for_tc <- chromatogram(raw_for, aggregationFun = 'sum')
plot(for_tc, col = sampl_colors)

for_tcb <- split(tic(raw_for), f = fromFile(raw_for))
boxplot(for_tcb, col = sampl_colors, ylab = "intensity", main = "Total ion current")

#plot raw data aquisition
filterFile(raw_for, file = c(55:63)) %>% 
  filterRt(rt = 8.45*60 + c(-120,120)) %>% 
  filterMz(mz = 416.979844600279 + c(-0.005, 0.005)) %>%
  plot(type = "XIC")

#plot internal standard EIC
mzrt_is <- read.csv('for_is.csv', row.names = 1)
for ( n in rownames(mzrt_is[c(13:16),])){
  mzr_1 <- mzrt_is[n,2] + c(-0.005, 0.005)
  rtr_1 <- mzrt_is[n,1] + c(-60, 60)
  chr_raw <- chromatogram(raw_for, mz = mzr_1, rt = rtr_1)
  png(file= paste(n, "_for_raw_001.png",sep="_"), width = 960, height = 960)
  plot(chr_raw, col=sampl_colors)
  dev.off()
}

cwp_for <- CentWaveParam(peakwidth  = c(5, 40), #Min/Max Peak Width
                         ppm             = 5, #Maximum tolerated m/z deviation
                         noise           = 1000, #Minimum intensity of true peaks
                         snthresh        = 15, #Signal to noise ratio cutoff
                         prefilter       = c(3, 2000),
                         fitgauss        = TRUE,
) #Prefilter for ROI, c(peaks, intensity)
xsn_for <- findChromPeaks(raw_for, param = cwp_for)

boxplot(split((chromPeaks(xsn_for)[,'rtmax'] - chromPeaks(xsn_for)[,'rtmin']), 
              f = chromPeaks(xsn_for)[,'sample']), col = sampl_colors, 
        ylab = 'Peak Width', xlab = 'File', main = 'Distibution of Peak Width')

mpp <- MergeNeighboringPeaksParam(expandRt = 2, ppm = 2)
xsn_for <- refineChromPeaks(xsn_for, mpp)
xsn_for <- refineChromPeaks(xsn_for, param = CleanPeaksParam(maxPeakwidth = 40))

boxplot(split((chromPeaks(xsn_for)[,'rtmax'] - chromPeaks(xsn_for)[,'rtmin']), 
              f = chromPeaks(xsn_for)[,'sample']), col = sampl_colors, 
        ylab = 'Peak Width', xlab = 'File', main = 'Distibution of Peak Width')

pdp_for <- PeakDensityParam(sampleGroups = xsn_for$sample_group,
                            bw      = 30, #Bandwidth of the smoothing kernel
                            binSize   = 0.005, #Binsize to slice up m/z dimension into
                            minFraction = 0.1, #Minimum fraction of samples a group must be present in
                            minSamples = 1, #Minimum number of samples a group must be present in 
                            maxFeatures = 20) #maximum number of features per group

xgroup_for <- groupChromPeaks(xsn_for, param = pdp_for)

peaks_for <- featureDefinitions(xgroup_for)
peaks_stats_for <- featureSummary(xgroup_for)
df_for <- featureValues(xgroup_for, method = "maxint", 
                        value = 'intb', missing = 0)
saveRDS(raw_for, 'for_XCMS_rawdata_031023.rds')
saveRDS(xgroup_for, 'for_XCMS_result_031023.rds')
write.csv(df_for, 'forensics_usgs1_all_031023.csv')
write.csv(peaks_for, 'forensics_usgs1_all_peaks_031023.csv')
