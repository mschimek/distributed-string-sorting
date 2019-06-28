Rscript singleTableMemory.r JSON/jsonTable_DToNWeak DToNWeak
Rscript singleTableMemory.r JSON/jsonTable_DToNSkewedWeak DToNSkewedWeak
#speedup
Rscript doubleTableSpeedup.r JSON/jsonTable_DToNWeak_speedup DToNWeak_speedup



Rscript singleTableMemory.r JSON/jsonTable_DToN DToN 200000000
Rscript singleTableMemory.r JSON/jsonTable_DToNSkewed DToNSkewed 200000000
#speedup
Rscript doubleTableSpeedup.r JSON/jsonTable_DToN_speedup DToN_speedup



Rscript multipleTablesCommunication.r JSON/jsonTable_Files Files
#speedup
Rscript multipleTablesSpeedup.r JSON/jsonTable_Files Files_speedup



#Rscript singleTableMemory.r JSON/jsonTable_CommonCrawl CommonCrawl 2074964317
#Rscript singleTableMemory.r JSON/jsonTable_CommonCrawlReduced CommonCrawlReduced 1790850388
#
#Rscript singleTableMemory.r JSON/jsonTable_Wiki Wikipedia 1092633438
#Rscript singleTableMemory.r JSON/jsonTable_WikiReduced WikipediaReduced 875259655
#
#Rscript singleTableMemory.r JSON/jsonTable_Suffix100GB Suffix100GB 447723
