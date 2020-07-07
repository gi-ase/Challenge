##first things first, check the working directory
##make sure that the VCF is pointed correctly to the correct directory with the correct file name is line 12
getwd()

############################################load all the packages needed###############################################################
library(VariantAnnotation)
library(devtools)
library(httr)
library(dplyr)
#######################################################################################################################################
#!! STOP!! Make sure that this part matches where you have the data saved, and if it is named the same!
v_total = readVcf('./Challenge_data.vcf')

#expand the file
v_expanded = expand(v_total, row.names=T)

#turn the file into a dataframe object
df3 = as.data.frame(info(v_expanded))

#Only pull the columns that are needed for the challenge
#v_total <--- only run this if want to see all the metadata, no need to run for the program to work.

# 1. Type of variation (Substitution, Insertion, Silent, Intergenic, etc.) 
# A. TYPE = The type of allele, either snp, mnp, ins, del, or complex.

# 2. Depth of sequence coverage at the site of variation
# A.DP = Total read depth at the locus. 

# 3. Number of reads supporting the variant
# A. AO = Alternate allele observations, with partial observations recorded fractionally. 

# 4. Percentage of reads supporting the variant vs. those supporting the reference
# A. PAIRED = Proportion of observed alternate alleles which are supported by properly paired read fragments.
# A. PAIREDR = Proportion of observed reference alleles which are supported by properly paired read fragments.

df4=subset(df3, select=c(TYPE,DP,AO,PAIRED,PAIREDR))

#Next, want to pull out the specific names for the chromosome variants to format it in the way ExAC likes
df4$CHROM = as.vector(seqnames(v_expanded)) #chromosome
c = fixed(v_expanded)[c("REF","ALT")]
d = data.frame(c)
df4$REF = d$REF #reference
df4$ALT = d$ALT #variant alt
df4$POS = start(v_expanded) #position

#### To finish question 1, need to check the data for duplicates, to see if there is a need for a ranking system.
dupcheck = any(duplicated(subset(df4, select = c(CHROM, POS, REF))))
#dupcheck comes back true, so there is potential for multiple possbilities for type of variant.

df4$VarID = paste0(df4$CHROM,"-",df4$POS,"-",df4$REF)
df4$TypeRanked = df4$TYPE
df4$TypeRanked = case_when(df4$TypeRanked == "snp" ~ 0, df4$TypeRanked == "mnp" ~ 1, df4$TypeRanked == "complex" ~ 2, df4$TypeRanked == "ins" ~ 3, df4$TypeRanked == "del" ~ 4)
df5 = df4 %>% group_by(VarID) %>% top_n(1,TypeRanked)
View(df5)


#build the variant name in the format that can use to search the ExAC database
VarNam = paste0(df5$CHROM,"-",df5$POS,"-",df5$REF,"-",df5$ALT)
df5$VariantName = VarNam


#Now for hte ExAC pull!! in bulk!! formatted to json 
ExACBulk = httr::POST(url="http://exac.hms.harvard.edu/rest/bulk/variant", body=jsonlite::toJSON(as.character(df5$VariantName)), encode = "json")
jcont = content(ExACBulk)

# pull that data from the ExAC bulk pull and merge it back into our data table
exacData = vector()
for (n in 1:nrow(df5))
{
  ExACinfo = jcont[[as.character(df5$VariantName[n])]]
  VariantName = as.character(df5$VariantName[n])
  alleleFreq = if (!(is.null(ExACinfo$variant$allele_freq))) paste(unlist(ExACinfo$variant$allele_freq),collapse = ",") else "NA"
  exacData = rbind(exacData, cbind(VariantName, alleleFreq))
}
colnames(exacData) = c("VariantName", "Allele_Frequency")



#AT THE FINAL COUNTDOWN *kazoo solo*
FinalTable = merge(df5, exacData, by="VariantName")


#Final table all cleaned up!
Final = subset(FinalTable, select = -c(VarID,TypeRanked) )
View(Final)
