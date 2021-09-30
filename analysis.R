library(data.table)

#Generating list of Essential Gene IDS
commonessentials <- read.csv("~/Achilles_common_essentials.csv")
View(commonessentials)
CE_mod1<-gsub(" \\(.*$","",commonessentials[,"gene"])
CE_GeneID<-gsub("\\).*$","",CE_mod1)
View(CE_mod1)


cn <- fread("C:/users/no/experimentsbyid/cp120/ccle_gene_cn.csv")
View(head(cn))
pickegenes <- read.csv("C:/users/no/experimentsbyid/cp120/tsss/CP120_Promoter_Project-Primer_Table_with_Sequences_Restriction_Sites_and_Final_Primers.csv")
View(pickegenes)
sub1 <- gsub(" \\(.*$", "", colnames(cn))
View(sub1)
cn_working <- cn
colnames(cn_working) <- sub1
colnames(cn) <- sub1
View(head(cn_working))
pickegenes$gene_name
picked <- which(colnames(cn_working) %in% pickegenes$gene_name)
cn_working <- cn_working[,c(1,  5391, 11192,
                            17473, 18126,
                            18933, 19790,
                            20969, 21042,
                            23419, 23781)]
View(cn_working)
Range <- range(cn_working[,2:length(colnames(cn_working))])
FullRange <- range(cn[,2:length(colnames(cn))])
FullRange
plot((cn_working$THAP1), breaks = 100)
CN_Subset_Range <- c(0.70, 1.30)

file <- fread("C:/users/no/experimentsbyid/cp120/Expression_input.csv")
file_working <- file
filecolnames <- gsub(" \\(.*$","", colnames(file_working))
colnames(file_working) <- filecolnames
View(file_working)
View(cn)
expression_lines <- file_working$V1
cn_lines <- cn$V1

#Which cell lines in the expression file do we have copy number data for?
expression_cell_lines_subset.idx <- which(file_working$V1 %in% cn$V1)
cell_lines_keep <- file_working$V1[expression_cell_lines_subset.idx]
View(cell_lines_keep)
#View(expression_subset.idx)
expressionfile_with_cn_data_lines <- file[expression_cell_lines_subset.idx,]
View(expressionfile_with_cn_data_lines)
expressionfile_with_cn_data_lines <- as.data.frame(expressionfile_with_cn_data_lines)

#Which genes in the expression file do we have copy number data for?
expression_genes_subset.idx <- which(colnames(file_working) %in% colnames(cn))
View(expression_genes_subset.idx)
genes_keep <- colnames(file_working)[expression_genes_subset.idx]
View(genes_keep)
expressionfile_with_cn_data_lines_and_Genes <- expressionfile_with_cn_data_lines[,which(colnames(file_working) %in% colnames(cn))]
View(expressionfile_with_cn_data_lines_and_Genes)


#Subsetting and sorting dataframes.
cn_working <- as.data.frame(cn_working)
cn_working_subsetted <- cn_working[which(cn_working$V1 %in% cell_lines_keep),which(colnames(cn_working) %in% genes_keep)]
View(cn_working_subsetted)

file_working <- as.data.frame(file_working)
View(file_working)
file_working_subsetted <- file_working[which(file_working$V1 %in% cell_lines_keep),which(colnames(file_working) %in% genes_keep)]
View(file_working_subsetted)

cn_working_subsetted_sort1 <- cn_working_subsetted[,order(names(cn_working_subsetted))]
View(cn_working_subsetted_sort1)
cn_working_fully_subsetted <- cn_working_subsetted_sort1[order(cn_working_subsetted$V1), ]
View(cn_working_fully_subsetted)
reorder_cn <- cn_working_fully_subsetted[,c(17562, 1:17561, 17563:length(colnames(cn_working_fully_subsetted)))]
View(reorder_cn)

file_working_subsetted_sort1 <- file_working_subsetted[,order(names(file_working_subsetted))]
View(file_working_subsetted_sort1)
file_working_subsetted_sort2 <- file_working_subsetted_sort1[order(file_working_subsetted_sort1$V1),]
View(file_working_subsetted_sort2)
match("V1", colnames(file_working_subsetted_sort2))
reorder_exp <- file_working_subsetted_sort2[,c(17562, 1:17561, 17563:length(colnames(file_working_subsetted_sort2)))]
which(colnames(reorder_cn) == colnames(reorder_exp))

#check for equality
aaa<- data.frame(colnames(reorder_cn), colnames(reorder_exp))
aaa$isequal <- (aaa$colnames.reorder_cn. == aaa$colnames.reorder_exp.)
bbb <- data.frame(reorder_cn$V1, reorder_exp$V1)
bbb$isequal <- (bbb$reorder_cn.V1 == bbb$reorder_exp.V1)
View(bbb)

#Write
cnout <- write.csv(reorder_cn, "C:/users/no/experimentsbyid/cp120/CN_Table_With_Corresponding_Expression_Data.csv", row.names = FALSE)
check <- fread("C:/users/no/experimentsbyid/cp120/CN_Table_With_Corresponding_Expression_Data.csv")
View(check)

expout <- write.csv(reorder_exp, "C:/users/no/experimentsbyid/cp120/Expression_Table_With_Corresponding_CN_Data.csv", row.names = FALSE)
check2 <- fread("C:/users/no/experimentsbyid/cp120/Expression_Table_With_Corresponding_CN_Data.csv")
View(check2)

check <- as.matrix(check)

#CN_range_mask
library(data.table)
isinrange <- between(check[,2:length(colnames(check))], 0.7, 1.30)
length(isinrange)
length(check)
View(isinrange)
check3 <- check2[,2:length(colnames(check2))]
check3[!isinrange] <- ""
check3[,1]<-check$V1
View(head(check3))
check4 <- data.matrix(check3)
View(head(check4))

check4 <- cbind(check2$V1,check4)
View(head(check4))
length(which(isinrange))
length(isinrange)

out <- write.csv(check4, "C:/users/no/experimentsbyid/cp120/Expression_with_CN_between_0.7_and_1.3.csv", row.names = FALSE)
inn <- fread("C:/users/no/experimentsbyid/cp120/Expression_with_CN_between_0.7_and_1.3.csv")
View(inn)

#Binning by CV

inn <- read.csv("~/ExperimentsbyID/CP120/Expression_with_CN_between_0.7_and_1.3.csv")

library(matrixStats)
library(dplyr)
library(ggplot2)
library(Matrix)

top  <- inn

#Removing genes which are not expressed in all cell lines.

top[is.na(top)] <- "CN_err"
View(top)
top[top == 0] <- NA
NoExpression<- grepl("NA", top)
OnlyExpressed<-top[!NoExpression]
View(OnlyExpressed)

Save <- write.csv(OnlyExpressed, "~/ExperimentsbyID/CP120/ExpressionTableWithCN_Errs_Marked_and_Only_Expressed_Genes.csv", row.names = FALSE)
check <- read.csv("~/ExperimentsbyID/CP120/ExpressionTableWithCN_Errs_Marked_and_Only_Expressed_Genes.csv")
View(check)

#Change CN_Errs to NAs
OnlyExpressed[OnlyExpressed == "CN_err"] <- NA
View(OnlyExpressed)

#Getting range for each gene across all cell lines.

Expression_working <- OnlyExpressed[,2:length(colnames(OnlyExpressed))]
View(Expression_working)
Ranges<- colRanges(data.matrix(Expression_working), na.rm = TRUE)
dfrange<-as.data.frame(Ranges)
Range<-dfrange$V2-dfrange$V1
View(Range)

#Getting average expression for each gene across all cell lines.
Mean_Expression<-colMeans(data.matrix(Expression_working), na.rm = TRUE)
df<-as.data.frame(Mean_Expression)
View(df)

#Quick plot and histogram to get a sense of what's going on. 24 breaks seemed to offer the best in terms of ease of viewing and also resolution.
plot(Mean_Expression)
hist(Mean_Expression, breaks=24)

#Breaking down the histogram into 4 bins. My thinking was that there seems to be a high population of genes with an average expression level between 3.5 and 5, then two populations flanking it which were less common, and then finally the extreme outlier population of above ~7. I thought I could start by breaking these 4 bins apart and sorting them by genes with lowest the standard deviation across all cell lines, and choose a few genes from especially the higher groups. I wasn't sure if it would be best to focus on the most highly expressed genes across all cell lines or get a nice range of them, but low STDEV seemed like a safe starting place.
x <- cut(Mean_Expression, c(-Inf, 3.5, 5, 7, 8,9,11,12, Inf), labels = c("<3.5", "3.5-5", "5-7", "7-8",'8-9','9-11','11-12','>12'))
y <- as.data.frame(x)

#Labelling the rows and columns of the new data.frame
colname<- c("Expression Bin")
rownames(y) <- rownames(df)
colnames(y) <- colname

#Making standard deviation function and running it on the genes which are expressed universally
colSd <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=sd, na.rm=na.rm)
Standard_Deviation<-colSd(data.matrix(Expression_working), na.rm = TRUE)

#Joins together the average expression level, standard deviation, range and binning ID of each gene.
joinedi<-cbind(Mean_Expression, Standard_Deviation, Range, y)
GeneNames<-rownames(joinedi)
joinedi<-cbind(GeneNames,joinedi)

View(joinedi)

##GENE ID MAKER
JoinedSub1<-gsub("^.*?\\.{2}","",rownames(joinedi))
EntrezIDs<-gsub("\\.","",JoinedSub1)
joinedj <- cbind(joinedi, EntrezIDs)

#rownames(joinedj)<-EntrezIDs
View(joinedj)

#IS GENE ESSENTIAL?
Is_Essential <- joinedj$EntrezIDs %in% CE_mod1
joined<-cbind(joinedj,Is_Essential)
rownames(joined)<-c()
joined<- merge.data.frame(x=joined,y=bed1.3, 'EntrezIDs')
View(joined)

#Making CV Column
joined$CV <- joined$Standard_Deviation/joined$Mean_Expression
View(joined)

hist(joinedinn$Mean_Expression, breaks = 100)

JoinedOut <- write.csv(joined, "~/ExperimentsbyID/CP120/CP120_Expression_Variance_Table_With_CN_Adjustment.csv", row.names = FALSE)
Firstrun <- read.csv("~/ExperimentsbyID/CP120/AnnotatedCompleteBEDListSortedByCV.csv")
joinedinn <- read.csv("~/ExperimentsbyID/CP120/CP120_Expression_Variance_Table_With_CN_Adjustment.csv")


View(joinedinn)
View(Firstrun)
hgA <- hist(Firstrun$CV, breaks = 100, plot = FALSE)
hgB <- hist(joinedinn$CV, breaks = 1000, plot = FALSE)

plot(hgB, col = "red", ylim = c(0, 1000))
plot(hgA, col = "blue", add = TRUE)
#BINMAKER
stdvhighest <- grepl(">12", joinedinn$Expression.Bin)
Highest <- joinedinn[stdvhighest,]
HighestSort<- Highest[order(Highest$CV),]
View(HighestSort)

stdvsecondhighest<- grepl("11-12", joinedinn$Expression.Bin)
SecondHighest <- joinedinn[stdvsecondhighest,]
SecondHighestSort<- SecondHighest[order(SecondHighest$CV),]
#View(SecondHighestSort)

stdvthirdhighest<- grepl("9-11", joinedinn$Expression.Bin)
ThirdHighest <- joinedinn[stdvthirdhighest,]
ThirdHighestSort<- ThirdHighest[order(ThirdHighest$CV),]
#View(ThirdHighestSort)

stdvfourthhighest<- grepl("8-9", joinedinn$Expression.Bin)
FourthHighest <- joinedinn[stdvfourthhighest,]
FourthHighestSort<- FourthHighest[order(FourthHighest$CV),]
#View(FourthHighestSort)

stdvfifthhighest<- grepl("7-8", joinedinn$Expression.Bin)
FifthHighest <- joinedinn[stdvfifthhighest,]
FifthHighestSort<- FifthHighest[order(FifthHighest$CV),]
#View(FifthHighestSort)

stdvSixthhighest<- grepl("5-7", joinedinn$Expression.Bin)
SixthHighest <- joinedinn[stdvSixthhighest,]
SixthHighestSort<- SixthHighest[order(SixthHighest$CV),]
#View(SixthHighestSort)

stdvSeventhhighest<- grepl("3.5-5", joinedinn$Expression.Bin)
SeventhHighest <- joinedinn[stdvSeventhhighest,]
SeventhHighestSort<- SeventhHighest[order(SeventhHighest$CV),]
#View(SeventhHighestSort)

EigthHighest <- joinedinn[stdvEigthhighest,]
EigthHighestSort<- EigthHighest[order(EigthHighest$CV),]
#View(EigthHighestSort)


Out<- rbind(HighestSort, SecondHighestSort, ThirdHighestSort, FourthHighestSort, FifthHighestSort, SixthHighestSort, SeventhHighestSort, EigthHighestSort)
rownames(Out) <- c()
View(Out)

secondpasswithcn_adjustments <- write.csv(Out, "~/ExperimentsbyID/CP120/CP120_Complete_First_Pass_On_CN_Adjusted_Expression_Table.csv", row.names = FALSE)
