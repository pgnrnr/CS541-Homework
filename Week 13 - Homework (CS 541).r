install.packages('versions')
install.packages('ggplot2')
library(versions)
#available.versions("sommer")
install.versions("sommer", "3.1")
library(sommer)
library(ggplot2)

geno = read.csv("mdp_num_geno.csv")
geno[1:10,1:5] #show a slice of the data frame

MarkerInfo = read.csv("MarkerInfo.csv")
head(MarkerInfo)

pheno = read.csv("mdp_pheno.csv")
head(pheno)

gmat = as.matrix(geno[,-1])
rownames(gmat) = geno[,1]
gmat[is.na(gmat)] = 0

#find some marker columns that are all zeroes and drop them
dropcols = which(colMeans(gmat) == 0)
gmat = gmat[,-dropcols]
str(gmat)

pcs = prcomp(gmat, center = T, scale = T, retx = T)
str(pcs)

PC_vars = pcs$sdev^2
total_var = sum(PC_vars)
percent_var = PC_vars/total_var*100
percent_var[1:5]

PCs_lines = data.frame(pcs$x[,1:2])
PCs_lines$Line = geno$Line
PCs_lines$group = ifelse(substr(PCs_lines$Line,1,3) == "CML", "CML",
                         ifelse(PCs_lines$Line %in% c("B37", "B73", "B104","A632"), "Stiff Stalk",
                                ifelse(PCs_lines$Line %in% c("MO17", "OH43", "B97", "OH7B"), "Non-Stiff Stalk", "NA")))

pc_plot = ggplot(PCs_lines, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = group))

pc_plot

PCs_pheno = merge(PCs_lines, pheno, by = "Line")
PCs_pheno$PC3 = pcs$x[,3]
head(PCs_pheno)

EHT_pc1 = lm(formula = EarHT ~ PC1 + PC2, data = PCs_pheno)
summary(EHT_pc1)

selector = c('zagl1.2', 'zagl1.6','PZA00236.7','PZA01228.2',
'PZA03351.1')
geno_sub = gmat[, selector]
geno_sub[geno_sub == 1] = 0
head(geno_sub)

Line = rownames(geno_sub)
geno_sub = cbind(Line,geno_sub)
rownames(geno_sub) = 1:nrow(geno_sub)
All_pheno = merge(PCs_pheno, geno_sub, by="Line")
head(All_pheno)

# Add code here

# Add code here

K = A.mat(X = gmat)
rownames(K) = rownames(gmat)
print(K[1:5,1:5])
heatmap(K, symm = T)

rand_mod = mmer2(fixed = dpoll ~ 1, random = ~g(Line), rcov = ~units, G = list(Line = K), data = All_pheno, silent = T)
summary(rand_mod)

Vg = rand_mod$var.comp$`g(Line)`
Verr = rand_mod$var.comp$units
Vp = Vg + Verr
h2 = Vg/Vp
print(paste("Vg:", Vg))
print(paste("Verr:", Verr))
print(paste("Vp:", Vp))
print(paste("Heritability:", h2))

mlm_ex = mmer2(fixed = dpoll ~ PC1, random = ~g(Line), rcov = ~units, G = list(Line = K), data = All_pheno, silent = T)
summary(mlm_ex)

# Add code here

ordered_lines = PCs_lines[order(PCs_lines$PC1), "Line"]
last = length(ordered_lines)
first = ceiling(0.8*last)
print(paste("First and last indices of lines in test set:", first, ",", last))
test_set = ordered_lines[first:last]
head(test_set)

pheno$dpoll_train = pheno$dpoll
pheno[pheno$Line %in% test_set, "dpoll_train"] = NA

gs_mod = mmer2(fixed = dpoll_train ~ 1, random = ~g(Line), rcov = ~units, G = list(Line = K), data = pheno, silent = T)
summary(gs_mod)

pheno$dpoll_pred = gs_mod$fitted.y

write.csv(pheno, ile = 'predictions.csv', row.names = FALSE)
