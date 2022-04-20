############################################################################################
############################################################################################
#### Copyright (c) 2017, Broad Institute
#### Redistribution and use in source and binary forms, with or without
#### modification, are permitted provided that the following conditions are
#### met:
####     Redistributions of source code must retain the above copyright
####     notice, this list of conditions and the following disclaimer.
####     Redistributions in binary form must reproduce the above copyright
####     notice, this list of conditions and the following disclaimer in
####     the documentation and/or other materials provided with the
####     distribution.
####     Neither the name of the Broad Institute nor the names of its
####     contributors may be used to endorse or promote products derived
####     from this software without specific prior written permission.
#### THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#### "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#### LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#### A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#### HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#### SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#### LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#### DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#### THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#### (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#### OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
############################################################################################
############################################################################################

#################################################################################################
###### Classifier for TCGA2017 expression subtypes
###### Reference: Kim et al, The Cancer Genome Atlas Expression Subtypes Stratify Response to Checkpoint Inhibition in Advanced Urothelial Cancer 
######            and Identify a Subset of Patients with High Survival Probability, European Urology, volume 75, 961-964 (2019)
#################################################################################################

library(RColorBrewer)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(reshape2)
library(grid)
library(gridExtra)
library(readxl)

get.SSEC.fold.short <- function(expr,expr.fold,expr.fold.up,cohort,marker0,W1) {
        n.sample <- ncol(expr)
        comm <- intersect(rownames(expr),marker0)
        W0.tmp <- W1[match(comm,rownames(W1),nomatch=0),]
        W0.tmp.norm <- t(apply(W0.tmp,1,function(x) x/sum(x)))
        H.tmp <- array(0,dim=c(5,n.sample))
        for (i in 1:n.sample) {
                X0.tmp <- as.matrix(expr.fold.up[match(comm,rownames(expr.fold.up),nomatch=0),i])
                x <- get.single.NMF(X0.tmp,W0.tmp,1.e-07,5)
                H.tmp[,i] <- x[[1]]
        }
        rownames(H.tmp) <- colnames(W0.tmp)
        colnames(H.tmp) <- colnames(expr.fold)
        H.tmp.norm <- apply(H.tmp,2,function(x) x/sum(x))
        g.tmp <- apply(H.tmp.norm,2,function(x) which.max(x))
        rownames(H.tmp.norm) <- c("Luminal","Luminal-infiltrated","Basal-squamous","Neuronal","Luminal-papiilary")
	return(list(g.tmp,H.tmp))
}

get.single.NMF <- function(X1,W0,tol,K) {
	X1[is.na(X1)] <- 0
        res <- NMF.W(X1,W0,tol,K)
        H1 <- res[[2]]
        g1 <- apply(H1,2,function(x) which.max(x))
	return(list(H1,g1))
}

NMF.W <- function(X,W,tol,K) {
        n.run <- 1
        n.iter <- 1000000
        eps <- 1.e-50
        N <- dim(X)[1]
        M <- dim(X)[2]
        meanX <- mean(X,na.rm=T)
        for (j in 1:n.run) {
                H <- matrix(runif(K * M)*meanX,ncol=M)
                X.ap <- W %*% H
                error.EU <- sum((X-X.ap)^2)
                del <- 1
                count <- 1
                while (del >= tol & count < n.iter) {
                        H <- H * (t(W) %*% X) / (t(W)%*%(W%*%H) + eps)
                        X.ap <- W %*% H
                        del <- abs(error.EU-sum((X-X.ap)^2))
                        error.EU <- sum((X-X.ap)^2)
                        if (count %% 100 == 0) cat(count,error.EU,del,'\n')
                        count <- count+1
                }
        }
        return(list(W,H))
}

gene.CSC <- c("CD44","KRT5","RPSA","ALDH1A1")
gene.basal <- c("CD44","CDH3","KRT1","KRT14","KRT16","KRT5","KRT6A","KRT6B","KRT6C")
gene.luminal <- c("CYP2J2","ERBB2","ERBB3","FGFR3","FOXA1","GATA3","GPX2","KRT18","KRT19","KRT20","KRT7","KRT8","PPARG","XBP1","UPK1A","UPK2")
gene.TP53 <- c("ACTG2","CNN1","MYH11","MFAP4","PGM5","FLNC","ACTC1","DES","PCP4")
gene.CLDN.down <- c("CLDN3","CLDN7","CLDN4")
gene.CLDN.up <- c("CDH1","VIM","SNAI2","TWIST1","ZEB1","ZEB2")
gene.NE <- c("CHGA","CHGB","SCG2","ENO2","SYP","NCAM1")
gene.CIS <- c("SH3BP1","CIC","MINK1","GAPVD1","SHOC2","FBXL5","MBD4","PPP2R5C","ERBB2IP","ARL5A","IL13RA1","SDCBP","BIRC2","SPOP","CDK19")
gene.CIS.down <- c("CRTAC1","CTSE","ANXA10","EEF1A2","IVL","PLA2G2A","FABP4","LAMB3","UPK3B","PADI3","PLEC1","ENTPD3","LOC56901","ZNF144","UPK2","TMPRSS4","LU",
                "FLJ20151","BBC3","BG1","SOX15","TNNI2","HOXA1","CYP11B2","INA","ITGB4","HOXB2","SOX15","SIAT4C","TRIM29","LAD1","LOC93408","CLCA4",
                "BMP7","LTBP3","BST2","C44A","KCNQ1","GRB7","FGFR3","MST1R","CYP2J2","LAD1","CA12","MAPRE3")
gene.CIS.up <- c("AKR1B10","FLNA","CXCR4","LYZ","PRG1","FLJ10134","CDH11","HOXA9","PRG1","IGHG3","RARRES1","IGKC","KYNU","COL15A1","SGCE","PDGFC","IGKC","EFEMP1","SPARC",
                "UAP1","DPYSL2","DCN","TUBB","LUM","IGLJ3","RARRES1","HLA-DQB1","CALD1","DF","HLA-DQA1","IGLJ3","MSN","ITM2A","S100A8","CLIC4","MAN1C1","LHFP","KPNA2",
                "TOP2A","COL3A1","CLECSF2","LOC115207","IGKC","NR3C1","KRAS2")
gene.CIS.up <- unique(gene.CIS.up)
gene.CIS.down <- unique(gene.CIS.down)
gene.squamous <- c("DSC1","DSC2","DSC3","DSG1","DSG2","DSG3","S100A7","S100A8")
gene.EMT <- c("ZEB1","ZEB2","VIM","SNAIL","TWIST1","FOXC2","CDH2")
gene.CCP <- c("FOXM1","CDC20","CDKN3","CDC2","KIF11","KIAA0101","NUSAP1","CENPF","ASPM","BUB1B","RRM2","DLGAP5","BIRC5","KIF20A","PLK1","TOP2A","TK1",
                "PBK","ASF1B","C18orf24","RAD54L","PTTG1","CDCA3","MCM10","PRC1","DTL","CEP55","RAD51","CENPM","CDCA8","ORC6L")
get.signature.score.per.sample <- function(expr.fold) {
        sig.CSC <- colMeans(expr.fold[rownames(expr.fold)%in%gene.CSC,],na.rm=T)
        sig.basal <- colMeans(expr.fold[rownames(expr.fold)%in%gene.basal,],na.rm=T)
        sig.luminal <- colMeans(expr.fold[rownames(expr.fold)%in%gene.luminal,],na.rm=T)
        sig.TP53 <- colMeans(expr.fold[rownames(expr.fold)%in%gene.TP53,],na.rm=T)
        sig.CLDN.up <- colMeans(expr.fold[rownames(expr.fold)%in%gene.CLDN.up,],na.rm=T)
        sig.CLDN.dn <- colMeans(expr.fold[rownames(expr.fold)%in%gene.CLDN.down,],na.rm=T)
        sig.NE <- colMeans(expr.fold[rownames(expr.fold)%in%gene.NE,],na.rm=T)
        sig.CIS.up <- colMeans(expr.fold[rownames(expr.fold)%in%gene.CIS.up,],na.rm=T)
        sig.CIS.dn <- colMeans(expr.fold[rownames(expr.fold)%in%gene.CIS.down,],na.rm=T)
        sig.squamous <- colMeans(expr.fold[rownames(expr.fold)%in%gene.squamous,],na.rm=T)
        sig.EMT <- colMeans(expr.fold[rownames(expr.fold)%in%gene.EMT,],na.rm=T)
        sig.CCP <- colMeans(expr.fold[rownames(expr.fold)%in%gene.CCP,],na.rm=T)
        sig.expr1 <- rbind(sig.luminal,sig.basal,sig.TP53,sig.CLDN.up,sig.CLDN.dn,sig.NE,sig.CIS.up,sig.CIS.dn,sig.CSC,sig.squamous,sig.EMT,sig.CCP)
        rownames(sig.expr1) <- c("Luminal","Basal","p53-like","Claudin.up","Claudin.dn","NE","CIS.up","CIS.dn","Bladder CSC","Squamous","EMT","Cell Cycle")
        return(list(sig.expr1))
}

plot.expr.fold.heatmap <- function(mat,scale0,cut.fold,g.Bayes) {
	mat[is.na(mat)] <- 0
        scale <- 1
        color.axis <- 'black'
        .theme_ss <- theme_bw(base_size=12) +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=12*scale, family="mono",face='bold',color=color.axis),
                axis.text.y = element_text(hjust = 0.5,size=8*scale0, family="mono",face='bold',color=color.axis),
                axis.text = element_text(size = 12*scale, family = "mono",color=color.axis),
                axis.title=element_text(face="bold",size=12*scale,color="black"),
                plot.title=element_text(face="bold",size=12*scale))
        #hc <- hclust(dist(mat,method="euclidean"),method="ward.D")
        #gene.ordering <- hc$labels[hc$order]
	gene.ordering <- rownames(mat)
	sample.ordering <- colnames(mat)
        x <- mat
        x[x > +cut.fold] <- +cut.fold
        x[x < -cut.fold] <- -cut.fold
        y <- data.frame(rownames(x),x)
        colnames(y) <- c("gene",colnames(x))
        df <- melt(y,id="gene")
        colnames(df) <- c("gene","sample","activity")
        df$gene <- factor(df$gene,levels=gene.ordering)
        df$sample <- factor(df$sample,levels=sample.ordering)
        p = ggplot(df,aes(x=sample,y=gene,fill=activity))+geom_tile() #geom_tile(colour="yellow")
        p = p + scale_fill_gradient2(low="green",mid="black",high ="red",name=paste("Log2(Fold",'\n',"Changes)",sep=""))
        p = p + .theme_ss
        p = p + ggtitle("Log2(Fold Changes)")
        p = p + xlab("Sample") + ylab("Differentially Expressed Genes")
        p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=12*scale))
        p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=12*scale))
        p = p + theme(legend.position="right")
	n.subtype <- table(g.Bayes)
	p = p + geom_vline(xintercept=sum(g.Bayes==1)+0.5,col='yellow',size=1)
	p = p + geom_vline(xintercept=sum(g.Bayes%in%c(1,2))+0.5,col='yellow',size=1)
	p = p + geom_vline(xintercept=sum(g.Bayes%in%c(1,2,3))+0.5,col='yellow',size=1)
	p = p + geom_vline(xintercept=sum(g.Bayes%in%c(1,2,3,4))+0.5,col='yellow',size=1)
	p = p + geom_hline(yintercept=nrow(mat)-sum(rownames(mat)%in%marker1)+0.5,col='yellow',size=1)
	p = p + geom_hline(yintercept=nrow(mat)-sum(rownames(mat)%in%c(marker1,marker2))+0.5,col='yellow',size=1)
	p = p + geom_hline(yintercept=nrow(mat)-sum(rownames(mat)%in%c(marker1,marker2,marker3))+0.5,col='yellow',size=1)
	p = p + geom_hline(yintercept=nrow(mat)-sum(rownames(mat)%in%c(marker1,marker2,marker3,marker4))+0.5,col='yellow',size=1)
        return(p)
}

############################
############################
CURRENT <- paste(getwd(),"/",sep="")
OUTPUT <- paste(CURRENT,"OUTPUT/",sep="")
system(paste("mkdir",OUTPUT,sep=" "))

############################
### Classifier
############################

#### read TCGA2017 data
load(file=paste("marker.classifier_for_TCGA2017.RData",sep=""))
load(file=paste("W1.classifier_for_TCGA2017.RData",sep=""))
marker0.original <- marker0
W1.original <- W1
W1.marker0 <- data.frame(W1[match(marker0.original,rownames(W1),nomatch=0),])
W1.marker0[,"g.Bayes"] <- apply(W1.marker0,1,function(x) which.max(x))
marker1 <- rownames(W1.marker0)[W1.marker0$g.Bayes==1]
marker2 <- rownames(W1.marker0)[W1.marker0$g.Bayes==2]
marker3 <- rownames(W1.marker0)[W1.marker0$g.Bayes==3]
marker4 <- rownames(W1.marker0)[W1.marker0$g.Bayes==4]
marker5 <- rownames(W1.marker0)[W1.marker0$g.Bayes==5]
marker.TCGA <- c(marker1,marker2,marker3,marker4,marker5)

################ please cutomize this part ####################
##### read expression data file - expression should be log2-transformed
load(file=paste("example.expression.RData",sep=""))
expr[expr==0] <- NA  ##### zero expression is replanced by NA before computing fold changes
expr.fold <- t(apply(expr,1,function(x) x-median(x,na.rm=T)))
expr.fold.up <- apply(expr.fold,2,function(x) ifelse(x>0,x,0))
###############################################################

cohort <- "TEST"
marker0 <- rownames(expr.fold)[rownames(expr.fold)%in%marker.TCGA] 
expr.fold.marker <- expr.fold[match(marker0,rownames(expr.fold),nomatch=0),]
x <- get.SSEC.fold.short(expr,expr.fold,expr.fold.up,cohort,marker.TCGA,W1)
g.Bayes <- x[[1]]
H.Bayes <- x[[2]]
rownames(H.Bayes) <- c("Lum","Lum_Inf","BS","Neuronal","Lum_Pap")
save(g.Bayes,file=paste(OUTPUT,paste(cohort,"g.Bayes.TCGA2017.RData",sep="."),sep=""))
save(H.Bayes,file=paste(OUTPUT,paste(cohort,"H.Bayes.TCGA2017.RData",sep="."),sep=""))

##########################################################
### df contains subtyle label information
### "TCGA2017" column refers to subtype labels
### First five columns in df represents a contribution from five TCGA subtypes and the subtype was determined by the maximum contribution.
### If you do a columnwise normalization of H.Bayes it will be a normalized association of each sample to five TCGA subtypes.
##########################################################
df <- data.frame(t(H.Bayes),g.Bayes)
df[,"sample"] <- rownames(df)
write.table(df,file=paste(OUTPUT,paste(cohort,"df.TCGA2017.txt",sep="."),sep=""), append = F, quote = F, sep = "\t",
       	eol = "\n", na = "NaN", dec = ".", row.names = F,
        col.names = T, qmethod = c("escape", "double"))

#########################
####### Gene expression signature
#########################
df <- read.delim(paste(OUTPUT,paste(cohort,"df.TCGA2017.txt",sep="."),sep=""),header=T,sep='\t',as.is=T)
df[,"TCGA2017"] <- df$g.Bayes
df$TCGA2017 <- gsub("1","Lum",df$TCGA2017)
df$TCGA2017 <- gsub("2","Lum_Inf",df$TCGA2017)
df$TCGA2017 <- gsub("3","BS",df$TCGA2017)
df$TCGA2017 <- gsub("4","Neuronal",df$TCGA2017)
df$TCGA2017 <- gsub("5","Lum_Pap",df$TCGA2017)

load(file=paste(OUTPUT,paste(cohort,"H.Bayes.TCGA2017.RData",sep="."),sep=""))
g.Bayes <- apply(H.Bayes,2,function(x) which.max(x))
g.tmp <- g.Bayes
H.tmp <- H.Bayes
rownames(H.tmp) <- paste(c("Lum","Lum_Inf","BS","Neuronal","Lum_Pap"),sep="")
H.tmp.norm <- apply(H.tmp,2,function(x) x/sum(x))
g1 <- names(g.tmp)[g.tmp==1]
g2 <- names(g.tmp)[g.tmp==2]
g3 <- names(g.tmp)[g.tmp==3]
g4 <- names(g.tmp)[g.tmp==4]
g5 <- names(g.tmp)[g.tmp==5]
order1 <- g1[order(H.tmp.norm[1,g.tmp==1],decreasing=T)]
order2 <- g2[order(H.tmp.norm[2,g.tmp==2],decreasing=T)]
order3 <- g3[order(H.tmp.norm[3,g.tmp==3],decreasing=T)]
order4 <- g4[order(H.tmp.norm[4,g.tmp==4],decreasing=T)]
order5 <- g5[order(H.tmp.norm[5,g.tmp==5],decreasing=T)]
sample0 <- c(order1,order2,order3,order4,order5)
gene0 <- rev(marker.TCGA[marker.TCGA%in%rownames(expr.fold.marker)])
fold.marker0 <- expr.fold[match(gene0,rownames(expr.fold),nomatch=0),match(sample0,colnames(expr.fold),nomatch=0)]
pdf(file=paste(OUTPUT,paste(cohort,"heatmap.fold.TCGA2017.pdf",sep="."),sep=""),width=12,height=0.1*nrow(fold.marker0))
        p <- plot.expr.fold.heatmap(fold.marker0,1.0,3.0,g.Bayes)
        plot(p)
dev.off()

sig.expr <- get.signature.score.per.sample(expr.fold)[[1]]
sig.expr <- sig.expr[,match(sample0,colnames(sig.expr),nomatch=0)]
sig.expr <- data.frame(t(sig.expr))
sig.expr[,"TCGA2017"] <- df$TCGA2017[match(rownames(sig.expr),df$sample,nomatch=0)]
sig.expr[,"sample"] <- rownames(sig.expr)
df.tmp <- melt(sig.expr,id.var=c("sample","TCGA2017"))
df.tmp$variable <- factor(df.tmp$variable,levels=colnames(sig.expr)[1:12])
df.tmp$TCGA2017 <- factor(df.tmp$TCGA2017,levels=c("Lum","Lum_Inf","BS","Neuronal","Lum_Pap"))
scale <- 0.8
alpha <- 0.5
p = ggplot(df.tmp,aes(x=TCGA2017,y=value))
p = p + geom_boxplot(aes(fill=TCGA2017),alpha=alpha)
p = p + geom_point(aes(fill=TCGA2017),size=1,shape=21,position=position_jitterdodge(jitter.width=0.5),alpha=alpha)
p = p + facet_wrap(~variable, scale = "free_y",ncol=6)
p = p + ggtitle(paste("Expression signature scores",sep=""))
p = p + theme(legend.position="none")
p = p + theme(strip.text.x = element_text(size=9,face='bold',angle=0))
p = p + theme(axis.text.x = element_text(angle=45,size=12*scale*1.0,hjust=1.0,vjust=1.0))
p = p + theme(axis.text.y = element_text(angle=0,size=12*scale*1.0))
p = p + xlab("") + ylab("Signature Scores")
pdf(file=paste(OUTPUT,paste(cohort,"boxplot.TCGA2017.pdf",sep="."),sep=""),width=8,height=6)
        plot(p)
dev.off()
write.table(sig.expr,file=paste(OUTPUT,paste(cohort,"sig_score.TCGA2017.txt",sep="."),sep=""), append = F, quote = F, sep = "\t",
       	eol = "\n", na = "NaN", dec = ".", row.names = F,
        col.names = T, qmethod = c("escape", "double"))
