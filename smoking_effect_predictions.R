##################### Instala��o de packages ##################################

#instala��o do Bioconductor e package GEOquery
source("http://bioconductor.org/biocLite.R")
biocLite() 

#livraria para carregar dataset
biocLite("GEOquery")
library(GEOquery)

#livraria de anota��o
biocLite("hgu133a")
library(hgu133a.db)

#livraria para efectuar filtragem
biocLite("genefilter")
library(genefilter) 

#livraria para modelos lineares(eBayes)
library(limma)

#livraria MLInterfaces para an�lise preditiva
biocLite("MLInterfaces")
library(MLInterfaces)

#enrichment analysis
biocLite('GOstats')
library(GOstats)

#instala��o do dataset
gds3257 = getGEO('GDS3257')
gds3257


#anota��o
#biocLite("annotationTools")


##################### Pr�-Processamento de dados #################################
#verificar o dataset e vari�veis do expressionset
expset = GDS2eSet(gds3257) #passar a expressionset
expset
dim(expset) #dimens�es

exp = exprs(expset) #dados de express�o
exp
exp[1:5,1:3] #amostra do dataset (probes nas linhas e amostras nas colunas)

sampleNames(expset) #nomes das amostras
featureNames(expset) #nomes das features

Meta(gds3257) #descri��o geral do dataset
Meta(gds3257)$description #descri��o do dataset
Meta(gds3257)$sample_organism #Organismo em estudo
Meta(gds3257)$sample_type #tipo de amostra em estudo
experimentData(expset)#informa��o sobre a experi�ncia

colnames(pData(expset)) #nomes das colunas/variaveis
levels(pData(expset)$tissue); levels(pData(expset)$individual); levels(pData(expset)$gender)
#valores que as variaveis podem tomar

#anota��o
expset@annotation='hgu133a'
annotation(expset)

#verificar dados omissos
sum(which(is.na(exp))) #nao possui dados omissos



#Pr�-processamento propriamente dito

#selecionar genes com desvio padr�o dos valores igual ao maior 
#do que tr�s vezes a mediana dos desvios padr�o 
biocLite("genefilter")
library(genefilter) 

sds = rowSds(exp)
sds
med = median(sds)
med
expset.f = expset[sds >= 3*med]
expset.f
dim(exprs(expset.f))#dimensao do dataset filtrado

sum(which(is.na(exprs(expset.f)))) #confirmar que nao tem NAs

#histograma do desvio padrao dos genes
hist(sds, breaks = 50, col = "red", xlab = 'Desvio padr�o', ylab = 'Frequ�ncia', main = NULL)
abline(v = med*3, col = 'green', lwd = 2)

#filtra genes cujo r�cio de maximo valor sobre o minimo valor de express�o seja 2
maximo = apply(exp, 1, max)
minimo = apply(exp, 1, min)
filtro = maximo/minimo >2  
GEO = expset[filtro,]
GEO
dim(GEO)

exp2 = exprs(GEO)
exprs(GEO) = scale(exp2) #normaliza��o dos dados
exprs(GEO)


################### Express�o diferencial ##################################
biocLite("genefilter")
library(genefilter)

#testes estatisticos
#Hip�tese nula: As m�dias dos n�veis de transcri��o das duas situa��es s�o id�nticas

#row t tests para tissue
levels(GEO$tissue) #pode ser tumor ou normal
#H0: a express�o em tecido normal e tumoral � id�ntica
tt = rowttests(GEO, "tissue") #Realizar os t-tests e verificar os p-values
names(tt)
tt$p.value
rank = order(tt$p.value)
p20 = rank[1:20] #20 genes com menor p-values (maior evid�ncia de express�o diferencial
tt$p.value[p20]
g=featureNames(GEO[p20])
annotation(expset)
biocLite("hgu133a")
library(hgu133a.db)
install.packages(org.Hs.eg.db)
unlist(mget(g, hgu133aSYMBOL)) #recuperar os nomes dos genes que correspondem �s probes identificadas
unlist(mget(g, hgu133aGENENAME))


#row t tests para individual 
levels(GEO$individual) #tissue pode ser tumor ou normal
#H0: a express�o de tecidos tumorais e normais � id�ntica
tt = rowttests(GEO, "individual") #Realizar os t-tests e verificar os p-values
names(tt)
tt$p.value
rank = order(tt$p.value)
p20 = rank[1:20] #20 genes com menor p-values (maior evid�ncia de express�o diferencial
tt$p.value[p20]
g=featureNames(GEO[p20])

annotation(expset)
biocLite("hgu133a")
library(hgu133a.db)
install.packages(org.Hs.eg.db)
unlist(mget(g, hgu133aSYMBOL)) #recuperar os simbolos dos genes que correspondem �s probes identificadas
unlist(mget(g, hgu133aGENENAME)) #recuperar os nomes dos genes


#GENERAL LINEAR MODELS
library(limma)#usado ebayes do package limma

#para o tecido
design.tissue = model.matrix(~GEO$tissue)
fit = lmFit(GEO, design.tissue) #especifica a hipotese que queremos testar
fit2 = eBayes(fit)
diff.tissue = topTable(fit2, coef = 2, adjust.method = "BH", 10) #verificar os 10 genes mais diferencialmente expressos
diff.tissue
gen_tis=unlist(mget(rownames(diff.tissue), hgu133aSYMBOL)) #obter nomes dos genes
gen_tis
verificar_na=is.na(gen_tis) #verificar apenas os que tem anota��o
gen_annot_tis=gen_tis[!is.na(gen_tis)]

gen_annot_tis

#verificar nome dos genes diferencialmente expressos
unlist(mget(rownames(diff.tissue), hgu133aGENENAME)) 

#para o individual (condi��o quanto a fumador, se �,era ou n�o �)
design.individual = model.matrix(~GEO$individual)
fit.individual = lmFit(GEO,design.individual)
fit2.individual = eBayes(fit.individual)

diff.individual = topTable(fit2.individual, coef = 2, adjust.method = "BH", 10)

diff.individual
gen_ind=unlist(mget(rownames(diff.individual), hgu133aSYMBOL)) #obter nomes dos genes

verificar_na=is.na(gen_ind) #verificar apenas os que tem anota��o
gen_annot_ind=gen_ind[!is.na(gen_ind)]
gen_annot_ind

#nomes dos genes diferencialmente expressos
unlist(mget(rownames(diff.individual), hgu133aGENENAME)) 




######################    CLUSTERING     ##########################
#para o tecido

#d� a express�o dos genes com express�o diferencial para o tecido

exp.tissue = as.matrix(exprs(GEO[diff.tissue$ID])) 
exp.tissue
matriz.distancias.tissue = dist(exp.tissue)  #cria a matriz de distancias

#clustering hierarquico segundo o metodo average
clustering.hierarquico.tissue = hclust(matriz.distancias.tissue, method = "average")  

plot(clustering.hierarquico.tissue) #grafico para visualizar o dendograma gerado pelo clustering
rc.tissue = rainbow(nrow(exp.tissue), start = 0, end= .3)
cc.tissue = rainbow(ncol(exp.tissue), start = 0, end= .3)
hv.tissue = heatmap(exp.tissue, RowSideColors = rc.tissue, ColSideColors = cc.tissue) #perfil de express�o para o tecido

kmeans.tissue = kmeans(exp.tissue, 3)  #kmeans para o tecido
names(kmeans.tissue)
kmeans.tissue$cluster


#para o individual

exp.individual = as.matrix(exprs(GEO[diff.individual$ID]))

matriz.dist.individual = dist(exp.individual)
clustering.individual = hclust(matriz.dist.individual, method = "average")

plot(clustering.individual)
rc.individual = rainbow(nrow(exp.individual), start = 0, end= .3)
cc.individual = rainbow(ncol(exp.individual), start = 0, end = .3)


hv.individual = heatmap(exp.individual, RowSideColors = rc.individual, ColSideColors = cc.individual) #perfil de express�o para o g�nero

kmeans.individual = kmeans(exp.individual, 3) #kmeans para o individual
names(kmeans.individual)
kmeans.individual$cluster


#######################CONSTRUCAO DE MODELOS DE PREDICAO############################################
##################### APENAS USANDO GENES DIFERENCIALMENTE EXPRESSOS ######################
#os dados que ser�o usado para os dados teste e treino foram coletados a partir dos dados gerados pelo eBayes
biocLite("MLInterfaces")
library(MLInterfaces)

#################Tissue
index.tissue = which(rownames(GEO) %in% rownames(diff.tissue))
GEO.machine.tissue = GEO[index.tissue,]

######### knn
knn.result.tissue = MLearn(tissue~., GEO.machine.tissue, knnI(k=1), 1:40)
confuMat(knn.result.tissue)   #matriz de confus�o
acc(confuMat(knn.result.tissue))  #acuracy
precision(confuMat(knn.result.tissue))  #precis�o
recall(confuMat(knn.result.tissue))  #sensibilidade
F1(confuMat(knn.result.tissue))  #especifidade

########## redes neuronais
#com leave-on-out
nnet.tissue = MLearn(tissue~., GEO.machine.tissue, nnetI, xvalSpec("LOO"), size= 3, decay =0.01)  #leave-on-out
confuMat(nnet.tissue)   #matriz de confus�o
acc(confuMat(nnet.tissue))    #acuracy
precision(confuMat(nnet.tissue))   #precis�o
recall(confuMat(nnet.tissue))  #sensibilidade
F1(confuMat(nnet.tissue))   #especificidade


#com cross-validation
nnet.tissue = MLearn(tissue~., GEO.machine.tissue, nnetI, xvalSpec("LOG", 5, balKfold.xvspec(5)),size= 3, decay =0.01)  #cross-validation
confuMat(nnet.tissue)   #matriz de confus�o
acc(confuMat(nnet.tissue))    #acuracy
precision(confuMat(nnet.tissue))   #precis�o
recall(confuMat(nnet.tissue))  #sensibilidade
F1(confuMat(nnet.tissue))   #especificidade

########### Arvores de decis�o
#cross-validation
arv.tissue  = MLearn(tissue~., GEO.machine.tissue, rpartI, xvalSpec("LOG", 5, balKfold.xvspec(5)))  #cross validation
confuMat(arv.tissue)    #matriz de confus�o
acc(confuMat(arv.tissue))   #acuracy
precision(confuMat(arv.tissue))  #precis�o
recall(confuMat(arv.tissue))   #sensibilidade
F1(confuMat(arv.tissue))   #especificidade

#leave-on.out
arv.tissue  = MLearn(tissue~., GEO.machine.tissue, rpartI, xvalSpec("LOO", 5, balKfold.xvspec(5)))  #leave-on-out
confuMat(arv.tissue)    #matriz de confus�o
acc(confuMat(arv.tissue))   #acuracy
precision(confuMat(arv.tissue))  #precis�o
recall(confuMat(arv.tissue))   #sensibilidade
F1(confuMat(arv.tissue))   #especificidade


#################Individual

index.individual = which(rownames(GEO) %in% rownames(diff.individual)) 
GEO.machine.individual = GEO[index.individual,]
GEO.machine.individual

####knn
knn.result.individual = MLearn(individual~., GEO.machine.individual, knnI(k=1), 1:40)
confuMat(knn.result.individual)   #matriz de confus�o
acc(confuMat(knn.result.individual)) #acuracy
precision(confuMat(knn.result.individual))  #precis�o
recall(confuMat(knn.result.individual))  #sensibilidade
F1(confuMat(knn.result.individual))  #especificiadade


#####redes neuronais 

#com leave-on-out
nnet.individual = MLearn(individual~., GEO.machine.individual, nnetI, xvalSpec("LOO"), size = 3, decay = 0.01)  #leave-on-out
confuMat(nnet.individual)  #matriz de confus�o
acc(confuMat(nnet.individual))  #acuracy
precision(confuMat(nnet.individual))  #precis�o
recall(confuMat(nnet.individual))    #sensibilidade
F1(confuMat(nnet.individual))    #especificidade

# cross-validation
nnet.individual = MLearn(individual~., GEO.machine.individual, nnetI, xvalSpec("LOG", 5, balKfold.xvspec(5)), size = 3, decay = 0.01)  #cross-validation
confuMat(nnet.individual)  #matriz de confus�o
acc(confuMat(nnet.individual))  #acuracy
precision(confuMat(nnet.individual))  #precis�o
recall(confuMat(nnet.individual))    #sensibilidade
F1(confuMat(nnet.individual))    #especificidade

###### Arvores de decis�o 
# leave-on-out
arv.individual = MLearn(individual~., GEO.machine.individual, rpartI, xvalSpec("LOO", 5, balKfold.xvspec(5)))  #leave-on-out
confuMat(arv.individual) #matriz de confus�o
acc(confuMat(arv.individual))  #acuracy
precision(confuMat(arv.individual))  #precis�o
recall(confuMat(arv.individual))   #sensibilidade
F1(confuMat(arv.individual))   #especificidade

# cross validation
arv.individual = MLearn(individual~., GEO.machine.individual, rpartI, xvalSpec("LOG", 5, balKfold.xvspec(5)))  #cross validation
confuMat(arv.individual) #matriz de confus�o
acc(confuMat(arv.individual))  #acuracy
precision(confuMat(arv.individual))  #precis�o
recall(confuMat(arv.individual))   #sensibilidade
F1(confuMat(arv.individual))   #especificidade





#######################CONSTRUCAO DE MODELOS DE PREDICAO############################################
########################### USANDO TODOS OS GENES ##########################
#os dados que ser�o usado para os dados teste e treino foram coletados a partir dos dados gerados pelo eBayes
biocLite("MLInterfaces")
library(MLInterfaces)

#################Tissue



########## redes neuronais
#com leave-on-out
nnet.tissue = MLearn(tissue~., GEO, nnetI, xvalSpec("LOO"), size= 3, decay =0.01)  #leave-on-out
confuMat(nnet.tissue)   #matriz de confus�o
acc(confuMat(nnet.tissue))    #acuracy
precision(confuMat(nnet.tissue))   #precis�o
recall(confuMat(nnet.tissue))  #sensibilidade
F1(confuMat(nnet.tissue))   #especificidade


#com cross-validation
nnet.tissue = MLearn(tissue~., GEO, nnetI, xvalSpec("LOG", 5, balKfold.xvspec(5)),size= 3, decay =0.01)  #cross-validation
confuMat(nnet.tissue)   #matriz de confus�o
acc(confuMat(nnet.tissue))    #acuracy
precision(confuMat(nnet.tissue))   #precis�o
recall(confuMat(nnet.tissue))  #sensibilidade
F1(confuMat(nnet.tissue))   #especificidade

########### Arvores de decis�o
#cross-validation
arv.tissue  = MLearn(tissue~., GEO, rpartI, xvalSpec("LOG", 5, balKfold.xvspec(5)))  #cross validation
confuMat(arv.tissue)    #matriz de confus�o
acc(confuMat(arv.tissue))   #acuracy
precision(confuMat(arv.tissue))  #precis�o
recall(confuMat(arv.tissue))   #sensibilidade
F1(confuMat(arv.tissue))   #especificidade

#leave-on.out
arv.tissue  = MLearn(tissue~., GEO, rpartI, xvalSpec("LOO", 5, balKfold.xvspec(5)))  #leave-on-out
confuMat(arv.tissue)    #matriz de confus�o
acc(confuMat(arv.tissue))   #acuracy
precision(confuMat(arv.tissue))  #precis�o
recall(confuMat(arv.tissue))   #sensibilidade
F1(confuMat(arv.tissue))   #especificidade


#################Individual

index.individual = which(rownames(GEO) %in% rownames(diff.individual)) 
GEO.machine.individual = GEO[index.individual,]
GEO.machine.individual

####knn
knn.result.individual = MLearn(individual~., GEO.machine.individual, knnI(k=1), 1:40)
confuMat(knn.result.individual)   #matriz de confus�o
acc(confuMat(knn.result.individual)) #acuracy
precision(confuMat(knn.result.individual))  #precis�o
recall(confuMat(knn.result.individual))  #sensibilidade
F1(confuMat(knn.result.individual))  #especificiadade


#####redes neuronais 

#com leave-on-out
nnet.individual = MLearn(individual~., GEO, nnetI, xvalSpec("LOO"), size = 3, decay = 0.01)  #leave-on-out
confuMat(nnet.individual)  #matriz de confus�o
acc(confuMat(nnet.individual))  #acuracy
precision(confuMat(nnet.individual))  #precis�o
recall(confuMat(nnet.individual))    #sensibilidade
F1(confuMat(nnet.individual))    #especificidade

# cross-validation
nnet.individual = MLearn(individual~., GEO, nnetI, xvalSpec("LOG", 5, balKfold.xvspec(5)), size = 3, decay = 0.01)  #cross-validation
confuMat(nnet.individual)  #matriz de confus�o
acc(confuMat(nnet.individual))  #acuracy
precision(confuMat(nnet.individual))  #precis�o
recall(confuMat(nnet.individual))    #sensibilidade
F1(confuMat(nnet.individual))    #especificidade

###### Arvores de decis�o 
# leave-on-out
arv.individual = MLearn(individual~., GEO, rpartI, xvalSpec("LOO", 5, balKfold.xvspec(5)))  #leave-on-out
confuMat(arv.individual) #matriz de confus�o
acc(confuMat(arv.individual))  #acuracy
precision(confuMat(arv.individual))  #precis�o
recall(confuMat(arv.individual))   #sensibilidade
F1(confuMat(arv.individual))   #especificidade

# cross validation
arv.individual = MLearn(individual~., GEO, rpartI, xvalSpec("LOG", 5, balKfold.xvspec(5)))  #cross validation
confuMat(arv.individual) #matriz de confus�o
acc(confuMat(arv.individual))  #acuracy
precision(confuMat(arv.individual))  #precis�o
recall(confuMat(arv.individual))   #sensibilidade
F1(confuMat(arv.individual))   #especificidade



########################    ENRICHMENT ANALYSIS    ###########################################

source("http://bioconductor.org/biocLite.R")
biocLite('GOstats')
library(GOstats)


#para o tissue
#overexpression
filt = nsFilter(GEO, require.entrez=T,
                  remove.dupEntrez=T, var.func=IQR, var.cutoff=0.5) #remover genes sem entrada e excluir genes ocm pouca variabilidade
ALLf = filt$eset
affyUniverse = featureNames(ALLf)
entrezUniverse = unlist(mget(affyUniverse, hgu133aENTREZID))#Recolhe conjunto de todos os genes para a anota��o dos dados
ttests = rowttests(ALLf, 'tissue') #para tissue
smPV = ttests$p.value < 0.05 #Testes t para determinar genes diferencialmente expressos e seus IDs com a=0.05
pvalFiltered = ALLf[smPV, ]
selectedEntrezIds =  unlist(mget(featureNames(pvalFiltered),
                                 hgu133aENTREZID))
#Criar par�metros e correr os testes estat�sticos hipergeom�tricos: grupos de genes do Gene
#ontology, genes sobre expressos
params = new("GOHyperGParams",
             geneIds=selectedEntrezIds, universeGeneIds=entrezUniverse,
             annotation="hgu133a.db", ontology="MF", pvalueCutoff=
               0.1, testDirection="over")
hgOver = hyperGTest(params)

#Lista de grupos com menores p- values (ordem crescente; d� contagens de genes no grupo
#alvo e total de genes no grupo)
hgOver
summary(hgOver) 

#para o tissue
#underexpression
filt = nsFilter(GEO, require.entrez=T,
                remove.dupEntrez=T, var.func=IQR, var.cutoff=0.5) #remover genes sem entrada e excluir genes ocm pouca variabilidade
ALLf = filt$eset
affyUniverse = featureNames(ALLf)
entrezUniverse = unlist(mget(affyUniverse, hgu133aENTREZID))#Recolhe conjunto de todos os genes para a anota��o dos dados
ttests = rowttests(ALLf, 'individual') #para tissue
smPV = ttests$p.value < 0.05 #Testes t para determinar genes diferencialmente expressos e seus IDs com a=0.05
pvalFiltered = ALLf[smPV, ]
selectedEntrezIds =  unlist(mget(featureNames(pvalFiltered),
                                 hgu133aENTREZID))
#Criar par�metros e correr os testes estat�sticos hipergeom�tricos: grupos de genes do Gene
#ontology, genes sobre expressos
params = new("GOHyperGParams",
             geneIds=selectedEntrezIds, universeGeneIds=entrezUniverse,
             annotation="hgu133a.db", ontology="MF", pvalueCutoff=
               0.05, testDirection="under")
hgUnder = hyperGTest(params)

#Lista de grupos com menores p- values (ordem crescente; d� contagens de genes no grupo
#alvo e total de genes no grupo)
hgUnder
summary(hgUnder) 












