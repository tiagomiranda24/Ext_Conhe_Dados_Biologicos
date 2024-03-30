---
title: "Projeto"
author: "Christian Neitzel, Diana Silva, Tiago Miranda"
date: "r format(Sys.time(), '%Y-%m-%d')"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

Objetivo do trabalho:

Este trabalho contém uma análise demonstrativa de expressão genética diferencial, utilizando amostras de sequenciamento direcionado de 526 tumores de bexiga e seus normais correspondentes através do MSK-IMPACT. O conjunto de dados selecionado inclui 526 amostras de mais de 400 pacientes, nas quais a expressão genética foi quantificada.

Conjunto de Dados selecionados:
https://portal.gdc.cancer.gov/projects/TCGA-ACC 
https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000178.v11.p8&phv=354519&phd=&pha=&pht=7516&phvf=&phdf=&phaf=&phtf=&dssp=1&consent=&temp=1

Em primeiro lugar, procedemos à instalação dos packages necessários para realizar a análise da expressão diferencial.


```{r}
# Instalação e carregamento dos pacotes necessários
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install(c("edgeR", "limma", "Glimma", "gplots", "org.Mm.eg.db", "RColorBrewer", "TCGAbiolinks", "DESeq2"))

# Carrega os pacotes instalados
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)
library(TCGAbiolinks)
library(DESeq2)

# Verifica se o BiocManager está instalado e instala se necessário
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Verifica se o pacote cBioPortalData está instalado e instala se necessário
if (!requireNamespace("cBioPortalData", quietly = TRUE))
  BiocManager::install("cBioPortalData")

# Carregamento do pacote cBioPortalData
library(cBioPortalData)

# Inicialização da API do cBioPortal
cbio <- cBioPortal()
```

## Including Plots

Nesta fase é realizado o carregamento dos dados clinicos do estudo,  explorados os ensaios do estudo disponiveis e organizados em dataframes:

```{r pressure, echo=FALSE}
# Obter dados clínicos
bladder_clin <- clinicalData(api = cbio, studyId = "bladder_msk_2023")
bladder_clin
dim(bladder_clin)  # Deve retornar 526 linhas e 29 colunas

# Carrega os dados do cBioPortal
bladder = cBioDataPack("bladder_msk_2023", ask = FALSE)
bladder

# Lista os nomes dos ensaios disponíveis no pacote dos dados do CBioPortal
names(assays(bladder))
# names(rowData(bladder)) # tipos de metadados associados a cada gene
# names(colData(bladder)) # tipos de metadados associados a cada amostra

# Summary dos dados
summary(bladder)

# Acessando os dados de cna e convertendo para um dataframe
dados_cna = assays(bladder)$cna
dados_cna_df = as.data.frame(dados_cna)

# Acessando os dados de cna_hg19 e convertendo para um dataframe
dados_cna_hg19 = assays(bladder)$cna_hg19.seg
dados_cna_hg19_df = as.data.frame(dados_cna_hg19)

# Acessando os dados de mutação e convertendo para um dataframe
dados_mutacao = assays(bladder)$mutations
dados_mutacao_df = as.data.frame(dados_mutacao)

```
```{r}

```
Análise Estatística de Dados de Cancro da Bexiga
1. Definição das Variáveis com as.factor'''

```{r}

# Idade em que a Sequenciação foi Reportada (Anos)
age = as.factor(bladder_clin$AGE_AT_SEQ_REPORTED_YEARS)

# Tipo de Cancro
cancer_type = as.factor(bladder_clin$CANCER_TYPE_DETAILED)

# Estado da Doença
disease_state = as.factor(bladder_clin$DISEASE_STATE)

# Tratamento com Erdafitinib
erdafitinib_treatment = as.factor(bladder_clin$ERDAFITINIB_TREATED)

# Estado de Sobrevivência Global desde o Tx com Erdafitinib
overall_status = as.factor(bladder_clin$ERDAFITINIB_TX_OS_STATUS)

# Estado Sem Progressão desde o Tx com Erdafitinib
progression_status = as.factor(bladder_clin$ERDAFITINIB_TX_PFS_STATUS)

# Painel de Genes
gene_panel = as.factor(bladder_clin$GENE_PANEL)

# Tipo MSI
msi_type = as.factor(bladder_clin$MSI_TYPE)

# Código Oncotree
oncotree_code = as.factor(bladder_clin$ONCOTREE_CODE)

# Raça
race = as.factor(bladder_clin$RACE)

# Classe da Amostra
sample_class = as.factor(bladder_clin$SAMPLE_CLASS)

# Sexo
sexo = as.factor(bladder_clin$SEX)

# Estado de Fumador
smoking_status = as.factor(bladder_clin$SMOKIMG_STATUS)

# Estado Somático
somatic_status = as.factor(bladder_clin$SOMATIC_STATUS)
```
2. Variáveis Numéricas

2.1. Idade em que a Sequenciação foi Reportada (Anos)
```{r}
# Remover valores ausentes e calcular estatísticas descritivas
summary(na.omit(bladder_clin$AGE_AT_SEQ_REPORTED_YEARS))

# Desvio padrão
sd(na.omit(bladder_clin$AGE_AT_SEQ_REPORTED_YEARS))

# Intervalo Interquartílico
IQR(na.omit(bladder_clin$AGE_AT_SEQ_REPORTED_YEARS))

# Barplot e boxplot
barplot(table(bladder_clin$AGE_AT_SEQ_REPORTED_YEARS))
boxplot(as.numeric(bladder_clin$AGE_AT_SEQ_REPORTED_YEARS),
        main = "Idade em que a Sequência foi Reportada (Anos)",
        xlab = "Idade", horizontal = TRUE)

```
2.2. Sobrevivência Global (Meses desde o Tx com Erdafitinib)
```{r}
# Remover valores ausentes e calcular estatísticas descritivas
summary(na.omit(bladder_clin$ERDAFITINIB_TX_OS_MONTHS))

# Desvio padrão
sd(na.omit(bladder_clin$ERDAFITINIB_TX_OS_MONTHS))

# Intervalo Interquartílico
IQR(na.omit(bladder_clin$ERDAFITINIB_TX_OS_MONTHS))

# Barplot e boxplot
barplot(table(bladder_clin$ERDAFITINIB_TX_OS_MONTHS))
boxplot(as.numeric(bladder_clin$ERDAFITINIB_TX_OS_MONTHS),
        main = "Sobrevivência Global (Meses desde o Tx com Erdafitinib)",
        xlab = "Meses", horizontal = TRUE)


```
2.3. Sobrevivência Sem Progressão (Meses desde o Tx com Erdafitinib)
```{r}
# Remover valores ausentes e calcular estatísticas descritivas
summary(na.omit(bladder_clin$ERDAFITINIB_TX_PFS_MONTHS))

# Desvio padrão
sd(na.omit(bladder_clin$ERDAFITINIB_TX_PFS_MONTHS))

# Intervalo Interquartílico
IQR(na.omit(bladder_clin$ERDAFITINIB_TX_PFS_MONTHS))

# Barplot e boxplot
barplot(table(bladder_clin$ERDAFITINIB_TX_PFS_MONTHS))
boxplot(as.numeric(bladder_clin$ERDAFITINIB_TX_PFS_MONTHS),

```
2.4 Sobrevivência Global em Meses desde o Sequenciamento (SEQUENCING_OS_MONTHS)
```{r}
# Remover valores missing
bladder_clin_filtered <- na.omit(bladder_clin)

# Resumo das medidas de tendencia central e dispersão
summary(bladder_clin_filtered$SEQUENCING_OS_MONTHS)

# Calcular o desvio padrão
sd(bladder_clin_filtered$SEQUENCING_OS_MONTHS)

# Calcular o Intervalo Interquartílico (IQR)
IQR(bladder_clin_filtered$SEQUENCING_OS_MONTHS)

# Barplot
barplot(table(bladder_clin_filtered$SEQUENCING_OS_MONTHS))

# Boxplot
boxplot(as.numeric(bladder_clin_filtered$SEQUENCING_OS_MONTHS),
        horizontal = TRUE)
```
