# Verifica se o BiocManager está instalado e instala se necessário
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("cBioPortalData", quietly = TRUE))
  BiocManager::install("cBioPortalData")

# Carregamento dos packages
library(cBioPortalData)

# Inicialização do API cBioPortal
cbio <- cBioPortal()

# Obter dados clínicos
bladder_clin <- clinicalData(api = cbio, studyId = "bladder_msk_2023")
#bladder_clin
dim(bladder_clin)  # Deve retornar 526 linhas e 29 colunas

# Carrega os dados do cBioPortal
bladder = cBioDataPack("bladder_msk_2023", ask = FALSE)
bladder

################################################################################
# Lista os nomes dos ensaios disponíveis no pacote dos dados do CBioPortal
names(assays(bladder))
# names(rowData(bladder)) # tipos de metadados associados a cada gene
# names(colData(bladder)) # tipos de metadados associados a cada amostra

# Summary dos dados
summary(bladder)

# Acessando os dados de cna e convertendo para um dataframe
dados_cna = assays(bladder)$cna
dados_cna_df = as.data.frame(dados_cna)

## Aceder aos ficheiros dos metadados (formato .txt) presentes dentro da pasta "bladder_msk_2023"
setwd("bladder_msk_2023")
meta_cna_hg19_seg <- read.delim("meta_cna_hg19_seg.txt", sep = ":", 
                                header = FALSE, 
                                row.names = 1, 
                                col.names = c("", "data")
                                )
meta_cna_hg19_seg

meta_cna <- read.delim("meta_cna.txt", sep = ":", 
                       header = FALSE, 
                       row.names = 1, 
                       col.names = c("", "data")
)
meta_cna

# Metadados do estudo
meta_study <- read.delim("meta_study.txt", sep = ":", 
                         header = FALSE, 
                         row.names = 1, 
                         col.names = c("", "data")
                         )
meta_study["citation", "data"] # Citação


# Acessando os dados de cna_hg19 e convertendo para um dataframe
dados_cna_hg19 = assays(bladder)$cna_hg19.seg
dados_cna_hg19_df = as.data.frame(dados_cna_hg19)

# Acessando os dados de mutação e convertendo para um dataframe
dados_mutacao = assays(bladder)$mutations
dados_mutacao_df = as.data.frame(dados_mutacao)

# Verificando a classe e dimensões dos dados de mutação
dim(dados_mutacao_df)

################################################################################

## Definir as variaveis com as.factor 
# age = as.factor(bladder_clin$AGE_AT_SEQ_REPORTED_YEARS)
# cancer_type = as.factor(bladder_clin$CANCER_TYPE_DETAILED)
# disease_state = as.factor(bladder_clin$DISEASE_STATE)
# erdafitinib_treatment = as.factor(bladder_clin$ERDAFITINIB_TREATED)
# overall_status = as.factor(bladder_clin$ERDAFITINIB_TX_OS_STATUS)
# progression_status = as.factor(bladder_clin$ERDAFITINIB_TX_PFS_STATUS)
# gene_panel = as.factor(bladder_clin$GENE_PANEL)
# msi_type = as.factor(bladder_clin$MSI_TYPE)
# oncotree_code = as.factor(bladder_clin$ONCOTREE_CODE)
# race = as.factor(bladder_clin$RACE)
# sample_class = as.factor(bladder_clin$SAMPLE_CLASS)
# overall_seq = as.factor(bladder_clin$ERDAFITINIB_TX_OS_MONTHS)
# sexo = as.factor(bladder_clin$SEX)
# smoking_status = as.factor(bladder_clin$SMOKIMG_STATUS)
# somatic_status = as.factor(bladder_clin$SOMATIC_STATUS)

"
1. Variáveis Numéricas:
   - Age at Which Sequencing was Reported (Years)
   - Overall Survival Months since Erdafitinib Tx
   - Progress Free Survival Months since Erdafitinib Tx
   - Overall Survival Months since Sequencing
     - **Análise Estatística Descritiva:** Calcule medidas de tendência central 
     (média, mediana) e dispersão (desvio padrão, intervalo) para entender a 
     distribuição da idade dos pacientes.
     - **Exploração com Gráficos:** Um histograma pode ser útil para visualizar 
     a distribuição da idade dos pacientes.
"

## Age at Which Sequencing was Reported (Years) (AGE_AT_SEQ_REPORTED_YEARS)
# Remover valores ausentes e calcular o valor mínimo, 1º quartil, mediana, média, 3º quartil e valor máximo da variável.
summary(na.omit(bladder_clin$AGE_AT_SEQ_REPORTED_YEARS))

# Remover valores ausentes e calcular o desvio padrão da variável.
sd(na.omit(bladder_clin$AGE_AT_SEQ_REPORTED_YEARS))

# Remover valores ausentes e calcular o Intervalo interquartílico.
IQR(na.omit(bladder_clin$AGE_AT_SEQ_REPORTED_YEARS))

# Barplot e boxplot da distribuição de idades
barplot(table(bladder_clin$AGE_AT_SEQ_REPORTED_YEARS))
boxplot(as.numeric(bladder_clin$AGE_AT_SEQ_REPORTED_YEARS, 
        main = "Idade em que a Sequência foi Reportada (Anos)", 
        xlab = "Idade", horizontal = TRUE))



## Overall Survival Months since Erdafitinib Tx (ERDAFITINIB_TX_OS_MONTHS)
# Remover valores ausentes e calcular o valor mínimo, 1º quartil, mediana, média, 3º quartil e valor máximo da variável.
summary(na.omit(bladder_clin$ERDAFITINIB_TX_OS_MONTHS))

# Remover valores ausentes e calcular o desvio padrão da variável.
sd(na.omit(bladder_clin$ERDAFITINIB_TX_OS_MONTHS))

# Remover valores ausentes e calcular o Intervalo interquartílico.
IQR(na.omit(bladder_clin$ERDAFITINIB_TX_OS_MONTHS))

# Barplot e boxplot
barplot(table(bladder_clin$ERDAFITINIB_TX_OS_MONTHS))
boxplot(as.numeric(bladder_clin$ERDAFITINIB_TX_OS_MONTHS,
        main = "Sobrevivência Geral (Meses desde o Tratamento com Erdafitinib)", 
        xlab = "Meses", 
        horizontal = TRUE))



## Progress Free Survival Months since Erdafitinib Tx (ERDAFITINIB_TX_PFS_MONTHS)
# Remover valores ausentes e calcular o valor mínimo, 1º quartil, mediana, média, 3º quartil e valor máximo da variável.
summary(na.omit(bladder_clin$ERDAFITINIB_TX_PFS_MONTHS))

# Remover valores ausentes e calcular o desvio padrão da variável
sd(na.omit(bladder_clin$ERDAFITINIB_TX_PFS_MONTHS))

# Remover valores ausentes e calcular o Intervalo interquartílico
IQR(na.omit(bladder_clin$ERDAFITINIB_TX_PFS_MONTHS))

# Barplot e boxplot
barplot(table(bladder_clin$ERDAFITINIB_TX_PFS_MONTHS))
boxplot(as.numeric(bladder_clin$ERDAFITINIB_TX_PFS_MONTHS, 
        horizontal = TRUE))



## Overall Survival Months since Sequencing (SEQUENCING_OS_MONTHS)
# Remover valores ausentes e calcular o valor mínimo, 1º quartil, mediana, média, 3º quartil e valor máximo da variável.
summary(na.omit(bladder_clin$SEQUENCING_OS_MONTHS))

# Remover valores ausentes e calcular o desvio padrão da variável
sd(na.omit(bladder_clin$SEQUENCING_OS_MONTHS))

# Remover valores ausentes e calcular o Intervalo interquartílico.
IQR(na.omit(bladder_clin$SEQUENCING_OS_MONTHS))

# Barplot e boxplot
barplot(table(bladder_clin$SEQUENCING_OS_MONTHS))
boxplot(as.numeric(bladder_clin$SEQUENCING_OS_MONTHS, 
        horizontal = TRUE))


"
2. **Variáveis Categóricas:**
   - Cancer Type
   - Cancer Type Detailed
   - Disease State
   - Erdafitinib Treated
   - Overall Survival Status since Erdafitinib Tx
   - Progression Free Status since Erdafitinib Tx
   - Race
   - Sample Class
   - Sex
   - Smoking Status
   - Somatic Status
     - **Análise Estatística Descritiva:** Calcule as frequências absolutas e 
     relativas para cada categoria.
     - **Exploração com Gráficos:** Gráficos de barras ou gráficos de pizza podem ser úteis para visualizar a distribuição das diferentes categorias.
   - MSI Type
   - Oncotree Code
     - **Análise Estatística Descritiva:** Calcule as frequências absolutas e 
     relativas para cada tipo de MSI (instabilidade de microssatélite) e para 
     cada código Oncotree.
     - **Exploração com Gráficos:** Gráficos de barras podem ser úteis para 
     visualizar a distribuição dos tipos de MSI e dos códigos Oncotree.
"

## Cancer Type - é semple "Blader Type"
## Cancer Type Detailed 

# Gráfico de barras para visualizar as frequências absolutas
barplot(table(bladder_clin$CANCER_TYPE_DETAILED), 
        main = "Frequência Absoluta por Tipo de Cancro",
        xlab = "Tipo de Câncer", ylab = "Frequência Absoluta")

# Calcular as frequências absolutas e relativas
freq_abs <- table(bladder_clin$CANCER_TYPE_DETAILED)
freq_rel <- prop.table(freq_abs)

# Criar um data frame com as frequências absolutas e relativas
freq_data <- data.frame(Categoria = names(freq_abs),
                        Frequencia_Absoluta = as.numeric(freq_abs),
                        Frequencia_Relativa = as.numeric(freq_rel))

# Gráfico de barras para as frequências absolutas
barplot(freq_abs, main = "Frequência Absoluta por Tipo de Câncer",
        xlab = "Tipo de Câncer", ylab = "Frequência Absoluta")

# Gráfico de barras para as frequências relativas
barplot(freq_rel, main = "Distribuição Relativa por Tipo de Câncer",
        xlab = "Tipo de Câncer", ylab = "Frequência Relativa")

# Criar o gráfico de pizza
pie(freq_rel, main = "Distribuição Relativa por Tipo de Cancro",
    labels = paste(names(freq_rel), ": ", round(freq_rel * 100, 2), "%"), 
    col = rainbow(length(freq_rel)),
    cex = 0.8)



## Disease State

# Gráfico de barras para visualizar as frequências absolutas
barplot(table(bladder_clin$DISEASE_STATE), 
        main = "Frequência Absoluta por Categoria",
        xlab = "Categoria", ylab = "Frequência Absoluta")

# Calcular as frequências absolutas e relativas
freq_abs <- table(bladder_clin$DISEASE_STATE)
freq_rel <- prop.table(freq_abs)

# Criar um data frame com as frequências absolutas e relativas
freq_data <- data.frame(Categoria = names(freq_abs),
                        Frequencia_Absoluta = as.numeric(freq_abs),
                        Frequencia_Relativa = as.numeric(freq_rel))

# Gráfico de barras para as frequências absolutas
barplot(freq_abs, main = "Frequência Absoluta por Categoria",
        xlab = "Categoria", ylab = "Frequência Absoluta")

# Gráfico de barras para as frequências relativas
barplot(freq_rel, main = "Distribuição Relativa por Categoria",
        xlab = "Categoria", ylab = "Frequência Relativa")

# Criar o gráfico de pizza
pie(freq_rel, main = "Distribuição Relativa por Categoria",
    labels = paste(names(freq_rel), ": ", round(freq_rel * 100, 2), "%"), 
    col = rainbow(length(freq_rel)),
    cex = 0.8)

## Erdafitinib Treated

# Gráfico de barras para visualizar as frequências absolutas
barplot(table(bladder_msk_2023_clinical_data$Erdafitinib.Treated), 
        main = "Frequência Absoluta por Categoria",
        xlab = "Categoria", ylab = "Frequência Absoluta")

# Calcular as frequências absolutas e relativas
freq_abs <- table(bladder_msk_2023_clinical_data$Erdafitinib.Treated)
freq_rel <- prop.table(freq_abs)

# Criar um data frame com as frequências absolutas e relativas
freq_data <- data.frame(Categoria = names(freq_abs),
                        Frequencia_Absoluta = as.numeric(freq_abs),
                        Frequencia_Relativa = as.numeric(freq_rel))

# Gráfico de barras para as frequências absolutas
barplot(freq_abs, main = "Frequência Absoluta por Categoria",
        xlab = "Categoria", ylab = "Frequência Absoluta")

# Gráfico de barras para as frequências relativas
barplot(freq_rel, main = "Distribuição Relativa por Categoria",
        xlab = "Categoria", ylab = "Frequência Relativa")

# Criar o gráfico de pizza
pie(freq_rel, main = "Distribuição Relativa por Categoria",
    labels = paste(names(freq_rel), ": ", round(freq_rel * 100, 2), "%"), 
    col = rainbow(length(freq_rel)),
    cex = 0.8)

## Overall Survival Status since Erdafitinib Tx

# Gráfico de barras para visualizar as frequências absolutas
barplot(table(bladder_msk_2023_clinical_data$Overall.Survival.Status.since.Erdafitinib.Tx), 
        main = "Frequência Absoluta por Categoria",
        xlab = "Categoria", ylab = "Frequência Absoluta")

# Calcular as frequências absolutas e relativas
freq_abs <- table(bladder_msk_2023_clinical_data$Overall.Survival.Status.since.Erdafitinib.Tx)
freq_rel <- prop.table(freq_abs)

# Criar um data frame com as frequências absolutas e relativas
freq_data <- data.frame(Categoria = names(freq_abs),
                        Frequencia_Absoluta = as.numeric(freq_abs),
                        Frequencia_Relativa = as.numeric(freq_rel))

# Gráfico de barras para as frequências absolutas
barplot(freq_abs, main = "Frequência Absoluta por Categoria",
        xlab = "Categoria", ylab = "Frequência Absoluta")

# Gráfico de barras para as frequências relativas
barplot(freq_rel, main = "Distribuição Relativa por Categoria",
        xlab = "Categoria", ylab = "Frequência Relativa")

# Criar o gráfico de pizza
pie(freq_rel, main = "Distribuição Relativa por Categoria",
    labels = paste(names(freq_rel), ": ", round(freq_rel * 100, 2), "%"), 
    col = rainbow(length(freq_rel)),
    cex = 0.8)



## Progression Free Status since Erdafitinib Tx

# Gráfico de barras para visualizar as frequências absolutas
barplot(table(bladder_msk_2023_clinical_data$Progression.Free.Status.since.Erdafitinib.Tx), 
        main = "Frequência Absoluta por Categoria",
        xlab = "Categoria", ylab = "Frequência Absoluta")

# Calcular as frequências absolutas e relativas
freq_abs <- table(bladder_msk_2023_clinical_data$Progression.Free.Status.since.Erdafitinib.Tx)
freq_rel <- prop.table(freq_abs)

# Criar um data frame com as frequências absolutas e relativas
freq_data <- data.frame(Categoria = names(freq_abs),
                        Frequencia_Absoluta = as.numeric(freq_abs),
                        Frequencia_Relativa = as.numeric(freq_rel))

# Gráfico de barras para as frequências absolutas
barplot(freq_abs, main = "Frequência Absoluta por Categoria",
        xlab = "Categoria", ylab = "Frequência Absoluta")

# Gráfico de barras para as frequências relativas
barplot(freq_rel, main = "Distribuição Relativa por Categoria",
        xlab = "Categoria", ylab = "Frequência Relativa")

# Criar o gráfico de pizza
pie(freq_rel, main = "Distribuição Relativa por Categoria",
    labels = paste(names(freq_rel), ": ", round(freq_rel * 100, 2), "%"), 
    col = rainbow(length(freq_rel)),
    cex = 0.8)


## Race

# Gráfico de barras para visualizar as frequências absolutas
barplot(table(bladder_msk_2023_clinical_data$Race), 
        main = "Frequência Absoluta por Categoria",
        xlab = "Categoria", ylab = "Frequência Absoluta")

# Calcular as frequências absolutas e relativas
freq_abs <- table(bladder_msk_2023_clinical_data$Race)
freq_rel <- prop.table(freq_abs)

# Criar um data frame com as frequências absolutas e relativas
freq_data <- data.frame(Categoria = names(freq_abs),
                        Frequencia_Absoluta = as.numeric(freq_abs),
                        Frequencia_Relativa = as.numeric(freq_rel))

# Gráfico de barras para as frequências absolutas
barplot(freq_abs, main = "Frequência Absoluta por Categoria",
        xlab = "Categoria", ylab = "Frequência Absoluta")

# Gráfico de barras para as frequências relativas
barplot(freq_rel, main = "Distribuição Relativa por Categoria",
        xlab = "Categoria", ylab = "Frequência Relativa")

# Criar o gráfico de pizza
pie(freq_rel, main = "Distribuição Relativa por Categoria",
    labels = paste(names(freq_rel), ": ", round(freq_rel * 100, 2), "%"), 
    col = rainbow(length(freq_rel)),
    cex = 0.8)

## Sample Class

# Gráfico de barras para visualizar as frequências absolutas
barplot(table(bladder_msk_2023_clinical_data$Sample.Class), 
        main = "Frequência Absoluta por Categoria",
        xlab = "Categoria", ylab = "Frequência Absoluta")

# Calcular as frequências absolutas e relativas
freq_abs <- table(bladder_msk_2023_clinical_data$Sample.Class)
freq_rel <- prop.table(freq_abs)

# Criar um data frame com as frequências absolutas e relativas
freq_data <- data.frame(Categoria = names(freq_abs),
                        Frequencia_Absoluta = as.numeric(freq_abs),
                        Frequencia_Relativa = as.numeric(freq_rel))

# Gráfico de barras para as frequências absolutas
barplot(freq_abs, main = "Frequência Absoluta por Categoria",
        xlab = "Categoria", ylab = "Frequência Absoluta")

# Gráfico de barras para as frequências relativas
barplot(freq_rel, main = "Distribuição Relativa por Categoria",
        xlab = "Categoria", ylab = "Frequência Relativa")

# Criar o gráfico de pizza
pie(freq_rel, main = "Distribuição Relativa por Categoria",
    labels = paste(names(freq_rel), ": ", round(freq_rel * 100, 2), "%"), 
    col = rainbow(length(freq_rel)),
    cex = 0.8)

## Sex

# Gráfico de barras para visualizar as frequências absolutas
barplot(table(bladder_msk_2023_clinical_data$Sex), 
        main = "Frequência Absoluta por Categoria",
        xlab = "Categoria", ylab = "Frequência Absoluta")

# Calcular as frequências absolutas e relativas
freq_abs <- table(bladder_msk_2023_clinical_data$Sex)
freq_rel <- prop.table(freq_abs)

# Criar um data frame com as frequências absolutas e relativas
freq_data <- data.frame(Categoria = names(freq_abs),
                        Frequencia_Absoluta = as.numeric(freq_abs),
                        Frequencia_Relativa = as.numeric(freq_rel))

# Gráfico de barras para as frequências absolutas
barplot(freq_abs, main = "Frequência Absoluta por Categoria",
        xlab = "Categoria", ylab = "Frequência Absoluta")

# Gráfico de barras para as frequências relativas
barplot(freq_rel, main = "Distribuição Relativa por Categoria",
        xlab = "Categoria", ylab = "Frequência Relativa")

# Criar o gráfico de pizza
pie(freq_rel, main = "Distribuição Relativa por Categoria",
    labels = paste(names(freq_rel), ": ", round(freq_rel * 100, 2), "%"), 
    col = rainbow(length(freq_rel)),
    cex = 0.8)


## Smoking Status

# Gráfico de barras para visualizar as frequências absolutas
barplot(table(bladder_msk_2023_clinical_data$Smoking.Status), 
        main = "Frequência Absoluta por Categoria",
        xlab = "Categoria", ylab = "Frequência Absoluta")

# Calcular as frequências absolutas e relativas
freq_abs <- table(bladder_msk_2023_clinical_data$Smoking.Status)
freq_rel <- prop.table(freq_abs)

# Criar um data frame com as frequências absolutas e relativas
freq_data <- data.frame(Categoria = names(freq_abs),
                        Frequencia_Absoluta = as.numeric(freq_abs),
                        Frequencia_Relativa = as.numeric(freq_rel))

# Gráfico de barras para as frequências absolutas
barplot(freq_abs, main = "Frequência Absoluta por Categoria",
        xlab = "Categoria", ylab = "Frequência Absoluta")

# Gráfico de barras para as frequências relativas
barplot(freq_rel, main = "Distribuição Relativa por Categoria",
        xlab = "Categoria", ylab = "Frequência Relativa")

# Criar o gráfico de pizza
pie(freq_rel, main = "Distribuição Relativa por Categoria",
    labels = paste(names(freq_rel), ": ", round(freq_rel * 100, 2), "%"), 
    col = rainbow(length(freq_rel)),
    cex = 0.8)

## Somatic Status
# Gráfico de barras para visualizar as frequências absolutas
barplot(table(bladder_msk_2023_clinical_data$Somatic.Status), 
        main = "Frequência Absoluta por Categoria",
        xlab = "Categoria", ylab = "Frequência Absoluta")

# Calcular as frequências absolutas e relativas
freq_abs <- table(bladder_msk_2023_clinical_data$Somatic.Status)
freq_rel <- prop.table(freq_abs)

# Criar um data frame com as frequências absolutas e relativas
freq_data <- data.frame(Categoria = names(freq_abs),
                        Frequencia_Absoluta = as.numeric(freq_abs),
                        Frequencia_Relativa = as.numeric(freq_rel))

# Gráfico de barras para as frequências absolutas
barplot(freq_abs, main = "Frequência Absoluta por Categoria",
        xlab = "Categoria", ylab = "Frequência Absoluta")

# Gráfico de barras para as frequências relativas
barplot(freq_rel, main = "Distribuição Relativa por Categoria",
        xlab = "Categoria", ylab = "Frequência Relativa")

# Criar o gráfico de pizza
pie(freq_rel, main = "Distribuição Relativa por Categoria",
    labels = paste(names(freq_rel), ": ", round(freq_rel * 100, 2), "%"), 
    col = rainbow(length(freq_rel)),
    cex = 0.8)


## MSI Type

# Gráfico de barras para visualizar as frequências absolutas
barplot(table(bladder_msk_2023_clinical_data$MSI.Type), 
        main = "Frequência Absoluta por Categoria",
        xlab = "Categoria", ylab = "Frequência Absoluta")

# Calcular as frequências absolutas e relativas
freq_abs <- table(bladder_msk_2023_clinical_data$MSI.Type)
freq_rel <- prop.table(freq_abs)

# Criar um data frame com as frequências absolutas e relativas
freq_data <- data.frame(Categoria = names(freq_abs),
                        Frequencia_Absoluta = as.numeric(freq_abs),
                        Frequencia_Relativa = as.numeric(freq_rel))

# Gráfico de barras para as frequências absolutas
barplot(freq_abs, main = "Frequência Absoluta por Categoria",
        xlab = "Categoria", ylab = "Frequência Absoluta")

# Gráfico de barras para as frequências relativas
barplot(freq_rel, main = "Distribuição Relativa por Categoria",
        xlab = "Categoria", ylab = "Frequência Relativa")

# Criar o gráfico de pizza
pie(freq_rel, main = "Distribuição Relativa por Categoria",
    labels = paste(names(freq_rel), ": ", round(freq_rel * 100, 2), "%"), 
    col = rainbow(length(freq_rel)),
    cex = 0.8)

## Oncotree Code

# Gráfico de barras para visualizar as frequências absolutas
barplot(table(bladder_msk_2023_clinical_data$Oncotree.Code), 
        main = "Frequência Absoluta por Categoria",
        xlab = "Categoria", ylab = "Frequência Absoluta")

# Calcular as frequências absolutas e relativas
freq_abs <- table(bladder_msk_2023_clinical_data$Oncotree.Code)
freq_rel <- prop.table(freq_abs)

# Criar um data frame com as frequências absolutas e relativas
freq_data <- data.frame(Categoria = names(freq_abs),
                        Frequencia_Absoluta = as.numeric(freq_abs),
                        Frequencia_Relativa = as.numeric(freq_rel))

# Gráfico de barras para as frequências absolutas
barplot(freq_abs, main = "Frequência Absoluta por Categoria",
        xlab = "Categoria", ylab = "Frequência Absoluta")

# Gráfico de barras para as frequências relativas
barplot(freq_rel, main = "Distribuição Relativa por Categoria",
        xlab = "Categoria", ylab = "Frequência Relativa")

# Criar o gráfico de pizza
pie(freq_rel, main = "Distribuição Relativa por Categoria",
    labels = paste(names(freq_rel), ": ", round(freq_rel * 100, 2), "%"), 
    col = rainbow(length(freq_rel)),
    cex = 0.8)

"""
3. **Variáveis Contínuas:**
   - Best Response to Erdafitinib
   - Fraction Genome Altered
   - Mutation Count
   - TMB (nonsynonymous)
   - Tumor Purity
     - **Análise Estatística Descritiva:** Calcule medidas de tendência central e dispersão para entender a distribuição dessas variáveis contínuas.
     - **Exploração com Gráficos:** Um boxplot pode ser útil para visualizar a distribuição dessas variáveis.
"""

## Best Response to Erdafitinib

# Remover valores ausentes e calcular o valor minimo, 1º quartil, mediana, média, 3º quartil e valor máximo da variavel.
summary_best_response <- summary(na.omit(bladder_msk_2023_clinical_data$Best.Response.to.Erdafitinib))
# Visualizar os valores
summary_best_response

# Remover valores ausentes e calcular o desvio padrão da variavel.
sd_best_response <- sd(na.omit(bladder_msk_2023_clinical_data$Best.Response.to.Erdafitinib))
# Visualizar o valor
sd_best_response

# Remover valores ausentes e calcular o Intervalo interquartilico.
IQR_best_response <- IQR(na.omit(bladder_msk_2023_clinical_data$Best.Response.to.Erdafitinib))
# Visualizar o valor
IQR_best_response

plot(table(bladder_msk_2023_clinical_data$Best.Response.to.Erdafitinib))

# Criar o gráfico de caixa
boxplot(bladder_msk_2023_clinical_data$Best.Response.to.Erdafitinib, 
        main = "Melhor resposta ao Erdafitinib", xlab = "", horizontal=T)


"""
3. Variáveis Contínuas:
   - Fraction Genome Altered
   - Mutation Count
   - TMB (nonsynonymous)
   - Tumor Purity
  
Análise Estatística Descritiva: Calcule medidas de tendência central e dispersão para entender a distribuição dessas variáveis contínuas.
- Exploração com Gráficos: Um boxplot pode ser útil para visualizar a distribuição dessas variáveis.

"""
## Fraction Genome Altered

# Remover valores ausentes e calcular o valor minimo, 1º quartil, mediana, média, 3º quartil e valor máximo da variavel.
summary_ <- summary(na.omit(bladder_msk_2023_clinical_data$"Fraction.Genome.Altered"))
# Visualizar os valores
summary_

# Remover valores ausentes e calcular o desvio padrão da variavel.
sd_ <- sd(na.omit(bladder_msk_2023_clinical_data$"Fraction.Genome.Altered"))
# Visualizar o valor
sd_

# Remover valores ausentes e calcular o Intervalo interquartilico.
IQR_ <- IQR(na.omit(bladder_msk_2023_clinical_data$"Fraction.Genome.Altered"))
# Visualizar o valor
IQR_

# Calcular a variância
var_ <- var(na.omit(bladder_msk_2023_clinical_data$"Fraction.Genome.Altered"))

# Visualizar o valor da variância
var_

# Criar o gráfico de caixa
boxplot(bladder_msk_2023_clinical_data$"Fraction.Genome.Altered", horizontal=T)



### Mutation Count
# Remover valores ausentes e calcular o valor minimo, 1º quartil, mediana, média, 3º quartil e valor máximo da variavel.
summary_ <- summary(na.omit(bladder_msk_2023_clinical_data$"Mutation.Count"))
# Visualizar os valores
summary_

# Remover valores ausentes e calcular o desvio padrão da variavel.
sd_ <- sd(na.omit(bladder_msk_2023_clinical_data$"Mutation.Count"))
# Visualizar o valor
sd_

# Remover valores ausentes e calcular o Intervalo interquartilico.
IQR_ <- IQR(na.omit(bladder_msk_2023_clinical_data$"Mutation.Count"))
# Visualizar o valor
IQR_

# Calcular a variância
var_ <- var(na.omit(bladder_msk_2023_clinical_data$"Mutation.Count"))

# Visualizar o valor da variância
var_

barplot(table(bladder_msk_2023_clinical_data$"Mutation.Count"))

# Criar o gráfico de caixa
boxplot(bladder_msk_2023_clinical_data$"Mutation.Count", horizontal=T)


### TMB (nonsynonymous)
# Remover valores ausentes e calcular o valor minimo, 1º quartil, mediana, média, 3º quartil e valor máximo da variavel.
summary_ <- summary(na.omit(bladder_msk_2023_clinical_data$"TMB..nonsynonymous."))
# Visualizar os valores
summary_

# Remover valores ausentes e calcular o desvio padrão da variavel.
sd_ <- sd(na.omit(bladder_msk_2023_clinical_data$"TMB..nonsynonymous."))
# Visualizar o valor
sd_

# Remover valores ausentes e calcular o Intervalo interquartilico.
IQR_ <- IQR(na.omit(bladder_msk_2023_clinical_data$"TMB..nonsynonymous."))
# Visualizar o valor
IQR_

barplot(table(bladder_msk_2023_clinical_data$"TMB..nonsynonymous."))
plot(table(bladder_msk_2023_clinical_data$"TMB..nonsynonymous."))

# Criar o gráfico de caixa
boxplot(bladder_msk_2023_clinical_data$"TMB..nonsynonymous.", horizontal=T)


#### Tumor Purity
# Remover valores ausentes e calcular o valor minimo, 1º quartil, mediana, média, 3º quartil e valor máximo da variavel.
summary_ <- summary(na.omit(bladder_msk_2023_clinical_data$"Tumor.Purity"))
# Visualizar os valores
summary_

# Remover valores ausentes e calcular o desvio padrão da variavel.
sd_ <- sd(na.omit(bladder_msk_2023_clinical_data$"Tumor.Purity"))
# Visualizar o valor
sd_

# Remover valores ausentes e calcular o Intervalo interquartilico.
IQR_ <- IQR(na.omit(bladder_msk_2023_clinical_data$"Tumor.Purity"))
# Visualizar o valor
IQR_

barplot(table(bladder_msk_2023_clinical_data$"Tumor.Purity"))
plot(table(bladder_msk_2023_clinical_data$"Tumor.Purity"))

# Criar o gráfico de caixa
boxplot(bladder_msk_2023_clinical_data$"Tumor.Purity", horizontal=T)

