library(cluster)
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(ape))

# Abrimos el archivo csv que generamos del alineamiento
InputData  <- read.table("/Users/carlosmichelmourradiaz/Documents/alignment_scores2.csv", header = TRUE, sep = ",", quote="")

# Generamos un dataframe unicamente con los datos que nos interesan
datos <- data.frame(Sequence_1 = InputData$Sequence_1, Sequence_2 = InputData$Sequence_2, Dissimilarity = InputData$Dissimilarity)

# Funcion que genera la matriz de scores
Generate_Matrix <- function(datos,rows,cols,value2){
  num_rows <- length(rows)
  num_cols <- length(cols)
  total <- num_cols*num_rows
  mat <- matrix(1:total,nrow = num_rows, ncol = num_cols)
  rownames(mat) <- rows
  colnames(mat) <- cols
  
  # Recorrer los valores y rellenar la matriz 
  for (row in rows){
    for (col in cols){
      index <- datos$Sequence_1 == row & datos$Sequence_2 == col
      if (length(datos$Dissimilarity[index]) > 0){
        mat[row,col] <- datos$Dissimilarity[index]
      }
      # Si no se encontro el alineamiento rellenar con un 0
      else{
        mat[row,col] <- value2
      }
    }
    
  }
  return(mat)
}

# Obtenemos los genes para las filas y columnas de la matriz
row_genes <- unique(datos$Sequence_1)
col_genes <- unique(datos$Sequence_2)

# Generamos la matriz
valores <- Generate_Matrix(datos,row_genes,col_genes,NA)

# hClust
ccom <- hclust(dist(valores), method = "complete")
csin <- hclust(dist(valores), method = "single")
cave <- hclust(dist(valores), method = "average")
cward <- hclust(dist(valores), method = "ward.D2")

# sacar los coeficientes de aglomeración de cada modelo de clustering
coeffcom <- coef(ccom)
coeffsin <- coef(csin)
coeffave <- coef(cave)
coeffward <- coef(cward)

# obtenemos los coeficientes de aglomeración
print("Complete")
coeffcom
print("Single")
coeffsin
print("Average")
coeffave
print("Ward")
coeffward

# Dendograma con complete
plot (ccom, hang = -1, main = "hclust Dendogram Complete")
cls3 <- cutree(ccom, k=3)

# plot de dendograma con single 
plot (csin, hang = -1, main = "hclust Dendogram Single")
cls4 <- cutree(cave, k=3)

# plot de dendograma con average
plot (cave, hang = -1, main = "hclust Dendogram Average")
cls5 <- cutree(cave, k=3)

# Plot de dendograma con ward 
plot (cward, hang = -1, main = "hclust Dendogram Ward")
cls6 <- cutree(cward, k=3)

# Guardar los dendogramas en formato tree
my_treecom <- as.phylo(ccom)
write.tree(phy=my_treecom, file="/Users/carlosmichelmourradiaz/Documents/complete_tree.tree")

my_treesin <- as.phylo(csin)
write.tree(phy=my_treesin, file="/Users/carlosmichelmourradiaz/Documents/single_tree.tree")

my_treecave <- as.phylo(cave)
write.tree(phy=my_treecave, file="/Users/carlosmichelmourradiaz/Documents/average_tree.tree")

my_treeward <- as.phylo(cward)
write.tree(phy=my_treeward, file="/Users/carlosmichelmourradiaz/Documents/ward_tree.tree")

# Agnes

# Generamos la matriz
valores2 <- Generate_Matrix(datos,row_genes,col_genes,0)

# Generamos los modelos de clustering con los diferentes metodos
aClustcom <- agnes(valores2, method = "complete")
aClustsin <- agnes(valores2, method = "single")
aClustave <- agnes(valores2, method = "average")
aClustward <- agnes(valores2, method = "ward")


# Observamos el dendograma generado con Agnes para el metodo complete
pltree(aClustcom, cex = 0.6, hang = -1, main = "agnes Dendrogram Complete")
aCoeffcom <- aClustcom$ac
rect.hclust(as.hclust(aClustcom), k = 6, border=2:4)
aCls3 <- cutree(as.hclust(aClustcom), k = 6)
fviz_cluster(list(data = valores2, cluster = aCls3))

# Observamos el dendograma generado con Agnes para el metodo single
pltree(aClustsin, cex = 0.6, hang = -1, main = "agnes Dendrogram Single")
aCoeffsin <- aClustsin$ac
rect.hclust(as.hclust(aClustsin), k=6, border=2:4)
aCls4 <- cutree(as.hclust(aClustsin), k = 6)
fviz_cluster(list(data = valores2, cluster = aCls4))

# Observamos el dendograma generado con Agnes para el metodo ward
pltree(aClustward, cex = 0.6, hang = -1, main = "agnes Dendrogram Ward")
aCoeffward <- aClustward$ac
rect.hclust(as.hclust(aClustward), k=6, border=2:4)
aCls6 <- cutree(as.hclust(aClustward), k = 6)
fviz_cluster(list(data = valores2, cluster = aCls6))

# Observamos los coeficientes de aglomeración
print("Complete")
aCoeffcom
print("Single")
aCoeffsin
print("Average")
aCoeffave
print("Ward")
aCoeffward


