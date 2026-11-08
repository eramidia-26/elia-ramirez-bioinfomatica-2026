##########################################################################################
#####################                                        			##################
#####################    TRABAJO ORDINARIA BIOINFORMÁTICA R             ##################
#####################                                       			##################
##########################################################################################

# Nombre: ________Elia Ramirez Díaz-Saavedra______________    

# Colores que se mantienen en todo el trabajo
col_control   <- "steelblue"
col_condicion <- "tomato"


### 1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente. 
### Recuerda incluir los nombres de las columnas como el encabezado.
### Examina el archivo con las funciones head(), tail() y dim(). (1 pto)

archivo <- "archivo-entrada.txt"
if (!file.exists(archivo)) archivo <- file.choose()

# header = TRUE: la primera fila contiene los nombres de las columnas
# sep = "\t": columnas separadas por tabulador
# check.names = FALSE: respeta nombres como "Gene Name" o "C1-TPM"
datos <- read.table(archivo, header = TRUE, sep = "\t", quote = "",
                    check.names = FALSE, stringsAsFactors = FALSE)

head(datos)
tail(datos)
dim(datos)     # filas = genes, columnas = gen + 6 muestras

# Nombres de gen como nombres de fila (make.unique por si hubiera repetidos)
rownames(datos) <- make.unique(as.character(datos[, 1]))

# Columnas de cada grupo: 2:4 = Control (C2, C1, C3); 5:7 = Condición (1W5-TPM, 1W6-TPM, 1W5-COV)
idx_control   <- 2:4
idx_condicion <- 5:7


### 2.- En nuestra tabla tenemos los valores de 619 genes para tres réplicas
### del Control y tres réplicas de una Condición. Calcula la media y la desviación
### estándar del control y del tratamiento para cada gen. Puedes emplear cualquier
### función que consideres como apply(). Añade esos valores como columnas nuevas a tu tabla original. 
### Serán 4 columnas nuevas en total, la media de Control, la DE de Control, la media de Condición y la DE de Condición. (1,5 pto)

# apply(..., 1, f) aplica la función f a cada fila (cada gen)
datos$Media_Control   <- apply(datos[, idx_control],   1, mean)
datos$DE_Control      <- apply(datos[, idx_control],   1, sd)
datos$Media_Condicion <- apply(datos[, idx_condicion], 1, mean)
datos$DE_Condicion    <- apply(datos[, idx_condicion], 1, sd)

head(datos)
dim(datos)     # deben aparecer 4 columnas más que al principio


### 3.- En una sola imagen, muestra un histograma para la media del Control y otro para la media de la Condicíon.
### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).

par(mfrow = c(1, 2))
hist(datos$Media_Control, breaks = 30, col = col_control,
     main = "Media del Control", xlab = "TPM medio", ylab = "Número de genes")
hist(datos$Media_Condicion, breaks = 30, col = col_condicion,
     main = "Media de la Condición", xlab = "TPM medio", ylab = "Número de genes")
par(mfrow = c(1, 1))


### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control. 
### Añade el título, colores y una una leyenda en la esquina superior izquierda. ( 1,5 pto)

# Los TPM tienen una distribución muy asimétrica (unos pocos genes con valores
# enormes). Para que el boxplot se pueda leer se representa log2(media + 1);
# el +1 evita problemas con las medias iguales a 0.
boxplot(log2(datos$Media_Control + 1), log2(datos$Media_Condicion + 1),
        names = c("Control", "Condición"),
        col   = c(col_control, col_condicion),
        main  = "Expresión media por gen: Control frente a Condición",
        ylab  = "log2(TPM medio + 1)")
legend("topleft", legend = c("Control", "Condición"),
       fill = c(col_control, col_condicion), bty = "n")


### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente: (2 pto)
### 		Realiza un test de comparativa de medias para cada gen entre condición y control. 
###			Guarda el p-value de la comparativa y anádelo como una columna extra a tu tabla original.

# Test de Shapiro-Wilk sobre todos los valores de cada grupo (619 x 3 = 1857 valores,
# por debajo del límite de 5000 que admite shapiro.test)
shapiro_control   <- shapiro.test(as.numeric(as.matrix(datos[, idx_control])))
shapiro_condicion <- shapiro.test(as.numeric(as.matrix(datos[, idx_condicion])))
shapiro_control
shapiro_condicion

# Regla: si los dos grupos son normales (p > 0.05) -> t-test;
#        si alguno no lo es (p <= 0.05)           -> Wilcoxon (no paramétrico)
usar_ttest <- shapiro_control$p.value > 0.05 && shapiro_condicion$p.value > 0.05
if (usar_ttest) {
  cat("\nLos datos son normales: se usa t-test.\n")
} else {
  cat("\nLos datos NO son normales (Shapiro-Wilk, p < 0.05): se usa test de Wilcoxon.\n")
}

# Un test por gen (control frente a condición). tryCatch devuelve NA si un gen
# da error (por ejemplo, valores constantes).
pvalor <- numeric(nrow(datos))
for (i in seq_len(nrow(datos))) {
  x <- as.numeric(datos[i, idx_control])
  y <- as.numeric(datos[i, idx_condicion])
  pvalor[i] <- tryCatch({
    if (usar_ttest) t.test(x, y)$p.value
    else wilcox.test(x, y, exact = FALSE)$p.value
  }, error = function(e) NA)
}

datos$P_value <- pvalor
head(datos)
sum(is.na(datos$P_value))     # genes sin p-valor calculable

# RESPUESTA: el test de Shapiro-Wilk indica que los datos no siguen una distribución
# normal (p < 0.05), por lo que se usó el test no paramétrico de Wilcoxon.


### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Sino te sale ninguno, justifica por qué. En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)

genes_DE <- datos[!is.na(datos$P_value) & datos$P_value < 0.05, ]
n_DE <- nrow(genes_DE)

# Sobreexpresado: media de la condición mayor que la del control; reprimido: menor
sobreexpresados <- sum(genes_DE$Media_Condicion > genes_DE$Media_Control)
reprimidos      <- sum(genes_DE$Media_Condicion < genes_DE$Media_Control)

cat("\nGenes diferencialmente expresados (p < 0.05):", n_DE, "\n")
cat("Sobreexpresados en la condición:", sobreexpresados, "\n")
cat("Reprimidos en la condición:", reprimidos, "\n")

if (n_DE == 0) {
  cat("\nNo sale ningún gen. P-valor mínimo obtenido:",
      signif(min(datos$P_value, na.rm = TRUE), 3), "\n")
  cat("Justificación: con solo 3 réplicas por grupo, el test de Wilcoxon no puede\n",
      "dar p-valores por debajo de ~0.08-0.1, ni aunque los 3 valores del control\n",
      "sean todos mayores (o menores) que los 3 de la condición. Con tan pocas\n",
      "réplicas el test no tiene potencia suficiente para alcanzar p < 0.05.\n")
}

if (n_DE > 0) write.table(genes_DE, "genes_DE.txt", sep = "\t", quote = FALSE)

# RESPUESTA: ningún gen tiene p < 0.05 (p-valor mínimo = 0.0636), así que hay 0 genes
# diferencialmente expresados, 0 sobreexpresados y 0 reprimidos. Con solo 3 réplicas
# por grupo, el test de Wilcoxon no tiene potencia suficiente para llegar a p < 0.05.


### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)

if (!requireNamespace("gplots", quietly = TRUE)) install.packages("gplots")
library(gplots)

# Si hay genes DE se usan esos. Si no hay ninguno (o solo uno), se representan los 30
# genes con menor p-valor, avisando en el título de que no son significativos.
if (n_DE >= 2) {
  genes_heatmap <- rownames(genes_DE)
  titulo <- "Heatmap de genes diferencialmente expresados"
} else {
  genes_heatmap <- rownames(datos)[order(datos$P_value)][1:30]
  titulo <- "30 genes con menor p-valor (ninguno significativo)"
}

# log2(TPM + 1) para que los genes muy expresados no tapen al resto, y se quitan
# los genes con valor constante, que no se pueden escalar.
matriz <- log2(as.matrix(datos[genes_heatmap, c(idx_control, idx_condicion)]) + 1)
matriz <- matriz[apply(matriz, 1, sd) > 0, , drop = FALSE]

heatmap.2(matriz,
          scale = "row", trace = "none",
          col = rev(redblue(256)),
          ColSideColors = c(rep(col_control, 3), rep(col_condicion, 3)),
          labRow = if (nrow(matriz) > 60) FALSE else rownames(matriz),
          cexCol = 0.8, margins = c(9, 9),
          main = titulo)

# RESPUESTA: como no hay genes diferencialmente expresados, el heatmap representa los
# 30 genes con menor p-valor; se observa la separación entre Control y Condición.


### EXTRA: guardar los gráficos como imágenes (no depende del panel Plots de RStudio)

png("histogramas.png", width = 1000, height = 500)
par(mfrow = c(1, 2))
hist(datos$Media_Control, breaks = 30, col = col_control,
     main = "Media del Control", xlab = "TPM medio", ylab = "Número de genes")
hist(datos$Media_Condicion, breaks = 30, col = col_condicion,
     main = "Media de la Condición", xlab = "TPM medio", ylab = "Número de genes")
par(mfrow = c(1, 1))
dev.off()

png("boxplot.png", width = 800, height = 600)
boxplot(log2(datos$Media_Control + 1), log2(datos$Media_Condicion + 1),
        names = c("Control", "Condición"),
        col   = c(col_control, col_condicion),
        main  = "Expresión media por gen: Control frente a Condición",
        ylab  = "log2(TPM medio + 1)")
legend("topleft", legend = c("Control", "Condición"),
       fill = c(col_control, col_condicion), bty = "n")
dev.off()

png("heatmap.png", width = 1000, height = 900)
heatmap.2(matriz,
          scale = "row", trace = "none",
          col = rev(redblue(256)),
          ColSideColors = c(rep(col_control, 3), rep(col_condicion, 3)),
          labRow = if (nrow(matriz) > 60) FALSE else rownames(matriz),
          cexCol = 0.8, margins = c(9, 9),
          main = titulo)
dev.off()

cat("\nImágenes guardadas en:", getwd(), "\n")
cat("Archivos: histogramas.png, boxplot.png, heatmap.png\n")
