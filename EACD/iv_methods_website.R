## ============================================================
## Variables Instrumentales en la Practica
## Econometria Aplicada y Ciencia de Datos - CIDE Fall 2026
## Carlos Brito
##
## Codigo R extraido de las diapositivas de la presentacion.
## Aplicacion: Retornos a la Educacion (Card, 1995)
## datos: (https://davidcard.berkeley.edu/data_sets.html)
## bajar zip.file y usar: nls.dat
## ============================================================

## ---- Paquetes ----
install.packages(c("readr", "dplyr", "AER", "gmm"))
library(readr)   # read_csv()
library(dplyr)   # select()
library(AER)     # ivreg()
library(gmm)     # gmm()

## ============================================================
## 1. Los Datos
## ============================================================
# Datos de Card (1995) sobre retornos a la educacion, con la
# proximidad a una universidad como instrumento de los anios
# de escolaridad acumulados.
#
# Modelo con cinco regresores mas la constante: educ, exper,
# expersq, black, south. El instrumento excluido es nearc4
# (universidad cercana).

# Column positions and names, taken directly from the codebook
layout <- data.frame(
  start = c(1,7,10,12,15,18,21,24,27,33,35,41,43,56,58,60,62,64,66,68,70,72,
            74,76,78,80,82,84,86,99,112,114,116,118,120,122,124,126,128,134,
            140,146,148,150,152,154,156,159,163,165,167,169),
  end   = c(5,7,10,13,16,19,22,25,31,33,39,41,54,56,58,60,62,64,66,68,70,72,
            74,76,78,80,82,84,97,110,112,114,116,118,120,122,124,126,132,138,
            144,146,148,150,152,154,157,161,163,165,167,169),
  name  = c("id","nearc2","nearc4","nearc4a","nearc4b","ed76","ed66","age76",
            "daded","nodaded","momed","nomomed","weight","momdad14","sinmom14",
            "step14","reg661","reg662","reg663","reg664","reg665","reg666",
            "reg667","reg668","reg669","south66","work76","work78","lwage76",
            "lwage78","famed","black","smsa76r","smsa78r","reg76r","reg78r",
            "reg80r","smsa66r","wage76","wage78","wage80","noint78","noint80",
            "enroll76","enroll78","enroll80","kww","iq","marsta76","marsta78",
            "marsta80","libcrd14"),
  stringsAsFactors = FALSE
)

# Build the widths vector for read.fwf (negative = gap to skip)
widths <- c()
prev_end <- 0
for (i in seq_len(nrow(layout))) {
  gap <- layout$start[i] - prev_end - 1
  if (gap > 0) widths <- c(widths, -gap)
  widths <- c(widths, layout$end[i] - layout$start[i] + 1)
  prev_end <- layout$end[i]
}

data.ingresos <- read.fwf(
  "add the file address here/nls.dat", 
  widths = widths,
  col.names = layout$name,
  strip.white = TRUE
)

data.ingresos <- data.ingresos %>%
  mutate(
    lwage   = lwage76,               # log wage (1976)
    educ    = ed76,                  # years of education (1976)
    exper   = age76 - ed76 - 6,      # potential experience (Mincer formula)
    expersq = exper^2,               # experience squared
    south   = reg76r                 # lived in the South in 1976
    # black and nearc4 already exist with those exact names — no change needed
  )

#final adjustments:
data.ingresos$constant <- 1
str(data.ingresos)
sapply(data.ingresos, class)
data.ingresos <- data.ingresos %>%
  mutate(across(c(lwage, educ, exper, expersq, black, south, nearc4, nearc2), as.numeric))

## ============================================================
## 2. Modelo Exactamente Identificado: Referencia con ivreg (IV)
## ============================================================

iv_ei <- ivreg(lwage ~ educ + exper + expersq + black + south |
                 . - educ + nearc4, data = data.ingresos)

summary(iv_ei)

## ============================================================
## 3. Replicando con Algebra de Matrices (IV)
## ============================================================
#   beta_hat = (Z'X)^(-1) Z'Y

# Keep only rows with no missing values in the variables we need
vars_needed <- c("lwage", "constant", "educ", "exper", "expersq",
                 "black", "south", "nearc4")

data.cc <- data.ingresos %>%
  filter(if_all(all_of(vars_needed), ~ !is.na(.)))

X <- data.matrix(select(data.cc, constant, educ, exper,
                        expersq, black, south))
Y <- data.matrix(select(data.cc, lwage))
Z <- data.matrix(select(data.cc, constant, nearc4, exper,
                        expersq, black, south))

N <- nrow(X)
k <- ncol(X)  # incluyendo la constante

b <- solve(t(Z) %*% X) %*% t(Z) %*% Y
b

## ============================================================
## 4. 2SLS Matriz de Varianzas (Homocedastica)
## ============================================================
# Con la matriz de proyeccion P_Z = Z(Z'Z)^(-1)Z' y
# sigma2_hat = N^(-1) u_hat'u_hat:
#
#   V_hat(beta_hat) = sigma2_hat (X'P_Z X)^(-1) * N/(N-k)
#
# (el factor N/(N-k) es el ajuste por grados de libertad que usa
# R por defecto)

u_hat <- Y - X %*% b
sigma2 <- as.numeric((1 / N) * t(u_hat) %*% u_hat)
P <- Z %*% solve(t(Z) %*% Z) %*% t(Z)

V <- sigma2 * solve(t(X) %*% P %*% X) * (N / (N - k))
sqrt(diag(V))

# Los errores estandar replicados coinciden con los reportados
# por ivreg

## ============================================================
## 5. 2SLS Matriz de Varianzas (Heterocedasticidad)
## ============================================================
#   S_hat = (1/N) Z'DZ,   D = diag[u_hat_i^2]

D <- diag(as.vector((Y - X %*% b)^2))
S_hat <- (1 / N) * t(Z) %*% D %*% Z
Vr <- N * solve(t(X) %*% Z %*% solve(t(Z) %*% Z) %*% t(Z) %*% X) %*%
  (t(X) %*% Z %*% solve(t(Z) %*% Z) %*% S_hat %*% solve(t(Z) %*% Z) %*% t(Z) %*% X) %*%
  solve(t(X) %*% Z %*% solve(t(Z) %*% Z) %*% t(Z) %*% X)

sqrt(diag(Vr))

#small sample correction:
dfc <- N / (N - k)
Vr_hc1 <- Vr * dfc
sqrt(diag(Vr_hc1))

## ============================================================
## 6. Modelo Sobreidentificado (r > q): Dos Instrumentos - ivreg (2SLS)
## ============================================================
# Agregamos un segundo instrumento excluido, nearc2 (universidad
# de 2 anios cercana):
#
#   beta_hat_2SLS = (X'P_Z X)^(-1) X'P_Z Y
iv_si <- ivreg(lwage ~ educ + exper + expersq + black + south |
                 . - educ + nearc4 + nearc2, data = data.ingresos)

summary(iv_si)

#2SLS by hand:
# Keep only rows with no missing values in the variables we need
vars_needed <- c("lwage", "constant", "educ", "exper", "expersq",
                 "black", "south", "nearc4", "nearc2")

data.cc <- data.ingresos %>%
  filter(if_all(all_of(vars_needed), ~ !is.na(.)))

X <- data.matrix(select(data.cc, constant, educ, exper,
                        expersq, black, south))
Y <- data.matrix(select(data.cc, lwage))
Z <- data.matrix(select(data.cc, constant, nearc4, nearc2, exper,
                        expersq, black, south))

N <- nrow(X)
k <- ncol(X)  # incluyendo la constante

P <- Z %*% solve(t(Z) %*% Z) %*% t(Z)
b <- solve(t(X) %*% P %*% X) %*% t(X) %*% P %*% Y
b

## ============================================================
## 7. Estimador Optimo de MGM con el Paquete gmm
## ============================================================
# Procedimiento en dos etapas:
#   1) Primera etapa con W = I para obtener beta_1_hat y construir S_hat.
#   2) Segunda etapa con W = S_hat^(-1) para obtener el estimador
#      optimo beta_GMM,O_hat.
#
# El coeficiente de educ cercano pero no
# identico al 2SLS, reflejando que ambos son validos
# pero difieren en eficiencia bajo heterocedasticidad.

gmm_opt <- gmm(lwage ~ educ + exper + expersq + black + south,
               ~ nearc4 + nearc2 + exper + expersq + black + south,
               vcov = "HAC", wmatrix = "optimal",
               type = "twoStep", data = data.ingresos)

summary(gmm_opt)

#Identity matrix as the weighting one:
gmm_ident <- gmm(lwage ~ educ + exper + expersq + black + south,
                 ~ nearc4 + nearc2 + exper + expersq + black + south,
                 vcov = "HAC", wmatrix = "ident",
                 data = data.ingresos)
summary(gmm_ident)

## ============================================================
## Test de Hausman (Durbin-Wu-Hausman): ¿es necesario usar VI?
##   H0: educ es exogena -> MCO es consistente (y mas eficiente)
##   H1: educ es endogena -> se necesita VI
## ============================================================
# Opcion rapida: ivreg calcula el test de Wu-Hausman automaticamente
# (junto con el de instrumentos debiles y el de Sargan) al pedir
# diagnosticos:
summary(iv_si, diagnostics = TRUE)

# Opcion "a mano" (regresion aumentada / control function):
# 1) Primera etapa: educ en funcion de los exogenos + instrumentos excluidos
first_stage <- lm(educ ~ exper + expersq + black + south + nearc4 + nearc2,
                  data = data.ingresos)
v_hat <- residuals(first_stage)

# 2) Se agrega v_hat como regresor adicional a la ecuacion estructural
aux_reg <- lm(lwage ~ educ + exper + expersq + black + south + v_hat,
              data = data.ingresos)

# 3) El t-test sobre v_hat ES el test de Hausman: si su coeficiente
#    es significativo, se rechaza H0
library(lmtest); library(sandwich)
coeftest(aux_reg, vcov = vcovHC(aux_reg, type = "HC1"))   # version robusta
summary(aux_reg)                                           # version clasica (no robusta)

## ============================================================
## Test de Hansen J (sobreidentificacion, robusto a heterocedasticidad)
##   H0: los instrumentos sobreidentificantes son validos (exogenos)
##   p-valor bajo -> se rechaza H0: al menos un instrumento
##   (nearc4 o nearc2) probablemente viola la exclusion restriction
## ============================================================
# El "Sargan" que reporta ivreg asume homocedasticidad. El analogo
# robusto (Hansen J) se obtiene estimando el mismo modelo por GMM
# con matriz de pesos optima:
library(gmm)
gmm_si <- gmm(lwage ~ educ + exper + expersq + black + south,
              ~ exper + expersq + black + south + nearc4 + nearc2,
              vcov = "HAC",        
              wmatrix = "optimal",
              type = "twoStep",
              data = data.ingresos)

summary(gmm_si)
# El bloque "J-Test" al final del summary es el estadistico de Hansen.
# Grados de libertad = (# instrumentos) - (# parametros) = 7 - 6 = 1,
# que corresponde a tener 2 instrumentos excluidos (nearc4, nearc2)
# para 1 variable endogena (educ).
# =============================================================