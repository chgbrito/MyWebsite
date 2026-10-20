# ============================================================
# Datos en panel en R
# Econometria Aplicada y Ciencia de Datos - CIDE Fall 2026
# Carlos Brito
#
#   - comportamiento_wide.csv <- panel NLSY (comportamiento
#     antisocial / autoestima / pobreza, 1990-1994), incluido
#     como el dataset `nlsy` en el paquete `panelr`.
#     Fuente: National Longitudinal Survey of Youth (US Dept.
#     of Labor); distribuido por Paul Allison / Statistical
#     Horizons.
#
#   - mlbook1.csv <- datos de la prueba de lenguaje de Snijders
#     & Bosker (1999), 2287 alumnos en 131 escuelas, incluidos
#     como el dataset `bdf` en el paquete `nlme` (viene con
#     R base).
# ============================================================

## 0. Paquetes -------------------------------------------------

pkgs <- c("tidyverse", "plm", "clubSandwich", "modelsummary",
          "panelr", "lmtest", "sandwich")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

library(tidyverse)
library(plm)
library(clubSandwich)
library(modelsummary)
library(panelr)   # trae el dataset `nlsy` (base para comportamiento_wide.csv)
library(nlme)     # viene con R base; trae el dataset `bdf` (base para mlbook1.csv)
dir.create("files", showWarnings = FALSE)

## 0.1 Comportamiento_wide.csv ----------------------------

data(nlsy, package = "panelr")
nlsy_wide <- as.data.frame(nlsy)
if (!"id" %in% names(nlsy_wide)) nlsy_wide$id <- seq_len(nrow(nlsy_wide))

comportamiento_wide <- nlsy_wide %>%
  dplyr::select(id, gender, black, hispanic, married, momage, momwork,
                childage, anti90, anti92, anti94, self90, self92, self94,
                pov90, pov92, pov94)

readr::write_csv(comportamiento_wide, "files/comportamiento_wide.csv")

## 0.2 mlbook1.csv -----------------------------------------

data(bdf, package = "nlme")

mlbook1 <- bdf %>%
  transmute(
    schoolnr = as.integer(as.character(schoolNR)),
    pupilnr  = as.integer(as.character(pupilNR)),
    iq_verb  = IQ.verb,
    sex      = as.integer(sex) - 1,        # 0/1
    minority = as.integer(Minority) - 1,   # 0/1
    repeatgr = as.integer(repeatgr) - 1,   # 0/1/2
    langpost = langPOST
  )

print(colSums(is.na(mlbook1)))
readr::write_csv(mlbook1, "files/mlbook1.csv")

# ============================================================
# 1. Datos en panel 
# ============================================================
# Trabajamos con los datos comportamiento_wide.csv: 
data.comp <- read_csv("files/comportamiento_wide.csv",
                       locale = locale(encoding = "latin1"))
colnames(data.comp)

# La base esta en formato wide; la pasamos a formato long:
data.comp <- data.comp %>%
  pivot_longer(c(anti90:anti94, self90:self94, pov90:pov94),
               names_to = c("measure", "year"),
               names_pattern = "(.*)(..)") %>%
  pivot_wider(names_from = measure,
              values_from = value)

colnames(data.comp)

# ============================================================
# 2. Estimadores para modelos de datos en panel
# ============================================================

## 2.1 Estimador de MCO --------------------------

m.mco <- plm(anti ~ self + pov,
             data = data.comp,
             model = "pooling",
             index = c("id", "year"))

modelsummary(list("MCO" = m.mco), stars = TRUE)

## 2.2 Estimador de MCO y errores agrupados ------

vcov <- list(NULL, clubSandwich::vcovCR(m.mco, type = 'CR1',
                                         cluster = data.comp$id))

modelsummary(
  list("MCO" = m.mco, "MCO, errores agrupados" = m.mco),
  vcov = vcov,
  stars = TRUE
)

## 2.3 Estimador de efectos fijos -----------------

m.fe <- plm(anti ~ self + pov,
            data = data.comp,
            model = "within",
            index = c("id", "year"))

modelsummary(
  list("MCO" = m.mco, "MCO, err. agrup." = m.mco, "Efectos fijos" = m.fe),
  vcov = list(NULL,
              clubSandwich::vcovCR(m.mco, type = 'CR1', cluster = data.comp$id),
              NULL),
  stars = TRUE
)

## 2.4 Estimador de efectos aleatorios ------------

m.re <- plm(anti ~ self + pov,
            data = data.comp,
            model = "random",
            index = c("id", "year"))

modelsummary(
  list("MCO" = m.mco, "MCO, err. agrup." = m.mco,
       "Ef. fijos" = m.fe, "Ef. aleatorios" = m.re),
  vcov = list(NULL,
              clubSandwich::vcovCR(m.mco, type = 'CR1', cluster = data.comp$id),
              NULL, NULL),
  stars = TRUE
)

## 2.5 Inclusion de caracteristicas invariantes en el tiempo

# Con efectos fijos no se puede estimar el efecto de genero (es
# invariante en el tiempo), asi que ese modelo queda descartado.
# La alternativa es MCO pooled o un modelo de efectos aleatorios.

m.sex.f <- plm(anti ~ self + pov + gender,
               data = data.comp, model = "within",
               index = c("id", "year"))

m.sex.a <- plm(anti ~ self + pov + gender,
               data = data.comp, model = "random",
               index = c("id", "year"))

modelsummary(
  list("Ef. fijos" = m.sex.f, "Ef. aleatorios (genero)" = m.sex.a),
  stars = TRUE
)

## 2.6 Prueba de Hausman --------------------------

phtest(m.fe, m.re)

# Se rechaza H0 de que los coeficientes estimados son iguales.
# Hay evidencia de que se prefiere un modelo de efectos fijos

# ============================================================
# 3. Equivalencia de estimadores
# ============================================================

## 3.1 Efectos fijos == MCO con dummies de individuos

m.fe <- plm(anti ~ self + pov,
            data = data.comp, model = "within",
            index = c("id", "year"))

m.dummy <- lm(anti ~ self + pov + factor(id),
              data = data.comp)

modelsummary(list("MCO" = m.fe, "MCO con dummies" = m.dummy), stars = TRUE)

## 3.2 Caracteristicas invariantes en el tiempo no identificables

summary(plm(anti ~ self + pov + black,
            data = data.comp, model = "within",
            index = c("id", "year")))
# (black omitido por ser invariante en el tiempo)

## 3.3 Efectos fijos == MCO en diferencias con respecto a la media

# Conservamos dos periodos consecutivos y observaciones completas.

data.comp.sub <- data.comp %>%
  dplyr::select(id, year, anti, self, pov) %>%
  filter(year == 90 | year == 92)

data.comp.sub <- data.comp.sub[complete.cases(data.comp.sub), ]

data.comp.sub <- data.comp.sub %>%
  group_by(id) %>%
  mutate(m.anti = mean(anti), m.self = mean(self),
         m.pov = mean(pov)) %>%
  mutate(dm.anti = anti - m.anti, dm.self = self - m.self,
         dm.pov = pov - m.pov)

m.fe.sub <- plm(anti ~ self + pov, data = data.comp.sub,
                model = "within", index = c("id", "year"))

m.demean <- lm(dm.anti ~ dm.self + dm.pov, data.comp.sub)

modelsummary(
  list("Efectos fijos" = m.fe.sub, "MCO demean" = m.demean),
  stars = TRUE
)

## 3.4 Efectos fijos == MCO en primeras diferencias

data.comp.sub <- data.comp.sub %>%
  group_by(id) %>%
  mutate(d.anti = anti - dplyr::lag(anti, order_by = year),
         d.self = self - dplyr::lag(self, order_by = year),
         d.pov  = pov  - dplyr::lag(pov,  order_by = year)) %>%
  ungroup()

m.difs <- lm(d.anti ~ -1 + d.self + d.pov, data = data.comp.sub)

modelsummary(
  list("Ef. fijos" = m.fe.sub, "MCO demean" = m.demean,
       "MCO 1ras dif." = m.difs),
  stars = TRUE
)

# ============================================================
# 4. Errores estandar
# ============================================================

## 4.1 Estructura de los datos --------------------
# mlbook1.csv: 2287 estudiantes en 131 escuelas:
data.examen <- read_csv("files/mlbook1.csv",
                         locale = locale(encoding = "latin1"))
colnames(data.examen)

## 4.2 Errores clasicos ---------------------------

m.mco <- lm(langpost ~ iq_verb + sex + minority + repeatgr,
            data = data.examen)

modelsummary(list("MCO" = m.mco), stars = TRUE)

## 4.3 Errores robustos ---------------------------

vcov <- list(NULL, "HC1")

modelsummary(
  list("MCO" = m.mco, "MCO, err. robustos" = m.mco),
  vcov = vcov,
  stars = TRUE
)

## 4.4 Estimador de efectos fijos de escuela ------

m.mco.ef <- lm(langpost ~ iq_verb + sex + minority + repeatgr +
                 factor(schoolnr), data = data.examen)

modelsummary(
  list("MCO" = m.mco, "MCO, err. rob." = m.mco, "Ef. fijos" = m.mco.ef),
  vcov = list(NULL, "HC1", "HC1"),
  stars = TRUE
)

## 4.5 Errores agrupados a nivel escuela ----------

clubSandwich::vcovCR(m.mco, type = 'CR1',
                      cluster = data.examen$schoolnr)

modelsummary(
  list("MCO" = m.mco, "MCO, err. rob." = m.mco, "Ef. fijos" = m.mco.ef,
       "Err. agrupados" = m.mco),
  vcov = list(NULL, "HC1", "HC1",
              clubSandwich::vcovCR(m.mco, type = 'CR1',
                                   cluster = data.examen$schoolnr)),
  stars = TRUE
)

## 4.6 Estimador de efectos fijos y errores agrupados
modelsummary(
  list("MCO" = m.mco, "Err. rob." = m.mco, "Ef. fijos" = m.mco.ef,
       "Err. agrup." = m.mco, "Ef. fijos + agrup." = m.mco.ef),
  vcov = list(NULL, "HC1", "HC1",
              clubSandwich::vcovCR(m.mco, type = 'CR1',
                                   cluster = data.examen$schoolnr),
              clubSandwich::vcovCR(m.mco.ef, type = 'CR1',
                                   cluster = data.examen$schoolnr)),
  stars = TRUE
)
