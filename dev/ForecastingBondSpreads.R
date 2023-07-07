### Bibliotecas Importadas --------------------
library(tidyverse)
library(rugarch)
library(fGarch)
library(fpp3)
library(gridExtra)

### Importando os dados --------------------------------------------------------
rawData <- readxl::read_xlsx("forecastbonds.xlsx", sheet = "bonds")

#verificando se todos os ativos têm a mesma qtd de dias
rawData |> glimpse()
rawData |> group_by(rawData$Ativo) |> summarise(min = min(Data),
                                                max = max(Data),
                                                observacoes = n())


# Iremos prever o yield dos bonds
bondsYield <- rawData |> select(Data, Ativo, YAS_BOND_YLD) |> 
                         rename(bondYield = YAS_BOND_YLD)
bondsYield <- bondsYield |> as_tsibble(key = Ativo, index = Data, regular = FALSE)


### Análise Exploratória e definição da base para modelo -----------------------
bondsYield |> autoplot(bondYield)

#Percebemos que há uma descontinuidade nos dados dos bonds de hidrovias
#Iremos considersar, neste primeiro momento, o intervalo que faz sentido.

bondsYield |> filter(Ativo %in% c('CMIGBZ24','HIDRVS25')) |> 
              mutate(bondYield_0 = lag(bondYield),
                     rDiario = bondYield/bondYield_0) |> 
              arrange(rDiario)
# a data em que acontece a descontinuidade é 2018-01-24. Utilizaremos dados
#posteriores a esta data

yields <- bondsYield |> filter(Ativo %in% c('CMIGBZ24','HIDRVS25'), Data > "2018-01-24") 


### estudo da série ---------------------------------------------

#Adicionando o retorno dos yields
yields <- yields |> mutate(rDiario = bondYield/lag(bondYield)-1)
#retirando o primeiro dia (retorno = NA)
yieldsClean <- yields |> group_by(Ativo) |> 
               filter(is.na(rDiario) == FALSE)

#plotando os graficos de yield e retornos
p1 <- yieldsClean |> autoplot(bondYield)
p2 <- yieldsClean |> autoplot(rDiario) +
          facet_wrap(~Ativo, ncol = 1)

grid.arrange(p1, p2)

## Plotando os graficos de ACF e PACF dos log-retornos
#ACF
yieldsClean |> ACF(rDiario) |> autoplot()

#PACF
yieldsClean |> PACF(rDiario) |> autoplot()

#analisando os gráficos de ACF e PACF, percebemos a presença de alguma autocorrelação
#em da série, o que nos faz inferir que a série não pode diretamente ser classificada
#como "estacionária". Para verificarmos se a série se trata de um ruído branco, utilizaremos
#os testes a seguir.

## Teste de Ljung-Box:
# para realizarmos o teste, precisamos ter em mente a definição das hipóteses nulas.

#H0: a série se comporta como um ruído branco
#H1: a série não se comporta como um ruído branco

#extraindo somente os retornos da série
yieldsTibble <- yieldsClean |> as_tibble()

rCMIGBZ24 <- yieldsTibble |> filter(Ativo == "CMIGBZ24") |> select(rDiario)
rHIDRVS25 <- yieldsTibble |> as_tibble() |> filter(Ativo == "HIDRVS25") |> select(rDiario)

#exibindo os resultados do teste
Box.test(rCMIGBZ24, lag=12, fitdf=1, type="Ljung-Box")
Box.test(rHIDRVS25, lag=12, fitdf=1, type="Ljung-Box")

#como o p-valor é muito baixo, rejeitamos a H0 e consideramos que a série não se
#comporta como um ruído branco.

### Modelos de Volatilidade
#Para testar os diferentes modelos de volatilidade, iremos utilizar a seguinte tabela
mtr <- crossing(m=c(0:2), n=c(0:2), p=c(0:3), q=c(0:3), dist=c("std"))

#com base nos diferentes paramentros da tabela, rodaremos um grid search para escolher
#o modelo que melhor perfroma através do Critério de Akaike

#modelo individual

garch_model <- ugarchspec(
  variance.model = list(model="fGARCH", submodel="TGARCH", garchOrder=c(2, 2)),
  mean.model = list(armaOrder = c(2, 3), include.mean = TRUE),
  distribution.model = 'std')

rDiario2 <- yieldsTibble %>% filter(Ativo=='CMIGBZ24') %>% select(rDiario) %>% pull()
rDiario2 <- rDiario2^2
rDiarioCEMIG <- yieldsTibble %>% filter(Ativo=='CMIGBZ24') %>% select(Data,rDiario) 
rDiarioCEMIG <- rDiarioCEMIG |> column_to_rownames(var = 'Data')
rDiarioCEMIG |> glimpse()

fit <- ugarchfit(spec=garch_model, 
                 data=rDiarioCEMIG, 
                 solver='solnp', 
                 solver.control=list(tol = 5e-8),
                 out.sample = 10)

resi <- residuals(fit, standardize=T)
ts.plot(resi)
acf(resi)
Box.test(resi,lag=12,type="Ljung")

fit |> infocriteria()

#### backtest

garchroll <- ugarchroll(spec=garch_model, 
             data = rDiarioCEMIG, 
             n.start = 1000,
             refit.window = c("recursive", "moving"),
             refit.every = 100)

backtest <- garchroll |> as.data.frame()

# Prediction error for the mean
e  <- backtest$Realized - backtest$Mu  
# Prediction error for the variance
d  <- e^2 - backtest$Sigma^2 
# Mean of prediction error
mean(d^2)




#### predict
pred <- ugarchboot(fit, method = c("Partial"),
                n.ahead = 30, 
                n.bootpred = 500)
show(pred)
plot(pred)


simulation <- ugarchsim(fit, 
                        n.sim = 60, 
                        m.sim = 100)
simulation |> plot()
resultSimula <- simulation |> fitted() |> as.data.frame()
resultSimula |> glimpse()
resultSimula <- resultSimula |> mutate(rMedio = rowMeans(resultSimula))

simulacaoCEMIG <- resultSimula |> select(rMedio)
simulacaoCEMIG |> glimpse()
simulacaoCEMIG |> write.csv(file = 'result.csv')

simulacaoCEMIG |> head()

resultado <- yieldsTibble |> 
             filter(Ativo == 'CMIGBZ24') |> 
             select(Data, bondYield, rDiario)


resultadof <- resultado |> tail(60)
resultadof |> mutate(lagyield = lag(bondYield, -59)) |> tail()

resultadof |> head()











#não funcionou bem
forecast <-ugarchforecast(fit,
           n.ahead = 100,
           n.roll = 9)

forecast |> plot()









#RODANDO BUSCA EM GRID
gridSearch <- bind_rows(cbind(tibble(model="fGARCH", submodel="TGARCH"), mtr),
              cbind(tibble(model="fGARCH", submodel="GARCH"), mtr)) %>%
              tibble::rownames_to_column("id")

#rodando para encontrar o melhor GARCH

Ativos <- c('CMIGBZ24', 'HIDRVS25')
GModels <- function(parms, series, prog = NULL) {
  
  if (!is.null(prog)) prog()
  
  if(parms$model=="iGARCH"){
    
    ## Configuring iGARCH model
    garch_model = ugarchspec(
      variance.model = list(model=parms$model, garchOrder=c(parms$m, parms$n)),
      mean.model = list(armaOrder = c(parms$p, parms$q), include.mean = TRUE),
      distribution.model = parms$dist)
  }
  
  else {
    
    ## Configuring other models
    garch_model = ugarchspec(
      variance.model = list(model=parms$model, submodel=parms$submodel, garchOrder=c(parms$m, parms$n)),
      mean.model = list(armaOrder = c(parms$p, parms$q), include.mean = TRUE),
      distribution.model = parms$dist
    )
  }
  
  
  ## Suppressing warning when the model doesnt converge
  suppressWarnings({fit <- ugarchfit(spec=garch_model, data=series, solver='solnp', solver.control=list(tol = 5e-8))})
  
  ## Return the fit model
  fit
}
find_best_garch <- function(asset, grid, df) {
  
  ## Selects the asset and it squared log return
  ret2 <- yieldsTibble %>% filter(Ativo==asset) %>% select(rDiario) %>% pull()
  
  ## Print which asset is being adjusted
  usethis::ui_info("Adjusting models for {asset}...")
  
  ## Show progress
  progressr::with_progress({
    prog <- progressr::progressor(nrow(grid))
    models <- grid %>%
      group_split(id) %>%
      purrr::map(GModels, series=ret2, prog)
  })
  
  safe_info <- purrr::possibly(infocriteria, tibble::tibble())
  
  ## Get model information
  suppressWarnings({
    info <- purrr::map(models, safe_info) %>%
      purrr::map(tibble::as_tibble, rownames = "criteria") %>%
      dplyr::bind_rows(.id = "id")
  })
  
  return(info)
  
  ## Collecting models info
  best <- info %>%
    dplyr::inner_join(grid, "id") %>%
    tidyr::pivot_wider(names_from = criteria, values_from = V1) %>%
    janitor::clean_names() %>%
    dplyr::arrange(akaike)
  
  ## Selecting the best parameters
  usethis::ui_info(c(
    "Best model:",
    "p <- {best$p[1]}",
    "q <- {best$q[1]}",
    "m <- {best$m[1]}",
    "n <- {best$n[1]}"
  ))
  
  best
}

Ativos |> set_names() |> 
  map(find_best_garch, grid=gridSearch) |> 
  bind_rows(.id="ativo")

fit.ga





          
      
    





