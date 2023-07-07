### Bibliotecas Importadas --------------------
library(tidyverse)
library(rugarch)
library(fGarch)
library(fpp3)
library(gridExtra)
library(purrr)

### Definição de algmas funções----------------


#função para encontrar o melhor modelo garch:
find_best_garch <- function(ativo, grid, df){
  #seleciona o retorno do ativo
  retDiario <- yieldsTibble |> filter(Ativo==ativo) |> select(rDiario) |> pull()
  
  #deixa visível ao usuário o ativo que está sendo escolhido
  usethis::ui_info("Adjusting models for {ativo}...")
  #mostra o progresso
  progressr::with_progress({
    prog <- progressr::progressor(nrow(grid))
    models <- grid %>%
      group_split(id) %>%
      purrr::map(GModels, series=retDiario, prog)
  })
  safe_info <- purrr::possibly(infocriteria, tibble::tibble())
  
  #pega as informações do modelo
  suppressWarnings({
    info <- purrr::map(models, safe_info) %>%
      purrr::map(tibble::as_tibble, rownames = "criteria") %>%
      dplyr::bind_rows(.id = "id")
  })
  return(info)
  
  #agrega as informações do modelo
  best <- info |> 
    inner_join(grid, "id") |> 
    pivot_wider(names_from = criteria, values_from = V1) |> 
    janitor::clean_names() |> 
    arrange(akaike)
  
  
  #seleciona os melhores paramestros
  usethis::ui_info(c(
    "Best model:",
    "p <- {best$p[1]}",
    "q <- {best$q[1]}",
    "m <- {best$m[1]}",
    "n <- {best$n[1]}"
  ))
  
  best
  
}

#função que roda o modelo GARCH
GModels <- function(parms, series, prog = NULL){
  
  if (!is.null(prog)) prog()
  
  #configurando o modelo iGARCH
  if(parms$model=="iGARCH"){
    garch_model = ugarchspec(
      variance.model = list(model=parms$model, garchOrder=c(parms$m, parms$n)),
      mean.model = list(armaOrder = c(parms$p, parms$q), include.mean = TRUE),
      distribution.model = parms$dist)
  }
  
  #configurando os outros modelos
  else{
    garch_model = ugarchspec(
      variance.model = list(model=parms$model, submodel=parms$submodel, garchOrder=c(parms$m, parms$n)),
      mean.model = list(armaOrder = c(parms$p, parms$q), include.mean = TRUE),
      distribution.model = parms$dist
    )
  }
  
  # ocultando os avisos quando o modelo não convergir
  suppressWarnings({fit <- ugarchfit(spec=garch_model, data=series, solver='solnp', solver.control=list(tol = 5e-8))})
  
  
  # modelo ajustado
  fit
}


### Importando os dados --------------------------------------------------------
rawDataFinal <- readxl::read_xlsx("forecastbondsfinal.xlsx", 
                                   sheet = "import",
                                   skip = 3)


rawData <- readxl::read_xlsx("forecastbonds_treated.xlsx", sheet = "bonds")

#verificando se todos os ativos têm a mesma qtd de dias:
rawData |> glimpse()
rawData |> group_by(rawData$Ativo) |> summarise(min = min(Data),
                                                max = max(Data),
                                                observacoes = n())
#Ativos que estão sendo analisados:
Ativos <- rawData$Ativo |> unique() 


# Criando o df Data - Ativo - Yield:
bondsYield <- rawData |> select(Data, Ativo, YAS_BOND_YLD) |> 
                         rename(bondYield = YAS_BOND_YLD)
bondsYield <- bondsYield |> as_tsibble(key = Ativo, index = Data, regular = FALSE)


### Análise Exploratória e definição da base para modelo -----------------------
bondsYield |> autoplot(bondYield)

#Percebemos que há uma descontinuidade nos dados dos bonds de hidrovias
yields <- bondsYield |> filter(Data > "2018-01-24") 


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
yieldsClean  |>  ACF(rDiario) |> autoplot()

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

#criando o df que será utilizado no modelo
yieldsTibble <- yieldsClean |> as_tibble()
yieldsSplitted <- yieldsTibble |> split(yieldsTibble$Ativo)

#realizando o teste de LJung-Box
whiteNoiseTest <- tibble(Ativo=as.character(), pValor=as.double())
for (j in Ativos){
  rAtivo <- yieldsTibble |> filter(Ativo == j) |> select(rDiario)
  Ljung <-Box.test(rAtivo, lag=12, fitdf=1, type="Ljung-Box")
  whiteNoiseTest <- whiteNoiseTest |> bind_rows(tibble(Ativo=j,pValor=Ljung$p.value)) |> as.data.frame()
}

whiteNoiseTest

#como o p-valor é muito baixo, rejeitamos a H0 e consideramos que a série não se
#comporta como um ruído branco.

### Modelos de Volatilidade ----------------------------------------------------

#Para testar os diferentes modelos de volatilidade, iremos utilizar a seguinte tabela----
mtr <- crossing(m=c(0:2), n=c(0:2), p=c(0:3), q=c(0:3), dist=c("std"))
GridSearch <- bind_rows(cbind(tibble(model="fGARCH", submodel="TGARCH"), mtr),
                        cbind(tibble(model="fGARCH", submodel="GARCH"), mtr)) %>%
  tibble::rownames_to_column("id")
#com base nos diferentes paramentros da tabela, rodaremos um grid search para escolher
#o modelo que melhor perfroma através do Critério de Akaike

#funções para rodar a busca em grid --------------------------------------------


          
     
#rodando a bsuca em grid (esta função demora para rodar) -----------------------
bestModels <- Ativos |> 
  set_names() |> 
  map(find_best_garch, grid=GridSearch) |> 
  bind_rows(.id="Ativo")

#plotando o melhor modelo por ativo 
modeloFinal <-  bestModels %>%
  filter(criteria == 'Akaike') |> 
  group_by(Ativo) |> 
  slice(which.min(V1)) |> 
  merge(GridSearch, by.x="id", by.y="id", all.x=FALSE, all.y=FALSE) |> 
  select(-id) |> 
  rename(criteriaValue=V1)
  
modeloFinal


#rodando a simulação do melhor modelo escolhido --------------------------------

finaldf <- yieldsTibble %>% filter(Ativo=='HIDRVS25') %>% select(Data,rDiario) |> 
  column_to_rownames(var = 'Data')


finalSpec <- ugarchspec(
  variance.model = list(model="fGARCH", submodel="GARCH", garchOrder=c(2, 1)),
  mean.model = list(armaOrder = c(2, 1), include.mean = TRUE),
  distribution.model = 'std')

finalFit <- ugarchfit(spec=finalSpec, 
                 data=finaldf, 
                 solver='solnp', 
                 solver.control=list(tol = 5e-8),
                 out.sample = 10)

finalsimulation <- ugarchsim(finalFit, 
                        n.sim = 60, 
                        m.sim = 100,
                        rseed=42)

finalsimulation |> plot()
finalsimulation |> summary()
finalsimulation |> fitted()

#-----------------------------------------------------------------------------


#df para salvar os dados do modelo
saveResult <- tibble(dataFutura=as.character(),
                     Ativo=as.character(), 
                     rPrevisto=as.double())

#df para salvar os erros
error <- tibble(Ativo=as.character(), 
                Realized=as.double(), 
                Mu=as.double(),
                Sigma=as.double())

for (i in c(1,nrow(modeloFinal))){

#o modelo roda somente no df com o retorno do ativo:
rAtivo <- yieldsTibble |>
  filter(Ativo == modeloFinal[i,]$Ativo) |>
  select(Data, rDiario) |>
  column_to_rownames(var = 'Data')

#especificação do modelo conforme o grid search: 
specModel <- ugarchspec(
              variance.model = list(model=modeloFinal[i,]$model, submodel=modeloFinal[i,]$submodel, garchOrder=c(modeloFinal[i,]$m, modeloFinal[i,]$n)),
              mean.model = list(armaOrder = c(modeloFinal[i,]$p, modeloFinal[i,]$q), include.mean = TRUE),
              distribution.model = modeloFinal[i,]$dist)

#fitando o modelo especificado 
fitModel <- ugarchfit(
            spec = specModel,
            data = rAtivo,
            solver = 'solnp',
            solver.control = list(tol = 5e-8),
            out.sample = 10)

#simulando diversos cenários 
simuModel <- ugarchsim(fitModel,
                       n.sim = 60,
                       m.sim = 100,
                       rseed = 42)

#o retorno previsto dos ativos é a média de todos os cenários simulados
resultSimu <- simuModel |> fitted() |> as.data.frame() 
resultSimu <- resultSimu |> mutate(rMedio = rowMeans(resultSimu)) |> rownames_to_column("dataFutura")

resultModel <- resultSimu |> select(dataFutura, rMedio)

saveResult <- saveResult |> bind_rows(tibble(dataFutura=resultModel$dataFutura,Ativo=modeloFinal[i,]$Ativo, rPrevisto=resultModel$rMedio))

## rodando o backtest do modelo especificado

garchroll <- ugarchroll(spec=specModel, 
                        data = rAtivo, 
                        n.start = 1000,
                        refit.window = c("recursive", "moving"),
                        refit.every = 100)

backtest <- garchroll |> as.data.frame()

dfError <- backtest |> select(Realized, Mu, Sigma)

error <- error |> bind_rows(tibble(Ativo=modeloFinal[i,]$Ativo,
                                        Realized=dfError$Realized, 
                                        Mu=dfError$Mu,
                                        Sigma=dfError$Sigma))

#erro de predição em relaçãoa a media
#e  <- backtest$Realized - backtest$Mu

#erro de predição para a variancia
#d  <- e^2 - backtest$Sigma^2
#erro de predição médio é dado por mean(d^2)



}

error1 <- error |>
          group_by(Ativo) |>
          reframe(e = Realized - Mu,
                  d = e ^ 2 - Sigma ^ 2)

errorFinal <- error1 |> 
              group_by(Ativo) |> 
              reframe(erro=mean(d^2))


saveResult |> write.csv("resultfinal.csv") 



