### Bibliotecas Importadas --------------------
library(tidyverse)
library(rugarch)
library(fGarch)
library(timeSeries)
library(fpp3)
library(gridExtra)
library(purrr)
library(bizdays)

### Definição de algmas funções----------------


#função para encontrar o melhor modelo garch:
find_best_garch <- function(ativo, grid, df){
  #seleciona o retorno do ativo
  retDiario <- bondsTibble |> dplyr::filter(Ativo==ativo) |> dplyr::select(rDiario) |> pull()
  
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

tidyData <- rawDataFinal |> 
            mutate_if(is.double, as.character) |> 
            pivot_longer(cols = !Dates, names_to = "Ativo")

tidyData <- tidyData |> dplyr::filter(!grepl("N/A", value))
tidyData <- tidyData |> mutate(Dates=as.Date(Dates),
                               value = as.double(value))

#Ativos que estão sendo analisados:
ativosAnalisados <- tidyData$Ativo |> unique() 


# Criando o df de time series
bondsTS <- tidyData |> as_tsibble(key = Ativo, index = Dates, regular = FALSE)


### Análise Exploratória e definição da base para modelo -----------------------
bondsTS |> autoplot()

### estudo da série ---------------------------------------------

#Adicionando o retorno dos yields
bondsTS <- bondsTS |>group_by(Ativo) |> mutate(rDiario = dplyr::lag(value)/value-1)

#retirando o primeiro dia (retorno = NA)
bondsTSClean <- bondsTS |> group_by(Ativo) |> 
                dplyr::filter(is.na(rDiario) == FALSE)

#plotando os graficos de yield e retornos (verficar se  retorno está correto)
p1 <- bondsTSClean |> autoplot(value)
p2 <- bondsTSClean |> autoplot(rDiario) +
          facet_wrap(~Ativo,ncol = 2)



## Plotando os graficos de ACF e PACF dos log-retornos
#ACF
bondsTSClean  |>  ACF(rDiario) |> autoplot()

#PACF
bondsTSClean |> PACF(rDiario) |> autoplot()

#analisando os gráficos de ACF e PACF, percebemos a presença de alguma autocorrelação
#em da série, o que nos faz inferir que a série não pode diretamente ser classificada
#como "estacionária". Para verificarmos se a série se trata de um ruído branco, utilizaremos
#os testes a seguir.

## Teste de Ljung-Box:
# para realizarmos o teste, precisamos ter em mente a definição das hipóteses nulas.

#H0: a série se comporta como um ruído branco
#H1: a série não se comporta como um ruído branco

#criando o df que será utilizado no modelo
bondsTibble <- bondsTSClean |> as_tibble()

#realizando o teste de LJung-Box
whiteNoiseTest <- tibble(Ativo=as.character(), pValor=as.double())
for (j in ativosAnalisados){
  rAtivo <- bondsTibble |> dplyr::filter(Ativo == j) |> dplyr::select(rDiario)
  Ljung <-Box.test(rAtivo, lag=12, fitdf=1, type="Ljung-Box")
  whiteNoiseTest <- whiteNoiseTest |> dplyr::bind_rows(tibble(Ativo=j,pValor=Ljung$p.value)) |> as.data.frame()
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

     
#rodando a bsuca em grid (esta função demora para rodar) -----------------------
bestModels <- ativosAnalisados |> 
  set_names() |> 
  map(find_best_garch, grid=GridSearch) |> 
  bind_rows(.id="Ativo")

#plotando o melhor modelo por ativo 
modeloFinal <-  bestModels %>%
  dplyr::filter(criteria == 'Akaike') |> 
  group_by(Ativo) |> 
  slice(which.min(V1)) |> 
  merge(GridSearch, by.x="id", by.y="id", all.x=FALSE, all.y=FALSE) |> 
  select(-id) |> 
  rename(criteriaValue=V1)
  
modeloFinal


#rodando a simulação do melhor modelo escolhido individualmente  --------------------------------

finaldf <- bondsTibble %>% dplyr::filter(Ativo=='ITAUBZ18') %>% select(Dates,rDiario) |> 
  column_to_rownames(var = 'Dates')


finalSpec <- ugarchspec(
  variance.model = list(model="fGARCH", submodel="GARCH", garchOrder=c(2, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = 'std')

finalFit <- ugarchfit(spec=finalSpec, 
                 data=finaldf, 
                 solver = 'hybrid') 
                 #solver.control=list(tol = 5e-8),
                 #out.sample = 10)

finalsimulation <- ugarchsim(finalFit, 
                        n.sim = 60, 
                        m.sim = 100,
                        rseed=42)

finalsimulation |> plot()
finalsimulation |> summary()
finalsimulation |> fitted()

#Rodando as especificações finais para cada ativo ------------------------------

#df para salvar os resultados
saveResult <- tibble(dataFutura=as.character(),
                     Ativo=as.character(), 
                     rPrevisto=as.double())

for (i in seq(1,nrow(modeloFinal))){

#o modelo roda somente no df com o retorno do ativo:
rAtivo <- bondsTibble |>
  dplyr::filter(Ativo == modeloFinal[i,]$Ativo) |>
  select(Dates, rDiario) |>
  column_to_rownames(var = 'Dates')

#especificação do modelo conforme o grid search: 
specModel <- ugarchspec(
              variance.model = list(model=modeloFinal[i,]$model, submodel=modeloFinal[i,]$submodel, garchOrder=c(modeloFinal[i,]$m, modeloFinal[i,]$n)),
              mean.model = list(armaOrder = c(modeloFinal[i,]$p, modeloFinal[i,]$q), include.mean = TRUE),
              distribution.model = modeloFinal[i,]$dist)

#fitando o modelo especificado 
fitModel <- ugarchfit(
            spec = specModel,
            data = rAtivo,
            solver = 'hybrid')
            #solver.control = list(tol = 5e-8))
            #out.sample = 10)

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

}


## rodando o backtest do modelo especificado----------

#df para salvar os erros do backtest
error <- tibble(Ativo=as.character(), 
                Realized=as.double(), 
                Mu=as.double(),
                Sigma=as.double())

for (i in seq(1,nrow(modeloFinal))){
  
#o modelo roda somente no df com o retorno do ativo:
rAtivo <- bondsTibble |>
  dplyr::filter(Ativo == modeloFinal[i,]$Ativo) |>
  select(Dates, rDiario) |>
  column_to_rownames(var = 'Dates')
  
  #especificação do modelo conforme o grid search: 
specModel <- ugarchspec(
    variance.model = list(model=modeloFinal[i,]$model, submodel=modeloFinal[i,]$submodel, garchOrder=c(modeloFinal[i,]$m, modeloFinal[i,]$n)),
    mean.model = list(armaOrder = c(modeloFinal[i,]$p, modeloFinal[i,]$q), include.mean = TRUE),
    distribution.model = modeloFinal[i,]$dist)


garchroll <- ugarchroll(spec=specModel,
                        data = rAtivo,
                        n.start = 100,
                        refit.window = c("recursive", "moving"),
                        solver = "hybrid")
                        #refit.every = 10)

backtest <- garchroll |> as.data.frame()

dfError <- backtest |> select(Realized, Mu, Sigma)

error <- error |> bind_rows(tibble(Ativo=modeloFinal[i,]$Ativo,
                                        Realized=dfError$Realized, 
                                        Mu=dfError$Mu,
                                        Sigma=dfError$Sigma))

}


errorAjustado <- error |>
          group_by(Ativo) |>
          reframe(e = Realized - Mu,
                  d = e ^ 2 - Sigma ^ 2)

errorFinal <- errorAjustado |> 
              group_by(Ativo) |> 
              reframe(erro=mean(d^2))


saveResult |> write.csv("resultfinal.csv") 


dataAtual <- bondsTibble$Dates |> max()
anoAtual <-  bondsTibble$Dates |> year() |> max()

teste <- saveResult |> group_by(Ativo) |> mutate(diaFuturo=row_number())
teste <- teste |> mutate(Dates=add.bizdays(dataAtual,diaFuturo))
teste <- teste |> rename(rDiario=rPrevisto)

passado <- bondsTibble |> select(Dates, Ativo, rDiario) |> group_by(Ativo)
futuro <-teste |> select(Dates, Ativo, rDiario)
consolidado<-bind_rows(passado,futuro)

consolidado |> dplyr::filter(year(Dates)>= anoAtual-1) |> 
               ggplot(aes(x=Dates, y=rDiario, color=Ativo))+
               geom_line()+
               geom_vline(xintercept=dataAtual)+
               facet_wrap(~Ativo,ncol = 2)
              
              


#----------
g<-7
rAtivoT <- bondsTibble |>
  dplyr::filter(Ativo == modeloFinal[g,]$Ativo) |>
  select(Dates, rDiario) |>
  column_to_rownames(var = 'Dates')

#especificação do modelo conforme o grid search: 
specModelT <- ugarchspec(
  variance.model = list(model=modeloFinal[g,]$model, submodel=modeloFinal[g,]$submodel, garchOrder=c(modeloFinal[g,]$m, modeloFinal[g,]$n)),
  mean.model = list(armaOrder = c(modeloFinal[g,]$p, modeloFinal[g,]$q), include.mean = TRUE),
  distribution.model = modeloFinal[g,]$dist)


garchrollT <- ugarchroll(spec=specModelT,
                        data = rAtivoT,
                        n.start = 200,
                        refit.window = c("recursive","moving"),
                        refit.every = 50)
                        


resume(garchrollT, solver="gosolnp",solver.control=list(tol=1e-4,delta=1e-5))


rAtivoT |> glimpse()

#usando o forecast ---------
g<-2
rAtivoT <- bondsTibble |>
  dplyr::filter(Ativo == modeloFinal[g,]$Ativo) |>
  select(Dates, rDiario) |>
  column_to_rownames(var = 'Dates')

#especificação do modelo conforme o grid search: 
specModelT <- ugarchspec(
  variance.model = list(model=modeloFinal[g,]$model, submodel=modeloFinal[g,]$submodel, garchOrder=c(modeloFinal[g,]$m, modeloFinal[g,]$n)),
  mean.model = list(armaOrder = c(modeloFinal[g,]$p, modeloFinal[g,]$q), include.mean = TRUE),
  distribution.model = modeloFinal[g,]$dist)

fitModelT <- ugarchfit(
  spec = specModelT,
  data = rAtivoT,
  solver = 'hybrid',
  out.sample = 100)
#solver.control = list(tol = 5e-8))
#out.sample = 10)

forecastsT <- ugarchforecast(fitModelT, 
                             n.ahead=60,
                             out.sample = 100,
                             n.roll = 90)
forecastsT |> fitted()
forecastsT |> plot()

simuModelT <- ugarchsim(fitModelT,
                       n.sim = 60,
                       m.sim = 100,
                       rseed = 42)
