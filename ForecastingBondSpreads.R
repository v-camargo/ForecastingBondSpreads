### Bibliotecas Importadas --------------------
library(tidyverse)

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
yields |> autoplot(bondYield)              




