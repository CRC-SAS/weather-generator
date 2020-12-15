# Weather Generator GAMWGEN

## Overview

Este paquete contiene la implementación de un generador estocástico diario y multisitio de series meteorológicas sintéticas. La generación de series sintéticas es un insumo básico para el análisis probabilista de las sequías y sus impactos en el sector agrícola. El generador desarrollado es muy flexible y capaz de generar secuencias de valores diarios de precipitación y temperaturas máxima y mínima. A partir de estas últimas se pueden derivar otras variables como la radiación solar y la evapotranspiración. 

El generador tiene dos variantes, una permite generar series para localidades puntuales, generador local, y otra que permite la generación en grillas regulares o en localidades arbitrarias, generador espacial. Esto la de la posibilidad al usuario de generar productos a medida para distintas aplicaciones. Además, el generador es capaz de utilizar variables auxiliares que producen simulaciones condicionadas para su uso en conjunto con modelos de cambio climático o de pronósticos estacionales para evaluar impactos en el largo o corto plazo, respectivamente. 

Esta flexibilidad se sustenta en el uso modelos generalizados aditivos (GAM) que permiten modelar con mucha precisión el comportamiento de las variables meteorológicas y capturar las propiedades estadísticas de los datos observados. La modelación de las variables meteorológicas se divide en dos: por un lado, la ocurrencia y monto de precipitación y por otro, las temperaturas máxima y mínima. La ocurrencia de precipitación se modela a través de un modelo probit mientras que al monto se lo hace a través de una distribución aleatoria Gamma. Para la temperatura se utiliza un modelo autorregresivo condicionado por la ocurrencia de lluvia. A su vez, estos modelos pueden ser espacialmente correlacionados con campos aleatorios Gaussianos que contemplan la variabilidad espacial y temporal regional. 

Además del generador, se incluyen diagnósticos estadísticos y gráficos para verificar la bondad de ajuste estadística de los GAM y validar que las series sintéticas sean consistentes con los registros históricos. Los diagnósticos son exhaustivos e incluyen todas las propiedades de las series que podrían afectar el desempeño de las mismas durante el análisis probabilista. 


## Installation
The *Weather Generator* is under active development and the only way to install it 
now is through the installation of the development version from GitHub. The easiest 
way to do this is via the `install_github()` function from package *devtools*. Make 
sure you have *devtools* installed, then run this command to install the package:

```r
devtools::install_github("CRC-SAS/weather-generator", ref = 'gamwgen')
```

If you are using a version of R prior to 3.6, before install the *Weather Generator* 
you must install a version of the R package called *pbkrtest* that is compatible with 
your R installation. For example, if your R version is 3.5, you must install *pbkrtest*
before install the *Weather Generator* executing this command:

```r
devtools::install_version("pbkrtest", version = "0.4-7")
```

Also note that, depending on your OS, R packages that will be installed as dependencies 
of the *Weather Generator* may require the installation of additional OS packages. 
Generally the installation process indicates which OS packages are required. For example, 
in a Debian 10 OS, you must run this commands to install all required OS packages:

```shell
sudo apt install libssl-dev
sudo apt install libcurl4-openssl-dev
sudo apt install libxml2-dev
sudo apt install libnlopt-dev
sudo apt install libnetcdf-dev
sudo apt install libudunits2-dev
sudo apt install libgdal-dev
```
