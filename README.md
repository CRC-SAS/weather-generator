# Weather Generator (GAM version)

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
devtools::install_version("pbkrtest", version = "0.4.7")
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
