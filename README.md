appSAT: Sea Level Anomaly in the Indian Ocean
================

appSat: using Sea Level Anomaly data to anticipate sea level drops in the Reunion Island reef lagoons
-----------------------------------------------------------------------------------------------------

-   **AIM** Display movies of the Sea Level anomaly for the World, Indian Ocean, Reunion Island.

-   **DATA**:
    -   Sea Level anomaly from satellite data: We used [Delayed Time Sea Level Anomaly data stem from Marine Copernicus data](http://marine.copernicus.eu/services-portfolio/access-to-products/?option=com_csw&view=details&product_id=SEALEVEL_GLO_SLA_MAP_L4_REP_OBSERVATIONS_008_027) over 1993/01/01 - 2016/01/11 and [Near-Real Time Sea Level Anomaly data stem from Marine Copernicus data](http://marine.copernicus.eu/services-portfolio/access-to-products/?option=com_csw&view=details&product_id=SEALEVEL_GLO_SLA_MAP_L4_REP_OBSERVATIONS_008_027) over the period 2014/04/08 - 2016/08/19. The NRT data are updated each year for this website. These products are fully described [here](http://marine.copernicus.eu/documents/PUM/CMEMS-SL-PUM-008-017-033.pdf).
    -   Sea level from a tide gauge: we used the tide gauge data located in Pointe des Galets, Reunion Island which are collected since 1967/03/01. The data and documentation are available [here](http://www.psmsl.org/data/obtaining/stations/1501.php) and [here](http://www.sonel.org/spip.php?page=maregraphe&idStation=1827.php).
-   **METHOD**:
    -   We extracted the SLA values (NRT & DT data) for the West part of Reunion Island where the tide gaude and the reef lagoons are located.
    -   We compare the times series of this mean SLA and the tide gauge time series. Autocorrelation in the time series is accounted for in the p-value estimates using [Pyper & Peterman, 1998](http://www.fishsciences.net/reports/CJFAS_u-d_6-28/CJ_55_p2127-40_Comparison_methods_to_account_for_autocorrelation_in_correlation_anal_fish_data.pdf).
    -   We extracted the time series of SLA for the whole Southwestern Indian Ocean. We then calculated the correlation between the Reunionese SLA time series and the SLA time series of all other pixels in this area. We did this by lagging the time series to see if we could predict future drop in sea level from previous drop in other area.
-   **RESULTS**:
    -   There is a very good correlation between SLA times series and sea level from the tide gauge
    -   The SLA in Reunion can be predicted from Eastern SLA, e.g. from the SLA Mauritius region 1 month ahead and from the SLA 500 km east from Reunion 2 month earlier
-   **DISCUSSION**
    -   SLA has a strong impact on the sea level in the Reunionese coral lagoons
    -   Sea level drops in these lagoons can be anticipated from SLA values East of Reunion

All these results can be seen at the [appSat website](https://scientific-contributions.shinyapps.io/appSat). You can get and run it easily:

``` r
library(shiny)

# Easiest way is to use runGitHub
runGitHub("appSat", "sylvainBonhommeau")

# Run a tar or zip file directly
runUrl("https://github.com/sylvainBonhommeau/appSat/archive/master.tar.gz")
runUrl("https://github.com/sylvainBonhommeau/appSat/archive/master.zip")
```
