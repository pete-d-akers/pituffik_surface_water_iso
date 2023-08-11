# pituffik_surface_water_iso
Data and R code for analyzing the surface water isotopes of the Pituffik Peninsula, Greenland

This repository contains the data and code needed to perform analyses on the variability of water isotopic composition (d18O, d2H, dxs) across the Pituffik Peninsula, Greenland. Water samples were collected from June through August 2018, in November 2018, and in July 2019.

The file pituffik_water_iso_analyses_script.r contains the R code to perform all analyses and plotting. It calls upon datafiles included here. Packages needed to run to code are all listed within the provided code.

The primary datafile is pituffik_H2O_iso_2018_2019.csv. This data contains every water sample collected and isotopically analyzed. Data columns for this file are:
sample_id: The unique sample identification for each water sample.
date: The date that the water was sampled on.
latitude: The latitude of the water sample site, in decimal degrees.
longitude: The longitude of the water sample site, in decimal degrees.
elevation: The elevation of the water sample site, in meters above sea level, extracted from ArcticDEM based on geographic coordinate.
d18O: The oxygen stable isotopic ratio of the water sample in delta notation, measured on a Picarro L2130-i. Multiply by 1000 to express in permil (‰).
d2H: The hydrogen stable isotopic ratio of the water sample in delta notation, measured on a Picarro L2130-i. Multiply by 1000 to express in permil (‰).
dxs: The deuterium excess of the water sample, calculated as dxs = 8*d2H - d18O. Multiply by 1000 to express in permil (‰).
site_name: The common name for the sampling site.
type: The type category of water source sampled.
type_number: The type category of water source sampled, expressed as numeric categories where 1 = lake, 2 = pool, 3 = stream, 4 = surface flow, 5 = snow or ice, 6 = precipitation rain event, and 7 = precipitation snow event.
laketype: The type category of lake sampled, if applicable.
laketype_number: The type category of lake sampled, expressed as numeric categories where 1 = endorheic, 2 = headwater, 3 = downstream, 4 = vale, 5 = proglacial, and 6 = altered.
lake_surface_area: The surface area of lake sampled, in m².
basin_name: The name of the watershed basin that the sample is located within.
alt_basin_name: The alternate or secondary name of the watershed basin that the sample is located within.
lakeshed_area: The surface area of the watershed of lake sampled, if applicable, in m².
dist_gris: The distance from the sample site to the nearest margin of the Greenland ice sheet, in m.
dist_ocean: The distance from the sample site to the nearest coast of the ocean, in m.
timing: A marker indicating samples that fell during one of three specific sampling periods where 1 = early summer 2018, 2 = late summer 2018, and 3 = mid-summer 2019.
main_lakes: A marker indicating whether a lake is a member of the main lakes region (1), if applicable.
qc_flag: A marker indicating samples that were flagged (1) during quality checking for vial cracking and/or evaporative water loss.

The file ptk_iso_wx_day.csv contains daily weather and water vapor isotope data taken at Pituffik Space Base from August 2017 through May 2020. Data previously published in Akers et al., 2020 (https://doi.org/10.5194/acp-20-13929-2020), and more information on data sourcing and collection can be found there. Data columns for this file are:
daybreak: The day that the rowdata covers.
wind_az: Mean azimuth of wind, in degrees.
wind_sp: Mean wind speed, in m/s.
H2O: Specific humidity of ambient air, measured on a Picarro L2130-i, in ppmv.
d18O: Mean oxygen stable isotopic ratio of ambient air in delta notation, measured on a Picarro L2130-i. Multiply by 1000 to express in permil (‰).
d2H: Mean hydrogen stable isotopic ratio of ambient air in delta notation, measured on a Picarro L2130-i. Multiply by 1000 to express in permil (‰).
dxs: Mean deuterium excess of ambient air, calculated as dxs = 8*d2H - d18O. Multiply by 1000 to express in permil (‰).
tavg: Mean of 10 min ambient air temperature readings at SMT site, in °C.
tdew: Mean of 10 min ambient air dew point temperature readings at SMT site, in °C.
pres: Mean of 10 min station air pressure readings at SMT site, in mb.
rh_ice: Mean of 10 min relative humidity with respect to ice readings at SMT site, in percent.
nao: Daily NAO index.
ao: Daily AO index.
aao: Daily AAO index.
pna: Daily PNA index.
seaice: Baffin Bay sea ice extent, in km².
Eness: Katabatic deviation of wind, defined as |100° - wind azimuth|, in degrees.
mo: The month of the rowdata.
sn: The season of the rowdata.
tmax_af: Maximum air temperature measured by USAF at THU airport, in °C.
tmin_af: Minimum air temperature measured by USAF at THU airport, in °C.
tavg_af: Mean air temperature measured by USAF at THU airport, in °C.
tdew_af: Mean dew point temperature measured by USAF at THU airport, in °C.
wind_pk_az_af: Azimuth of maximum wind speed measured by USAF at THU airport, in degrees.
wind_pk_sp_af: Maximum wind speed measured by USAF at THU airport, in m/s.
prcp_snow: Precipitation amount of snow measured by USAF at THU airport, in mm.
prcp_H2O: Liquid precipitation amount measured by USAF at THU airport, in mm.
snow_depth: Depth of snow on ground measured by USAF at THU airport, in cm.
prcp_occur: Precipitation recorded during day by USAF at THU airport, including trace amounts.

The file pfk_gnip.csv contains monthly precipitation isotope data for Pituffik from 1966–1971, collected as part of GNIP (https://nucleus.iaea.org/wiser). Data columns for this file are:
site: The site identification.
date_start: The starting date of precipitation collection.
date_end: The ending data of preciptiation collection.
d18O: Oxygen stable isotopic ratio of collected precipitation in delta notation. Multiply by 1000 to express in permil (‰).
d2H: Hydrogen stable isotopic ratio of collected precipitation in delta notation. Multiply by 1000 to express in permil (‰).
dxs: Deuterium excess of collected precipitation, calculated as dxs = 8*d2H - d18O. Multiply by 1000 to express in permil (‰).
precip_amt: Amount of collected precipitation, in mm.
air_temp: Mean air temperature during collection period, in °C.
vapor_pres: The mean vapor pressure during collection period, in kPa.
type: Indicating that samples are GNIP for later comparison with Pituffik surface water data.
qc_flag: Marker for values flagged during quality checking.

The two files of 2018_daily_pet_pituffik.nc and 2019_daily_pet_pituffik.nc are potential evapotranspiration (PET) spatial data extracted for the Pituffik region from the global database provided by Singer et al., 2021 (https://doi.org/10.1038/s41597-021-01003-9). They are used by the R script to extract time series of PET for Pituffik.

