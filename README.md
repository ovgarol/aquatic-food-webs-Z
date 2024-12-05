# z
databases and scripts for a simple food web model

Code and data available in https://github.com/ovgarol/aquatic-food-webs-Z

## Requirements
### Software
  R, including libraries: scales, latex2exp, igraph, mclust, stringr

### Input data
  database.csv: database of optimal prey size of aquatic predators
  
  minimal_model.csv: database of feeding guilds in aquatic ecosystems

  physio_data: limit1 to limit6 # energetic and mechanical limits described in Portalier et al. (2019). 

  NameYYYY.csv  # data collected and available in original databases https://doi.org/10.4319/lo.1994.39.2.0395,  https://doi.org/10.3354/meps08716,  http://dx.doi.org/10.3354/meps09502, https://doi.org/10.1007/s00227-022-04102-2, https://doi.org/10.1098/rspb.2014.2103, https://doi.org/10.1890/07-1551.1, https://doi.org/10.1890/0012-9658(2006)87[2411:CBRINF]2.0.CO;2 and https://idata.idiv.de/ddm/Data/ShowData/283?version=3

### External data
  283_2_FoodWebDataBase_2018_12_10.csv: GATEWAy database available in https://idata.idiv.de/ddm/Data/ShowData/283?version=3

## Execution
Run the following scripts in the enumerated sequence:

Data-processing-1.R # calculates body size and optimal prey size of provided databases in database.csv and minimal_model.csv 

Data-tester.R # checks data source files

Figure-script-1.R   # calculates the Z-model using the data collated in database.csv and minimal_model.csv

Figure-script-2.R   # applies to generated model in the script above and compare it to the sites of the database GATEWAy database  

Figure-script-3.R   # calculates the horizontal bands

Figure-script-4.R   # plot map of studied ecosystems

Figure-script-5.R   # physical limits for preu acquisition

Figure-script-6.R   # Accuracy and complexity of artificial food webs (supplement B)

Figure-script-7.R   # food web topology (supplement C)

Figure-script-8.R   # comparison of food webs reconstructions sub-samples

Figure-script-9.R   # parsimony analysis using AIC (supplement D part 2)

Figure-script-10.R   # applies the size-only model and compare it to the sites of the database GATEWAy database  

Figure-script-11.R   # residual analysis for optimal number of feeding guilds


## Additional data
WORMS_names_taxonomy.csv: taxonomy of the predator species included in database.csv

## License
If not stated otherwise, the entire analysis software is licensed under
  the GNU Public License version 3 or later.
  See <http://www.gnu.org/licenses/gpl-3.0.txt> for the complete terms.

## Documentation
  see text and equations in Material and Methods of 
