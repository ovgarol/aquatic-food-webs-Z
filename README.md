# z
databases and scripts for a simple food web model

Code and data available in https://github.com/ovgarol/z

## Requirenments
### Software
  R, including libraries: scales, latex2exp

### Input data
  database.csv: database of optimal prey size of aquatic predators
  
  minimal_model.csv: database of feeding guilds in aquatic ecosystems

### External data
  283_2_FoodWebDataBase_2018_12_10.csv: GATEWAy database available in https://idata.idiv.de/ddm/Data/ShowData/283?version=3

## Execution
Run the following two scripts in sequence

Figure-script-1.R   # calculates the Z-model using the data collated in database.csv and minimal_model.csv

Figure-script-2.R   # applies to generated model in the script above and compare it to the sites of the database GATEWAy database  

Figure-script-3.R   # calculates the horizontal bands

Figure-script-4.R   # plot map of studied ecosystems

Figure-script-5.R   # physical limits for preu acquisition

Figure-script-6.R   # Accuracy and complexity of artificial food webs (supplement B)

Figure-script-7.R   # ood web topology (supplement C)

Figure-script-8.R   # comparison of food webs reconstructions sub-samples


## Additional data
WORMS_names_taxonomy.csv: taxonomy of the predator species included in database.csv

## License
If not stated otherwise, the entire analysis software is licensed under
  the GNU Public License version 3 or later.
  See <http://www.gnu.org/licenses/gpl-3.0.txt> for the complete terms.

## Documentation
  see text and equations in Material and Methods of 
