This README_DB.txt file was generated on 2022-05-03 by Coralie Picoche & Frédéric Barraquand.

The dataset contained in DBWithMonthlyPhotoTeich.csv is analysed in Looking for compensation at multiple scales in a wetland bird community (2022), Ecology and Evolution, by Frédéric Barraquand, Coralie Picoche, Christelle Aluome, Laure Carassou & Claude Feigné. https://doi.org/10.1002/ece3.8876 

It contains monthly abundances for birds at the Teich reserve in Arcachon Bay, France. The species list is made mostly of wader and waterfowl species, although some other groups are present as well. This dataset is dubbed the "monthly photograph", as we attempt to create a snapshot of monthly abundance to homogeneize data over the long-term duration of this study (1973-2016, with most data covering 1981-2016). 


GENERAL INFORMATION

1. Title of Dataset: Data from: Looking for compensation at multiple scales in a wetland bird community

2. Author Information

	A. First author, and corresponding author

		Name: Frédéric Barraquand

		Institution: Institute of Mathematics of Bordeaux, CNRS, France

		Email: frederic.barraquand@u-bordeaux.fr


	B. Co-author

		Name: Coralie Picoche

		Institution: Institute of Mathematics of Bordeaux, University of Bordeaux, France

		Email: coralie.picoche@u-bordeaux.fr

	C. Co-author

		Name: Christelle Aluome

		Institution: ISPA, Bordeaux Science Agro & INRAE, France

		Email: christelle.aluome@inrae.fr

	D. Co-author

		Name: Laure Carassou

		Institution: EABX, INRAE, France

		Email: laure.carassou@inrae.fr



	E. Co-author

		Name: Claude Feigné

		Institution: Teich Ornithological Reserve, PNR Landes Gascogne, France



3. Date of data collection: 1973 - 2016

4. Geographic location of data collection:  Teich Ornithological Reserve, Arcachon Bay, France (44.64°N / -1.02°E)

5. Funding source: Landes Gascogne regional park and Teich municipality for data collection, and LabEx COTE (ANR‐10‐LABX‐45) for data formatting


SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data:  Creative Commons Attribution 4.0 International License

2. Publications that cite or use the data: Barraquand, F. et al. (2022) Looking for compensation at multiple scales in a wetland bird community (2022), Ecology and Evolution. DOI: 10.1002/ece3.8876

3. Recommended citation for this dataset: Barraquand, F. et al. (2022) Data from: Looking for compensation at multiple scales in a wetland bird community (2022) DOI: 10.5061/dryad.zpc866t9v


METHODOLOGICAL INFORMATION

1. Description of methods used for collection/generation of data: The staff of the Teich reserve have patrolled the 120 ha of wetlands that constitute the reserve multiple times every month to count and identify the birds using binoculars. Birdwatchers visiting the reserve were also encouraged to report their observations to the staff. The reserve comprises 18 different sectors in which observations were made. Before 2007, observations corresponding to the maximum for the month and the whole reserve were only given for a single day (the 15th of the month). This corresponds to Protocol 1 in the dataset. From 2007 onwards, observations were provided by sector. 

2. Methods for processing the data: In order to homogeneize the data, we reconstructed a "monthly snapshot". After 2007, we needed first to compute the sum of abundances over all sectors throughout the reserve for a given date, then use the maximum observed abundance over a month.  See https://doi.org/10.5281/zenodo.6367146 for the corresponding code, in R (run with 3.6.3).


DATA-SPECIFIC INFORMATION FOR: DBWithMonthlyPhotoTeich.csv

1. Number of variables: 13

2. Number of rows: 118699

3. Variable List: 
	
	Ref: identifier of the data
	
	Nom_espece: French name of the bird species
	
	Nom_latin: Latin name of the bird species
	
	Date: Complete date for the observation. 
		
	Jour: Day of the observation
	
	Mois: Month of the observation
	
	Annee: Year of the observation
	
	Jour_de_l_année: Day in the year, i.e. number of days since the beginning of the year
	
	ID.Lieu.dit: Identifier for the place where the observation was made
	
	Lieu_dit: Name of the place where the observation was made
	
	Nombre: Number of observed birds
	
	Protocol: "0" corresponds to raw data (not a monthly observation, or the maximum over a month), "1" corresponds to fully standardized reserve-scale "monthly snapshot" (i.e. the maximum abundance observed throughout a month, with a large spatial coverage, corresponding to >10 sectors surveyed) while "2" corresponds to milder spatial coverage or more concentrated observations (<10 sectors with data).
	
	number_of_LDs_surveyed: Number of sectors over which observations were summed when needed.


