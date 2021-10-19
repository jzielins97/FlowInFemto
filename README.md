# FlowInFemto
Program for calculating background correlation function with including eliptic flow. Heavy Ion Physics, High Energy Physics.

## Usage
This program should be used for calculating background from anizotropic flow for the femtoscopic correlation funtions. It uses a database that stores information about flow coefficients. This data was measured in experiments (like ALICE). The database provides name of the experiment as well as the exact reference to the table in HepData where the values can be found. The connection with the database is handled by the program and the user should only provide information about the analysis. The main program is called cfylmFLOW.cxx. Inside is defined entire analysis and creation of the correlation function (background). After compiling the program (simple make command) a CfBackgroundFlow executable is created.

To run program a user needs:
- pid of first particle;
- pT disrtibution of first particle (in root file as a TH1D object called "hpt");
- pid of second particle;
- pT distribution of second particle (in root file as TH1D object called "hpt");
- N - number of pairs to create
//optional
- centrality of the collision (format "low-high%", i.e. "0-10%") Note: make sure that the database has the same format. If you want to calculate for "0-10%" but database has points "0-5%" and "5-10%" it will return no values to the program;
- eta (i.e. "<0.8", "0.5<||<1.2");
- energy (in GeV);
- experiment name (i.e. "ALICE")

Note: there is only one table availbale for now called 'PbPb.' There is no data for other collisions.

### Expanding the database 
To add data to the database, user can use flowSQL.C macro. To do so load macro in root interface by running:
````C
.L flowSQL.C
````
Create TSQLServer object to connect to the database:
````C
TSQLServer* srv = connectSQL('username','password');
````
You need an user in the database that has privlidges to edit tables.
There are two main functions for filling the database: readCSV() and readCSVnotHist(). The first should be used for the data in CSV file that is storred in a historam format (there are pT_LOW and pT_HIGH values). The CSV file should only have the first row with column names and rows with data (as it appears in HEPdata file) so all the informational lines should be removed. The second function can be used for the data saved as graphs (so only pT has no bin limits).
  
Both functions have two versions, one for all charged particles (so not identified particles) and one for identified particles (with given PDG code).

## Prerequisition
- ROOT 6
- GSL

## Instalation
make 
