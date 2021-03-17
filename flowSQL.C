#include <TMySQLServer.h>

/*
* function to connect with MySQL server
* returns pointer to the server
*/
TSQLServer* connectSQL(const char* userName, const char* password){
  TSQLServer* server = TMySQLServer::Connect("mysql://localhost/flow",userName,password);
  return server;
}

/*
* function for running a statement in in MySQL to look through the data
* (statements like SHOW or SELECT)
*   INPUT:
*     server - pointer to the databse
*     statement - statement to run in MySQL
*   OUTPUT:
*     printing out MySQL response
*   RETURN (void):
*
*/
void runSentence(TSQLServer* server, const char* statement){
  TSQLStatement* stmt = server->Statement(statement, 100);
  if(stmt->Process()){
    stmt->StoreResult();
    for(int i=0; i<stmt->GetNumFields();i++){
      TString result = stmt->GetFieldName(i);
      if(result == "") std::cout<<"NULL";
      else std::cout<<result;
      std::cout<<" | ";
    }
    std::cout<<std::endl;

    while(stmt->NextResultRow()){
      std::cout<<"| ";
      for(int i=0; i<stmt->GetNumFields();i++){
        TString result = stmt->GetString(i);
        if(result == "") std::cout<<"NULL";
        else std::cout<<result;
        std::cout<<" | ";
      }
      std::cout<<std::endl;
    }
  }
}

/*
* function for drwaing TGraph of V2 vs pT to check data sets
* INPUT:
*    server     - pointer to MySQL server with data
*    tName      - name of the table with data
*    v_param    - which parameter to draw
*    centrality - centrality percentage for which the graph is drawn
*
*/
void draw(TSQLServer* server, TString tName, TString v_param, TString centrality){
  TSQLStatement* stmt = server->Statement(Form("SELECT pT,%s FROM %s WHERE (centrality = \'%s\') AND %s IS NOT NULL",v_param.Data(),tName.Data(),centrality.Data(), v_param.Data()),100);

  std::vector<double>* val1 = new std::vector<double>();
  std::vector<double>* val2 = new std::vector<double>();

  int n = stmt->Process();
  std::cout<<n<<std::endl;
  if(stmt->Process()){
    stmt->StoreResult();
    int i=0;
    while(stmt->NextResultRow()){
      val1->push_back(stmt->GetDouble(0));
      val2->push_back(stmt->GetDouble(1));

      printf("%5f GeV | %lf \n",val1->at(i), val2->at(i));
      i++;
    }
  }else{
    std::cout<<"ERROR: there was a problem receiving data from the database"<<std::endl;
    return;
  }

  double pT[val1->size()];
  double v[val2->size()];
  for(int i=0; i<val1->size();i++){
    pT[i] = val1->at(i);
    v[i] = val2->at(i);
  }

  TCanvas* c = new TCanvas("c",Form("%s vs p_{T} (centrality = %s)",v_param.Data(), centrality.Data()), 10, 10, 900, 800);
  TGraph* gr = new TGraph(val1->size(), pT, v);
  gr->SetTitle(Form("%s vs p_{T} (centrality = %s);p_{T} (GeV); %s",v_param.Data(), centrality.Data(), v_param.Data()));
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.2);
  gr->Draw("AP");
}


/*
* function for drwaing TGraph of V2 vs pT to check data sets (for a givien particle)
* INPUT:
*    server     - pointer to MySQL server with data
*    tName      - name of the table with data
*    v_param    - which parameter to draw
*    centrality - centrality percentage for which the graph is drawn
*    pdg        - pdg code of the particle
*
*/
void drawPDG(TSQLServer* server, TString tName, TString v_param, TString centrality, int pdg){
  TSQLStatement* stmt = server->Statement(Form("SELECT pT,%s FROM %s WHERE (centrality = \'%s\') AND %s IS NOT NULL AND pdg = %d",v_param.Data(),tName.Data(),centrality.Data(), v_param.Data(), pdg),100);

  std::vector<double>* val1 = new std::vector<double>();
  std::vector<double>* val2 = new std::vector<double>();

  int n = stmt->Process();
  std::cout<<n<<std::endl;
  if(stmt->Process()){
    stmt->StoreResult();
    int i=0;
    while(stmt->NextResultRow()){
      val1->push_back(stmt->GetDouble(0));
      val2->push_back(stmt->GetDouble(1));

      printf("%5f GeV | %lf \n",val1->at(i), val2->at(i));
      i++;
    }
  }else{
    std::cout<<"ERROR: there was a problem receiving data from the database"<<std::endl;
    return;
  }

  double pT[val1->size()];
  double v[val2->size()];
  for(int i=0; i<val1->size();i++){
    pT[i] = val1->at(i);
    v[i] = val2->at(i);
  }

  TCanvas* c = new TCanvas("c",Form("%s vs p_{T} (centrality = %s)",v_param.Data(), centrality.Data()), 10, 10, 900, 800);
  TGraph* gr = new TGraph(val1->size(), pT, v);
  gr->SetTitle(Form("%s vs p_{T} (centrality = %s);p_{T} (GeV); %s",v_param.Data(), centrality.Data(), v_param.Data()));
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.2);
  gr->Draw("AP");
}


/*
* funtion to show all the data available in the table
*   INPUT:
*      server - pointer to MySQL server
*      tName  - name of the table which should be written out
*   RETURN:
*      number of entries in the table that have been written out.
*      if there was a problem with table, function returns -1
*/
int getData(TSQLServer* server, TString tName){

  TSQLStatement* stmt = server->Statement(Form("SELECT * FROM %s",tName.Data()),100);

  int entries = 0;
  printf("E (GeV) | centrality | pT (GeV) | v2 | v3 \n");
  if(stmt->Process()){
    stmt->StoreResult();
    while(stmt->NextResultRow()){
      Double_t energy = stmt->GetDouble(1);
      TString centrality = stmt->GetString(2);
      Double_t pT = stmt->GetDouble(3);
      Double_t u_pT = stmt->GetDouble(4);
      Double_t v2 = stmt->GetDouble(5);
      Double_t u_v2 = stmt->GetDouble(6);
      Double_t v3 = stmt->GetDouble(7);
      Double_t u_v3 = stmt->GetDouble(8);
      entries++;

      printf("%5f GeV | %10s | %.3f GeV | %.6f | %.6f\n",energy,centrality.Data(),pT, v2, v3);
    }
  }else{
    std::cout<<"ERROR: there was a problem receiving data from the database"<<std::endl;
    return -1;
  }

  return entries;
}

/*
* Function for creating new MySQL tables in the "flow" database
* for different collision systems (for example AuAu or pp)
*   INPUT:
*     server - pointer to MySQL server
*     tName  - name for the table you want to create
*   RETURN:
*     int - 1 if the table has been created successfully and -1 if there was
*           a problem. 0 if a table with that name already exists
*/
int createTable(TSQLServer* server, const char* tName ){
  TSQLStatement* stmt = server->Statement("SHOW TABLES;",100);
  if(stmt->Process()){
    stmt->StoreResult();
    while(stmt->NextResultRow()){
        TString result = stmt->GetString(0);
        if(result == tName){
          std::cout<<"INFO: table named "<<tName<<" already exists!"<<std::endl;
          return 0;
        }
      }
  }else{
    std::cout<<"ERROR: there was a problem receiving data from the database!"<<std::endl;
    return -1;
  }

  TString statement = Form("CREATE TABLE %s(",tName);
  statement += "id int NOT NULL AUTO_INCREMENT,";
  statement += "pdg int,";
  statement += "energy double,";
  statement += "centrality varchar(8),";
  statement += "eta varchar(10),";
  statement += "pT double,";
  statement += "pT_LOW double,";
  statement += "pT_HIGH double,";
  statement += "v2 double,";
  statement += "v2_statP double,";
  statement += "v2_statM double,";
  statement += "v2_sysP double,";
  statement += "v2_sysM double,";
  statement += "v3 double,";
  statement += "v3_statP double,";
  statement += "v3_statM double,";
  statement += "v3_sysP double,";
  statement += "v3_sysM double,";
  statement += "reference varchar(60) NOT NULL,";
  statement += "experiment varchar(8),";
  statement += "PRIMARY KEY (id) );";
  runSentence(server, statement.Data());
  return 1;
}


/*
* function to insert data into MySQL table names tTable
*    INPUT:
*       server     - MySQL server at which table is lovated
*       tName      - name of the table
*       v_param    - name of the parameter
*       energy     - value for energy ( sqrt(s_NN) ) in GeV
*       centrality - centrality of the collision
*                    in format for example "0-5%"
*       pT         - value of pT in GeV
*       v          - value of v parameter
*    RETURN:
*       1 if data was succesfully inserted and -1 if something failed
*/
int insertData(TSQLServer* server,
               TString tName,
               TString v_param,
               double energy,
               TString centrality,
               TString eta,
               double pT,
               double pT_low,
               double pT_high,
               double v,
               double statP,
               double statM,
               double sysP,
               double sysM,
               TString reference,
               TString experiment)
{

  TString v_statP = v_param;
  v_statP += "_statP";
  TString v_statM = v_param;
  v_statM += "_statM";
  TString v_sysP = v_param;
  v_sysP += "_sysP";
  TString v_sysM = v_param;
  v_sysM += "_sysM";
  TSQLStatement* stmt = server->Statement(Form("INSERT INTO %s (energy, centrality, eta, pT, pT_LOW, pT_HIGH, %s, %s, %s, %s, %s, reference, experiment) VALUES (%lf, \"%s\", \"%s\", %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, \"%s\",\"%s\")",
                                                tName.Data(),
                                                v_param.Data(),
                                                v_statP.Data(),
                                                v_statM.Data(),
                                                v_sysP.Data(),
                                                v_sysM.Data(),
                                                energy,
                                                centrality.Data(),
                                                eta.Data(),
                                                pT,
                                                pT_low,
                                                pT_high,
                                                v,
                                                statP,
                                                statM,
                                                sysP,
                                                sysM,
                                                reference.Data(),
                                                experiment.Data() ), 1000);

  if(stmt->Process())
    return 1;
  else
    return -1;
}

/*
* function for filling MySQL table with data from CSV file. CSV file is a tmp
* file that is created by removing first info lines from HEPdata file
* tmp format: pT, pT low, pT high, v2, stat +, stat-, sys +, sys -
*    INPUT:
*      fileName   - name of the CSV file (tmp format)
*      server     - pointer to the MySQL server (from connectSQL() function)
*      tName      - name of a table to which data should be added
*      v_param    - name of the parameter (for example "v2")
*      energy     - energy of the collision (in GeV)
*      centrality - centrality of the collision (format "low-high%")
*      eta        - pseudorapidity
*      reference  - reference to the table from HEPdata or directly to the publication
*      experiment - name of the experiment that did the measurements
*   RETURN:
*      number of entries added to the table.
*         If there was an error it will return -1
*         If this data already is in the table function returns 0
*/

int readCSV(const char* fileName, TSQLServer* server, TString tName, TString v_param, double energy, TString centrality, TString eta, TString reference, TString experiment){
  int entries = 0;

  std::ifstream file(fileName);

  if(!file.is_open()){
    std::cout<<"ERROR: Cannot open file "<<fileName<<std::endl;
    return -1;
  }

  TSQLStatement* stmt = server->Statement(Form("SELECT COUNT(%s) FROM %s WHERE reference == %s AND pdg == NULL AND %s IS NOT NULL;", v_param.Data(), tName.Data(), reference.Data(), v_param.Data()),100);
  stmt->StoreResult();
  if(stmt->Process()){
    stmt->NextResultRow();
    if(stmt->GetInt(0) != 0){
      std::cout<<"INFO: values from a table with this reference ("<<reference<<") has been already inserted into the "<<tName<<" table!"<<std::endl;
      return 0;
    }
  }

  std::string line, colname;
  double value, pT, pT_low, pT_high, v, sysP, sysM, statP, statM;

  if(file.good()){
    std::getline(file,line);

    while(std::getline(file,line)){
      std::stringstream ss(line); //line with titles of columns
      //csv file is saved as: pT, pT low, pT high, v2, stat +, stat-, sys +, sys -
      int colID = 0;
      while(ss>>value){
      	if(colID == 0) pT = value;
        if(colID == 1) pT_low = value;
        if(colID == 2) pT_high = value;
      	if(colID == 3) v = value;
        if(colID == 4) statP = value;
        if(colID == 5) statM = value;
        if(colID == 6) sysP = value;
        if(colID == 7) sysM = value;

      	if(ss.peek() == ',') ss.ignore();
      	colID++;
      }
      insertData(server, tName, v_param, energy, centrality, eta, pT, pT_low, pT_high, v, statP, statM, sysP, sysM, reference, experiment);
      entries++;
      //std::cout<<pT<<" "<<v2<<std::endl;
    }
  }

  return entries;
}

/*
* function for filling MySQL table with data from CSV file. The data format is 
* different but function works as the function above. The data in CSV file was
* not saved as a histogram in HEPdata.
*   tmp file format: pT, v2, stat +, stat-, sys +, sys -
*/
int readCSVnotHist(const char* fileName, TSQLServer* server, TString tName, TString v_param, double energy, TString centrality, TString eta, TString reference, TString experiment){
  int entries = 0;

  std::ifstream file(fileName);

  if(!file.is_open()){
    std::cout<<"ERROR: Cannot open file "<<fileName<<std::endl;
    return -1;
  }

  TSQLStatement* stmt = server->Statement(Form("SELECT COUNT(%s) FROM %s WHERE reference == %s AND pdg == NULL AND %s IS NOT NULL;", v_param.Data(), tName.Data(), reference.Data(), v_param.Data()),100);
  stmt->StoreResult();
  if(stmt->Process()){
    stmt->NextResultRow();
    if(stmt->GetInt(0) != 0){
      std::cout<<"INFO: values from a table with this reference ("<<reference<<") has been already inserted into the "<<tName<<" table!"<<std::endl;
      return 0;
    }
  }

  std::string line, colname;
  double value, pT, pT_low, pT_high, v, sysP, sysM, statP, statM;
  double pT_last = 0.0;
  if(file.good()){
    std::getline(file,line);
    while(std::getline(file,line)){
      std::stringstream ss(line); //line with titles of columns
      //csv file is saved as: pT, v2, stat +, stat-, sys +, sys -
      int colID = 0;
      while(ss>>value){
      	if(colID == 0){
           pT = value;
           pT_low = pT_last;
           pT_high = pT;
         }
      	if(colID == 1) v = value;
        if(colID == 2) statP = value;
        if(colID == 3) statM = value;
        if(colID == 4) sysP = value;
        if(colID == 5) sysM = value;

      	if(ss.peek() == ',') ss.ignore();
      	colID++;
      }
      insertData(server, tName, v_param, energy, centrality, eta, pT, pT_low, pT_high, v, statP, statM, sysP, sysM, reference, experiment);
      entries++;
      //std::cout<<pT<<" "<<v2<<std::endl;
    }
  }

  return entries;
}


/********** alteranative function for specified particle with pdg ************/
/********* same functions as above with additional argument for pdg **********/
/*
* function to insert data into MySQL table names tTable
*/
int insertData(TSQLServer* server,
               TString tName,
               TString v_param,
               int pdg,
               double energy,
               TString centrality,
               TString eta,
               double pT,
               double pT_low,
               double pT_high,
               double v,
               double statP,
               double statM,
               double sysP,
               double sysM,
               TString reference,
               TString experiment)
{

  TString v_statP = v_param;
  v_statP += "_statP";
  TString v_statM = v_param;
  v_statM += "_statM";
  TString v_sysP = v_param;
  v_sysP += "_sysP";
  TString v_sysM = v_param;
  v_sysM += "_sysM";
  TSQLStatement* stmt = server->Statement(Form("INSERT INTO %s (pdg, energy, centrality, eta, pT, pT_LOW, pT_HIGH, %s, %s, %s, %s, %s, reference, experiment) VALUES (%d, %lf, \"%s\", \"%s\", %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, \"%s\",\"%s\")",
                                                tName.Data(),
                                                v_param.Data(),
                                                v_statP.Data(),
                                                v_statM.Data(),
                                                v_sysP.Data(),
                                                v_sysM.Data(),
                                                pdg,
                                                energy,
                                                centrality.Data(),
                                                eta.Data(),
                                                pT,
                                                pT_low,
                                                pT_high,
                                                v,
                                                statP,
                                                statM,
                                                sysP,
                                                sysM,
                                                reference.Data(),
                                                experiment.Data() ), 1000);

  if(stmt->Process())
    return 1;
  else
    return -1;
}

/*
* function for filling MySQL table with data from CSV file
*   tmp file format: pT, pT low, pT high, v2, stat +, stat-, sys +, sys -
*/

int readCSV(const char* fileName, TSQLServer* server, TString tName, TString v_param, int pdg, double energy, TString centrality, TString eta, TString reference, TString experiment){
  int entries = 0;

  std::ifstream file(fileName);

  if(!file.is_open()){
    std::cout<<"ERROR: Cannot open file "<<fileName<<std::endl;
    return -1;
  }

  TSQLStatement* stmt = server->Statement(Form("SELECT COUNT(%s) FROM %s WHERE reference == %s AND pdg == %d AND %s IS NOT NULL;", v_param.Data(), tName.Data(), reference.Data(), pdg, v_param.Data()),100);
  stmt->StoreResult();
  if(stmt->Process()){
    stmt->NextResultRow();
    if(stmt->GetInt(0) != 0){
      std::cout<<"INFO: values from a table with this reference ("<<reference<<") has been already inserted into the "<<tName<<" table!"<<std::endl;
      return 0;
    }
  }

  std::string line, colname;
  double value, pT, pT_low, pT_high, v, sysP, sysM, statP, statM;

  if(file.good()){
    std::getline(file,line);

    while(std::getline(file,line)){
      std::stringstream ss(line); //line with titles of columns
      //csv file is saved as: pT, pT low, pT high, v2, stat +, stat-, sys +, sys -
      int colID = 0;
      while(ss>>value){
      	if(colID == 0) pT = value;
        if(colID == 1) pT_low = value;
        if(colID == 2) pT_high = value;
      	if(colID == 3) v = value;
        if(colID == 4) statP = value;
        if(colID == 5) statM = value;
        if(colID == 6) sysP = value;
        if(colID == 7) sysM = value;

      	if(ss.peek() == ',') ss.ignore();
      	colID++;
      }
      insertData(server, tName, v_param, pdg, energy, centrality, eta, pT, pT_low, pT_high, v, statP, statM, sysP, sysM, reference, experiment);
      entries++;
      //std::cout<<pT<<" "<<v2<<std::endl;
    }
  }

  return entries;
}

/*
* function to fill data that was not saved as a histogram in the HEPdata
*   tmp file format: pT, v2, stat +, stat-, sys +, sys -
*/

int readCSVnotHist(const char* fileName, TSQLServer* server, TString tName, TString v_param, int pdg, double energy, TString centrality, TString eta, TString reference, TString experiment){
  int entries = 0;

  std::ifstream file(fileName);

  if(!file.is_open()){
    std::cout<<"ERROR: Cannot open file "<<fileName<<std::endl;
    return -1;
  }

  TSQLStatement* stmt = server->Statement(Form("SELECT COUNT(%s) FROM %s WHERE reference == %s AND pdg == %d AND %s IS NOT NULL;", v_param.Data(), tName.Data(), reference.Data(), pdg, v_param.Data()),100);
  stmt->StoreResult();
  if(stmt->Process()){
    stmt->NextResultRow();
    if(stmt->GetInt(0) != 0){
      std::cout<<"INFO: values from a table with this reference ("<<reference<<") has been already inserted into the "<<tName<<" table!"<<std::endl;
      return 0;
    }
  }

  std::string line, colname;
  double value, pT, pT_low, pT_high, v, sysP, sysM, statP, statM;
  double pT_last = 0.0;
  if(file.good()){
    std::getline(file,line);
    while(std::getline(file,line)){
      std::stringstream ss(line); //line with titles of columns
      //csv file is saved as: pT, v2, stat +, stat-, sys +, sys -
      int colID = 0;
      while(ss>>value){
        if(colID == 0){
           pT = value;
           pT_low = pT_last;
           pT_high = pT;
           pT_last=pT;
         }
        if(colID == 1) v = value;
        if(colID == 2) statP = value;
        if(colID == 3) statM = value;
        if(colID == 4) sysP = value;
        if(colID == 5) sysM = value;

        if(ss.peek() == ',') ss.ignore();
        colID++;
      }
      insertData(server, tName, v_param, pdg, energy, centrality, eta, pT, pT_low, pT_high, v, statP, statM, sysP, sysM, reference, experiment);
      entries++;
      //std::cout<<pT<<" "<<v2<<std::endl;
    }
  }

  return entries;
}
