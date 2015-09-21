#include "double_cone.hh"
#include <fstream>

using namespace std;
int main()
{
  ifstream f("ww_14TeV.dat", ios::in);
  int npars;
  int nevents = 1;
  for (int ievent = 0; ievent < nevents; ievent ++) {
     f >> npars; //number of particles
     int id;
     double px, py, pz, e;
     vector<fastjet::PseudoJet> inputPars;
     inputPars.resize(0);
     for (int i = 0; i < npars; i ++) { //read particles
        f >> id >> px >> py >> pz >> e;
        e = sqrt(px*px+py*py+pz*pz);
        fastjet::PseudoJet tempPar(px, py, pz, e);
        //if (fabs(tempPar.eta()) < 5)
        inputPars.push_back(tempPar);
     }
     //if (ievent == 0) continue;
     vector<fastjet::PseudoJet> jets; //jets to be found
     vector<vector<fastjet::PseudoJet> > jetsPars; //particles in each jet
     vector<fastjet::PseudoJet> restPars; //unused particles
     double beta = 11, gamma = 9.5;
     int nMax = 5; //max number of jets to be found, -1 = find all jets
     jets = double_cone::find_jets(inputPars, beta, gamma, jetsPars, restPars, nMax);

     //output some information
     cout << "Event " << ievent << endl;
     cout << jets.size() << " jets found." << endl;

     for (unsigned i = 0; i < jets.size(); i ++) {
        cout << "jet " << i << " contains " << jetsPars[i].size() << " particles" << endl;
        cout << "  momentum (px, py, pz, e, m) = " << jets[i].px() << "\t" << jets[i].py() << "\t" << jets[i].pz() << "\t" << jets[i].e() << "\t" << jets[i].m() << endl;
	for (unsigned j = 0; j < jetsPars[i].size(); j ++) {
  	    cout << "    particle " << j << " momentum = " << jetsPars[i][j].px() << "\t" << jetsPars[i][j].py() << "\t" << jetsPars[i][j].pz() << "\t" << jetsPars[i][j].e() << "\t" << jetsPars[i][j].m() << endl;
        }
     }

     cout << restPars.size() << " particles unused." << endl << endl;

  }
  return 0;
}



