/**********************************************************************
 JET-II jet finding algorithm in the double anti-KT cone approximation
 Authors:     Yang Bai, Zhenyu Han, Ran Lu
 Version:     1.00
 Reference:   arXiv
 Last update: Sep 20, 2015
**********************************************************************/

#ifndef DOUBLECONE
#define DOUBLECONE

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
using namespace std;
namespace double_cone{

#define RSUBMAX 0.7
#define NSTEP 14
/*
Usage: Prepare the particles you want to find jets from in fastjet::PseudoJet format ("inputPars"). The function returns the jets found, as well as particles in each jet ("jetsPars") and unused particles ("restPars"). The function terminates when "njets" jets are found, or, if njets = -1, when all particles are used.
*/
vector <fastjet::PseudoJet> find_jets(vector<fastjet::PseudoJet> &inputPars, double beta, double gamma, vector<vector<fastjet::PseudoJet> > &jetsPars, vector<fastjet::PseudoJet> &restPars, int njets = -1, int ET2 = 1);


//cluster
vector<fastjet::PseudoJet> resolve(vector<fastjet::PseudoJet> &pars, fastjet::ClusterSequence &resolved_seq, double R, fastjet::JetAlgorithm algorithm = fastjet::antikt_algorithm)
{
   fastjet::Strategy               strategy = fastjet::Best;
   fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
   fastjet::JetDefinition          jetDef(algorithm, R, recombScheme, strategy);

   fastjet::ClusterSequence tempseq(pars, jetDef);
   resolved_seq = tempseq;
   vector<fastjet::PseudoJet> inclusiveJets = resolved_seq.inclusive_jets();
   vector<fastjet::PseudoJet> sortedJets    = sorted_by_pt(inclusiveJets);

   return sortedJets;
}

fastjet::PseudoJet sumMomenta(vector<fastjet::PseudoJet> &pj)
{

   double ptotal[4];
   unsigned npars = pj.size();

   for (int i = 0; i < 4; i ++) ptotal[i] = 0;

   for (unsigned i = 0; i < npars; i ++) {
       ptotal[0] += pj[i].e();
       ptotal[1] += pj[i].px();
       ptotal[2] += pj[i].py();
       ptotal[3] += pj[i].pz();
   }

   fastjet::PseudoJet totalJ(ptotal[1], ptotal[2], ptotal[3], ptotal[0]);
   totalJ.set_user_index(npars);
   return totalJ;
}

double Jbeta_had(fastjet::PseudoJet j, double beta)
{
   double pt = j.perp();
   double m  = j.m();
   double Et = sqrt(pt*pt + m*m);
   return Et - beta*m*m/Et;
}

double dot(fastjet::PseudoJet &p1, fastjet::PseudoJet &p2)
{
   return p1.e()*p2.e() - (p1.px()*p2.px() +  p1.py()*p2.py() +  p1.pz()*p2.pz());
}

//calculte jet function
double Jbeta_gamma_restframe(vector<fastjet::PseudoJet> &p, double beta, double gamma, int ET2 = 1)
{

    fastjet::PseudoJet pjet = sumMomenta(p);
    double pt = pjet.perp();
    double ET = sqrt(pt*pt + pjet.m()*pjet.m());

    double tempJ = Jbeta_had(pjet, beta);

    double mJ = pjet.m();
    double mJsq = mJ*mJ;

    if (p.size() < 2) {
      if (ET2 == 1) tempJ = tempJ*ET;
      return tempJ;
    }
    for (unsigned i = 0; i < p.size(); i ++)
    for (unsigned j = 0; j < p.size(); j ++) {
        double pJdotpi = dot(p[i], pjet);
        double pJdotpj = dot(p[j], pjet);
	double pidotpj = dot(p[i], p[j]);
        //double temp    = mJsq*pidotpj - pJdotpi*pJdotpj;
	//temp           = temp*temp/(mJsq*pJdotpi*pJdotpj);
	double temp = mJsq*(pidotpj*pidotpj/pJdotpi/pJdotpj);
        tempJ         += gamma*temp/ET;
    }
    tempJ -= gamma*mJsq/ET;
    if (ET2 == 1) tempJ = tempJ*ET;
    return tempJ;
}



/*find one jet from particles in 'inputPars', return the jet as a fastjet::PseudoJet, return the particles in the jet by 'jetPars', return the remaining particles by 'restPars'

ET2 = 1: use ET2 as the jet function
ET2 = 0: use ET as the jet function*/

fastjet::PseudoJet find_one_jet(vector<fastjet::PseudoJet> &inputPars, double beta, double gamma, vector<fastjet::PseudoJet> &jetPars, vector<fastjet::PseudoJet> &restPars, int ET2 = 1)
{

    vector<fastjet::PseudoJet> pars = inputPars; //make a local copy of the input particles
    for (unsigned i = 0; i < pars.size(); i ++) pars[i].set_user_index(i);

    //prepair cones by running anti_kt
    double RsubMax = RSUBMAX;
    double RsubMin = RsubMax/NSTEP;
    double Rstep = RsubMin;
    fastjet::ClusterSequence clustSeq[NSTEP];
    vector<fastjet::PseudoJet> jets;
    vector<vector<int> > cones;
    //vector<double> cone_sizes;

    for (int istep = 0; istep < NSTEP; istep ++)
    {
        double Rsub =  RsubMin + istep*Rstep;
        vector<fastjet::PseudoJet> subjets =  resolve(pars, clustSeq[istep], Rsub, fastjet::antikt_algorithm);

	for (unsigned isub = 0; isub < subjets.size(); isub ++)
	{
	    vector<fastjet::PseudoJet> cons = clustSeq[istep].constituents(subjets[isub]);
	    vector<int> cons_index;
	    cons_index.resize(0);
	    for (unsigned ipar = 0; ipar < cons.size(); ipar ++)
	       cons_index.push_back(cons[ipar].user_index());
	    cones.push_back(cons_index);
	}
   }


   double maxJ = -1, maxcone1, maxcone2;
   jetPars.resize(0);
   vector<fastjet::PseudoJet> subjet1;
   vector<fastjet::PseudoJet> subjet2;

   for (int icone1 = cones.size() - 1; icone1 >= 0; icone1 --) //loop over cones, from bigger to smaller
   for (int icone2 = icone1; icone2 >= 0; icone2 --) {
      vector<fastjet::PseudoJet> conePars;
      vector<fastjet::PseudoJet> conePars1;
      vector<fastjet::PseudoJet> conePars2;
      conePars.resize(0);
      conePars1.resize(0);
      conePars2.resize(0);

      for (unsigned i = 0; i < cones[icone1].size(); i ++) {
 	  conePars.push_back(pars[cones[icone1][i]]);
	  conePars1.push_back(pars[cones[icone1][i]]);
      }

      for (unsigned i = 0; i < cones[icone2].size(); i ++)
	  conePars2.push_back(pars[cones[icone2][i]]);

      if (icone2 != icone1) {
	  for (unsigned i = 0; i < cones[icone2].size(); i ++) {
      	        int itemp = cones[icone2][i];
		bool had = false;
		for (unsigned j = 0; j < cones[icone1].size(); j ++)
		    if (itemp == cones[icone1][j]) {

		       had = true;
		       break;
                    }
		if (!had) conePars.push_back(pars[itemp]);
	  }
       }

       double jf = Jbeta_gamma_restframe(conePars, beta, gamma, ET2);

       if (jf > maxJ) {
	   maxJ     = jf;
	   maxcone1 = icone1;
	   maxcone2 = icone2;
	   jetPars = conePars;
       }

    }

    fastjet::PseudoJet jet = sumMomenta(jetPars);

    restPars.resize(0);
    for (unsigned i = 0; i < pars.size(); i ++) {
        bool used = false;
        for (unsigned j = 0; j < jetPars.size(); j ++) {
            if (jetPars[j].user_index() == pars[i].user_index()) {
	        used = true;
	        break;
	    }
        }

        if (!used) restPars.push_back(pars[i]);
    }

    //reset the user index to the original
    for (unsigned i = 0; i < jetPars.size(); i ++) jetPars[i].set_user_index(inputPars[jetPars[i].user_index()].user_index());
    for (unsigned i = 0; i < restPars.size(); i ++) restPars[i].set_user_index(inputPars[restPars[i].user_index()].user_index());

    return jet;
}

vector <fastjet::PseudoJet> find_jets(vector<fastjet::PseudoJet> &inputPars, double beta, double gamma, vector<vector<fastjet::PseudoJet> > &jetsPars, vector<fastjet::PseudoJet> &restPars, int njets, int ET2)
{

   vector<fastjet::PseudoJet> currentJets;
   vector<fastjet::PseudoJet> currentPars = inputPars;
   jetsPars.resize(0);
   while (currentPars.size() > 0) {
       if (njets > 0 && currentJets.size() >= njets) break;
       vector<fastjet::PseudoJet> jetPars;
       vector<fastjet::PseudoJet> currentRestPars;
       fastjet::PseudoJet jet = find_one_jet(currentPars, beta, gamma, jetPars, currentRestPars, ET2);
       currentJets.push_back(jet);
       jetsPars.push_back(jetPars);
       currentPars = currentRestPars;
   }
   restPars = currentPars;
   return currentJets;
}


//end namespace double_cone
}
#endif
