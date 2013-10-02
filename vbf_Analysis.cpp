/*
 c++ -o vbf_Analysis `root-config --glibs --cflags` vbf_Analysis.cpp


*/

#include "LHEF.h"
#include "TH1.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include <map>
#include <cmath>
#include <algorithm>


using namespace std ;


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


struct mySort: public std::binary_function<TLorentzVector *, TLorentzVector *, bool>
{
  bool operator() (TLorentzVector * x, TLorentzVector * y)
    {
      return x->Pt () > y->Pt () ;
    }
} ;

	double getDeltaPhi(double phi1, double phi2  )
	{
	const double PI = 3.14159265;
	double result = phi1 - phi2;
	if(result > PI)
	{result = result - 2 * PI;}
	if(result <= (-1 * PI))
	{result = result + 2 * PI;}
	result = TMath::Abs(result);
	return result;
	}

               float delta_phi(float phi1, float phi2)
                {
                        const float PI=2.0*acos(0.);
                        const float TWOPI=2.0*PI;
                        float PHI=fabs(phi1-phi2);
                        return (PHI<=PI)? PHI : TWOPI-PHI;
                }


	float delta_R(float eta1, float phi1, float eta2, float phi2)
	{	
	const float PI=2.0*acos(0.);
	const float TWOPI=2.0*PI;
	float deta = eta1-eta2;
	float dphi = delta_phi(phi1,phi2);
	return sqrt(deta*deta + dphi*dphi);
	}	


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


void fillHistos (ifstream & ifs) 
{
  LHEF::Reader reader (ifs) ;
  TH1F h_deta ("h_deta", "h_deta", 100, 0, 10) ;
  TH1F h_jeta ("h_jeta", "h_jeta", 100, -5, 5) ;
  TH1F h_j1eta ("h_j1eta", "h_j1eta", 100, -5, 5) ;
  TH1F h_j2eta ("h_j2eta", "h_j2eta", 100, -5, 5) ;
  TH1F h_j1pt ("h_j1pt", "h_j1pt", 100, 0, 200) ;
  TH1F h_j2pt ("h_j2pt", "h_j2pt", 100, 0, 200) ;


   //PG histograms
//method 1 hist
// tag Jets
   TH1F h_vbf_jj_pt ("h_vbf_jj_pt", "h_vbf_jj_pt", 100, 0, 500) ;
   TH1F h_vbf_jj_eta ("h_vbf_jj_eta", "h_vbf_jj_eta", 100, -5, 5) ;
   TH1F h_vbf_jj_e ("h_vbf_jj_e", "h_vbf_jj_e", 100, 0, 500) ;
   TH1F h_vbf_jj_phi ("h_vbf_jj_phi", "h_vbf_jj_phi", 100, -3, 3) ;
   TH1F h_vbf_jj_m ("h_vbf_jj_m", "h_vbf_jj_m", 100, 0, 2000) ;


   TH1F h_vbf_aj_pt ("h_vbf_aj_pt", "h_vbf_aj_pt", 100, 0, 500) ;
   TH1F h_vbf_aj_eta ("h_vbf_aj_eta", "h_vbf_aj_eta", 100, -5, 5) ;
   TH1F h_vbf_aj_e ("h_vbf_aj_e", "h_vbf_aj_e", 100, 0, 500) ;
   TH1F h_vbf_aj_phi ("h_vbf_aj_phi", "h_vbf_aj_phi", 100, -3, 3) ;

   TH1F h_vbf_bj_pt ("h_vbf_bj_pt", "h_vbf_bj_pt", 100, 0, 500) ;
   TH1F h_vbf_bj_eta ("h_vbf_bj_eta", "h_vbf_bj_eta", 100, -5, 5) ;
   TH1F h_vbf_bj_e ("h_vbf_bj_e", "h_vbf_bj_e", 100, 0, 500) ;
   TH1F h_vbf_bj_phi ("h_vbf_bj_phi", "h_vbf_bj_phi", 100, -3, 3) ;

   TH1F h_vbf_jj_deta ("h_vbf_jj_deta", "h_vbf_jj_deta", 100, -5, 5) ;
   TH1F h_vbf_jj_dphi ("h_vbf_jj_dphi", "h_vbf_jj_dphi", 100, -3, 3) ;
// W JETS
   TH1F h_vbf_wjj_pt ("h_vbf_wjj_pt", "h_vbf_wjj_pt", 100, 0, 500) ;
   TH1F h_vbf_wjj_eta ("h_vbf_wjj_eta", "h_vbf_wjj_eta", 100, -5, 5) ;
   TH1F h_vbf_wjj_e ("h_vbf_wjj_e", "h_vbf_wjj_e", 100, 0, 500) ;
   TH1F h_vbf_wjj_phi ("h_vbf_wjj_phi", "h_vbf_wjj_phi", 100, -3, 3) ;
   TH1F h_vbf_wjj_m ("h_vbf_wjj_m", "h_vbf_wjj_m", 100, 0, 500) ;

   TH1F h_vbf_waj_pt ("h_vbf_waj_pt", "h_vbf_waj_pt", 100, 0, 500) ;
   TH1F h_vbf_waj_eta ("h_vbf_waj_eta", "h_vbf_waj_eta", 100, -5, 5) ;
   TH1F h_vbf_waj_e ("h_vbf_waj_e", "h_vbf_waj_e", 100, 0, 500) ;
   TH1F h_vbf_waj_phi ("h_vbf_waj_phi", "h_vbf_waj_phi", 100, -3, 3) ;

   TH1F h_vbf_wbj_pt ("h_vbf_wbj_pt", "h_vbf_wbj_pt", 100, 0, 500) ;
   TH1F h_vbf_wbj_eta ("h_vbf_wbj_eta", "h_vbf_wbj_eta", 100, -5, 5) ;
   TH1F h_vbf_wbj_e ("h_vbf_wbj_e", "h_vbf_wbj_e", 100, 0, 500) ;
   TH1F h_vbf_wbj_phi ("h_vbf_wbj_phi", "h_vbf_wbj_phi", 100, -3, 3) ;
   
   TH1F h_vbf_lvjj_pt ("h_vbf_lvjj_pt", "h_vbf_lvjj_pt", 100, 0, 500) ;
   TH1F h_vbf_lvjj_eta ("h_vbf_lvjj_eta", "h_vbf_lvjj_eta", 100, -5, 5) ;
   TH1F h_vbf_lvjj_e ("h_vbf_lvjj_e", "h_vbf_lvjj_e", 100, 0, 500) ;
   TH1F h_vbf_lvj_m ("h_vbf_lvj_m", "h_vbf_lvj_m", 100, 0, 500) ;

   TH1F h_vbf_lvjj_phi ("h_vbf_lvjj_phi", "h_vbf_lvjj_phi", 100, -3, 3) ;
   TH1F h_vbf_lvjj_m ("h_vbf_lvjj_m", "h_vbf_lvjj_m", 100, 0, 1000) ;

   TH1F h_vbf_wjj_deta ("h_vbf_wjj_deta", "h_vbf_wjj_deta", 100, -5, 5) ;
   TH1F h_vbf_wjj_dphi ("h_vbf_wjj_dphi", "h_vbf_wjj_dphi", 100, -3, 3) ;

   TH1F h_vbf_l_pt ("h_vbf_l_pt", "h_vbf_l_pt", 100, 0, 500) ;
   TH1F h_vbf_l_eta ("h_vbf_l_eta", "h_vbf_l_eta", 100, -5, 5) ;
   TH1F h_vbf_l_e ("h_vbf_l_e", "h_vbf_l_e", 100, 0, 500) ;
   TH1F h_vbf_l_phi ("h_vbf_l_phi", "h_vbf_l_phi", 100, -3, 3) ;
   TH1F h_vbf_nu_e ("h_vbf_nu_e", "h_vbf_nu_e", 100, 0, 500) ;

   TH1F h_vbf_l_MET_deltaphi ("h_vbf_l_MET_deltaphi", "h_vbf_l_MET_deltaphi", 100, -3, 3) ;
   TH1F h_vbf_lW_hW_deltaphi ("h_vbf_lW_hW_deltaphi", "h_vbf_lW_hW_deltaphi", 100, -3, 3) ;

   TH1F h_vbf_l_tj1_dR ("vbf_l_tj1_dR", "vbf_l_tj1_dR", 100, 0, 6) ;
   TH1F h_vbf_l_tj2_dR ("vbf_l_tj2_dR", "vbf_l_tj2_dR", 100, 0, 6) ;

   TH1F h_vbf_l_wj1_dR ("vbf_l_wj1_dR", "vbf_l_wj1_dR", 100, 0, 6) ;
   TH1F h_vbf_l_wj2_dR ("vbf_l_wj2_dR", "vbf_l_wj2_dR", 100, 0, 6) ;

   TH1F h_vbf_l_wjj_dR ("vbf_l_wjj_dR", "vbf_l_wjj_dR", 100, 0, 6) ;
   TH1F h_vbf_l_tjj_dR ("vbf_l_tjj_dR", "vbf_l_tjj_dR", 100, 0, 6) ;




	
/*	tmva_t->Branch("vbf_event",&vbf_event, "vbf_event/I");
	tmva_t->Branch("vbf_aj_id",&vbf_aj_id, "vbf_aj_id/I");
	tmva_t->Branch("vbf_bj_id",&vbf_bj_id, "vbf_bj_id/I");
	tmva_t->Branch("vbf_waj_id",&vbf_waj_id, "vbf_waj_id/I");
	tmva_t->Branch("vbf_wbj_id",&vbf_wbj_id, "vbf_wbj_id/I");		
*/


  long ieve = 0 ;
  while ( reader.readEvent() ) 
    {
      if (ieve % 1000 == 0) std::cout << "event " << ieve << "\n" ;
      ++ieve;
      
      std::vector<int> bosons ;      
      vector<pair<int, TLorentzVector> > fs_particles ;  //PG final state particles
      //PG loop over particles in the event
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart)
        {
          //PG outgoing particle          
          if (reader.hepeup.ISTUP.at (iPart) == 1)
            {
              TLorentzVector particle 
                (
                  reader.hepeup.PUP.at (iPart).at (0), //PG px
                  reader.hepeup.PUP.at (iPart).at (1), //PG py
                  reader.hepeup.PUP.at (iPart).at (2), //PG pz
                  reader.hepeup.PUP.at (iPart).at (3) //PG E
                ) ;
//              fs_particles[reader.hepeup.IDUP.at (iPart)] = particle ;            
              fs_particles.push_back (pair<int, TLorentzVector> (reader.hepeup.IDUP.at (iPart), particle)) ;
            } //PG outgoing particle
        } //PG loop over particles in the event

      TLorentzVector Z_cand ;
      TLorentzVector H_cand ;
      vector <TLorentzVector *> Jets ;
      vector <TLorentzVector *> leptons ;
      vector <TLorentzVector *> neutrinos ;
      vector <TLorentzVector *> gluons ;

      for (int iPart = 0 ; iPart < fs_particles.size () ; ++iPart)
	{
//	std::cout<<"  Particle Id  "<<fs_particles.at (iPart).first<<" Pt   "<<(fs_particles.at (iPart).second).Pt()<<std::endl;
	}

      for (int iPart = 0 ; iPart < fs_particles.size () ; ++iPart)
        {
          if (abs (fs_particles.at (iPart).first) < 11) //PG tag Jets
            {
              Jets.push_back (&(fs_particles.at (iPart).second)) ;
              continue ;
            }
          if (abs (fs_particles.at (iPart).first) == 21) //PG gluons
            {
              gluons.push_back (&(fs_particles.at (iPart).second)) ;
              continue ;
            }

          if (abs (fs_particles.at (iPart).first) == 12 ||
              abs (fs_particles.at (iPart).first) == 14)
            {
               neutrinos.push_back (&(fs_particles.at (iPart).second)) ;
            }

          if (abs (fs_particles.at (iPart).first) == 11 ||
              abs (fs_particles.at (iPart).first) == 13)
            {
               leptons.push_back (&(fs_particles.at (iPart).second)) ;
               Z_cand += fs_particles.at (iPart).second ; //PG electron or muon
            }
          H_cand += fs_particles.at (iPart).second ;
        }

      //PG VBF cuts

      sort (Jets.begin (), Jets.end (), mySort ()) ;
      sort (leptons.begin (), leptons.end (), mySort ()) ;
      sort (neutrinos.begin (), neutrinos.end (), mySort ()) ;


	 TLorentzVector vbf_ajp(0,0,0,0), vbf_bjp(0,0,0,0);
         TLorentzVector wjj_ajp(0,0,0,0), wjj_bjp(0,0,0,0); 
         TLorentzVector lep(0,0,0,0), nu(0,0,0,0);

         float best_detatagjj = 0; // float best_mtagjj =0;
         float best_mjj = 0; // float best_mjj =0;

   float vbf_jj_e =-999,   vbf_jj_pt =-999,   vbf_jj_eta=-999,  vbf_jj_phi =-999, vbf_jj_m=-999;   
   float vbf_aj_e =-999,   vbf_aj_pt =-999,   vbf_aj_eta=-999,  vbf_aj_phi =-999;   
   float vbf_bj_e =-999,   vbf_bj_pt =-999,   vbf_bj_eta=-999,  vbf_bj_phi =-999;   
   float vbf_jj_deta=-999; Float_t vbf_jj_dphi=-999; 
    int   vbf_jj_type=0,   vbf_n_excj=0,   vbf_n_exfj=0,   vbf_n_gdjj=0;

   float vbf_wjj_e =-999,   vbf_wjj_pt =-999,   vbf_wjj_eta=-999,  vbf_wjj_phi =-999, vbf_wjj_m=-999;   
   float vbf_waj_e =-999,   vbf_waj_pt =-999,   vbf_waj_eta=-999,  vbf_waj_phi =-999;   
   float vbf_wbj_e =-999,   vbf_wbj_pt =-999,   vbf_wbj_eta=-999,  vbf_wbj_phi =-999;   
   float vbf_wjj_deta=-999; Float_t vbf_wjj_dphi=-999;
   float vbf_lvjj_e =-999,   vbf_lvjj_pt =-999,   vbf_lvjj_eta=-999,  vbf_lvjj_phi =-999, vbf_lvjj_m =-999, vbf_lvj_m=-999;
   float vbf_l_pt =-999,   vbf_l_e =-999,   vbf_l_eta=-999,  vbf_l_phi =-999, vbf_nu_e =-999;

   float vbf_l_MET_deltaphi=-999, vbf_lW_hW_deltaphi=-999,vbf_l_tj1_dR=-999,vbf_l_tj2_dR=-999,vbf_l_wj1_dR=-999,vbf_l_wj2_dR=-999,vbf_l_wjj_dR=-999,vbf_l_tjj_dR=-999;



// 	 float jess 1.0; 
         int   n_excj =0, n_exfj = 0, n_gdjj = 0, jj_type = 0, tag_i_id = -1, tag_j_id = -1, wjj_a_id = -1, wjj_b_id = -1;



         for ( size_t i=0; i < leptons.size(); ++i)
         {
        if ((fabs(leptons.at(i)->Eta()) > 2.5) && (leptons.at(i)->Pt() < 25) ) continue;
            TLorentzVector l;
            l.SetPxPyPzE (leptons.at(i)->Px(), leptons.at(i)->Py(), leptons.at(i)->Pz(), leptons.at(i)->E()
              );
	lep=l;
	}

         for ( size_t i=0; i < neutrinos.size(); ++i)
         {
        if (neutrinos.at(i)->Et() < 25) continue;
            TLorentzVector n;
            n.SetPxPyPzE (neutrinos.at(i)->Px(), neutrinos.at(i)->Py(), neutrinos.at(i)->Pz(), neutrinos.at(i)->E()
              );
	nu=n;
        }



       //method1a A two tag Jets from complete eta region
         for ( size_t i=0; i < Jets.size(); ++i)
         {
        if ((fabs(Jets.at(i)->Eta()) > 4.7 )&& (Jets.at(i)->Pt()<20)) continue;
            TLorentzVector i_p;
            i_p.SetPxPyPzE (Jets.at(i)->Px(), Jets.at(i)->Py(), Jets.at(i)->Pz(), Jets.at(i)->E()
              );
            for (size_t j=i+1; j <Jets.size(); ++j)
            {
        if ((fabs(Jets.at(j)->Eta()) > 4.7) && (Jets.at(i)->Pt()<20))continue;
               TLorentzVector j_p;
            j_p.SetPxPyPzE (Jets.at(j)->Px(), Jets.at(j)->Py(), Jets.at(j)->Pz(), Jets.at(j)->E()
              ) ;
               if ( (Jets.at(i)->Eta()*Jets.at(j)->Eta())>0 )  continue;     // 1.  have to be one forward, one backward
               if ( (fabs(Jets.at(i)->Eta()-Jets.at(j)->Eta())<3.5) || (sqrt ((i_p+j_p).M2()) <500)) continue;// 2.Tag pair delta eta>3.5, Mjj>500
                //cout<<"  "<<sqrt ((i_p+j_p).M2())<<endl;
               // if find more than one combinations
               //if ( (fas(i_Eta-j_Eta)>best_detatagjj) )      // 3   Select best combination with maximum deta Eta
               if ( sqrt((i_p+j_p).M2()) > best_mjj )
               {                          // 3   Select best combination with maximum Mjj because of the bad angular resolution in the HF
                  best_detatagjj = fabs(Jets.at(i)->Eta()-Jets.at(j)->Eta());
                  n_gdjj++;
                  best_mjj = sqrt ((i_p+j_p).M2());
                  tag_i_id = i;
                  tag_j_id = j;
                  vbf_ajp = i_p;
                  vbf_bjp = j_p;
               } //loop to find two tag Jets with highest mjj
           } // over over loop over Nmax-1 reco Jets inside 
       } //loop over Nmax reco Jets 
                if (tag_i_id !=-1 && tag_j_id != -1)
                {
                  vbf_jj_e      = (vbf_ajp+ vbf_bjp).E();
                  vbf_jj_pt     = (vbf_ajp+ vbf_bjp).Pt();
                  vbf_jj_eta    = (vbf_ajp+ vbf_bjp).Eta();
                  vbf_jj_phi    = (vbf_ajp+ vbf_bjp).Phi();
                  vbf_jj_m      = best_mjj;
                  vbf_aj_e      = (vbf_ajp).E();
                  vbf_aj_pt     = (vbf_ajp).Pt();
                  vbf_aj_eta    = (vbf_ajp).Eta();
                  vbf_aj_phi    = (vbf_ajp).Phi();
                  //vbf_aj_m      = (i_p).M();
                  vbf_bj_e      = (vbf_bjp).E();
                  vbf_bj_pt     = (vbf_bjp).Pt();
                  vbf_bj_eta    = (vbf_bjp).Eta();
                  vbf_bj_phi    = (vbf_bjp).Phi();
                 // vbf_bj_m      = (j_p).M();

                  vbf_jj_deta   =vbf_aj_eta-vbf_bj_eta;
                  vbf_jj_dphi   = vbf_aj_phi-vbf_bj_phi;
                  //cout<<"  "<<vbf_jj_dphi<<endl;
        } //loop  
        // method1a B

	if (tag_i_id!=-1&&tag_j_id!=-1) 
	{  
	for ( int k=0; k < (int) Jets.size(); ++k) 
	{ 
        if (fabs(Jets.at(k)->Eta()) > 4.7) continue;
	if ( k!=tag_i_id && k!= tag_j_id && wjj_ajp.Pt()!=0 && wjj_bjp.Pt()==0 ) 
	{
	int Bj = k;  
           wjj_bjp.SetPxPyPzE (Jets.at(Bj)->Px(), Jets.at(Bj)->Py(), Jets.at(Bj)->Pz(), Jets.at(Bj)->E()
              );
         wjj_b_id=Bj;                                                                                       
	}
	if ( k!=tag_i_id&&k!=tag_j_id&&wjj_ajp.Pt()==0 && wjj_bjp.Pt()==0 ) 
	{
	int Aj = k;  
           wjj_ajp.SetPxPyPzE (Jets.at(Aj)->Px(), Jets.at(Aj)->Py(), Jets.at(Aj)->Pz(), Jets.at(Aj)->E()
              );
         wjj_a_id=Aj;                                                                                       
	}
	}// loop over reco Jets upto Nmax                                                     
	} // loop of insuring having allready two tagJets

if ( wjj_a_id!=-1 && wjj_b_id!=-1)
{ 
                    //    two W Jets
                        vbf_wjj_e      = (wjj_ajp+wjj_bjp).E();
                        vbf_wjj_pt     = (wjj_ajp+wjj_bjp).Pt();
                        vbf_wjj_eta    = (wjj_ajp+wjj_bjp).Eta();
                        vbf_wjj_phi    = (wjj_ajp+wjj_bjp).Phi();
                        vbf_wjj_m      = sqrt((wjj_ajp+wjj_bjp).M2());
                        vbf_waj_e      = (wjj_ajp).E();
                        vbf_waj_pt     = (wjj_ajp).Pt();
                        vbf_waj_eta    = (wjj_ajp).Eta();
                        vbf_waj_phi    = (wjj_ajp).Phi();
                        //vbf_waj_m      = (wjj_ajp).M();
                        vbf_wbj_e      = (wjj_bjp).E();
                        vbf_wbj_pt     = (wjj_bjp).Pt();
                        vbf_wbj_eta    = (wjj_bjp).Eta();
                        vbf_wbj_phi    = (wjj_bjp).Phi();
                        vbf_wjj_deta=vbf_waj_eta-vbf_wbj_eta;
                        vbf_wjj_dphi= vbf_waj_phi-vbf_wbj_phi;
}
if (tag_i_id!=-1 && tag_j_id!=-1 && wjj_a_id!=-1 && wjj_b_id!=-1 && lep.Pt()>30.)
	{
        		    vbf_lvjj_e      = (lep+nu+wjj_ajp+wjj_bjp).E();
		            vbf_lvjj_pt     = (lep+nu+wjj_ajp+wjj_bjp).Pt();
        		    vbf_lvjj_eta    = (lep+nu+wjj_ajp+wjj_bjp).Eta();
		            vbf_lvjj_phi    = (lep+nu+wjj_ajp+wjj_bjp).Phi();
		            vbf_lvjj_m      = sqrt((lep+nu+wjj_ajp+wjj_bjp).M2());
                            vbf_lvj_m      = (lep+nu).Mt();
//                            vbf_lvj_m      = sqrt((lep+nu).M2());

                            vbf_l_pt = lep.Pt();
                            vbf_l_eta = lep.Eta();
                            vbf_l_phi = lep.Phi();
                            vbf_l_e = lep.E();

                            vbf_nu_e = nu.Et();

                            vbf_l_MET_deltaphi = getDeltaPhi(vbf_l_phi, nu.Phi());
                            vbf_lW_hW_deltaphi = getDeltaPhi((lep+nu).Phi(), vbf_wjj_phi);

                            vbf_l_tj1_dR = delta_R(vbf_l_eta,vbf_l_phi,vbf_aj_eta,vbf_aj_phi);
                            vbf_l_tj2_dR = delta_R(vbf_l_eta,vbf_l_phi,vbf_bj_eta,vbf_bj_phi);

                            vbf_l_wj1_dR = delta_R(vbf_l_eta,vbf_l_phi,vbf_waj_eta,vbf_waj_phi);
                            vbf_l_wj2_dR = delta_R(vbf_l_eta,vbf_l_phi,vbf_wbj_eta,vbf_wbj_phi);

                            vbf_l_tjj_dR = delta_R(vbf_aj_eta,vbf_aj_phi, vbf_bj_eta,vbf_bj_phi);
                            vbf_l_wjj_dR = delta_R(vbf_waj_eta,vbf_waj_phi, vbf_wbj_eta,vbf_wbj_phi);

	}

      double tj_dEta = fabs (Jets.at (0)->Eta () - Jets.at (1)->Eta ()) ; 
      h_deta.Fill (tj_dEta) ; 
      h_j1eta.Fill (Jets.at (0)->Eta ()) ;
      h_j2eta.Fill (Jets.at (1)->Eta ()) ;
      h_j1pt.Fill (Jets.at (0)->Pt ()) ;
      h_j2pt.Fill (Jets.at (1)->Pt ()) ;

        if (tag_i_id!=-1 && tag_j_id!=-1 && wjj_a_id!=-1 && wjj_b_id!=-1 && vbf_aj_pt>30. && vbf_bj_pt>30. && vbf_waj_pt>30. && vbf_wbj_pt>30 && vbf_wjj_m>30. && vbf_lvj_m>30. && vbf_l_pt>30. && vbf_nu_e >30); 
{
	h_vbf_jj_pt.Fill(vbf_jj_pt);
        h_vbf_jj_eta.Fill(vbf_jj_eta);
        h_vbf_jj_phi.Fill(vbf_jj_phi);
        h_vbf_jj_e.Fill(vbf_jj_e);
        h_vbf_jj_m.Fill(vbf_jj_m);


        h_vbf_aj_pt.Fill(vbf_aj_pt);
        h_vbf_aj_eta.Fill(vbf_aj_eta);
        h_vbf_aj_phi.Fill(vbf_aj_phi);
        h_vbf_aj_e.Fill(vbf_aj_e);

        h_vbf_bj_pt.Fill(vbf_bj_pt);
        h_vbf_bj_eta.Fill(vbf_bj_eta);
        h_vbf_bj_phi.Fill(vbf_bj_phi);
        h_vbf_bj_e.Fill(vbf_bj_e);

        h_vbf_jj_deta.Fill(vbf_jj_deta);
        h_vbf_jj_dphi.Fill(vbf_jj_dphi);

        h_vbf_wjj_pt.Fill(vbf_wjj_pt);
        h_vbf_wjj_eta.Fill(vbf_wjj_eta);
        h_vbf_wjj_phi.Fill(vbf_wjj_phi);
        h_vbf_wjj_e.Fill(vbf_wjj_e);
        h_vbf_wjj_m.Fill(vbf_wjj_m);

        h_vbf_waj_pt.Fill(vbf_waj_pt);
        h_vbf_waj_eta.Fill(vbf_waj_eta);
        h_vbf_waj_phi.Fill(vbf_waj_phi);
        h_vbf_waj_e.Fill(vbf_waj_e);

        h_vbf_wbj_pt.Fill(vbf_wbj_pt);
        h_vbf_wbj_eta.Fill(vbf_wbj_eta);
        h_vbf_wbj_phi.Fill(vbf_wbj_phi);
        h_vbf_wbj_e.Fill(vbf_wbj_e);

        h_vbf_wjj_deta.Fill(vbf_wjj_deta);
        h_vbf_wjj_dphi.Fill(vbf_wjj_dphi);

        h_vbf_lvjj_pt.Fill(vbf_lvjj_pt);
        h_vbf_lvjj_eta.Fill(vbf_lvjj_eta);
        h_vbf_lvjj_phi.Fill(vbf_lvjj_phi);
        h_vbf_lvjj_e.Fill(vbf_lvjj_e);
        h_vbf_lvjj_m.Fill(vbf_lvjj_m);
        h_vbf_lvj_m.Fill(vbf_lvj_m);

        h_vbf_l_pt.Fill(vbf_l_pt);
        h_vbf_l_eta.Fill(vbf_l_eta);
        h_vbf_l_phi.Fill(vbf_l_phi);
        h_vbf_l_e.Fill(vbf_l_e);
        h_vbf_nu_e.Fill(vbf_nu_e);

    h_vbf_l_MET_deltaphi.Fill(vbf_l_MET_deltaphi); 
    h_vbf_lW_hW_deltaphi.Fill(vbf_lW_hW_deltaphi); 
    h_vbf_l_tj1_dR.Fill(vbf_l_tj1_dR); 
    h_vbf_l_tj2_dR.Fill(vbf_l_tj2_dR); 
    h_vbf_l_wj1_dR.Fill(vbf_l_wj1_dR); 
    h_vbf_l_wj1_dR.Fill(vbf_l_wj1_dR); 
    h_vbf_l_wjj_dR.Fill(vbf_l_wjj_dR); 
    h_vbf_l_tjj_dR.Fill(vbf_l_tjj_dR); 


}

      for (int k = 0 ; k < Jets.size () ; ++k) h_jeta.Fill (Jets.at (k)->Eta ()) ;
      if (gluons.size () > 0) h_jeta.Fill (gluons.at (0)->Eta ()) ;
        

    } // Now loop over all events

//  TFile histosFile ("cuts_phantom_lheAnalysis_.root","recreate") ;
  TFile histosFile ("unweighted_500Kevents_VBFcuts.root","recreate") ;

//  TFile histosFile ("cuts_Madgraph_lheAnalysis_.root","recreate") ;
//  TFile histosFile ("cuts_vbfnlo_lheAnalysis_.root","recreate") ;


  histosFile.cd () ;
  h_deta.Write () ;
  h_j1eta.Write () ;
  h_j2eta.Write () ;
  h_jeta.Write () ;
  h_j1pt.Write () ;
  h_j2pt.Write () ;

  h_vbf_jj_pt.Write ();
  h_vbf_jj_eta.Write ();
  h_vbf_jj_phi.Write ();
  h_vbf_jj_e.Write ();
  h_vbf_jj_m.Write ();

  h_vbf_aj_pt.Write ();
  h_vbf_aj_eta.Write ();
  h_vbf_aj_phi.Write ();
  h_vbf_aj_e.Write ();

  h_vbf_bj_pt.Write ();
  h_vbf_bj_eta.Write ();
  h_vbf_bj_phi.Write ();
  h_vbf_bj_e.Write ();

  h_vbf_jj_deta.Write ();
  h_vbf_jj_dphi.Write ();

  h_vbf_wjj_pt.Write ();
  h_vbf_wjj_eta.Write ();
  h_vbf_wjj_phi.Write ();
  h_vbf_wjj_e.Write ();
  h_vbf_wjj_m.Write ();

  h_vbf_waj_pt.Write ();
  h_vbf_waj_eta.Write ();
  h_vbf_waj_phi.Write ();
  h_vbf_waj_e.Write ();

  h_vbf_wbj_pt.Write ();
  h_vbf_wbj_eta.Write ();
  h_vbf_wbj_phi.Write ();
  h_vbf_wbj_e.Write ();

  h_vbf_lvjj_pt.Write ();
  h_vbf_lvjj_eta.Write ();
  h_vbf_lvjj_phi.Write ();
  h_vbf_lvjj_e.Write ();
  h_vbf_lvjj_m.Write ();
  h_vbf_lvj_m.Write ();

  h_vbf_l_pt.Write();
  h_vbf_l_eta.Write();
  h_vbf_l_phi.Write();
  h_vbf_l_e.Write();
  h_vbf_nu_e.Write();
  h_vbf_l_MET_deltaphi.Write();
  h_vbf_lW_hW_deltaphi.Write();
  h_vbf_l_tj1_dR.Write();
  h_vbf_l_tj2_dR.Write();
  h_vbf_l_wj1_dR.Write();
  h_vbf_l_wj1_dR.Write();
  h_vbf_l_wjj_dR.Write();
  h_vbf_l_tjj_dR.Write();

  h_vbf_wjj_deta.Write ();
  h_vbf_wjj_dphi.Write ();


  histosFile.Close () ;

  return ;
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


int main (int argc, char **argv) 
{
  ifstream ifs ;

//  ifs.open ("/uscms_data/d3/ajay/VBF_WW_Prod/CMSSW_5_3_2_patch4/src/VBF@NLO/8TeV/Wplus/wwplus-sm-event-unweighted-seed0.lhe") ; 
//  ifs.open ("/uscms_data/d3/ajay/VBF_WW_Prod/CMSSW_5_3_2_patch4/src/total.lhe") ;

	ifs.open ("/uscms_data/d3/ajay/VBF_WW_Prod/CMSSW_5_3_2_patch4/src/LHEactions/Final_lhe/unweighted_500Kevents_VBFcuts.lhe") ;

//  ifs.open ("/uscms_data/d3/ajay/VBF_WW_Prod/CMSSW_5_3_2_patch4/src/LHEactions/Final_lhe/total.lhe") ;
//  ifs.open ("/uscms_data/d3/ajay/VBF_WW_Prod/CMSSW_5_3_2_patch4/src/LHEactions/Final_lhe/madgraph-sm-event-unweighted.lhe") ;
//    ifs.open ("/uscms_data/d3/ajay/VBF_WW_Prod/CMSSW_5_3_2_patch4/src/LHEactions/Final_lhe/wwminusplus-sm-event-unweighted.lhe") ;

	fillHistos (ifs) ; 
	ifs.close () ;
  
  // Now we are done.
  return 0 ;
}
