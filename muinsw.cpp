/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Colvar.h"
#include "ActionRegister.h"
#include "tools/Communicator.h"
#include "tools/Pbc.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <string>
#include <cmath>
#include <cstdlib>

using namespace std;

namespace PLMD{
  namespace colvar{

   
    //+PLUMEDOC COLVAR DELTAU 
    /*
      Calculate the solute and solvent concentration in a planar shell of the box
      !!! The derivative contains only the outer boundary terms, the resulting bias force is non-conservative !!! Only for CmuMD

      \par Examples


      \verbatim
      d: DIPOLE GROUP=1-10
       PRINT FILE=output STRIDE=5 ARG=5
      \endverbatim
      (see also \ref PRINT)

      \attention 
      If the total charge Q of the group in non zero, then a charge Q/N will be subtracted to every atom,
      where N is the number of atoms. This implies that the dipole (which for a charged system depends
      on the position) is computed on the geometric center of the group.


    */
    //+ENDPLUMEDOC
   
    class muinsw : public Colvar {
      //      SwitchingFunction switchingFunction; //instance sw.f
      int gsz,func, Natoms, id, gsz2, Nrot, Nconf;
      double l_func, u,  coff_Du, coff_switch, kT, zdist, eps_Wall, sigma_Wall, r_cut_Wall, wallv;
      double qOO, qHH, qOH, C0, sigmaOO, sigmaHH, sigmaOH, epsilonOO, epsilonHH, epsilonOH, AOO, AHH, AOH, BOO, BHH, BOH;  //Potential parameters
      double R_minOO, R_minHH, R_minOH, R_max, r_cut, r_skin, r0, r_switch;                                                       //Structural parameters
      double max_der, V_minOO, V_minHH, V_minOH, dV_minOO, dV_minHH, dV_minOH;                                                      //Potential parameters
      double rtest, rstep, test_ph, test_ph_qs, test_dph, test_dph_qs;                                                      //Potential testing tools
      double Z0;
      bool quad, nopbcz, gridfile, rigid;
      vector<double> TIP3Pcharge;
      vector<double> TIP3Psigma;
      vector<double> TIP3Pepsilon;
      vector<double> Wallpar;
      vector<double> LJr;
      vector<double> grid_w;
      vector<int> griddim;
      vector<AtomNumber> at_list;
      
      //CHANGE (inserted molecule attributes)
      
      int Nat;
      Vector rcom;
      double Mtot;
      double Maxri;
      vector<double> ins_m;
      vector<Vector> ins_r;
      vector<Vector> ins_s;
      vector<Vector> ins_rot;
      vector<string> sp_ins;
      // verlet list structure
      
      //WIDOM & RANDOM grid
      bool randomgrid,widom;      
      vector<Vector> sgrid;

      struct  dulist_s{
	
	double rcut;
	double rskin;
	int  step;
	vector<int> nn;
	vector<vector<int> > ni;
	vector<Vector> posit;
      } dulist;
      
    public:

      
      muinsw(const ActionOptions&);
      virtual void calculate();
      string to_string(int s);
      double swfon(double z, double Coff);
      double swfoff(double z, double Coff);
      double dswf(double z, double Coff);
      double phi_tip3p(double r, double qq, double Apar, double Bpar);
      double phi_tip3p_qs(double r, double rmin, double rmax, double qq, double Apar, double Bpar, double coff, double Vmin, double dmin, double r_switch);
      double dphi(double rq, double r, double rmin, double rmax, double sigma, double eps);
      double dphiq(double rq, double r, double rmin, double rmax, double sigma, double eps, double dmin);
      double dphi_tip3p(double r, double qq, double Apar, double Bpar);
      double dphi_tip3p_qs(double r, double rmin, double rmax, double qq, double Apar, double Bpar, double coff, double Vmax, double dmin, double r_switch);
      void newlist(vector<Vector> &grid, dulist_s &dulist);
      void checklist(vector<Vector> &grid, dulist_s &dulist);
      void rotate(Vector angle, Vector& vin, Vector& vout);
      
      static void registerKeywords(Keywords& keys);
      //debugfile DBG
      ofstream fdbg;
      ofstream fdbg2;
      ofstream rdbg;
      ifstream inmol;
      ifstream inangles;
    };

    PLUMED_REGISTER_ACTION(muinsw,"MUINSW")

    void muinsw::registerKeywords(Keywords& keys){
      Colvar::registerKeywords(keys);
      keys.add("atoms","GROUP","the group of atoms we are calculating the dipole moment for");
      keys.add("optional","CHARGES","TIP3P charges (default:   qO = -0.834 C\t qH = 0.417 C)");
      keys.add("optional","LJSIGMA","TIP3P LJ sigma parameters (default:   sigmaOO = 0.31507 nm\t sigmaHH = 0.04000 nm\t sigmaOH = 0.17753 nm)");
      keys.add("optional","LJEPSILON","TIP3P LJ epsilon parameters (default:   epsilonOO = 0.63681228 kJ/mol\t epsilonHH = 0.19259280 kJ/mol\t epsilonOH = 0.35001648 kJ/mol)");
      keys.add("optional","LJR","Lennard-Jones radii (default: R_min = 0.9*sigmaLJ \t R_max =  0.18 nm");
      keys.add("optional","R_CUT","Verlet lists, default r_cut = R_max");
      keys.add("optional","R_SWITCH","Distance at which the long-range switching function starts to act up to R_MAX (default: 0.9*R_MAX)");
      keys.add("optional","MAX_DER","Max gradient value determining soft-core threshold (default: max_der = -400 kJ/mol AA");
      keys.addFlag("QUAD",false,"Set to TRUE if for the quadratic soft-core (R_0 will be ignored)");
      keys.add("optional","R_SKIN","Verlet lists, default r_skin = 1.25 R_max");
      keys.add("compulsory","GRID","space grid");
      keys.addFlag("GRIDFILE",false,"the space grid is read from an input file grid.dat");
      keys.add("optional","KT","kT, default 0.73167 kJ/mol =>  88 K");
      keys.add("optional","COFF_DU","cutoff exp(-beta*Du)");
      keys.add("optional","COFF_SWITCH","cutoff switching function");

      //CHANGES
      keys.add("compulsory","INSERT","Molecule file containing the inserted molecule geometry and masses.");
      keys.add("optional","ROTATE","Angle file containing the rotations for the test insertion of the molecule.");
      keys.addFlag("RIGID",false,"The distances in inserted vector are rigid"); 
      keys.addFlag("RANDOMGRID",false,"random uniform grid.");
      keys.addFlag("WIDOM",false,"generates a random uniform grid.");
      keys.remove("NOPBC");
    }
    
    muinsw::muinsw(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao)
    {
      //Read atom groups
      
      parseAtomList("GROUP",at_list);
      
      Natoms=at_list.size();
      //CHANGE add var number of molecules 
      //
      //log check
      log.printf("Number of atoms:\t%d\n",Natoms);

      //GRID DEFINITION (Uniform grid)

      griddim.resize(3);
      parseVector("GRID",griddim);
      
      gsz=1;
      
      for(int i=0; i<3; i++) gsz=gsz*griddim[i];
      grid_w.resize(gsz);
      fill(grid_w.begin(), grid_w.end(),1.); //already divided by the box volume
      
      sgrid.resize(gsz); //normalized grid vector

      //log check
      
      log.printf("gridsize:\t%d\n",gsz);

      //if RANDOM, generate random grid
      parseFlag("RANDOMGRID",randomgrid);
      if(randomgrid){ //random insertion grid
	log.printf("Random distributed insertion points\n");
	double r1,r2,r3;
	for(i=0; i<gsz; i++){
	  r1=static_cast <float> (rand()) / static_cast <float> (RAND_MAX+1.);
	  r2=static_cast <float> (rand()) / static_cast <float> (RAND_MAX+1.);
	  r3=static_cast <float> (rand()) / static_cast <float> (RAND_MAX+1.);
	  sgrid[i]=Vector(r1,r2,r3);
	}
      }else{ //regular grid
	log.printf("Regularly distributed insertion points: %d X %d X %d\n",griddim[0],griddim[1],griddim[2]);
	for(i=0; i<griddim[0]; ++i){
	  for(j=0; j<griddim[1]; ++j){
	    for(k=0; k<griddim[2]; ++k){
	      index=griddim[1]*griddim[2]*i+griddim[2]*j+k; //index
	      sgrid[i]=Vector(di[0]*i,di[1]*j,di[2]*k);
	    }
	  }
	}
      }
      parseFlag("WIDOM",widom);
      if(widom){
	log.printf("Widom flag is active, random insertion grid generated at each step\n");
	log.printf("WARNING: No bias force calculation\n");
      }

      //INSERTION MOLECULE INPUT
      
      string instring;
      parse("INSERT",instring);
      inmol.open(instring.c_str());
      if(!inmol.is_open()){
	log.printf("Insert file not found, ERROR.\n");
      }else{
	Nat=0; //number of atoms equals number of lines
	Mtot=0; //Total mass
	rcom.zero();
	string buff;
	inmol >> Nat;
	getline(inmol,buff); //next line
	getline(inmol,buff); //skip comment line
	while(!inmol.eof()){ 
	  string Symbol;
	  double Mass;
	  Vector r;
	  inmol >> Symbol >> r[0] >> r[1] >> r[2] >> Mass;
	  getline(inmol,buff); 
          r = r*0.1;       //Read from Angstrom coordi
	  if(!inmol.eof()){//fill arrays and compute rcom
	    sp_ins.push_back(Symbol);
	    ins_m.push_back(Mass);
	    ins_r.push_back(r);
	    rcom += Mass*delta(ins_r[0],ins_r[Nat]) ; 
	    Mtot+=Mass;
	    //Nat++;
	  }
	}
	rcom = rcom/Mtot + ins_r[0];
	inmol.close();
      }
      
      //Refer input coordinates to CoM and print log
      log.printf("Inserting %d atom molecule.\n",Nat);
      log.printf("Com in %.4f,%.4f,%.4f and total mass %.4f (uma).\n",rcom[0],rcom[1],rcom[2],Mtot);
      Maxri=0.0;
      for(int i=0; i<Nat; i++){
	log.printf("Atom %d: %s  M = %.4f \tx = %.4f\ty = %.4f\tz = %.4f\n",i,sp_ins[i].c_str(),ins_m[i],ins_r[i][0],ins_r[i][1],ins_r[i][2]);
	ins_r[i] = delta(rcom,ins_r[i]);
	log.printf("In COM reference: \tx = %.4f\ty = %.4f\tz = %.4f\n",ins_r[i][0],ins_r[i][1],ins_r[i][2]);
	//log.printf("Scaled: \tx = %.4f\ty = %.4f\tz = %.4f\n",ins_s[i][0],ins_s[i][1],ins_s[i][2]);
	if(ins_r[i].modulo()>Maxri) Maxri = ins_r[i].modulo(); // update max distance from COM
      }

      log.printf("Max d from COM = %.4f\n",Maxri);

      //Rotations: generate N ins_r[i] vectors, to be used for the insertion in each gridpoint
      //The ins_r[i] vectors can be generated via input, or via random direction generation
      string inang;
      parse("ROTATE",inang);
      inangles.open(inang.c_str());
      Nrot=0;
      if(!inangles.is_open()){
	log.printf("Insertion angles file not found, molecule will not be rotated.\n");
      }else{
	string buff;
	log.printf("The following rotations are tested:\n");
	while(!inangles.eof()){ //count lines
	  Vector theta;
	  inangles >> theta[0] >> theta[1] >> theta[2];
	  getline(inangles,buff); 
	  if(!inangles.eof()){ //generate Nrot-th rotation
	    //collect rotation angle arrays
	    ins_rot.push_back(theta);
	    Nrot++;
	  }
	}
	inangles.close();
      }
      
      //CHANGE automatic generation of rotations & random generation (for WIDOM)
      Nconf=Nrot+1; //original configuration plus rotations
      //ROTATE ins_r
      if(Nrot>0){
	for(int i=0; i<Nrot; i++){
            log.printf("Rotation %d : R_X %.4f\t R_Z %.4f \t R_X %.4f [degrees]\n",i+1,ins_rot[i][0], ins_rot[i][1],ins_rot[i][2]);
            ins_rot[i] = ins_rot[i]*M_PI/180.0; //insertion angles in radiants;
            for(int j=0; j<Nat; j++){
	      Vector vout;
	      rotate(ins_rot[i],ins_r[j],vout);
	      ins_r.push_back(vout);
	      log.printf("\tAtom %d : x = %.4f\ty = %.4f\tz = %.4f\n",j,ins_r[(i+1)*Nat+j][0],ins_r[(i+1)*Nat+j][1],ins_r[(i+1)*Nat+j][2]);
	  }
	}
      }
      
      parseFlag("RIGID",rigid);
      if(rigid){
	log.printf("The inserted molecule geometry is rigid.\n");
      }else{
	log.printf("The inserted molecule geometry fluctuates with the box tensor.\n");
      }

            
      //Coulomb parameters
      TIP3Pcharge.resize(2);
      TIP3Pcharge[0] = -0.834;   //Oxygen charge (C)
      TIP3Pcharge[1] =  0.417;   //Hydrogen charge
      parseVector("CHARGES",TIP3Pcharge);
      C0  = 139.028429916;                  //Electrostatic constant x squared electron charge in kJ/mol*nm --> 332.0637 (kcal/mol*AA) * 4.1868 J->cal * 0.1 AA->nm
      qOO = C0*TIP3Pcharge[0]*TIP3Pcharge[0];  //Product of charges 
      qHH = C0*TIP3Pcharge[1]*TIP3Pcharge[1];  //
      qOH = C0*TIP3Pcharge[0]*TIP3Pcharge[1];  //

      //LJ parameters  
      
      //sigma
      TIP3Psigma.resize(3);
      TIP3Psigma[0] = 0.31507;    //O-O LJ sigma in nm
      TIP3Psigma[1] = 0.04000;    //H-H LJ sigma in nm
      TIP3Psigma[2] = 0.17753;    //O-H LJ sigma in nm
      parseVector("LJSIGMA",TIP3Psigma);
      sigmaOO = TIP3Psigma[0];
      sigmaHH = TIP3Psigma[1];
      sigmaOH = TIP3Psigma[2];

      //epsilon
      TIP3Pepsilon.resize(3);
      TIP3Pepsilon[0] = 0.63681228;  //O-O LJ epsilon in kJ/mol 
      TIP3Pepsilon[1] = 0.19259280;  //H-H LJ epsilon in kJ/mol
      TIP3Pepsilon[2] = 0.35001648;  //O-H LJ epsilon in kJ/mol
      parseVector("LJEPSILON",TIP3Pepsilon);
      epsilonOO = TIP3Pepsilon[0];
      epsilonHH = TIP3Pepsilon[1];
      epsilonOH = TIP3Pepsilon[2];

      //LJ AB parameters
      AOO = 4.0*epsilonOO*pow(sigmaOO,12.0);  // A = 4*eps*sig^12
      AHH = 4.0*epsilonHH*pow(sigmaHH,12.0);
      AOH = 4.0*epsilonOH*pow(sigmaOH,12.0);
   
      BOO = 4.0*epsilonOO*pow(sigmaOO,6.0);   // B = 4*eps*sig^6 
      BHH = 4.0*epsilonHH*pow(sigmaHH,6.0);
      BOH = 4.0*epsilonOH*pow(sigmaOH,6.0);

      R_max = 1.8;   //R_max defined the region for switching to zero the potential tail, default 18 AA
      //log check
      log.printf("Cut-off  radius= %.4f nm\n",R_max);

      //Calculate single atom-atom interaction gradient cut-off
      //Can be better done, this is just for testing
      max_der = -400.0;
      parse("MAX_DER",max_der);

      r_switch=0.9*R_max;
      parse("R_SWITCH",r_switch);
      if(r_switch>=R_max){
        log.printf("Error, R_SWITCH>R_CUT, this cannot be!");
      }


      rtest = R_max;
      rstep = 0.01;    //Just for test purpose
      for(;;){
         rtest -= rstep;
         test_dph = dphi_tip3p(rtest,qOO,AOO,BOO);
         test_ph  =  phi_tip3p(rtest,qOO,AOO,BOO);
         if(test_dph > max_der){
           continue;
         }else{
          R_minOO  = rtest;
          V_minOO  = test_ph;
          dV_minOO = test_dph;
          break;
         }
      }

      rtest = R_max;
      rstep = 0.01;    //Just for test purpose
      for(;;){
         rtest -= rstep;
         test_dph = dphi_tip3p(rtest,qHH,AHH,BHH);
         test_ph  =  phi_tip3p(rtest,qHH,AHH,BHH);
         if(test_dph > max_der){
           continue;
         }else{
          R_minHH  = rtest;
          V_minHH  = test_ph;
          dV_minHH = test_dph;
          break;
         }
      }

      rtest = R_max;
      rstep = 0.01;    //Just for test purpose
      for(;;){
         rtest -= rstep;
         test_dph = dphi_tip3p(rtest,qOH,AOH,BOH);
         test_ph  =  phi_tip3p(rtest,qOH,AOH,BOH);
         if(test_dph > max_der){
           continue;
         }else{
          R_minOH  = rtest;
          V_minOH  = test_ph;
          dV_minOH = test_dph;
          break;
         }
      }

      //exponential cutoff (default = 300.0)
      coff_Du=300.0;
      parse("COFF_DU",coff_Du);

      //log check
      log.printf("Exponential cut-off = %.4f kT\n",coff_Du);

      //Switch cutoff (default = 32.0)
      coff_switch=32.0;
      parse("COFF_SWITCH",coff_switch);
      
      //log check
      log.printf("Switch cut-off = %.4f kT\n",coff_switch);

      //////////////// CHECK POTENTIALS //////////////
      rtest = 0.00;
      rstep = 0.005;
      ofstream potoo;
      potoo.open("PotOO.dat"); //file 
      ofstream pothh;
      pothh.open("PotHH.dat"); //file 
      ofstream potoh;
      potoh.open("PotOH.dat"); //file 
      for(;;){
         rtest += rstep;
         test_ph     = phi_tip3p(rtest,qOO,AOO,BOO);
         test_ph_qs  = phi_tip3p_qs(rtest,R_minOO,R_max,qOO,AOO,BOO,coff_switch,V_minOO,dV_minOO,r_switch);
         test_dph    = dphi_tip3p(rtest,qOO,AOO,BOO);
         test_dph_qs = dphi_tip3p_qs(rtest,R_minOO,R_max,qOO,AOO,BOO,coff_switch,V_minOO,dV_minOO,r_switch);
         potoo << rtest << "\t" << test_ph << "\t" << test_dph << "\t" << test_ph_qs << "\t" << test_dph_qs << endl; //O-O potential

         test_ph     = phi_tip3p(rtest,qHH,AHH,BHH);
         test_ph_qs  = phi_tip3p_qs(rtest,R_minHH,R_max,qHH,AHH,BHH,coff_switch,V_minHH,dV_minHH,r_switch);
         test_dph    = dphi_tip3p(rtest,qHH,AHH,BHH);
         test_dph_qs = dphi_tip3p_qs(rtest,R_minHH,R_max,qHH,AHH,BHH,coff_switch,V_minHH,dV_minHH,r_switch);
         pothh << rtest << "\t" << test_ph << "\t" << test_dph << "\t" << test_ph_qs << "\t" << test_dph_qs << endl; //H-H potential

         test_ph     = phi_tip3p(rtest,qOH,AOH,BOH);
         test_ph_qs  = phi_tip3p_qs(rtest,R_minOH,R_max,qOH,AOH,BOH,coff_switch,V_minOH,dV_minOH,r_switch);
         test_dph    = dphi_tip3p(rtest,qOH,AOH,BOH);
         test_dph_qs = dphi_tip3p_qs(rtest,R_minOH,R_max,qOH,AOH,BOH,coff_switch,V_minOH,dV_minOH,r_switch);
         potoh << rtest << "\t" << test_ph << "\t" << test_dph << "\t" << test_ph_qs << "\t" << test_dph_qs << endl; //O-H potential

         if(rtest > R_max){
           break;
         }else{
          continue;
         }
      }
      potoo.close();
      pothh.close();
      potoh.close();

      /////////////// END CHECK POTENTIALS ////////////
       
      //log check
      log.printf("TIP3P Coulomb parameters:\tqO = %.4f\tqH = %.4f\n",TIP3Pcharge[0],TIP3Pcharge[1]);
      log.printf("TIP3P LJ sigma parameters:\tsigmaOO = %.4f\tsigmaHH = %.4f\tsigmaOH = %.4f\n",sigmaOO,sigmaHH,sigmaOH);
      log.printf("TIP3P LJ epsilon parameters:\tepsilonOO = %.4f\tepsilonHH = %.4f\tepsilonOH = %.4f\n",epsilonOO,epsilonHH,epsilonOH);
      log.printf("SOFT-CORE radii:\tR_min O-O = %.4f\tR_min H-H = %.4f\tR_min O-H = %.4f\n",R_minOO,R_minHH,R_minOH);
      log.printf("LONG RANGE switching:\tR_max = %.4f\n",R_max);
      log.printf("LONG RANGE SWITCH for potential sigmoidal switch:\tR_SWITCH = %.4f\n", r_switch);
      
      
      
      parseFlag("QUAD",quad);
      if(quad){
	log.printf("Quadratic soft-core\n");
      }

      //init Verlet lists
      //CHANGE check plumed 2.0 lists and atoms -> molecules 
      r_cut=R_max + Maxri; //molecule size is included 
      parse("R_CUT",r_cut);
      r_skin=1.25*r_cut; 
      parse("R_SKIN",r_skin);

      
      //log check
      log.printf("Verlet list:\tR_cut = %.4f\tR_skin = %.4f\n",r_cut,r_skin);


      dulist.rcut=r_cut;  
      dulist.rskin=r_skin;  
      dulist.step=0;
      id=0;
      //create dulist.pos, array containing atom positions X
      Vector ze;
      ze.zero();

      //init neigh lists 
            
      dulist.posit.resize(Natoms);
      fill(dulist.posit.begin(), dulist.posit.end(), ze);
      //create dulist.nn array containing the number of neighbors per gridpoint X
      dulist.nn.resize(gsz);
      fill(dulist.nn.begin(), dulist.nn.end(), 0);
      //create dulist.ni array containing the atom index of the neighbors X
      dulist.ni.resize(gsz);
      for (int i = 0; i < gsz; ++i){
	dulist.ni[i].resize(Natoms);
	fill(dulist.ni[i].begin(), dulist.ni[i].end(),0);
      }
      
      //kT [kJ/mol] (default = 2.494 => 300K)
      
      kT=2.494; 
      parse("KT",kT);
      
      checkRead();
      addValueWithDerivatives(); 
      setNotPeriodic();
      log.flush();
      //log atom list
      log.printf("Atom indexes\n");
      for(unsigned int i=0; i<at_list.size(); ++i){
	log.printf("  %d", at_list[i].serial());
      }
      log.printf("  \n");
      requestAtoms(at_list);
      log.printf("Atoms requested  \n");
      
      unsigned int rank; 
      
      rank=comm.Get_rank(); //Rank of pr
      
      //fdbg.open("dbg.dat");
      //fdbg2.open("dbg2.dat");
      //fdbg2/fdbg.flush(); //DBG
    }
    
    void muinsw::rotate(Vector angle, Vector& vin, Vector& vout){
      //sin and cos of proper euler angles
      double c[3];
      double s[3];
      for(int i=0; i<3; i++){
	c[i]=cos(angle[i]);
	s[i]=sin(angle[i]);
      }
      //XZX rotation proper euler angles
      Tensor R = Tensor(c[1],-c[2]*s[1],s[1]*s[2],c[0]*s[1],c[0]*c[1]*c[2]-s[0]*s[2],-c[2]*s[0]-c[0]*c[1]*s[2],s[0]*s[1],c[0]*s[2]+c[1]*c[2]*s[0],c[0]*c[2]-c[1]*s[0]*s[2]);
      vout = matmul(R,vin);
    }
    

    //fermi switchon (auto-width (rmin-r0)/20)
    double muinsw::swfon(double z, double Coff){
      double sig;
      if( z < -Coff){
	sig=0.0;
      }else if(z > Coff){
	sig=1.0;
      }else{
	sig=1.0/(exp(-z)+1.0);
      }
      return(sig);
    }
     
    //fermi switchoff (auto-width (rmin-r0)/20)
    double muinsw::swfoff(double z, double Coff){
      double sig;
      if( z < -Coff){
	sig=0.0;
      }else if(z > Coff){
	sig=1.0;
      }else{
	sig=1.0/(exp(z)+1.0);
      }
      return(sig);
    }
    
    //fermi switchon derivative  
    double muinsw::dswf(double z, double Coff){
      double dsig;
      if(fabs(z) > Coff){
	dsig=0.0;
      }else{
	dsig=0.5/(1.0+cosh(z));
      }
      return(dsig);
    }
    

    //TIP3P Potential w/ quadratic soft-core and sigmoidla long-range switching
    double muinsw::phi_tip3p_qs(double r, double rmin, double rmax, double qq, double Apar, double Bpar, double coff, double Vmin, double dmin, double r_switch){
      double phi, ph, r6, r12;
      if(r <= rmin){
        phi=0.5*dmin*(r*r/rmin-rmin)+Vmin; 
      }else if(r >= r_switch){
        r6  = pow(r,6.0);
        r12 = pow(r,12.0);
        ph  = qq/r + Apar/r12 - Bpar/r6;
        double Z0=0.5*(r_switch+rmax);
        double w=0.05*(rmax-r_switch);                     //One can add a key for the smearing paramter (0.05)
        double z=(r-Z0)/w;
        double soff=swfoff(z,coff);
        phi=ph*soff;
      }else{
        r6  = pow(r,6.0);
        r12 = pow(r,12.0);
        phi = qq/r + Apar/r12 - Bpar/r6;
      }
      return(phi);
    }

    //TIP3P Potential derivative w/ quadratic soft-core and sigmoidla long-range switching
    double muinsw::dphi_tip3p_qs(double r, double rmin, double rmax, double qq, double Apar, double Bpar, double coff, double Vmax, double dmin, double r_switch){
      double dphi, dph, ph, r2, r6, r7, r12, r13;
      if(r <= rmin){
        dphi=dmin*r/rmin;
      }else if(r >= r_switch){
        r2  = pow(r,2.0);
        r6  = pow(r,6.0);
        r7  = pow(r,7.0);
        r12 = pow(r,12.0);
        r13 = pow(r,13.0);
        ph  = qq/r + Apar/r12 - Bpar/r6;
        dph = -qq/r2 - 12.0*Apar/r13 + 6.0*Bpar/r7;
        double Z0=0.5*(r_switch+rmax);
        double w=0.05*(rmax-r_switch);
        double z=(r-Z0)/w;
        double soff=swfoff(z,coff);
        double ds=-dswf(z,coff)/w;
        dphi=ph*ds+dph*soff;
      }else{
        r2  = pow(r,2.0);
        r7  = pow(r,7.0);
        r13 = pow(r,13.0);
        dphi = -qq/r2 - 12.0*Apar/r13 + 6.0*Bpar/r7;
      }
      return(dphi);
    }


    //TIP3P test for cut-off calculation
    double muinsw::phi_tip3p(double r, double qq, double Apar, double Bpar){
      double phi, r6, r12;
      r6  = pow(r,6.0);
      r12 = pow(r,12.0);
      phi = qq/r + Apar/r12 - Bpar/r6;
      return(phi);
    }
                                      

    //TIP3P derivative test for cut-off calculation
    double muinsw::dphi_tip3p(double r, double qq, double Apar, double Bpar){
      double dphi, r2, r7, r13;
      r2  = pow(r,2.0);
      r7  = pow(r,7.0);
      r13 = pow(r,13.0);
      dphi = -qq/r2 - 12.0*Apar/r13 + 6.0*Bpar/r7;
      return(dphi);
    }
    
    
    // verlet list formation
    //CHANGE cycle on molecule number (???)    

    void muinsw::newlist(vector<Vector> &grid, dulist_s &dulist)
    {
      int i,j,gridsize,listsize;
      Vector rij, test;
      double mod_rij;
      
      gridsize=grid.size(); //get size of the grid
      listsize=at_list.size(); //get the size of the list

      //fill dulist.pos with new list positions

      for(j=0; j<listsize; ++j){
	dulist.posit[j] = getPosition(j);
      }
      
      unsigned int stride;
      unsigned int rank; 
      
      stride=comm.Get_size();  //Number of processes
      rank=comm.Get_rank(); //Rank of pr

      for(j=rank; j<gridsize; j+=stride) {     // sum over grid
	dulist.nn[j]=0; //reset nlistsize
	for(i=0; i<listsize; ++i) {                                           // sum over atoms
	  rij = pbcDistance(dulist.posit[i],grid[j]);
	  mod_rij=rij.modulo();
	  if (mod_rij < dulist.rskin){ //if distance < rskin
	    dulist.ni[j][dulist.nn[j]]=i; //index in neighlist
	    dulist.nn[j]++; //increment nn 
	  }	  
	}
      }
    }
    
    //CHANGE cycle on molecule number (???)
    void muinsw::checklist(vector<Vector> &grid, dulist_s &dulist)
    {
      unsigned int j; 
      double dr=(dulist.rskin-dulist.rcut)*0.5;
      Vector rij;
        
      for (j=0; j<at_list.size(); ++j) { 
	rij = pbcDistance(getPosition(j),dulist.posit[j]); //check position variations
       	if( fabs(rij[0])>dr || fabs(rij[1])>dr || fabs(rij[2])>dr ) { 
	  id++; //rebuild list counter
	  newlist(grid,dulist); 
	  break; 
	}
      }
    } 
    
    // calculator
    void muinsw::calculate()    
    {
      int i,j,k,c;
      //Parallel parameters
      
      unsigned int stride;
      unsigned int rank; 
      
      stride=comm.Get_size();  //Number of processes
      rank=comm.Get_rank(); //Rank of pr

      //Box size
      
      vector<double> LBC(3);
      for(i=0;i<3;++i) LBC[i]=getBox()[i][i];
            
      double Vbox=getBox().determinant(); //box volume
      //Space grid
      
      Vector di;

      for(i=0;i<3;++i){
	di[i]=1./(1.*griddim[i]); //grid2 pace
      }
      
      if(dulist.step==0 && !rigid){ //pbc on insertion vectors
	for(i=0;i<Nconf*Nat;i++){
	  ins_s.push_back(getPbc().realToScaled(ins_r[i])); //scaled points
	  //fdbg<< i << "\t" << ins_s[i][0] <<"\t"<< ins_s[i][1]<<"\t"<< ins_s[i][2]<<endl;
	  //for(j=0;j<3;j++) ins_s[i][j]=Tools::pbc(ins_s[i][j]); //within PBC
	  //fdbg<< i << "\t" << ins_s[i][0] <<"\t"<< ins_s[i][1]<<"\t"<< ins_s[i][2]<<endl;
	  //ins_r[i]=getPbc().scaledToReal(ins_s[i]); //Unscale
	  //fdbg<< i << "\t" << ins_r[i][0] <<"\t"<< ins_r[i][1]<<"\t"<< ins_r[i][2]<<endl;
	}
      }

      vector<Vector> xgrid(gsz);
      Vector gs;
      int index;

      
      //derivative vector
      vector<Vector> deriv(getNumberOfAtoms());
      
      //virial contribution
      Tensor zerot;
      zerot.zero();
      Tensor virial;
      virial.zero();
      Vector ze;
      ze.zero();
      //init derivatives
      fill(deriv.begin(), deriv.end(), ze);

      //gridpoints (random or regular)
      
      if(widom){ //generate random insertion grid every step
	double r1,r2,r3;
	for(i=0; i<gsz; i++){
	  r1=static_cast <float> (rand()) / static_cast <float> (RAND_MAX+1.);
	  r2=static_cast <float> (rand()) / static_cast <float> (RAND_MAX+1.);
	  r3=static_cast <float> (rand()) / static_cast <float> (RAND_MAX+1.);
	  sgrid[i]=Vector(r1,r2,r3);
	}
      }


      for(i=0; i<gsz; ++i) xgrid[i]=matmul(transpose(getBox()),sgrid0[i]);

     
      //create neigh lists  
      
      if(dulist.step==0){
	newlist(xgrid,dulist);
      }else{
	if(widom){ //no checklist for widom
	  newlist(xgrid,dulist);
	}else{
	  checklist(xgrid,dulist);
	}
      }
      ++dulist.step; 

      
      double Su=0;
      double Du,modrkj,modq,zu,fDu,dfDu, ph;
      Vector rkj;
      int atomid;
      string sp_real;
      double beta=1./kT;
      
      double invM=1./gsz/Nconf;  //Number of insertions

      //cycle on space grid

      for(i=rank; i<gsz; i+=stride){
	//auxiliary array of vectors for derivatives (the contribution for each ins configuration is separated)
	vector<Vector> phiprime(dulist.nn[i]);
	//auxiliary array of tensors for derivatives  (the contribution for each ins configuration is separated)
	vector<Tensor> phit(dulist.nn[i]);
	
	for(c=0; c<Nconf; c++){ //loop on configurations (different Du values)      
	  
	  //initialize arrays and Du
	  fill(phiprime.begin(), phiprime.end(), ze);
	  fill(phit.begin(), phit.end(), zerot);
	  Du=0; 	  

	  
	  //loop on atoms of the inserted test-molecule (contributing to the same Du)
	  for(k=0; k<Nat; ++k){
	    if(beta*Du>coff_Du)   break; //speed up calculation (it might introduce an error)
	    
	    Vector insx;
	    if(rigid) {
	      insx=xgrid[i]+ins_r[k+c*Nat]; //insertion vector (rigid)
	    }else{
	      insx=xgrid[i]+getPbc().scaledToReal(ins_s[k+c*Nat]); //insertion vector (scaled)
	    }
	    
	    //fdbg2<<setprecision(10);
	    //fdbg2<< insx[0]<<"\t"<< insx[1]<<"\t"<< insx[2]<<endl;
	    
	    
	    //loop over neighbor atoms
	    
	    for(j=0; j<dulist.nn[i]; ++j){
	      if(beta*Du>coff_Du)   break;//speed up calculation (it might introduce an error)
	      
	      atomid = dulist.ni[i][j]; //atomindex
	      //Delta U

	      rkj = pbcDistance(insx,getPosition(atomid)); 
	      //fdbg setprecision(8) << modrij << "\t" << scientific << ph << "\t" << dph << endl; //DBG
	      modrkj=rkj.modulo();
	      modq=modrkj*modrkj;
	      
	      //CHANGE once k and j are associated to atoms, set the sigma,epsilon and Coulombic parameter
	      double mass_real;
	      double ph,dph;
              mass_real = getMass(atomid); 
              if(mass_real >= 10.){
                sp_real = "O";
              }else{
                sp_real = "H";
              }

	      if(sp_real != sp_ins[k]){
                if(quad){
                  ph  =  phi_tip3p_qs(modrkj, R_minOH, R_max, qOH, AOH, BOH, coff_switch, V_minOH, dV_minOH, r_switch);
                  dph = dphi_tip3p_qs(modrkj, R_minOH, R_max, qOH, AOH, BOH, coff_switch, V_minOH, dV_minOH, r_switch);
                }else{
                  ph  =  phi_tip3p(modrkj, qOH, AOH, BOH);
                  dph = dphi_tip3p(modrkj, qOH, AOH, BOH);
                }
              }else{
                   if(sp_real == "O"){
                     if(quad){
                       ph  =  phi_tip3p_qs(modrkj, R_minOO, R_max, qOO, AOO, BOO, coff_switch, V_minOO, dV_minOO, r_switch);
                       dph = dphi_tip3p_qs(modrkj, R_minOO, R_max, qOO, AOO, BOO, coff_switch, V_minOO, dV_minOO, r_switch);
                     }else{
                       ph  =  phi_tip3p(modrkj, qOO, AOO, BOO);
                       dph = dphi_tip3p(modrkj, qOO, AOO, BOO);
                     }
                   }else{
                     if(quad){
                       ph  =  phi_tip3p_qs(modrkj, R_minHH, R_max, qHH, AHH, BHH, coff_switch, V_minHH, dV_minHH, r_switch);
                       dph = dphi_tip3p_qs(modrkj, R_minHH, R_max, qHH, AHH, BHH, coff_switch, V_minHH, dV_minHH, r_switch);
                     }else{
                       ph  =  phi_tip3p(modrkj, qHH, AHH, BHH);
                       dph = dphi_tip3p(modrkj, qHH, AHH, BHH);
                     }

                   }
              }
	      
	      Du+=ph; //update potential difference
	      Vector phiprimej=rkj*dph/modrkj;

	      phiprime[j] += phiprimej; //potential derivative (k-th atom contribution on j)

	      if(rigid){
		rkj = rkj+ins_r[k+c*Nat]; //refer to CoM
	      }
	      phit[j] += Tensor(rkj,phiprimej); //Virial component (k-th atom contribution on j)
	      //fdbg << i << "\t" << c << "\t" << k <<"\t"<<  j <<"\t"<< setprecision(10) << modrkj << "\t" << ph << "\t" << Du << endl;
	      //if(atomid==302) fdbg << i << "\t" << k << "\t"<<  j << "\t" << setprecision(10) << modrkj << "\t" << delta(insx,getPosition(atomid)).modulo() << "\t" <<  ph << "\t" << dph << endl;

	    }
	  }
	  
	  //exponential (with cut off at coff_Du)
	  if(beta*Du>coff_Du){
	    fDu=0.0;
	    dfDu=0.0;
	  }else{
	    fDu=exp(-beta*(Du));
	    dfDu=-beta*fDu;
	  }
	
	  Su+=fDu*invM*grid_w[i]; //evaluate integral
	  double dSoV=dfDu*invM*grid_w[i];
	
	  //derivatives -- cycle again over neighbors
	  //CHANGE cycle on all atoms of the neighborhood list
	  if(!widom){ //derivatives are not used
	    for(j=0; j < dulist.nn[i]; ++j){
	      atomid = dulist.ni[i][j]; //atomindex
	      deriv[atomid]+=dSoV*phiprime[j];
	      virial-=dSoV*phit[j]; //Tensor(phiprime[j],rkj); //Virial component
	    }
	  }
	}
	vector<Vector>().swap(phiprime);
	vector<Tensor>().swap(phit);
      }
      
      comm.Sum(Su);
      comm.Sum(deriv);
      comm.Sum(virial);
      
      //Logaritm of the sum
      if(Su==0){
	setValue(kT*coff_Du); //When the integral is 0 mu=coff_Du* kT
      }else{
	setValue(-kT*(std::log(Su)));
      }
      for(j=0; j < Natoms; ++j)	{
	setAtomsDerivatives(j, -kT/Su*deriv[j]);
      }
      
      //Add total volume derivative
      
      setBoxDerivatives(-kT/Su*virial);//+SuT);
    }
  }
}
