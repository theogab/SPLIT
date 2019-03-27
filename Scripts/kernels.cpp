#include <Rcpp.h>
using namespace Rcpp;

// Find points closer than eps from P
// [[Rcpp::export]]
IntegerVector epsN(NumericVector LP, double d1){
  IntegerVector epsP = seq(1, LP.size());
  return epsP[LP <= d1];
}

//Concatenate IntegerVectors
// [[Rcpp::export]]
IntegerVector conc(IntegerVector x, IntegerVector y){
  IntegerVector out=no_init(x.size()+y.size());
  std::merge( x.begin(), x.end(), y.begin(), y.end(), out.begin() ) ;
  return out;
}

//Any function
// [[Rcpp::export]]
bool is_any_f(LogicalVector x){
  return is_true(any(x == false));
}

// Expand cluster
// [[Rcpp::export]]
IntegerVector Tdbscan(NumericMatrix D, double d1, double minPts){
  int C = 0;
  int n = D.nrow();
  IntegerVector PC(n);
  LogicalVector PM(n);
  // Clustering procedure
  for(int i = 0; i < n; ++i){
    if(!PM[i]){
      PM[i] = true;
      NumericVector LP = D(i,_);
      IntegerVector NP = epsN(LP, d1) - 1;
      if(NP.size() < minPts){
        //Considered as noise
        PC[i] = 0;
      } else {
        //Define a new cluster
        C++;
        PC[i] = C;
        IntegerVector tovisit = NP;
        IntegerVector seeds = NP;
        //Select points that belong to that cluster
        while(is_any_f(PM[tovisit])){
          //Visiting neighbours
          int m = tovisit.size();
          for(int j = 0; j < m; ++j){
            int Pp = tovisit[j];
            if(!PM[Pp]){
              PM[Pp] = true;
              //Finding neighbours of Pp
              IntegerVector NPp = epsN(D(Pp,_), d1) - 1;
              NP = unique(conc(NP, NPp));
              if(NPp.size() >= minPts){
                seeds = unique(conc(seeds, NPp));
              }
            }
          }
          tovisit = setdiff(NP, tovisit);
        }
        //Rcout << "The cluster number " << C << " has " << NP.size() << " elements"<< std::endl;
        for(int k = 0; k < NP.size(); ++k){
          if(!(PC[NP[k]] > 0)){
            PC[NP[k]] = C;
          }
        }
      }
    }
  }
  return PC;
}


// Allopatry kernel
// [[Rcpp::export]]
Rcpp::IntegerVector AllopatryKernel(Rcpp::NumericMatrix dist1, double ds, Rcpp::IntegerVector Speciesj, int spr){
  if(ds > 0){
    //Selecting speciesj habitat
    int nrow = sum(Speciesj), ncol = sum(Speciesj);
    Rcpp::NumericMatrix dist1j(nrow,ncol);
    IntegerVector idx = seq(1,Speciesj.size()) -1;
    IntegerVector occ = idx[Speciesj ==1];
    for(int i = 0; i < occ.size(); ++i){
      for(int j = 0; j < occ.size(); ++j){
        dist1j(i,j) = dist1(occ[i],occ[j]);
      }
    }
    //Clustering
    IntegerVector clust = Tdbscan(dist1j, ds, 5);
    IntegerVector newsp = Speciesj;
    int j = 0;
    for(int i = 0; i < Speciesj.size(); ++i){
      if(Speciesj[i] == 1) {
        newsp[i] = clust[j];
        j++;
      }
    }
    return newsp;
  }
  else{
    return Speciesj;
  }
} 

// Sympatry kernel
// [[Rcpp::export]]
Rcpp::IntegerVector SympatryKernel(Rcpp::IntegerVector Speciesj, double p, int NbAllo, double spr, Rcpp::IntegerVector Hab){
  //test if there is a sympatric event in each cell of the species range
  IntegerVector newsp = Speciesj;
  int k = 1 + NbAllo;
  for(int i = 0; i < Speciesj.size(); ++i){
    if(Speciesj[i] == 1){
      if(Hab[i] == 3){
        if(rbinom(1, 1, p*spr)[0] == 1){
          newsp[i] = k;
          k++;
        }
      }
      if(Hab[i] == 30){
        if(rbinom(1, 1, p)[0] == 1){
          newsp[i] = k;
          k++;
        }
      }
    }
  }
  return newsp;
}


// Dispersion kernel
// [[Rcpp::export]]
Rcpp::IntegerMatrix DispersionKernel(Rcpp::NumericMatrix dist2, double d, Rcpp::IntegerMatrix MatrixSp, int SpMax, double K, Rcpp::IntegerVector Hab, Rcpp::IntegerVector NewHab, double ad){
  // Size of the p/a matrix
  int nsc = dist2.ncol();
  int nsp = MatrixSp.ncol();
  
  // Remove non occupied cells from input (optimisation)
  LogicalVector occ(MatrixSp.nrow());
  for(int i = 0; i < occ.size(); ++i){
    if(sum(MatrixSp(i,_)) > 0){
      occ[i] = true;
    }
  }
  NumericMatrix dist2j(sum(occ),nsc);
  IntegerMatrix MatrixSpj(sum(occ), nsp);
  IntegerVector Habj(sum(occ));
  int j = 0;
  for(int i = 0; i < occ.size(); ++i){
    if(occ[i]){
      dist2j(j,_) = dist2(i,_);
      MatrixSpj(j,_) = MatrixSp(i,_);
      Habj[j] = Hab[i];
      j++;
    }
  }
  
  // Create the new Matrix
  Rcpp::IntegerMatrix NewMatrixSp(nsc, nsp);
  
  //species dispersion capacity
  NumericVector Draw = rweibull(nsp,1,d);
  
  // evaluating each cell
  for(int i = 0; i < nsc; ++i){
    if(min(dist2j(_,i)) < max(Draw)){
      // Carrying capacity of the cell
      double Kc = SpMax;
      if(NewHab[i] == 3) {Kc = Kc*K;}
      LogicalVector spe(nsp);
      // Cells on the range of the focus cell depending on species dispersion capacity
      LogicalVector RCells = dist2j(_,i) <= max(Draw);
      // nb of evaluated cells
      int j = 0;
      // richness of the cell
      int r = 0;
      // looking at reachable cells
      while(j < sum(RCells) & r < Kc & r < nsp){
        // closest suitable cell
        int n = which_min(dist2j(_,i));
        // species occupying that cell
        LogicalVector speC = MatrixSpj(n,_) == 1;
        // can they colonize that cell
        speC = speC & (Draw >= dist2j(n,i));
        // Is there an adaptation problem?
        if(Habj[n] == 30 & NewHab[i] == 3){
          speC = speC & (rbinom(nsp, 1, ad) == 1);
        }
        spe = spe | speC; 
        // Preparing the next while step
        j++;
        dist2j(n,i) = 1000;
        r = sum(spe);
        // Is there to much species?
        while(r > Kc){
          NumericVector DrawspeC = Draw[speC];
          speC = speC & Draw != max(DrawspeC);
          spe = spe & Draw != max(DrawspeC);
          r = sum(spe);
        }
      }
      IntegerVector pa(nsp);
      pa[spe] = 1;
      NewMatrixSp(i,_) = pa;
    }
  }
  return NewMatrixSp;
}

// Test if the species i trop temp or both
// [[Rcpp::export]]
IntegerVector TestRange(IntegerMatrix MatrixSpecies, IntegerVector Hab){
  int NSp = MatrixSpecies.ncol();
  IntegerVector PHab(2);
  PHab[0] = 30;
  PHab[1] = 3;
  IntegerVector out(NSp);
  for(int i = 0; i < NSp; i++){
    IntegerVector Speciesj = MatrixSpecies(_,i);
    if(sum(Speciesj) != 0){
      IntegerVector Occ = Hab[Speciesj == 1];
      out[i] = sum(match(unique(Occ), PHab));
    }
  }
  return out;
}


// Dispersion and Speciation for all species at time step i
// [[Rcpp::export]]
IntegerMatrix SpDiv(double ds, double d, double p, IntegerMatrix MatrixSpecies, NumericMatrix dist1, NumericMatrix dist2, int SpMax, double K, double spr, double ad, Rcpp::IntegerVector Hab, Rcpp::IntegerVector NewHab){
  
  // Initialisation
  int NSp = MatrixSpecies.ncol();
  int NbNewSp = 0;
  IntegerVector NbAllo(NSp);
  
  // Speciation processes
  for(int i = 0; i < NSp; i++){
    IntegerVector Speciesj = MatrixSpecies(_,i);
    if(sum(Speciesj) != 0){
      // Allopatry
      IntegerVector ClustAllo = AllopatryKernel(dist1, ds, Speciesj, spr);
      //Rcout << "test = " << sum(ClustAllo) << std::endl;
      NbAllo[i] = max(ClustAllo);
      // Sympatry
      IntegerVector ClustSym;
      if(p > 0){
        ClustSym = SympatryKernel(Speciesj, p, NbAllo[i], spr, Hab);
      }
      else{
        ClustSym = ClustAllo;
      }
      MatrixSpecies(_,i) = ClustSym;
      if(sum(ClustSym) != 0){
        NbNewSp += max(ClustSym) - 1;
      }
    }
  }
  
  // Compute Matrix with new species
  IntegerMatrix NewMat(MatrixSpecies.nrow(),NSp + NbNewSp);
  CharacterVector Spec(NSp + NbNewSp);
  int l = NSp;
  
  for(int i = 0; i < NSp; i++){
    IntegerVector Speciesj = MatrixSpecies(_,i);
    if(sum(Speciesj) != 0){
      NewMat(_,i) = (Speciesj == 1);
      Spec(i) = "anc";
      if(max(Speciesj) > 1){
        for(int j = 2; j <= max(Speciesj); j++){
          NewMat(_,l) = (Speciesj == j);
          String ID = "a";
          ID += i + 1;
          ID += "_d";
          ID += j - 1;
          if(j <= NbAllo[i]){
            ID += "_allo";
          }
          else{
            ID += "_sym";
          }
          Spec(l) = ID;
          l++;
        }
      }
    }
  }
  
  // Disperse species in new habitat
  IntegerMatrix Out = DispersionKernel(dist2, d, NewMat, SpMax, K, Hab, NewHab, ad);
  colnames(Out) = Spec;
  return Out;
}


