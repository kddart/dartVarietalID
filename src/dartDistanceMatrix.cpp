#include <Rcpp.h>
#include <string.h>
using namespace Rcpp;
#include <stdint.h>
#ifdef CPU_INTRINSICS
#include <nmmintrin.h>
#endif
#include <algorithm>
#include <vector>
#include <thread>


#define TWO(c)   ((uint64_t)0x1u << (c))
#define MASK(c)  (((uint64_t)(-1)) / (TWO(TWO(c)) + 1u))
#define COUNT(x,c) ((x) & MASK(c)) + (((x) >> (TWO(c))) & MASK(c))

static inline uint32_t bitcountParallel( uint32_t nWord )
{
  nWord = COUNT(nWord, 0) ;
  nWord = COUNT(nWord, 1) ;
  nWord = COUNT(nWord, 2) ;
  nWord = COUNT(nWord, 3) ;
  nWord = COUNT(nWord, 4) ;
  //n = COUNT(n, 5) ;   // for 64-bit integers (DAVE: ???)
  return nWord;
}

static inline uint32_t bitcountParallel ( uint64_t nWord )
{
  nWord = COUNT(nWord, 0) ;
  nWord = COUNT(nWord, 1) ;
  nWord = COUNT(nWord, 2) ;
  nWord = COUNT(nWord, 3) ;
  nWord = COUNT(nWord, 4) ;
  nWord = COUNT(nWord, 5) ;   // for 64-bit integers
  return (uint32_t)nWord ;
}


template< class T_SUM, class T_WORDS >
static inline T_SUM sumBitcountAnd( T_WORDS* pn1, T_WORDS* pn2, T_SUM nWords )
{
  T_SUM nSum = 0;
  for( int nWord = 0; nWord < nWords; ++nWord )
  {
#ifndef CPU_INTRINSICS
    nSum += bitcountParallel( (*pn1) & (*pn2) );
#else
#if defined(_M_X64) || defined(__x86_64__) || defined(__x86_64)
    nSum += (uint32_t)_mm_popcnt_u64( (*pn1) & (*pn2) );
#else
    nSum += bitcountParallel( (*pn1) & (*pn2) );
#endif
#endif
    ++pn1;
    ++pn2;
  }
  return nSum;
}

static std::string getBitcountConfiguration(){
#ifndef CPU_INTRINSICS
  return "Not using bitcounting intrinsics.";
#elif defined(_M_X64) || defined(__x86_64__) || defined(__x86_64)
  return "Using bitcounting intrinsics (_mm_popcnt_u64).";
#else
  return "Bitcount intrinsics not supported.";
#endif
}


double subsetFunc(int p, int P, int t){
  double a = 1.0;
  double b = -2.0*t;
  double c = (double)p*t*t/(double)P;
  return ( -1.0 * b - std::sqrt( b*b - 4.0 * a * c ) ) / ( 2.0 * a );
}


std::vector<int> getMarker1SubsetIdx(int nRows,int parts){
  parts = std::min(nRows - 1,parts);
  std::vector<int> subsets;
  subsets.push_back(0);
  if(parts>1){
    for(int i=1;i<parts;i++){
      subsets.push_back((int)std::round(subsetFunc(i,parts,nRows)));
    }
  }
  subsets.push_back(nRows);
  return subsets;
}

#define FIDDLE_FACTOR 0.000001

class PackedGenotypes{
public:
  uint64_t** bA;
  uint64_t** b0;
  uint64_t** b1;
  uint64_t** b2;
  int nrow;
  int ncol;
  int ncolPacked;

  PackedGenotypes(const IntegerMatrix &genotypes, int nCores, bool isSampleWise){

    int lncol;
    int lncolPacked;
    std::function<int(int&,int&)> getGeno;
    if(isSampleWise){
      this->nrow = genotypes.ncol();
      lncol = this->ncol = genotypes.nrow();
      lncolPacked = this->ncolPacked = (int)std::ceil(genotypes.nrow() / 64.0);
      getGeno = [genotypes](int& r, int& c) { return genotypes.at(c,r); };
    }else{
      this->nrow = genotypes.nrow();
      lncol = this->ncol = genotypes.ncol();
      lncolPacked = this->ncolPacked = (int)std::ceil(genotypes.ncol() / 64.0);
      getGeno = [genotypes](int& r, int& c) { return genotypes.at(r,c); };
    }


    uint64_t** lb0 = this->b0 = new uint64_t*[this->nrow];
    uint64_t** lb1 = this->b1 = new uint64_t*[this->nrow];
    uint64_t** lb2 = this->b2 = new uint64_t*[this->nrow];
    uint64_t** lbA = this->bA = new uint64_t*[this->nrow];
    for(int r = 0; r < this->nrow; r++){
      lb0[r] = new uint64_t[lncolPacked];
      lb1[r] = new uint64_t[lncolPacked];
      lb2[r] = new uint64_t[lncolPacked];
      lbA[r] = new uint64_t[lncolPacked];
    }
    std::vector<int> primaryDimIdx = getMarker1SubsetIdx(this->nrow,nCores);
    std::vector<std::thread> workers;
    for (int i = 0; i < (int)primaryDimIdx.size() - 1; i++) {
      int lower = primaryDimIdx.at(i);
      int upper = primaryDimIdx.at(i+1);
      workers.push_back(std::thread([lower,upper, genotypes, lb0, lb1, lb2, lbA, lncolPacked, lncol, getGeno] () {
        for(int r = lower; r < upper; r++){
          uint64_t* b0p = lb0[r];
          uint64_t* b1p = lb1[r];
          uint64_t* b2p = lb2[r];
          uint64_t* bAp = lbA[r];
          for(int cp=0; cp < lncolPacked; cp ++){
            b0p[cp] = 0;
            b1p[cp] = 0;
            b2p[cp] = 0;
            bAp[cp] = 0;
          }
          for(int c = 0;c < lncol; c++){
            int g = getGeno(r,c);
            if(g!=NA_INTEGER){
              uint64_t shift = ((uint64_t)1 << (c % 64)) ;
              switch(g) {
              case 0:
                b0p[(int)(c/64)] += shift;
                break;
              case 1:
                b1p[(int)(c/64)] += shift;
                break;
              case 2:
                b2p[(int)(c/64)] += shift;
                break;
              default :
                std::string str = "Unknown genotype value " + std::to_string(g) + " on row " + std::to_string(r)+", column " + std::to_string(c);
              ::Rf_error(str.c_str());
              break;
              }
              bAp[(int)(c/64)] += shift * (uint64_t)1;
            }
          }
        }
      }));
    }
    std::for_each(workers.begin(), workers.end(), [](std::thread& t) {
      t.join();
    });
  }
};



// [[Rcpp::export]]
NumericMatrix dartDistanceMatrixCpp(
    const IntegerMatrix & snpData,
    const LogicalVector & sampleWiseAnalysis = false,
    const IntegerVector & nThreads = 1) {

  Rcpp::Rcout << getBitcountConfiguration() << std::endl;

// #if defined(__x86_64__) || defined(__x86_64)
  bool isSampleWise = sampleWiseAnalysis[0];

  int primaryDim = snpData.nrow();
  if(isSampleWise){
    primaryDim = snpData.ncol();
  }

  int nCores = std::min(nThreads.at(0),primaryDim - 1);
  PackedGenotypes snpDataPacked(snpData,nCores,isSampleWise );
  NumericMatrix distances(primaryDim,primaryDim);
  if(isSampleWise){
    rownames(distances) = colnames(distances) = colnames(snpData);
  }else{
    rownames(distances) = colnames(distances) = rownames(snpData);
  }

  // vector container stores threads
  std::vector<int> primaryDimIdx = getMarker1SubsetIdx(primaryDim,nCores);
  std::vector<std::thread> workers;
  for (int c = 0; c < (int)primaryDimIdx.size() - 1; c++) {
    int lower = primaryDimIdx.at(c);
    int upper = primaryDimIdx.at(c+1);
    workers.push_back(std::thread([lower,upper, &distances, snpDataPacked, primaryDim] () {
        uint64_t** bA = snpDataPacked.bA;
      uint64_t** b0 = snpDataPacked.b0;
      uint64_t** b1 = snpDataPacked.b1;
      uint64_t** b2 = snpDataPacked.b2;
      int ncolPacked = snpDataPacked.ncolPacked;
      for(int m=lower;m<upper;m++){
        distances(m,m) = 0.0f;
        for(int m2=m+1;m2<primaryDim;m2++){
          int presentMatchs = sumBitcountAnd(bA[m],bA[m2], ncolPacked);
          float score = NA_REAL;
          if(presentMatchs > 0){
            int refAndHetMatchs = sumBitcountAnd(b0[m],b2[m2], ncolPacked);
            int hetAndRefMatchs = sumBitcountAnd(b2[m],b0[m2], ncolPacked);

            int snpAndHetMatchs = sumBitcountAnd(b0[m],b2[m2], ncolPacked);
            int hetAndSnpMatchs = sumBitcountAnd(b2[m],b0[m2], ncolPacked);

            int refAndSnpMatches =  sumBitcountAnd(b0[m],b1[m2], ncolPacked);
            int snpAndRefMatches =  sumBitcountAnd(b1[m],b0[m2], ncolPacked);

            int fullScores = refAndSnpMatches + snpAndRefMatches;
            int halfScores = refAndHetMatchs + hetAndRefMatchs + snpAndHetMatchs + hetAndSnpMatchs;

            score = ( (float)fullScores + (float)halfScores / 2.0f ) / (float)presentMatchs ;
          }
          distances(m,m2) = score;
          distances(m2,m) = score;
        }
      }
    }));
  }
  std::for_each(workers.begin(), workers.end(), [](std::thread& t) {
    t.join();
  });
  return distances;
// #else
//   ::Rf_error("OS not supported");
// #endif
}
