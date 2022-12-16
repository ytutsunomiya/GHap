#include <Rcpp.h>
#include <bitset>
#include <iostream>
#include <fstream>
#include <math.h>

// Function: sliceCpp
// License: GPLv3 or later
// Modification date: 16 Dec 2022
// Written by: Yuri Tani Utsunomiya, Adam Taiti Harth Utsunomiya
// Contact: ytutsunomiya@gmail.com, adamtaiti@gmail.com
// Description: slice GHap object

using namespace Rcpp;
using namespace std;

// [[Rcpp::export(.sliceCpp)]]
IntegerVector sliceCpp(const char* binfile,
                       const int mode,
                       const int nvars,
                       const int nids,
                       const int phased,
                       const int imp,
                       IntegerVector iidx,
                       IntegerVector hidx,
                       IntegerVector vidx){
  
  // Control variables
  const int nids_in = iidx.length();
  const int nvars_in = vidx.length();
  
  // Output vector
  IntegerVector geno(phased*nids_in*nvars_in);
  
  // Open file connection
  ifstream binfilecon(binfile, ifstream::binary);
  
  // If connection is good
  
  if(binfilecon.is_open()){
    
    // Phase: variant by individual mode
    if(mode == 0 | mode == 1){
      
      int tbits = 2*nids % 8;
      int nbytes = (2*nids + tbits)/8;
      
      // Loop by variant
      for(int i = 0; i < nvars_in; i++){
        
        // Initialize haplotypes
        IntegerVector hap(nbytes*8);

        // Position pointer in binary file
        int start = (vidx[i]-1)*nbytes + mode;
        binfilecon.seekg(start, binfilecon.beg);
        
        // Read individual into memory
        char* buffer = new char [nbytes];
        binfilecon.read(buffer, nbytes);
        bitset<8> bits;
        int a = 0;
        for(int b = 0; b < nbytes; b++){
          bits = buffer[b];
          hap[a] = bits[7];
          hap[a + 1] = bits[6];
          hap[a + 2] = bits[5];
          hap[a + 3] = bits[4];
          hap[a + 4] = bits[3];
          hap[a + 5] = bits[2];
          hap[a + 6] = bits[1];
          hap[a + 7] = bits[0];
          a = a + 8;
        }
        hap.erase(2*nids,hap.length());
        hap = hap[(hidx-1)];
        geno[Rcpp::seq(i*2*nids_in, 2*nids_in*(i + 1) - 1)] = hap;
      }
    }
    
    // Phase: individual by variant mode
    if(mode == 2){
      
      int tbits = 2*nvars % 8;
      int nbytes = (2*nvars + tbits)/8;
      
      // Loop by individual
      for(int i = 0; i < nids_in; i++){
        
        // Initialize haplotypes
        IntegerVector hap1(nbytes*4);
        IntegerVector hap2(nbytes*4);
        
        // Position pointer in binary file
        int start = (iidx[i]-1)*nbytes + 1;
        binfilecon.seekg(start, binfilecon.beg);
        
        // Read individual into memory
        char* buffer = new char [nbytes];
        binfilecon.read(buffer, nbytes);
        bitset<8> bits;
        int a = 0;
        for(int b = 0; b < nbytes; b++){
          bits = buffer[b];
          hap1[a] = bits[7];
          hap2[a] = bits[6];
          hap1[a+1] = bits[5];
          hap2[a+1] = bits[4];
          hap1[a+2] = bits[3];
          hap2[a+2] = bits[2];
          hap1[a+3] = bits[1];
          hap2[a+3] = bits[0];
          a = a + 4;
        }
        hap1.erase(nvars,hap1.length());
        hap2.erase(nvars,hap2.length());
        hap1 = hap1[(vidx-1)];
        hap2 = hap2[(vidx-1)];
        geno[Rcpp::seq(2*i*nvars_in, nvars_in*(2*i + 1) - 1)] = hap1;
        geno[Rcpp::seq(nvars_in*(2*i + 1), 2*nvars_in*(i + 1) - 1)] = hap2;
      }
    }
    
    // Plink or Haplo
    if(mode == 3){
      
      int tbits = 2*nids % 8;
      int nbytes = (2*nids + tbits)/8;
      int miss = 0;
      if(imp == 0){
        miss = NA_INTEGER;
      }
      // IntegerVector code = IntegerVector::create(Named("11",0), Named("01",1),
      //                                            Named("00",2), Named("10",miss));
      unordered_map<std::string, int>code;
      code.insert({"11", 0});code.insert({"01", 1});
      code.insert({"00", 2});code.insert({"10", miss});
      
      // Loop by variant
      for(int i = 0; i < nvars_in; i++){
        
        // Initialize genotypes
        IntegerVector gen(nbytes*4);
        
        // Position pointer in binary file
        int start = (vidx[i]-1)*nbytes + mode;
        binfilecon.seekg(start, binfilecon.beg);
        
        // Read individual into memory
        char* buffer = new char [nbytes];
        binfilecon.read(buffer, nbytes);
        bitset<8> bits;
        int a = 0;
        for(int b = 0; b < nbytes; b++){
          bits = buffer[b];
          // gen[a] = code[to_string(bits[0]) + to_string(bits[1])];
          // gen[a + 1] = code[to_string(bits[2]) + to_string(bits[3])];
          // gen[a + 2] = code[to_string(bits[4]) + to_string(bits[5])];
          // gen[a + 3] = code[to_string(bits[6]) + to_string(bits[7])];
          gen[a] = code.at(to_string(bits[0]) + to_string(bits[1]));
          gen[a + 1] = code.at(to_string(bits[2]) + to_string(bits[3]));
          gen[a + 2] = code.at(to_string(bits[4]) + to_string(bits[5]));
          gen[a + 3] = code.at(to_string(bits[6]) + to_string(bits[7]));
          a = a + 4;
          // for(size_t d = 0; d < bits.size(); d = d + 2){
          //   gen[a] = code.at(to_string(bits[d]) + to_string(bits[d+1]));
          //   a++;
          // }
        }
        gen.erase(nids,gen.length());
        gen = gen[(hidx-1)];
        geno[Rcpp::seq(i*nids_in, nids_in*(i + 1) - 1)] = gen;
      }
    }
    
  }
  
  // Close file connection
  binfilecon.close();
  
  // Return output
  return geno;
}
