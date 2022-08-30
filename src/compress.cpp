#include <Rcpp.h>
#include <bitset>
#include <iostream>
#include <fstream>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector compress(const char* infile,
                       const char* outfile,
                       const int nbits,
                       const int tbits,
                       const int mode){
  
  // Expected number of characters
  int echar = 2*nbits;
  int nbytes = (nbits + tbits)/8;
  
  // Build trailing bits
  String trail;
  if(tbits > 0){
    for(int i=0; i<tbits; i++){
      trail += "0"; 
    }
  }
  
  // Open connection
  FILE* infilecon = fopen(infile, "r");
  ofstream outfilecon(outfile, ios::out | ios::binary);
  
  // Print magic number
  if(mode == 1){
    string magic = "00000001";
    bitset<8> bits(magic);
    outfilecon.write((const char*)&bits, 1);
  }else if(mode == 2){
    string magic = "00000010";
    bitset<8> bits(magic);
    outfilecon.write((const char*)&bits, 1);
  }
  
  // Read file one line at a time
  char* line = NULL;
  int nline = 0;
  size_t linesize = 0;
  while((getline(&line, &linesize, infilecon)) != -1){
    int nchar = strlen(line);
    nline++;
    if(nchar == echar | nchar == echar - 1){
      string sline(line);
      sline.erase(remove_if(sline.begin(), sline.end(), ::isspace), sline.end());
      sline += trail;
      for(int i=0; i < nbytes; i++){
        bitset<8> bits(sline, i*8, 8);
        outfilecon.write((const char*)&bits, 1);
      }
    }else{
      fclose(infilecon);
      outfilecon.close();
      cerr << "[ERROR] Expected " << echar << " characters in line " << nline << " but found" << nchar;
      return -1;
    }
    if(nline % 1000 == 0){
      cout << nline << " lines processed\r";
    }
  }
  
  // Close files and exit
  fclose(infilecon);
  outfilecon.close();
  return 0;
  
}
