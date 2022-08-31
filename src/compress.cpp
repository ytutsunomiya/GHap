#include <Rcpp.h>
#include <bitset>
#include <iostream>
#include <fstream>

// Function: compress
// License: GPLv3 or later
// Modification date: 31 Aug 2022
// Written by: Yuri Tani Utsunomiya, Adam Taiti Harth Utsunomiya
// Contact: ytutsunomiya@gmail.com, adamtaiti@gmail.com
// Description: Auxiliary function to compress phased data into GHap binary

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
int compress(const char* infile,
             const char* outfile,
             const int nunits,
             const int nbits,
             const int tbits,
             const int fmode){
  
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
  int exitcode = 0;
  
  // Print magic number
  string unit = "markers";
  if(fmode == 1){
    unit = "markers";
    string magic = "00000001";
    bitset<8> bits(magic);
    outfilecon.write((const char*)&bits, 1);
  }else if(fmode == 2){
    unit = "individuals";
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
      cerr << "\n[ERROR] Expected " << echar << " characters in line " << nline << " but found " << nchar << endl;
      exitcode = 1;
      goto exit;
    }
    if(nline % 1000 == 0){
      cout << nline << " " << unit << " written to file\r";
    }
    if(nline > nunits){
      cerr << "\n[ERROR] Expected " << nunits << " lines in the phase file but found more!" << endl;
      exitcode = 1;
      goto exit;
    }
  }
  
  // Close files and exit
  exit:
    fclose(infilecon);
    outfilecon.close();
    return exitcode;
  
}
