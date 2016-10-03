// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <algorithm>
#include <wordexp.h>


//' Counts the lines of a text file.
//' 
//' @param fName file name
//' @return The number of lines
//' @export
// [[Rcpp::export]]
long lineCount(std::string fName) { 
    
    unsigned long number_of_lines = 0;
    
    wordexp_t exp_result;
    wordexp(fName.c_str(), &exp_result, WRDE_NOCMD);
    char* sanitizedFilename = exp_result.we_wordv[0];
    
    
    std::ifstream fin(sanitizedFilename);

    std::string line;
    while (getline(fin, line))
        ++number_of_lines;
    //std::cout << "Number of lines in text file " << number_of_lines;
    
    fin.close();
    wordfree(&exp_result);
    return number_of_lines;
}




//' Sub-samples reads from a FASTQ-file.
//' 
//' @param fqR1Name file name of R1-FASTQ-file
//' @param fqR2Name file name of R2-FASTQ-file
//' @param nbrDesiredReads long number of reads to sample from FASTQ file
//' @param fOutName file name for output
//' @export
// [[Rcpp::export]]
int subSampleFastq(std::string fqR1Name, std::string fqR2Name, unsigned long nbrDesiredReads, std::string fOutName) {
  
  const unsigned long nbrLines = lineCount(fqR1Name);
  const unsigned long nbrReads = nbrLines / 4;
  
  if ( nbrLines <= 0 ) {
    std::cout << "No lines found in file!\n";
    return -1;
  }

  if (  nbrLines % 4 != 0 ) {
    std::cout << "Number of lines in file not multiple of 4!\n";
    return -2;
  }
  
  if ( nbrDesiredReads > nbrReads ) {
    std::cout << "Requested more reads than available in file!\n";
    return -3;
  }
  
  std::vector<long> readEntries(nbrReads);
  for (size_t i = 0; i < nbrReads; ++i){
    readEntries[i] = i+1;
  }
  
  std::shuffle(readEntries.begin(), readEntries.end(),  std::default_random_engine(std::time(NULL)));
  //showContainer(readEntries);
  
  
  // sort first nbrDesiredReads entries
  std::sort(readEntries.begin(), std::next(readEntries.begin(), nbrDesiredReads));
  for (auto v : readEntries) std::cout << v << "\n";

  
  std::vector<long> outp(nbrDesiredReads);
  std::adjacent_difference(readEntries.begin(), std::next(readEntries.begin(), nbrDesiredReads), outp.begin());
    
  for (auto v : outp) std::cout << v << "\n";
  
  
  // read in file
  std::ifstream finR1(fqR1Name);
  std::ifstream finR2(fqR2Name);
  std::ofstream fout(fOutName, std::ios::trunc);
  
  if (!finR1) {
    std::cout << fqR1Name << " can not be read!\n";
    return 7;
  }
  
  if (!finR2) {
    std::cout << fqR2Name << " can not be read!\n";
    return 7;
  }
  
  if (!fout) {
    std::cout << fOutName << " can not be opened!\n";
    return 7;
  }
  
  std::string line;
  
  
  // sub-sample the reads
  for (size_t i = 0; i < nbrDesiredReads; ++i){
    std::cout << readEntries[i] << std::endl;
    // fast-forward
    for (size_t j = 0; j < 4 * (outp[i]-1); ++j) {
      getline(finR1, line);
      getline(finR2, line);
    }
    
    // copy entry into output file: here interleaved mode!
    for (int k = 0; k < 4; ++k) {
      getline(finR1, line);
      fout << line << "\n";
    }
    for (int k = 0; k < 4; ++k) {
      getline(finR2, line);
      fout << line << "\n";
    }
  }
  
  finR1.close();
  finR2.close();
  fout.close();
  
  return 0;
} 


/*** R
subSampleFastq(fqR1Name = "/Users/kuhnmat/test1a.fq", fqR2Name = "/Users/kuhnmat/test2a.fq", 3, fOutName = "/Users/kuhnmat/testa.subfq")
*/
