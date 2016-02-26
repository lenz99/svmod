#include <Rcpp.h>
using namespace Rcpp;

#include <fstream>
#include <iostream>
#include <random>
#include <wordexp.h>

//' Counts the lines of a text file.
//' 
//' @param fName file name
//' @return The number of lines
//' @export
// [[Rcpp::export]]
long lineCount(std::string fName) { 
    
    long number_of_lines = 0;
    
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
//' @param fName file name of FASTQ-file
//' @export
// [[Rcpp::export]]
int subSampleFastq(std::string fName) {
  
  int nbrLines = lineCount(fName);
  
  if ( nbrLines <= 0 ) {
    std::cout << "No lines found in file!\n";
    return -1;
  }

  if (  nbrLines % 4 != 0 ) {
    std::cout << "Number of lines in file not multiple of 4!\n";
    return -2;
  }
  
  
  typedef std::minstd_rand G;
  typedef std::uniform_int_distribution<> D; //ZZZ C++0x?
  
  G g(time(NULL));
  D d(1, nbrLines/4);
  
  int c = 0;
  int N = 5;
  int nextNum = -1;
  
  std::ifstream fin(fName.c_str()); //ZZZ memory-mapped stuff?
  
  for (int i = 0; i < N; ++i){ 
      nextNum = d(g); // ZZZ duplicates!
      std::cout << "Next number is " << nextNum << ".\n";
      c += nextNum;
  }//rof
  
  return c;
  
  //ZZZ continue here! int c (sum the random numbers) is just a test
  
} 
