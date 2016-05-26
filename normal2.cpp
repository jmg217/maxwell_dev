#include <iostream>
#include <random>
#include <fstream>
int main()
{
  const int nrolls=10000;  // number of experiments
  const int nstars=100;    // maximum number of stars to distribute

  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0.0,1.0);
/*
  int p[11]={};

  for (int i=0; i<nrolls; ++i) {
    double number = distribution(generator);
    if ((number>=-3.0)&&(number<4.0)) ++p[int(number)];
  }

  std::cout << "normal_distribution (0.0,1.0):" << std::endl;

  for (int i=-3; i<4; ++i) {
    std::cout << i << "- " << (i+1) << ": ";
    std::cout << std::string(p[i]*nstars/nrolls,'*') << std::endl;
  }
*/

std::ofstream outFile("normal.txt");
 


for (int i=0;i<nrolls; ++i){
  double number = distribution(generator);
  outFile << number << std::endl;
 
}
 
outFile.close();



  return 0;
}

