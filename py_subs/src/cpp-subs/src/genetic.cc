#include <cstdio>
#include "trm/subs.h"
#include "trm/array1d.h"

// Note that most routines here need to know the precise number format.
// At the moment this is of the form +14.7e which gives a number of the form
// [+-]n.nnnnnnne[+-]nn

/** Translates a set of numbers into a string suitable for genetic programs
 * \param vals the vector of values to translate. Uses scientific notation of the
 * form [+-]n.nnnnnnne[+-]nn
 */
std::string Subs::Genetic::model_to_string(const Subs::Array1D<double>& vals) {
  std::string temp = "";
  char store[15];
  for(int i=0; i<vals.size(); i++){
    sprintf(store, "%#+14.7e", vals[i]);
    temp += store;
  }
  return temp;
}

/** Translates a string to a set of numbers
 * \param model a string consisting of a concatenation of all the values.
 * \param vals the vector of values extracted from the string
 * form [+-]n.nnnnnnne[+-]nn
 */
void Subs::Genetic::string_to_model(const std::string& model, Subs::Array1D<double>& vals) {
  std::string::size_type n = model.length();
  if(int(n / 14) != vals.size())
    throw Subs_Error("Subs::Genetic::string_to_model: std::string and std::vector have conflicting numbers");
  for(int i=0; i<vals.size(); i++){
    std::istringstream istr(model.substr(14*i,14));
    istr >> vals[i];
    if(!istr)
      throw Subs_Error("Subs::Genetic::string_to_model: failed to translate " + model.substr(14*i,14));
  }
}

/** Crosses two models to create a third. Randomly picks a number of characters from
 * 1 to N-1 where N is the length of the two input strings and takes these from the first
 * model, and fills out the rest from the second model
 * \param model1 the first model
 * \param model2 the second model, must have same number of characters as the first
 * \param seed   the seed integer for the random number generator
 * \return string representing the new model
 */
std::string Subs::Genetic::cross(const std::string& model1, const std::string& model2, INT4& seed) {
  if(model1.length() != model2.length())
    throw Subs_Error("Subs::Genetic::cross: two models have differing numbers of characters");
  std::string::size_type n = 1 + std::string::size_type(Subs::ran2(seed)*(model1.length()-2));
  return model1.substr(0,n) + model2.substr(n);
}


/** Mutates a model. This works by randomising the digits, excluding those in the exponent
 * \param model the input model, which will be mutated at the end
 * \param seed   the seed integer for the random number generator
 * \param rate the mutation rate in terms of the mean number of digits altered. 
 */
void Subs::Genetic::mutate(std::string& model, INT4& seed, double rate) {
  
  std::string temp = model;
  int nparam  = model.length()/14;
  int ndigit  = nparam*(14 - 6);
  // The 10/9 factor accounts for the chance that the mutation actually just returns
  // the same digit
  double prob = 10.*rate/ndigit/9.;
  int digit;
  for(int i=0; i<nparam; i++){

    // First digit
    if(ran2(seed) < prob){
      do{
	digit = int(10*Subs::ran2(seed));
      }while(digit > 9);
      model[14*i+1] = Subs::digit_to_char(digit);
    }

    // Digits after decimal point, excluding exponent
    for(int j=14*i+3; j<14*i+10; j++){
      if(ran2(seed) < prob){
	do{
	  digit = int(10*Subs::ran2(seed));
	}while(digit > 9);
	model[j] = Subs::digit_to_char(digit);
      }
    }
  }
}

