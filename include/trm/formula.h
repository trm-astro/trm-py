#ifndef TRM_FORMULA
#define TRM_FORMULA

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include "trm/subs.h"

//! Namespace for formula functions and classes
namespace Formula {

    //! Exception class for Formula
    class Formula_Error : public std::runtime_error {
    public:
	Formula_Error(const std::string& msg = "") : runtime_error(msg) {}
    };

    //! Basic structure to store tree formula data
    /** The structure store the variable name, number or operation,
     * the number of inputs needed by the node and a pointer to an array
     * of pointers to the input nodes.
     */

    struct Fnode {

	//! Default constructor to set correct null values   
	Fnode() : opvarnum(), ninput(0), input(NULL), value(0.) {}

	//! A string representing an operator, variable or number
	std::string opvarnum;
    
	//! Number of inputs to this node
	/** If ninput = 0, then the node should be the end of branch and opvarnum
	 * be a variable name to be looked up or a number. If ninput < 0 then
	 * the node is the end of a branch and value is the numerical value. If ninput > 0
	 * then opvarnum should be an operation such as sqrt. 
	 */
	int ninput;

	//! Array of ninput pointers to the input nodes
	Fnode **input;  

	//! The value if ninput < 0 
	double value;
    };

    //! Loads a tree given an expression
    void loader(std::string expression, Fnode *form_node);

    //! Strips redundant outer brackets from an expression
    void stripper(std::string& expression);

    //! Computes the value of an expression 
    double valuer(const Fnode* form_node, const std::map<std::string,double>& vars);

    //! Substitutes the values of known variables.
    void substitute(Fnode* form_node, const std::map<std::string, double>& vars);

    //! Checks whether all variables are specified
    void checker(const Fnode* form_node, const std::map<std::string, double>& vars);  

    //! Ouputs an expression tree in algebraic form
    void lister(const Fnode* form_node, int &level, std::ostream& ostr);

    //! Deletes memory associated with a tree
    void deleter(Fnode* form_node);

    //! Copies a tree
    void copier(const Fnode* old_formula, Fnode* new_formula);

    //! Forms a tree to calculate derivatives of an expression
    void derivative(const Fnode* form_node, Fnode* derv_node, const std::string& variable);
  
    //! Prunes a tree of dead wood
    bool pruner(Fnode* form_node);

    //! Formula class to represent mathematical formulae

    /** Formula is a class for representing mathematucal formulae. It does
     * so by creating a tree structure given an expression which can then
     * be used to evaluate the formula given the values of the variables.
     *
     * Formula can handle quite complex expressions such as a*sin(TWOPI*(x-x0)/P)/(1.e3+2.*t)
     * and it can compute their derivatives. While much slower than direct evaluation in a program
     * it allows the user to vary expression at run time rather than compile time.
     */
    class Formula {
    public:
	//! Standard constructor
	Formula() : head(NULL) {}
      
	//! Equation loading constructor
	/** From a string defining the expression this sets the tree
	 * \param expression An expression
	 */
	Formula(const std::string& expression); 
      
	//! Destructor
	~Formula(); 
      
	//! Copy constructor
	Formula(const Formula &eqn);
      
	//! Assignment operator
	Formula &operator=(const Formula &eqn); 
      
	//! Lists the formula tree
	void list(std::ostream& ostr) const; 
      
	//! Evaluates formula
	double value(const std::map<std::string,double>& vars) const;
      
	//! substitutes values of known variables
	void subst(const std::map<std::string,double>& vars);
      
	//! Check whether all variables have been specified
	void check(const std::map<std::string,double>& vars) const;
      
	//! Computes a derivative formula
	Formula deriv(const std::string& variable) const;
      
    private:
	Fnode *head; // head of the tree
    };
  
};

#endif
