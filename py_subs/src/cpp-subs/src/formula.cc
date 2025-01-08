#include <cstdlib>
#include <iostream>
#include <string>
#include <map>
#include "trm/formula.h"
#include "trm/constants.h"

/** Constructor which loads an expression as a tree. Most
 * of the work is done with the recursive routine 'loader';
 * 'pruner' is applied in an attempt to simply and hence speed
 * up the result.
 * \param expression the expression, e.g. a*(b+c/d)/(e+f)
 * \exception Throws Formula::Formula_Error exceptions if
 * the expression does not parse.
 */
Formula::Formula::Formula(const std::string& expression) {
    head = new Fnode;
    try {
	loader(expression, head);
	while(pruner(head));
    }
    catch(const Formula_Error& err){
	if(head){
	    deleter(head);
	    delete head;
	}
	throw Formula_Error(std::string("Formula::Formula(const std::string&) constructor failed, error = ") + err.what());
    }
}

/** Destructor. Deletes memory
 * associated with all nodes of a tree including the top-level one,
 */
Formula::Formula::~Formula(){
    if(head != NULL){
	deleter(head);
	delete head;
    }
}

/** Copy constructor. Copies entire tree
 * \param eqn the tree to copy
 */
Formula::Formula::Formula(const Formula& eqn){
    if(eqn.head != NULL){
	head = new Fnode;
	copier(eqn.head,head);
    }else{
	head = NULL;
    }
}

/** Assignment. wipe out tree to be assigned
 * whatever it is, and replcae with new version.
 * \param eqn the tree to be copied
 */
Formula::Formula& Formula::Formula::operator=(const Formula& eqn){
    if(this == &eqn) return *this;
    
    if(eqn.head != NULL){
	if(head == NULL) head = new Fnode;
	copier(eqn.head,head);
    }else{
	if(head != NULL){
	    delete head;
	    head = NULL;
	}
    }
    return *this;
}

/** Prints out an expression equivalent to a tree
 * \param ostr the output stream to print to.
 */
void Formula::Formula::list(std::ostream& ostr) const {
    int level = 0;
    if(head != NULL)
	lister(head, level, ostr);
    else
	ostr << "Empty Formula" << std::endl;
}

/** Tree evaluator. Computes the value of an expression given the values
 * of all of its variables.
 * \param vars map of variable names and values. Must cover all the variables needed.
 * \exception Throws a Formula::FOrmula_Error exception if all variables not present
 */
double Formula::Formula::value(const std::map<std::string,double>& vars) const {
    if(head != NULL)
	return valuer(head, vars);
    else
	return 0.;
}

/** Substitutes known variables by their value. This should I hope speed up
 * later evaluation of the formula.
 * \param vars map of variable names and values. Must cover all the variables needed.
 * \exception Throws a Formula::FOrmula_Error exception if all variables not present
 */
void Formula::Formula::subst(const std::map<std::string,double>& vars) {
    if(head != NULL){
	substitute(head, vars);
	while(pruner(head));
    }
}

/** Tree checker. Checks that all necessary variables are present without computing anything.
 * Also comes back with false if the tree is empty
 * \param vars map of variable names and values. Must cover all the variables needed.
 * \return true if all variables are present in vars. False if not or the tree is empty.
 */
void Formula::Formula::check(const std::map<std::string,double>& vars) const {
    if(head != NULL)
	checker(head, vars);
    else
	throw Formula_Error("formula empty");
}

/** Computes a tree which represents the derivative of the object tree
 *  As long as each possible operation has a derivative then this routine
 *  can work.
 *
 */
Formula::Formula Formula::Formula::deriv(const std::string& variable) const {
  
    Formula temp;
  
    if(head){
	temp.head = new Fnode;
	derivative(head, temp.head, variable);
	while(pruner(temp.head));

    }else{
	temp.head = NULL;
    }
  
    return temp;
}

/**
 * Formula::loader loads up a tree structure given an expression. 
 * It does so by recursively parsing the expression which is  
 * stored in a tree structure for later fast evaluation. Includes many
 * checks for valid expressions so it should be reasonably
 * difficult to get an invalid expression past it.
 *
 * Example expressions that loader can handle:
 * <ul>
 * <li>a+b
 * <li>2*a-b
 * <li>2.3e-5*a-b
 * <li>2.3e-5*(a-b)
 * <li>x/y/(a/b)
 * <li>sqrt(a+b)/sqrt(a-b)
 * </ul>
 * loader recognises the following functions: sqrt, sqr, cos, sin, exp, pow, ln, inv (=1/x).
 * loader also recognises various standard values, namely ZERO, UNIT, PI, TWOPI, MUNIT
 * (=-1), VLIGHT (km/s), DAY (secs), HALPHA, HBETA, HGAMMA, HDELTA (wavelengths in air, 
 * angstroms). 
 *
 * \param expression expression string to be parsed
 * \param form_node  top of the treee.
 */

void Formula::loader(std::string expression, Fnode* form_node){
  
    // List of recognised functions and their number of arguments. 
    std::map<std::string,int> functions;
    functions["sqrt"] = 1;
    functions["sqr"]  = 1;
    functions["cos"]  = 1;
    functions["sin"]  = 1;
    functions["exp"]  = 1;
    functions["pow"]  = 2;
    functions["ln"]   = 1;

    std::string lop, buff;
  
    std::map<std::string,int>::const_iterator fit;
    int i, posn = 0, nargs;
    int precedence = 0, arg_first = 0, arg_last = 0;

    // Start by stripping redundant brackets and white space; then wind through
    // expression to identify last operation; bit tricky.

    stripper(expression);
  
    bool start = false, found_op = false, is_a_number = false;
    bool is_a_func = false, num_var_func = false;
    bool exp_sign_next = false, exp_next = false;
    bool had_dot = false, had_an_exp = false;
    bool all_blank = true;
    int nc = 0, depth = 0;
  
    std::string::const_iterator sit = expression.begin();
    while(sit != expression.end()){
	if(*sit != ' '){              // skip blanks
      
	    all_blank = false;
      
	    // Don't check anything if we are inside brackets
	    if(depth == 0){
	
		// First check for leading + or - signs.
		// precedence -- larger equals stronger.
		// e.g. binary * beats binary +
	
		if(*sit == '-' && !start){
		    lop        = "u-"; // unary minus
		    posn       = nc;
		    found_op   = true;
		    precedence = 10;

		}else if(*sit == '+' && !start){
		    lop        = "u+"; // unary plus
		    posn       = nc;
		    found_op   = true;
		    precedence = 10;

		}else if((*sit == '*' || *sit == '/' ||
			  *sit == '.' || *sit == ',') && !start){
		    // expression should not start with any of the above
		    // although they are legal at other places.
	  
		    throw Formula_Error("Formula::loader: one of */., in illegal position in expresssion = " + expression);
	  
		}else{
	  
		    if(*sit != '(' && *sit != ')'){
	    
			// Multiplication and division
			if(*sit == '*' || *sit == '/'){
	      
			    if(num_var_func && is_a_number && (exp_sign_next || (exp_next && !had_an_exp)))
				throw Formula_Error("Formula::loader: character number " + Subs::str(nc+1) + " = " + *sit + " is invalid within a number");
	      
			    is_a_func    = false;
			    num_var_func = false;
			    if(!found_op || precedence >= 5){
				lop          = *sit;
				posn         = nc;
				found_op     = true;
				precedence   = 5;
			    }
	      
			    // addition and subtraction. have to guard against being
			    // in a number and just about to get the exponent
			    // as in 2.0e+01 or 3.45e-5 as indicated by exp_sign_next
		
			}else if(!exp_sign_next && (*sit == '+' || *sit == '-')){
			    is_a_func    = false;
			    num_var_func = false;
			    if(!found_op || precedence >= 1){
				lop          = *sit;
				posn         = nc;
				found_op     = true;
				precedence   = 1;
			    }
		
			}else{
		
			    if(num_var_func){
		  
				// We are in the process of building a number, a variable
				// or a function. Must past a few tests on the way.
				// Numbers are the tough ones; variables etc can be 
				// almost anything.
		  
				if(is_a_number){
		    
				    // A digit is always OK, so only check if not.
				    // checks are
				    //
				    // 1) + or - after an e (if not a digit)
				    // 2) only digits after e
				    // 3) only e or . if not after an exponent
				    // 4) no more than one .
		    
				    if(!isdigit(*sit) && ((exp_sign_next &&  *sit != '+' && *sit != '-') || exp_next ||
							  (!exp_sign_next && *sit != 'e' && *sit != '.') || (had_dot && *sit == '.')))
					throw Formula_Error(std::string("Formula::loader: invalid number (character = ") + *sit + 
							    std::string(") in expression = ") + expression);
		    
				    if(exp_next)
					had_an_exp = true;
				    if(exp_sign_next){
					exp_sign_next = false;
					exp_next      = true;
				    }
		    
				    exp_sign_next = (*sit == 'e');
				    had_dot       = (*sit == '.');
		    
				}
				buff += *sit;
		  
			    }else{
		  
				// First element which is not one of ()+-*/ and
				// so should be a number, variable or function.
				// numbers must start with a digit, variables or
				// functions must start with [a-z] or [A-Z].
		  
				if(!(is_a_number = isdigit(*sit)) && !isalpha(*sit))
				    throw Formula_Error("Formula::loader: character number " + Subs::str(nc+1) + " = " + *sit + " is invalid");
				num_var_func = true;
				buff = *sit;
				if(is_a_number){
				    exp_sign_next = exp_next = false;
				    had_an_exp = had_dot = false;
				}
			    }
			}
		    }else if(num_var_func && is_a_number){
			throw Formula_Error("Formula::loader: character number " + Subs::str(nc+1) + " = " + *sit + " is invalid");
		    }
		}
	    }
	
	    // Set depth 
	    if(*sit == '('){
		depth++;
		if(num_var_func){
	    
		    // We have a function. Need to check that it is one
		    // we know about.
		    if((fit = functions.find(buff)) == functions.end())
			throw Formula_Error("Formula::loader: function " + buff + " not recognised.");
	    
		    if(!found_op || precedence >= 20){
			lop        = buff;
			found_op   = true;
			precedence = 20;
			is_a_func  = true;
			arg_first  = nc+1;
			arg_last   = 0;
		    }
		    num_var_func = false;
		}
	  
	    }else if(*sit == ')'){
		depth--;
		if(depth < 0)
		    throw Formula_Error("Formula::loader: unmatched ) found in expression = " + expression);
	    }
	    if(is_a_func && depth == 0 && !arg_last)
		arg_last = nc - 1;
	}
	start = true;
	sit++;
	nc++;
    }

    if(depth != 0)
	throw Formula_Error("Formula::loader: unmatched ( found in expression = " + expression);
  
    if(all_blank)
	throw Formula_Error("Formula::loader: expression blank");

    // If no operation found, then the expression is either
    // a simple number or variable and the parsing stops.

    if(!found_op){

	// load up end node
	form_node->opvarnum = expression;
	if(isdigit(form_node->opvarnum[0])){

	    std::istringstream istr(form_node->opvarnum);
	    istr >> form_node->value;
	    if(!istr)
		throw Formula_Error("Formula::loader: failed to translate " + form_node->opvarnum + " as a number.");
	    form_node->ninput   = -1;

	}else{

	    if(form_node->opvarnum == "ZERO"){
		form_node->value  =  0.;
		form_node->ninput = -1;
	    }else if(form_node->opvarnum == "UNIT"){
		form_node->value  =  1.;
		form_node->ninput = -1;
	    }else if(form_node->opvarnum == "MUNIT"){
		form_node->value  = -1.;
		form_node->ninput = -1;
	    }else if(form_node->opvarnum == "PI"){
		form_node->value  =  Constants::PI;
		form_node->ninput = -1;
	    }else if(form_node->opvarnum == "TWOPI"){
		form_node->value  =  Constants::TWOPI;
		form_node->ninput = -1;
	    }else if(form_node->opvarnum == "VLIGHT"){
		form_node->value  =  Constants::C/1000.;
		form_node->ninput = -1;
	    }else if(form_node->opvarnum == "DAY"){
		form_node->value  =  Constants::DAY;
		form_node->ninput = -1;
	    }else if(form_node->opvarnum == "HALPHA"){
		form_node->value  =  Constants::HALPHA;
		form_node->ninput = -1;
	    }else if(form_node->opvarnum == "HBETA"){
		form_node->value  =  Constants::HBETA;
		form_node->ninput = -1;
	    }else if(form_node->opvarnum == "HGAMMA"){
		form_node->value  =  Constants::HGAMMA;
		form_node->ninput = -1;
	    }else if(form_node->opvarnum == "HDELTA"){
		form_node->value  =  Constants::HDELTA;
		form_node->ninput = -1;
	    }else{
		form_node->ninput = 0;
	    }
	}
	return;
    }

    // At this stage we have gone through whole string and should
    // have found last operation. We now need to extract its arguments
    // which may themselves be other expressions. If the last operation
    // is a function, its arguments should be comma separated as in
    // pow((x+y),3). We know the positions of the brackets from 
    // arg_first and arg_last, so we just need to locate commas
    // and check that number of arguments is OK. If the last operation
    // is */+- (binary) then arguments are everything to the left and right
    // etc.

    std::vector<std::string> args;
    if(lop == "+" || lop == "-" || lop == "*" || lop == "/"){

	// binary +, -, *, /
	nargs = 2;
	args.resize(nargs);
	args[0] = expression.substr(0,posn);
	args[1] = expression.substr(posn+1);
      
    }else if(lop == "u-" || lop == "u+"){
      
	// unary - or +
	nargs = 1;
	args.resize(nargs);
	args[0] = expression.substr(posn+1);
      
    }else{

	args.resize(fit->second);
	sit   = expression.begin();
	sit  += arg_first;
	nc    = arg_first;
	nargs = 0;
	while(nc <= arg_last && nargs < fit->second){
      
	    if(*sit == '('){
		depth++;
	    }else if(*sit == ')'){
		depth--;
	    }
	    if(depth == 0){
		if(*sit == ','){
		    args[nargs] = expression.substr(arg_first, nc-arg_first);
		    arg_first = nc + 1;
		    nargs++;
		}else if(nc == arg_last){
		    args[nargs] = expression.substr(arg_first, nc-arg_first+1);
		    nargs++;
		}
	    }    	  
	    sit++;
	    nc++;
      
	}
	if(nargs != fit->second){
	    if(nargs > fit->second)
		throw Formula_Error("Formula::loader: too many arguments for operation = " + lop + " in expression = " + expression);
	    else
		throw Formula_Error("Formula::loader: too few arguments for operation = " + lop + " in expression = " + expression);
	}
    }
  
    // Create new node
    form_node->opvarnum = lop;
    form_node->ninput   = nargs;
    form_node->input    = new Fnode* [nargs];      

    // First set pointers to null as it makes it easier to clear up
    // in case of disaster.
    for(i=0; i<nargs; i++)
	form_node->input[i] = NULL;
      
    // Now a bit of recursion
    for(i=0; i<nargs; i++){
	form_node->input[i] = new Fnode;      
	loader(args[i], form_node->input[i]);
    }

}

/** Formula::stripper recursively strips redundant pairs of brackets from an expression (including
 * any leading/trailing blanks).
 * \param expression input expression, modified on output
 * \exception Throws exceptions if brackets don't match up. 
 */

void Formula::stripper(std::string& expression){
  
    if(expression == "") return;

    int depth = 0;
    std::string::size_type first, length = 0, last, n = 0;
    first = expression.find_first_not_of(" ");
    last  = expression.find_last_not_of(" ");
    if(first != std::string::npos)
	expression = expression.substr(first,last-first+1);

    std::string::iterator sit = expression.begin(); 
    bool all_blank = true;
    while(sit != expression.end()){
	if(*sit != ' '){
	    all_blank = false;
	    if(*sit == '('){
		depth++;
		if(depth == 1) 
		    first = n + 1;
	    }else if(*sit == ')'){
		depth--;
		if(depth < 0)
		    throw Formula_Error("Formula::stripper: unmatched ) found in expression = " + expression);
		if(depth == 0) 
		    length = n - first;
	    }else if(depth == 0){
		return;
	    }
	}
	sit++;
	n++;
    }
    if(depth != 0)
	throw Formula_Error("Formula::stripper: unmatched ( found in expression = " + expression);

    if(all_blank) return;
  
    // If we have got here then there are a pair of enclosing brackets
    // Strip them off
    expression = expression.substr(first,length);
  
    // Now recursively look for any other such pairs.
    stripper(expression);

}

/** Formula::copier copies one tree to another
 * \param old_node the head of the tree to copy from
 * \param copy_node the head of the tree to copy to
 */
void Formula::copier(const Fnode* old_node, Fnode* copy_node){
  
    if(old_node == NULL || copy_node == NULL)
	throw Formula_Error("Formula::copier: one of two input pointers is null");
  
    if(copy_node == old_node)
	return;

    // Clear the memory of the node to be copied to
    deleter(copy_node);

    copy_node->ninput   = old_node->ninput;
    copy_node->opvarnum = old_node->opvarnum;
    copy_node->value    = old_node->value;

    if(old_node->ninput > 0){
	copy_node->input = new Fnode* [old_node->ninput];

	for(int i=0; i<old_node->ninput; i++){
	    copy_node->input[i] = new Fnode;      
	    copier(old_node->input[i], copy_node->input[i]);
	}
    }

}
    
/** deleter deallocates all the memory of a tree below a given node
 * \param form_node the top level of the tree which is not deleted
 */  
void Formula::deleter(Fnode* form_node){
    if(form_node != NULL){
	if(form_node->ninput > 0){
	    for(int i=0;i<form_node->ninput;i++){
		deleter(form_node->input[i]);
		delete form_node->input[i];
	    }
	    delete form_node->input;
	    form_node->ninput = 0;
	}
    }
}
  
/** valuer computes the value of a particular expression
 * given variable values
 * \param form_node the top of the tree
 * \param vars     the list of variable names and equivalent values
 * \return The value of the expression will be returned
 */
  
double Formula::valuer(const Fnode* form_node, const std::map<std::string, double>& vars){
    
    // Maximum number of arguments and their values
    const int MAX_ARG = 3;
    double arg[MAX_ARG];
    
    if(form_node->ninput > 0){
      
	// Operation node
	for(int i=0;i<form_node->ninput;i++)
	    arg[i] = valuer(form_node->input[i],vars);
      
	if(form_node->opvarnum == "*"){
	    return arg[0]*arg[1];
	}else if(form_node->opvarnum == "/"){
	    return arg[0]/arg[1];
	}else if(form_node->opvarnum == "+"){
	    return arg[0]+arg[1];
	}else if(form_node->opvarnum == "-"){
	    return arg[0]-arg[1];
	}else if(form_node->opvarnum == "u+"){
	    return arg[0];
	}else if(form_node->opvarnum == "u-"){
	    return -arg[0];
	}else if(form_node->opvarnum == "sqrt"){
	    return sqrt(arg[0]);
	}else if(form_node->opvarnum == "sqr"){
	    return arg[0]*arg[0];
	}else if(form_node->opvarnum == "cos"){
	    return cos(arg[0]);
	}else if(form_node->opvarnum == "sin"){
	    return sin(arg[0]);
	}else if(form_node->opvarnum == "exp"){
	    return exp(arg[0]);
	}else if(form_node->opvarnum == "pow"){
	    return pow(arg[0],arg[1]); 
	}else if(form_node->opvarnum == "ln"){
	    return log(arg[0]);
	}else{
	    throw Formula_Error("Formula::valuer: unrecognised operation = " + form_node->opvarnum);
	}

    }else if(form_node->ninput < 0){

	// Number node
	return form_node->value;
    
    }else{
      
	// Variable node
	std::map<std::string,double>::const_iterator fit = vars.find(form_node->opvarnum);
	if(fit == vars.end())
	    throw Formula_Error("Formula::valuer: could not recognize variable = " + form_node->opvarnum);
	else
	    return fit->second;
    }

}

/** substitute switches the values of all variables in a formula input via 
 * the vars variable for numbers. This is to speed up execution time.
 * \param form_node the top of the tree
 * \param vars     the list of variable names and equivalent values; does not
 * have to cover all variables, juts those which are known.
 * \return The value of the expression will be returned
 */
  
void Formula::substitute(Fnode* form_node, const std::map<std::string, double>& vars){
    
    if(form_node->ninput > 0){
      
	// Operation node
	for(int i=0;i<form_node->ninput;i++)
	    substitute(form_node->input[i],vars);
    
    }else if(form_node->ninput == 0){
      
	// Variable node
	std::map<std::string,double>::const_iterator fit = vars.find(form_node->opvarnum);
	if(fit != vars.end()){
	    form_node->ninput   = -1;
	    form_node->opvarnum = "NUMBER";
	    form_node->value    = fit->second;
	}
    }
}

/** checker checks that all variable values are specified.
 * \param form_node the top of the tree
 * \param vars     the list of variable names and equivalent values
 * \exception Throw a Forumla_Error if one of the variables have been specified
 */
  
void Formula::checker(const Fnode* form_node, const std::map<std::string, double>& vars){
    
    if(form_node->ninput > 0){
      
	// Operation node
	for(int i=0;i<form_node->ninput;i++)
	    checker(form_node->input[i],vars);

    }else if(form_node->ninput == 0){
      
	// Variable node
	std::map<std::string,double>::const_iterator fit = vars.find(form_node->opvarnum);
	if(fit == vars.end())
	    throw Formula_Error("Formula::checker: variable \'" + form_node->opvarnum + "\' not specified.");
  
    }
}


/** Formula::lister outputs an expression equivalent to a tree. It may not be
 * particularly elegant (it tends to overdo the brackets), but it should 
 * be understandable.
 * \param form_node the top node of the tree
 * \param level keep track of bracketting level. Set to 0 at the start
 * \param ostr output stream
 */

void Formula::lister(const Fnode* form_node, int &level, std::ostream& ostr){

    level++;
    if(form_node->ninput > 0){
    
	if(form_node->opvarnum == "*" || 
	   form_node->opvarnum == "/" ||
	   form_node->opvarnum == "+" ||
	   form_node->opvarnum == "-"){

	    if(level > 1) ostr << "(";
	    lister(form_node->input[0],level,ostr);
	    level--;
	    ostr << form_node->opvarnum;
	    lister(form_node->input[1],level,ostr);
	    level--;

	    if(level != 1)
		ostr << ")";

	}else{

	    if(form_node->opvarnum == "u+"){
		ostr << "+(";
	    }else if(form_node->opvarnum == "u-"){
		ostr << "-(";
	    }else{
		ostr << form_node->opvarnum << "(";
	    }

	    for(int i=0; i<form_node->ninput; i++){
		if(i > 0)
		    ostr << ",";
		lister(form_node->input[i],level,ostr);
		level--;
	    }
	    ostr << ")";
	}
    }else if(form_node->ninput < 0){
	ostr << form_node->value;
    }else{
	ostr << form_node->opvarnum;
    }
}

/** Formula::pruner removes redundant structures from a tree in order
 * to improve evaluation efficiency. After an application of pruner, 
 * if any pruning has happened, it may be possible to carry out more, 
 * so this routine should be applied to the point where nothing more happens.
 *
 * \param form_node the top level node of the tree
 * \return 'true' means that you can prune some more.
 */

bool Formula::pruner(Fnode* form_node){

    bool again = false;
    if(form_node){
	if(form_node->ninput > 0){
      
	    // First prune any inputs to this node
	    bool all_const = true;
	    for(int i=0; i<form_node->ninput; i++){
		if(pruner(form_node->input[i])) again = true;
		if(form_node->input[i]->ninput >= 0) all_const = false;
	    }

	    // Now see if the current node is prunable
	    if(all_const){
		// all inputs to the node are constants and so we can carry out the operation
		// and close out the node
		std::map<std::string, double> tvars;
		form_node->value  = valuer(form_node, tvars);
		form_node->ninput = -1;
		deleter(form_node);

	    }else if((form_node->opvarnum == "sqr"  && form_node->input[0]->opvarnum == "sqrt") ||
		     (form_node->opvarnum == "sqrt" && form_node->input[0]->opvarnum == "sqr") ||
		     (form_node->opvarnum == "u-"   && form_node->input[0]->opvarnum == "u-") ||
		     (form_node->opvarnum == "ln"   && form_node->input[0]->opvarnum == "exp") ||
		     (form_node->opvarnum == "exp"  && form_node->input[0]->opvarnum == "ln")){

		// a few straight simple inverse combinations: save argument two steps down the chain,
		// delete the chain, copy back in saved argument to cut out the two inverse steps
		Fnode* tnode = new Fnode;
		copier(form_node->input[0]->input[0], tnode);
		deleter(form_node);
		
		// Copy back relevant argument to top level
		copier(tnode,form_node);
		deleter(tnode);
		delete tnode;
		again = true;

	    }else if(form_node->opvarnum == "*"){

		if(form_node->input[0]->opvarnum == "ZERO" || form_node->input[1]->opvarnum == "ZERO" ||
		   (form_node->input[0]->ninput < 0 && form_node->input[0]->value == 0.) ||
		   (form_node->input[1]->ninput < 0 && form_node->input[1]->value == 0.)){

		    // If either input to a multiplication is zero then result will always
		    // be zero in which case we can destroy the input trees 
		    deleter(form_node);

		    form_node->opvarnum = "ZERO";
		    form_node->ninput = -1;
		    form_node->value  =  0.;
		    again = true;

		}else if(form_node->input[0]->opvarnum == "UNIT" ||
			 (form_node->input[0]->ninput < 0 && form_node->input[0]->value == 1.)){

		    // If input 0 = 1, save argument 1, and delete input trees
		    Fnode* tnode = new Fnode;
		    copier(form_node->input[1], tnode);
		    deleter(form_node);

		    // Copy argument 1 to top level, delete temporary copy
		    copier(tnode, form_node);
		    deleter(tnode);
		    delete tnode;
		    again = true;

		}else if(form_node->input[0]->opvarnum == "MUNIT" ||
			 (form_node->input[0]->ninput < 0 && form_node->input[0]->value == -1.)){

		    // If input 0 = -1, save -(argument 1)
		    Fnode* tnode = new Fnode;
		    copier(form_node->input[1], tnode);
		    deleter(form_node);

		    form_node->input    = new Fnode* [1];
		    form_node->input[0] = new Fnode;
		    form_node->ninput   = 1;
		    form_node->opvarnum = "u-";

		    // Copy argument 1 to top level, delete temporary copy
		    copier(tnode, form_node->input[0]);
		    deleter(tnode);
		    delete tnode;
		    again = true;
	  
		}else if(form_node->input[1]->opvarnum == "UNIT" ||
			 (form_node->input[1]->ninput < 0 && form_node->input[1]->value == 1.)){
	  
		    // If input 1 is unit then save argument 0
		    Fnode* tnode = new Fnode;
		    copier(form_node->input[0], tnode);
		    deleter(form_node);

		    // Copy temporary to top level then delete the old copy
		    copier(tnode, form_node);
		    deleter(tnode);
		    delete tnode;
		    again = true;
	  
		}else if(form_node->input[1]->opvarnum == "MUNIT" ||
			 (form_node->input[1]->ninput < 0 && form_node->input[1]->value == -1.)){
	  
		    // If input 1 = -1 then save argument -(argument 0)
		    Fnode* tnode = new Fnode;
		    copier(form_node->input[0], tnode);
		    deleter(form_node);

		    form_node->input    = new Fnode* [1];
		    form_node->input[0] = new Fnode;
		    form_node->ninput   = 1;
		    form_node->opvarnum = "u-";

		    // Copy temporary to top level then delete the old copy
		    copier(tnode, form_node->input[0]);
		    deleter(tnode);
		    delete tnode;
		    again = true;

		}

	    }else if(form_node->opvarnum == "/"){

		if(form_node->input[0]->opvarnum == "ZERO" ||
		   (form_node->input[0]->ninput < 0 && form_node->input[0]->value == 0.)){

		    // If first input is zero then result will
		    // always be zero
		    deleter(form_node);

		    form_node->opvarnum = "ZERO";
		    form_node->ninput = -1;
		    form_node->value  =  0.;
		    again = true;

		}else if(form_node->input[1]->opvarnum == "UNIT" ||
			 (form_node->input[1]->ninput < 0 && form_node->input[1]->value == 1.)){
	  
		    // If second input is unit then just pass through other argument
		    Fnode* tnode = new Fnode;
		    copier(form_node->input[0], tnode);
		    deleter(form_node);

		    // Copy argument 0 to top level, delete old copy
		    copier(tnode, form_node);
		    deleter(tnode);
		    delete tnode;
		    again = true;

		}else if(form_node->input[1]->opvarnum == "MUNIT" ||
			 (form_node->input[1]->ninput < 0 && form_node->input[1]->value == -1.)){
	  
		    // If second input is -1 then convert to u-
		    Fnode* tnode = new Fnode;
		    copier(form_node->input[0], tnode);
		    deleter(form_node);

		    form_node->input    = new Fnode* [1];
		    form_node->input[0] = new Fnode;
		    form_node->ninput   = 1;
		    form_node->opvarnum = "u-";

		    copier(tnode, form_node->input[0]);
		    deleter(tnode);
		    delete tnode;
		    again = true;

		}

	    }else if(form_node->opvarnum == "+"){

		if(form_node->input[0]->opvarnum == "ZERO" ||
		   (form_node->input[0]->ninput < 0 && form_node->input[0]->value == 0.)){


		    // If input is zero then just pass through other argument
		    Fnode* tnode = new Fnode;
		    copier(form_node->input[1], tnode);
		    deleter(form_node);

		    // Copy argument 0 to top level, delete old copy
		    copier(tnode, form_node);
		    deleter(tnode);
		    delete tnode;
		    again = true;

		}else if(form_node->input[1]->opvarnum == "ZERO" ||
			 (form_node->input[1]->ninput < 0 && form_node->input[1]->value == 0.)){

		    // If input is zero then just pass through other argument
		    Fnode* tnode = new Fnode;
		    copier(form_node->input[0], tnode);
		    deleter(form_node);

		    // Copy argument 0 to top level, delete old copy
		    copier(tnode, form_node);
		    deleter(tnode);
		    delete tnode;
		    again = true;

		}else if(form_node->input[1]->opvarnum == "u-"){

		    // Save argument
		    Fnode* tnode = new Fnode;
		    copier(form_node->input[1]->input[0], tnode);
		    deleter(form_node->input[1]);

		    // Copy back relevant argument up one level
		    copier(tnode,form_node->input[1]);
		    form_node->opvarnum = "-";
		    deleter(tnode);
		    delete tnode;
		    again = true;
		}

	    }else if(form_node->opvarnum == "-"){

		if(form_node->input[1]->opvarnum == "ZERO" ||
		   (form_node->input[1]->ninput < 0 && form_node->input[1]->value == 0.)){
	  
		    // If second input is zero then just pass through other argument
		    Fnode* tnode = new Fnode;
		    copier(form_node->input[0], tnode);
		    deleter(form_node);

		    // Copy argument 0 to top level, delete old copy
		    copier(tnode, form_node);
		    deleter(tnode);
		    delete tnode;
		    again = true;

		}else if(form_node->input[0]->opvarnum == "ZERO" ||
			 (form_node->input[0]->ninput < 0 && form_node->input[0]->value == 0.)){

		    // If first input is zero then switch to u-
		    Fnode* tnode = new Fnode;
		    copier(form_node->input[1], tnode);
		    deleter(form_node);

		    // Rejig the node
		    form_node->opvarnum = "u-";
		    form_node->input    = new Fnode* [1];
		    form_node->input[0] = tnode; 
		    form_node->ninput   = 1;
		    again = true;

		}else if(form_node->input[1]->opvarnum == "u-"){

		    // Save argument
		    Fnode* tnode = new Fnode;
		    copier(form_node->input[1]->input[0], tnode);
		    deleter(form_node->input[1]);

		    // Copy back relevant argument up one level
		    copier(tnode,form_node->input[1]);
		    form_node->opvarnum = "+";
		    deleter(tnode);
		    delete tnode;
		    again = true;
	  
		}

	    }else if(form_node->opvarnum == "u+"){

		// Save argument
		Fnode* tnode = new Fnode;
		copier(form_node->input[0], tnode);
		deleter(form_node);

		// Copy back saved argument to top level
		copier(tnode,form_node);
		deleter(tnode);
		delete tnode;
		again = true;

	    }else if(form_node->opvarnum == "pow"){

		if(form_node->input[1]->opvarnum == "ZERO" ||
		   (form_node->input[1]->ninput < 0 && form_node->input[1]->value == 0.)){

		    // If second input is zero then result will just equal 1.
		    deleter(form_node);

		    form_node->opvarnum = "UNIT";
		    form_node->ninput   = -1;
		    form_node->value    =  1.;
		    again = true;

		}else if(form_node->input[1]->opvarnum == "UNIT" ||
			 (form_node->input[1]->ninput < 0 && form_node->input[1]->value == 1.)){
	  
		    // If second input = 1 then just pass through first argument
		    Fnode* tnode = new Fnode;
		    copier(form_node->input[0], tnode);
		    deleter(form_node);

		    // Copy argument 0 to top level, delete old copy
		    copier(tnode, form_node);
		    deleter(tnode);
		    delete tnode;
		    again = true;

		}else if(form_node->input[1]->opvarnum == "MUNIT" ||
			 (form_node->input[1]->ninput < 0 && form_node->input[1]->value == -1.)){
	  
		    // If second input is -1 then convert to a division on 1 by the first input
		    Fnode* tnode = new Fnode;
		    copier(form_node->input[0], tnode);
		    deleter(form_node);

		    form_node->input    = new Fnode* [2];
		    form_node->input[0] = new Fnode;
		    form_node->input[1] = new Fnode;
		    form_node->ninput   = 2;
		    form_node->opvarnum = "/";

		    form_node->input[0]->opvarnum = "UNIT";
		    form_node->input[0]->ninput   = -1;
		    form_node->input[0]->value    =  1.;

		    // Copy old first argument to second input, delete temporary copy
		    copier(tnode, form_node->input[1]);
		    deleter(tnode);
		    delete tnode;
		    again = true;

		}else if(form_node->input[1]->ninput < 0 && form_node->input[1]->value == 2.){
	  
		    // If second input == 2. then change to sqr
		    Fnode* tnode = new Fnode;
		    copier(form_node->input[0], tnode);
		    deleter(form_node);

		    form_node->input    = new Fnode* [1];
		    form_node->input[0] = new Fnode;
		    form_node->ninput   = 1;
		    form_node->opvarnum = "sqr";

		    // Copy old first argument to top level, delete temporary copy
		    copier(tnode, form_node->input[0]);
		    deleter(tnode);
		    delete tnode;
		    again = true;

		}else if(form_node->input[1]->ninput < 0 && form_node->input[1]->value == 0.5){
	  
		    // If second input == 0.5 then change to sqrt
		    Fnode* tnode = new Fnode;
		    copier(form_node->input[0], tnode);
		    deleter(form_node);

		    form_node->input    = new Fnode* [1];
		    form_node->input[0] = new Fnode;
		    form_node->ninput   = 1;
		    form_node->opvarnum = "sqrt";

		    // Copy old first argument to top level, delete temporary copy
		    copier(tnode, form_node->input[0]);
		    deleter(tnode);
		    delete tnode;
		    again = true;

		}

	    }
	}
    }
    return again;
}

/** derivative constructs a tree representing the partial derivative of
 * an expression with respect to one of its variables. It is automatic
 * although it can produce rather clumsy structures which typically you may want 
 * to attempt to improve with 'pruner'.
 * \param form_node the head of the tree of the expression whose derivative is to be taken
 * \param derv_node the head of the derivative tree
 * \param variable  the variable to take the partial derivative with respect to
 */
void Formula::derivative(const Fnode* form_node, Fnode* derv_node, const std::string& variable){
  
    if(form_node == NULL)
	throw Formula_Error("Formula::derivative expression tree is null");

    if(derv_node == 0)
	throw Formula_Error("Formula::derivative derivative tree is null");

    Fnode *tnode, *tnode1;

    if(form_node->ninput > 0){

	// Operation node. 
	if(form_node->opvarnum == "*"){

	    // Derivative is sum of two terms
	    derv_node->opvarnum = "+";
	    derv_node->ninput   = 2;
	    derv_node->input    = new Fnode* [2];

	    // First term is da1*a2
	    tnode = derv_node->input[0] = new Fnode;
	    tnode->opvarnum = "*";
	    tnode->ninput   =  2;
	    tnode->input    = new Fnode* [2];

	    tnode->input[0] = new Fnode;
	    derivative(form_node->input[0],tnode->input[0],variable);

	    tnode->input[1] = new Fnode;
	    copier(form_node->input[1],tnode->input[1]);

	    // Second term is da2*a1
	    tnode = derv_node->input[1] = new Fnode;
	    tnode->opvarnum = "*";
	    tnode->ninput   =  2;
	    tnode->input    = new Fnode* [2];

	    tnode->input[0] = new Fnode;
	    derivative(form_node->input[1],tnode->input[0],variable);

	    tnode->input[1] = new Fnode;
	    copier(form_node->input[0],tnode->input[1]);


	}else if(form_node->opvarnum == "/"){

	    // Derivative is difference of two terms
	    derv_node->opvarnum = "-";
	    derv_node->ninput   = 2;
	    derv_node->input    = new Fnode* [2];

	    // First term is da1/a2
	    tnode = derv_node->input[0] = new Fnode;
	    tnode->opvarnum = "/";
	    tnode->ninput   =  2;
	    tnode->input    = new Fnode* [2];

	    tnode->input[0] = new Fnode;
	    derivative(form_node->input[0],tnode->input[0],variable);

	    tnode->input[1] = new Fnode;
	    copier(form_node->input[1],tnode->input[1]);

	    // Second term is (da2*(a1/sqr(a2)))
	    tnode = derv_node->input[1] = new Fnode;
	    tnode->opvarnum = "*";
	    tnode->ninput   =  2;
	    tnode->input    = new Fnode* [2];

	    tnode->input[0] = new Fnode;
	    derivative(form_node->input[1],tnode->input[0],variable);

	    // a1/sqr(a2)
	    tnode = tnode->input[1] = new Fnode;
	    tnode->opvarnum = "/";
	    tnode->ninput   =  2;
	    tnode->input    = new Fnode* [2];

	    tnode->input[0] = new Fnode;
	    copier(form_node->input[0],tnode->input[0]);

	    tnode = tnode->input[1] = new Fnode;
	    tnode->opvarnum = "sqr";
	    tnode->ninput   =  1;
	    tnode->input    = new Fnode* [1];

	    tnode->input[0] = new Fnode;
	    copier(form_node->input[1],tnode->input[0]);

	}else if(form_node->opvarnum == "+"){

	    // Derivative is sum of two terms
	    derv_node->opvarnum = "+";
	    derv_node->ninput   = 2;
	    derv_node->input    = new Fnode* [2];

	    // First term is da1
	    derv_node->input[0] = new Fnode;
	    derivative(form_node->input[0],derv_node->input[0],variable);

	    // Second term is da2
	    derv_node->input[1] = new Fnode;
	    derivative(form_node->input[1],derv_node->input[1],variable);

	}else if(form_node->opvarnum == "-"){

	    // Derivative is difference of two terms
	    derv_node->opvarnum = "-";
	    derv_node->ninput   = 2;
	    derv_node->input    = new Fnode* [2];

	    // First term is da1
	    derv_node->input[0] = new Fnode;
	    derivative(form_node->input[0],derv_node->input[0],variable);

	    // Second term is da2
	    derv_node->input[1] = new Fnode;
	    derivative(form_node->input[1],derv_node->input[1],variable);

	}else if(form_node->opvarnum == "u+"){

	    // Derivative is one term
	    derv_node->opvarnum = "u+";
	    derv_node->ninput   = 1;
	    derv_node->input    = new Fnode* [1];

	    // da1
	    derv_node->input[0] = new Fnode;
	    derivative(form_node->input[0],derv_node->input[0],variable);

	}else if(form_node->opvarnum == "u-"){

	    // Derivative is one term
	    derv_node->opvarnum = "u-";
	    derv_node->ninput   = 1;
	    derv_node->input    = new Fnode* [1];

	    // da1
	    derv_node->input[0] = new Fnode;
	    derivative(form_node->input[0],derv_node->input[0],variable);

	}else if(form_node->opvarnum == "sqrt"){

	    // da/(2*sqrt(a))
	    derv_node->opvarnum = "/";
	    derv_node->ninput   = 2;
	    derv_node->input    = new Fnode* [2];

	    derv_node->input[0] = new Fnode;
	    derivative(form_node->input[0],derv_node->input[0],variable);

	    // 2*sqrt(a)
	    tnode1 = tnode = derv_node->input[1] = new Fnode;
	    tnode->opvarnum = "*";
	    tnode->ninput   =  2;
	    tnode->input    = new Fnode* [2];

	    tnode = tnode->input[0] = new Fnode;
	    tnode->opvarnum = "2";
	    tnode->value    =  2.;
	    tnode->ninput   =  -1;

	    // sqrt(a)
	    tnode = tnode1->input[1] = new Fnode;
	    tnode->opvarnum = "sqrt";
	    tnode->ninput   =  1;
	    tnode->input    = new Fnode* [1];

	    tnode->input[0] = new Fnode;
	    copier(form_node->input[0],tnode->input[0]);

	}else if(form_node->opvarnum == "sqr"){

	    // da*(2*a)
	    derv_node->opvarnum = "*";
	    derv_node->ninput   = 2;
	    derv_node->input    = new Fnode* [2];

	    derv_node->input[0] = new Fnode;
	    derivative(form_node->input[0],derv_node->input[0],variable);

	    // 2*a
	    tnode1 = tnode = derv_node->input[1] = new Fnode;
	    tnode->opvarnum = "*";
	    tnode->ninput   =  2;
	    tnode->input    = new Fnode* [2];

	    tnode = tnode->input[0] = new Fnode;
	    tnode->opvarnum = "2";
	    tnode->value    =  2.;
	    tnode->ninput   = -1;

	    tnode1->input[1] = new Fnode;
	    copier(form_node->input[0],tnode1->input[1]);

	}else if(form_node->opvarnum == "cos"){

	    // da*(-sin(a))
	    derv_node->opvarnum = "*";
	    derv_node->ninput   = 2;
	    derv_node->input    = new Fnode* [2];

	    derv_node->input[0] = new Fnode;
	    derivative(form_node->input[0],derv_node->input[0],variable);

	    // -sin(a)
	    tnode = derv_node->input[1] = new Fnode;
	    tnode->opvarnum = "u-";
	    tnode->ninput   =  1;
	    tnode->input    = new Fnode* [1];

	    // sin(a)
	    tnode = tnode->input[0] = new Fnode;
	    tnode->opvarnum = "sin";
	    tnode->ninput   =  1;
	    tnode->input    = new Fnode* [1];

	    tnode->input[0] = new Fnode;
	    copier(form_node->input[0],tnode->input[0]);

	}else if(form_node->opvarnum == "sin"){

	    // da*(cos(a))
	    derv_node->opvarnum = "*";
	    derv_node->ninput   = 2;
	    derv_node->input    = new Fnode* [2];

	    derv_node->input[0] = new Fnode;
	    derivative(form_node->input[0],derv_node->input[0],variable);

	    tnode = derv_node->input[1] = new Fnode;
	    tnode->opvarnum = "cos";
	    tnode->ninput   =  1;
	    tnode->input    = new Fnode* [1];

	    tnode->input[0] = new Fnode;
	    copier(form_node->input[0],tnode->input[0]);

	}else if(form_node->opvarnum == "exp"){

	    // da*(exp(a))
	    derv_node->opvarnum = "*";
	    derv_node->ninput   = 2;
	    derv_node->input    = new Fnode* [2];

	    derv_node->input[0] = new Fnode;
	    derivative(form_node->input[0],derv_node->input[0],variable);

	    tnode = derv_node->input[1] = new Fnode;
	    tnode->opvarnum = "exp";
	    tnode->ninput   =  1;
	    tnode->input    = new Fnode* [1];

	    tnode->input[0] = new Fnode;
	    copier(form_node->input[0],tnode->input[0]);

	}else if(form_node->opvarnum == "pow"){

	    // Derivative is sum of two terms
	    derv_node->opvarnum = "+";
	    derv_node->ninput   = 2;
	    derv_node->input    = new Fnode* [2];

	    // First term is da1*(a2*pow(a1,(a2-1)))
	    tnode = derv_node->input[0] = new Fnode;
	    tnode->opvarnum = "*";
	    tnode->ninput   =  2;
	    tnode->input    = new Fnode* [2];

	    tnode->input[0] = new Fnode;
	    derivative(form_node->input[0],tnode->input[0],variable);

	    // a2*pow(a1,(a2-1)))
	    tnode = tnode->input[1] = new Fnode;
	    tnode->opvarnum = "*";
	    tnode->ninput   =  2;
	    tnode->input    = new Fnode* [2];

	    tnode->input[0] = new Fnode;
	    copier(form_node->input[1],tnode->input[0]);

	    // pow(a1,(a2-1))) 
	    tnode = tnode->input[1] = new Fnode;
	    tnode->opvarnum = "pow";
	    tnode->ninput   =  2;
	    tnode->input    = new Fnode* [2];

	    tnode->input[0] = new Fnode;
	    copier(form_node->input[0],tnode->input[0]);

	    // a2-1 
	    tnode = tnode->input[1] = new Fnode;
	    tnode->opvarnum = "-";
	    tnode->ninput   =  2;
	    tnode->input    = new Fnode* [2];

	    tnode->input[0] = new Fnode;
	    copier(form_node->input[1],tnode->input[0]);

	    tnode = tnode->input[1] = new Fnode;
	    tnode->opvarnum = "UNIT";
	    tnode->value    =  1.;
	    tnode->ninput   = -1;

	    // Second term is da2*(ln(a1)*pow(a1,a2))
	    tnode = derv_node->input[1] = new Fnode;
	    tnode->opvarnum = "*";
	    tnode->ninput   =  2;
	    tnode->input    = new Fnode* [2];

	    tnode->input[0] = new Fnode;
	    derivative(form_node->input[1],tnode->input[0],variable);

	    // ln(a1)*pow(a1,a2))
	    tnode1 = tnode = tnode->input[1] = new Fnode;
	    tnode->opvarnum = "*";
	    tnode->ninput   =  2;
	    tnode->input    = new Fnode* [2];

	    tnode = tnode->input[0] = new Fnode;
	    tnode->opvarnum = "ln";
	    tnode->ninput   =  1;
	    tnode->input    = new Fnode* [1];

	    tnode->input[0] = new Fnode;
	    copier(form_node->input[0],tnode->input[0]);

	    // pow(a1,a2)
	    tnode = tnode1->input[1] = new Fnode;
	    tnode->opvarnum = "pow";
	    tnode->ninput   =  2;
	    tnode->input    = new Fnode* [2];

	    tnode->input[0] = new Fnode;
	    copier(form_node->input[0],tnode->input[0]);

	    tnode->input[1] = new Fnode;
	    copier(form_node->input[1],tnode->input[1]);

	}else if(form_node->opvarnum == "ln"){

	    // da/a
	    derv_node->opvarnum = "/";
	    derv_node->ninput   = 2;
	    derv_node->input    = new Fnode* [2];

	    derv_node->input[0] = new Fnode;
	    derivative(form_node->input[0],derv_node->input[0],variable);

	    derv_node->input[1] = new Fnode;
	    copier(form_node->input[0],derv_node->input[1]);

	}else{
	    throw Formula_Error("Formula::derivative: unrecognised operation " + form_node->opvarnum);
	}

    }else{

	// Variable/number node. Either the variable matches
	// in which case the derivative = 1, or it does not in which case
	// the derivative = 0
	derv_node->opvarnum = new char [5];

	derv_node->ninput   = -1;
	if(form_node->opvarnum == variable){
	    derv_node->opvarnum = "UNIT";
	    derv_node->value    = 1.;
	}else{
	    derv_node->opvarnum = "ZERO";
	    derv_node->value    = 0.;
	}
    }
}
    






