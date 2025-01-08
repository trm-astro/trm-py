#include <map>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip>
#include "trm/subs.h"
#include "trm/header.h"

namespace Subs {

    class Varprop;

    //! A class for handling command input from the user.

    /**
     * Input handles the tricky task of getting user input. Its chief feature
     * is that it preserves values of inputs from one call to another of a function
     * by storage inside files. This can save a great deal of effort. It also can check ranges
     * of input variables. It uses two types of files, one associated with whatever the
     * particular command or program is called and another 'global' file which all commands
     * & programs (within a given set) access.
     *
     * There are many examples but basically one must first construct an Input, then 'sign in'
     * all the variable names one will be prompting for, and then finally input the values.
     *
     * One the Input goes out of scope and is destructed, the new values are stored within
     * the default files. 
     *
     * Input also supports inter-program communication through the global defaults file. Programs
     * and commands can set and access variables within this file. 
     */
    class Input : public std::map<std::string,Varprop> {
    public:

	//! Sets a variable to be local or global
	/** The defaults for parameters are stored in either a file particular to the program or
	 * a global file accessed by several programs. This allows different programs to access
	 * the same parameter. This enum specifies which of the two cases applies.
	 */
	enum Var_control1 {
	    LOCAL, /**< Indicates the parameter is stored in a file local to the program */  
	    GLOBAL /**< Indicates the parameter is stored in a global file accessed by several programs */  
	};

	//! Sets prompting characteristic of a variable
	enum Var_control2 {PROMPT,  NOPROMPT};

	//! Constructor given specific file names and the argument list 
	Input(const std::string& global_file, const std::string& local_file, int argc, char* argv[]);

	//! Constructor given the argument list as char* array and directories 
	Input(int argc, char* argv[], const std::string& directory_environment_var, const std::string& default_directory);

	//! Constructor given the argument list as strings and directories 
	Input(const std::vector<std::string>& args, const std::string& directory_environment_var, const std::string& default_directory);

	//! Destructor
	~Input();

	//! Saves inputs to the default files.
	void save() const;

	//! Sign in a variable.
	void sign_in(const std::string& varname, Var_control1 store_where, Var_control2 prompt);

	//! Get a range checked value
	template <class T>    
	void get_value(std::string name, T& value, const T& initial_value, const T& low, const T& high, 
		       const std::string& prompt_string);

	//! Get a range checked vector of doubles
	void get_value(std::string name, std::vector<double>& value, const double& initial_value, 
		       const double& low, const double& high, const std::vector<double>::size_type nvec,
		       const std::string& prompt_string);

	//! Get a string
	void get_value(std::string name, std::string& value, const std::string& initial_value, const std::string& prompt_string);
    
	//! Get a boolean
	void get_value(std::string name, bool& value, const bool& initial_value, const std::string& prompt_string);

	//! Get a character
	void get_value(std::string name, char& value, const char& initial_value, const std::string& allowed_values,
		       const std::string& prompt_string);

	//! Set default for char* strings
	void set_default(std::string name, const char* value);

	//! Set default (over-rides those from files)
	template <class T>    
	void set_default(std::string name, const T& value);

	//! Store a value in the global default file (inter program communication)
	template <class T>    
	void add_to_global(const std::string& name, const T& value);

	//! Copies an item in the global defaults
	void global_copy(const std::string& name, const std::string& copy);

	//! Retrieves a value from the global default file (inter program communication)
	template <class T>    
	void get_from_global(const std::string& name, T& value) const;

	//! Retrieves a value from the global default file (inter program communication)
	Subs::Hitem* get_from_global(const std::string& name) const;

	//! Checks whether a value is present in the global defaults
	bool is_a_global(const std::string& name) const;

	//! Prints out the global defaults 
	void print_global(std::ostream& ostr) const;

	//! Exception class
	/**
	 * Input objects throw Input_Error exceptions
	 */
	class Input_Error : public Subs_Error {
	public:
	    //! Default constructor.
	    Input_Error() : Subs_Error() {};
	    //! Constructor storing a message
	    Input_Error(const std::string& str) : Subs_Error(str) {}
	};
    
    private:

	// File names for storage & retrieval of global and local defaults
	std::string global_file_, local_file_; 

	// Argument counter, initialised to zero, incremented for each 'get_value' call
	size_t narg;

	// Control whether to automatically take default value (or initial value
	// if default not set), to list every value, to prompt in all cases (over-riding
	// prompt attribute of variables) and whether to interact with default files at all
	bool set_to_default, list_values, prompt_all, def_files;

	// List of command line arguments. These are initially stored from those sent to the program,
	// trimmed of the special keywords 'list', 'prompt' and 'initialise' and then trimmed of all the
	// 'variable=value' ones which are stored separately. This trimming happens as the variables
	// are signed in.
	std::vector<std::string> simple_args;

	// This is where the 'variable=value' arguments are stored.
	std::vector<std::string> named_args;

	// Whether the global and local default files were found
	bool global_exist, local_exist;

	// The global and local defaults stored in memory
	Header global_defs, local_defs;

	// The command name
	std::string command;

	// Routines needed for storage of defaults by template functions
	void local_store(const std::string& name, const bool& value) {local_defs.set(name, new Hbool(value, "Last set by: " + command));}
	void local_store(const std::string& name, const char& value) {local_defs.set(name, new Hchar(value, "Last set by: " + command));}
	void local_store(const std::string& name, const int& value) {local_defs.set(name, new Hint(value, "Last set by: " + command));}
	void local_store(const std::string& name, const unsigned int& value){local_defs.set(name, new Huint(value, "Last set by: " + command));}
	void local_store(const std::string& name, const long int& value){local_defs.set(name, new Hlint(value, "Last set by: " + command));}
	void local_store(const std::string& name, const unsigned long int& value){local_defs.set(name, new Hulint(value, "Last set by: " + command));}
	void local_store(const std::string& name, const double& value){local_defs.set(name, new Hdouble(value, "Last set by: " + command));}
	void local_store(const std::string& name, const float& value){local_defs.set(name, new Hfloat(value, "Last set by: " + command));}
	void local_store(const std::string& name, const std::string& value){local_defs.set(name, new Hstring(value, "Last set by: " + command));}
	void local_store(const std::string& name, const std::vector<double>& value){local_defs.set(name, new Hdvector(value, "Last set by: " + command));}

	void global_store(const std::string& name, const bool& value) {global_defs.set(name, new Hbool(value, "Last set by: " + command));}
	void global_store(const std::string& name, const char& value) {global_defs.set(name, new Hchar(value, "Last set by: " + command));}
	void global_store(const std::string& name, const int& value){global_defs.set(name, new Hint(value, "Last set by: " + command));}
	void global_store(const std::string& name, const unsigned int& value){global_defs.set(name, new Huint(value, "Last set by: " + command));}
	void global_store(const std::string& name, const long int& value){global_defs.set(name, new Hlint(value, "Last set by: " + command));}
	void global_store(const std::string& name, const unsigned long int& value){global_defs.set(name, new Hulint(value, "Last set by: " + command));}
	void global_store(const std::string& name, const double& value){global_defs.set(name, new Hdouble(value, "Last set by: " + command));}
	void global_store(const std::string& name, const float& value){global_defs.set(name, new Hfloat(value, "Last set by: " + command));}
	void global_store(const std::string& name, const std::string& value){global_defs.set(name, new Hstring(value, "Last set by: " + command));}
	void global_store(const std::string& name, const std::vector<double>& value){global_defs.set(name, new Hdvector(value, "Last set by: " + command));}

    }; 

    /** Function object used for testing for a variable name at the start of a string as in
     * 'name=7' given just 'name'. Case insensitive
     */
    class Check_start_arg{
    public:
	//! Constructor which stores the test string in upper case
	Check_start_arg(const std::string& chk_string) : check_string(Subs::toupper(chk_string)) {}
	//! Tests whether the internal string can be found in 'test_string'
	bool operator()(const std::string& test_string) {
	    std::string::size_type n = Subs::toupper(test_string).find(check_string);
	    return (n == 0);
	}
    private:
	std::string check_string;
    };

    /** Function object for checking the argument list. This is like Subs::Input::Check_start_arg 
     * but it tests for a complete (case-insensitive) match to the string
     */
    class Check_all_arg{
    public:
	//! Constructor which stores the test string
	Check_all_arg(const std::string& chk_string) : check_string(Subs::toupper(chk_string)) {}
	//! Checks whether 'test_string' is the same as the internal string
	bool operator()(const std::string& test_string) {
	    return (Subs::toupper(test_string) == check_string);
	}
    private:
	std::string check_string;
    };

    //! Structure storing variable's properties
    struct Varprop {
	//! Local or global default
	Input::Var_control1 store_where;
	//! Prompt characteristics
	Input::Var_control2 prompt;
    };

    //! Gets the local and global default file names
    void get_default_files(const char* ENV_DIR, const char* DEF_DIR, const std::string& command, std::string& global, std::string& local);

    //! Support function to get over problem of too low a precision for doubles and floats
    template <class T>    
    void print_value(const T& value, std::ostream& ostr){
	ostr << value;
    }

    //! Specialisation for floats
    template <> 
    inline void print_value<float>(const float& value, std::ostream& ostr){
	ostr << std::setprecision(7) << value << std::setprecision(0);
    }

    //! Specialisation for doubles
    template <> 
    inline void print_value<double>(const double& value, std::ostream& ostr){
	ostr << std::setprecision(17) << value << std::setprecision(0);
    }
  

};

/** This function supports range-checked user input. That is, input of doubles, ints etc with
 * a minimum and maximum that will be tested.
 * \param name the name of the variable (input)
 * \param value its value (returned)
 * \param initial_value its initial value
 * \param low the lowest value it can have
 * \param high the highest value it can have
 * \param prompt_string the string to use in any prompt
 */
template <class T>
void Subs::Input::get_value(std::string name, T& value, const T& initial_value, const T& low, const T& high,  
			    const std::string& prompt_string){

    bool check = true;

    // Convert variable to upper-case
    name = Subs::toupper(name);
    
    // Check to see if the variable has been correctly signed in
    std::map<std::string,Varprop>::const_iterator mp = find(name);
    if(mp == end()) throw Input_Error("Variable = " + name + " has not been signed in");
    
    // No point carrying on if zero range
    if(low == high){

	value = low;

	// Store new value in defaults
	if(mp->second.store_where == LOCAL){
	    local_store(name, value);
	}else{
	    global_store(name, value);
	}
    
	// Report value if wanted
	if(list_values){
	    std::cout << name << " = ";
	    print_value(value, std::cout);
	    std::cout << std::endl;
	}
	return;
    }
  
    // Scan named_args to see if a variable of this name has been specified
    std::vector<std::string>::iterator vsi = find_if(named_args.begin(), named_args.end(), Check_start_arg(name + "="));
    if(vsi == named_args.end()){
    
	// No named value, see if there is a value by position
	if(mp->second.prompt == PROMPT && narg < simple_args.size()){

	    // A '\' is interpreted meaning take the default for this and all subsequent arguments (unless
	    // it is a named argument
      
	    if(simple_args[narg] == "\\"){
		value = initial_value;
		Header::Hnode* node;
		if(mp->second.store_where == LOCAL)
		    node = local_defs.find(name);
		else
		    node = global_defs.find(name);
		if(node->has_data()) node->value->get_value(value);
		set_to_default = true;
		value = value < low ? low : (value > high ? high : value);

	    }else if(Subs::toupper(simple_args[narg]) == "MIN"){
		value = low;
	
	    }else if(Subs::toupper(simple_args[narg]) == "MAX"){
		value = high;

	    }else if(simple_args[narg].find('@') == 0){
		simple_args[narg].erase(0,1);
		get_from_global(simple_args[narg], value);	

	    }else{

		std::istringstream ist(simple_args[narg]);
		ist >> value;
		if(!ist) throw Input_Error("Error translating " + simple_args[narg] + " for input = " + name);
		ist.clear();

	    }

	}else{

	    // No value specified on the command line, or set to default. Determine the default
	    T default_value = initial_value;
	    Header::Hnode* node;
	    if(mp->second.store_where == LOCAL)
		node = local_defs.find(name);
	    else
		node = global_defs.find(name);
	    if(node->has_data()) node->value->get_value(default_value);
	   
	    default_value = default_value < low ? low : (default_value > high ? high : default_value);

	    if(mp->second.prompt == PROMPT && !set_to_default){
		std::string entry = "?";
		std::cout << name << " - " << prompt_string << " [";
		print_value(default_value,std::cout);
		std::cout << "]: " << std::flush;
		while(entry == "?"){
		    getline(std::cin,entry);

		    if(entry == ""){
			value = default_value;
		    }else if(entry == "\\"){
			value = default_value;
			set_to_default = true;
		    }else if(entry == "?"){
			std::cout << "Parameter = " << name << ", range: ";
			print_value(low, std::cout);
			std::cout << " to ";
			print_value(high, std::cout);
			std::cout  << std::endl;
		    }else if(Subs::toupper(entry) == "MIN"){
			value = low;
		    }else if(Subs::toupper(entry) == "MAX"){
			value = high;
		    }else	if(entry.find('@') == 0){
			entry.erase(0,1);
			get_from_global(entry, value);
		    }else{

			std::string buff;
			std::istringstream ist(buff);
			if(entry.find(">>") == 0){
			    ist.str(entry.substr(2));
			    check = false;
			}else{
			    ist.str(entry);
			}
			ist >> value;
			if(!ist) throw Input_Error("Error translating " + entry + " for input = " + name);
			ist.clear();
		    }
		}
	    }else{
		value = default_value;
	    }
	}

	// Increment argument counter
	if(mp->second.prompt == PROMPT) narg++;
      
    }else{

	// Value named on command line. Translate. Argument counter not
	// incremented so that named arguments do not count for argument counting.

	std::string entry = vsi->substr(name.length()+1);
	std::string buff;
	std::istringstream ist(buff);
	if(entry.find(">>") == 0){
	    ist.str(entry.substr(2));
	    ist >> value;
	    if(!ist) 
		throw Input_Error("Error translating " + vsi->substr(name.length()+1) + " for input = " + name);
	    ist.clear();
	    check = false;
	}else if(Subs::toupper(entry) == "MIN"){
	    value = low;
	}else if(Subs::toupper(entry) == "MAX"){
	    value = high;
	}else if(entry.find('@') == 0){
	    entry.erase(0,1);
	    get_from_global(entry, value);    
	}else{
	    ist.str(entry);
	    ist >> value;
	    if(!ist) 
		throw Input_Error("Error translating " + vsi->substr(name.length()+1) + " for input = " + name);
	    ist.clear();
	}

    }

    // Range check
    if(check && (value < low || value > high)) 
	throw Input_Error("Value = " + Subs::str(value) + " out of range " + Subs::str(low) + " to " + Subs::str(high) + " of input = " + name);

    // Store new value in defaults
    if(mp->second.store_where == LOCAL){
	local_store(name, value);
    }else{
	global_store(name, value);
    }

    // Report value if wanted
    if(list_values){
	std::cout << name << " = ";
	print_value(value, std::cout);
	std::cout << std::endl;
    }

}

/** This function allows you to override any default loaded from disk files
 * which can be useful if there is a preferable run-time default. 
 * \param name the name of the variable (input)
 * \param value the new default
 */

template <class T>
void Subs::Input::set_default(std::string name, const T& value){

    // Convert variable to upper-case  
    name = Subs::toupper(name);

    // Check to see if the variable has been correctly signed in
    std::map<std::string,Varprop>::const_iterator mp = find(name);
    if(mp == end()) throw Input_Error("Variable = " + name + " has not been signed in");
  
    // Set default
    if(mp->second.store_where == LOCAL){
	local_store(name, value);
    }else{
	global_store(name, value);
    }
}

/** This function allows you to add a value to the global defaults even if
 * it is not an officially 'signed in' variable. The name is stored with its
 * case preserved.
 * \param name the name of the variable
 * \param value the new default
 */

template <class T>
void Subs::Input::add_to_global(const std::string& name, const T& value){

    global_store(name, value);

}

/** This function allows you to retrieve a value to the global defaults even if
 * it is not an officially 'signed in' variable. The name is case sensitive. This
 * version is for when you (think that you) know the type of the value to be returned
 * \param name the name of the variable (input)
 * \param value the new default
 * \exception Throws a Subs_Error if the name cannot be found or if the value is not
 * retrievable for the type of 'value'
 */
template <class T>
void Subs::Input::get_from_global(const std::string& name, T& value) const {
    try{
	global_defs[name]->get_value(value);
    }
    catch(const Subs_Error& err){
	throw Subs_Error("void Subs::Input::get_from_global(const std::string&, T&): " + err);
    }
}

