#include <fstream>
#include <sstream>
#include <algorithm>
#include "trm/subs.h"
#include "trm/input.h"
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/**
 * Constructor -- reads in two files specifying the default values of the variables and stores
 * command-line arguments
 * \param global_file the name of the file where the global defaults are stored.
 * \param local_file  the name of the file where the local defaults are stored.
 * \param argc the number of command line arguments
 * \param argv the command line arguments
 * \exception Throws a Subs::Input::Input_Error if the defaults cannot be loaded or
 * if special keywords such as 'prompt' are specified more than once.
 */
Subs::Input::Input(const std::string& global_file, const std::string& local_file, int argc, char *argv[]) 
    : global_file_(global_file), local_file_(local_file), narg(0), set_to_default(false), 
      list_values(false), prompt_all(false), def_files(true) {

    // store the command name, strip prefix up to final '/'
    command = argv[0];
    std::string::size_type slash = command.find_last_of('/');
    if(slash != std::string::npos)
	command = command.substr(slash+1);

    // Store the rest of the arguments
    for(int i=1; i<argc; i++) simple_args.push_back(argv[i]);

    // Check for special arguments, 'list', 'prompt' and 'initialise'. If found, set respective
    // flags and remove from argument list. Guard against repetition
    std::vector<std::string>::iterator vsi;

    // 'list' means list all input values
    vsi = find_if(simple_args.begin(), simple_args.end(), Check_all_arg("list"));
    if(vsi != simple_args.end()){
	list_values = true;
	simple_args.erase(vsi);
    }
    vsi = find_if(simple_args.begin(), simple_args.end(), Check_all_arg("list"));
    if(vsi != simple_args.end())
	throw Input_Error("Cannot repeat keyword 'list' in argument list");

    // 'prompt' forces prompting of all values 
    vsi = find_if(simple_args.begin(), simple_args.end(), Check_all_arg("prompt"));
    if(vsi != simple_args.end()){
	prompt_all = true;
	simple_args.erase(vsi);
    }
    vsi = find_if(simple_args.begin(), simple_args.end(), Check_all_arg("prompt"));
    if(vsi != simple_args.end())
	throw Input_Error("Cannot repeat keyword 'prompt' in argument list");

    // 'nodefs' means that the default files will not be accessed
    vsi = find_if(simple_args.begin(), simple_args.end(), Check_all_arg("nodefs"));
    if(vsi != simple_args.end()){
	def_files = false;
	simple_args.erase(vsi);
    }
    vsi = find_if(simple_args.begin(), simple_args.end(), Check_all_arg("nodefs"));
    if(vsi != simple_args.end())
	throw Input_Error("Cannot repeat keyword 'nodefs' in argument list");

    // Read the defaults
    if(def_files){
	std::ifstream istr(global_file_.c_str(), std::ios::binary);
	if(istr){
	    try{
		global_defs.read(istr,false);
	    }
	    catch(...){
		throw Input_Error(
		    std::string("Error occurred trying to load global defaults from \"") + global_file_ +
		    std::string("\".\nIt may be corrupted and you might want to delete it."));
	    }
	}
	istr.close();
    
	istr.open(local_file_.c_str(), std::ios::binary);
	if(istr){
	    try{
		local_defs.read(istr,false);
	    }
	    catch(...){
		throw Input_Error(
		    std::string("Error occurred trying to load local defaults from \"") + local_file_ +
		    std::string("\".\nIt may be corrupted and you might want to delete it."));
	    }
	}
	istr.close();
    }
}

/*
 * This constructor works out the standard file names of the local and global default files 
 * from an environment variable and a default if the environment variable does not exist.
 * It is called with reversed arguments to the other constructor to give a unique signature.
 * \param argc the number of command line arguments
 * \param argv the command line arguments
 * \param directory_environment_var the name of the environment variable specifying the directory where the
 * default files are to be found
 * default_directory the name of the directory to search if the environment variable is not set
 */
Subs::Input::Input(int argc, char *argv[], const std::string& directory_environment_var, const std::string& default_directory)
    : narg(0), set_to_default(false), list_values(false), prompt_all(false), def_files(true) {

    // get the file names
    get_default_files(directory_environment_var.c_str(), default_directory.c_str(), argv[0], global_file_, local_file_);

    // Same code as other constructor from now. More efficient to repeat it than to try to use other constructor
    // and copying since that involves unecessary destruction and re-loading.

    // store the command name
    command = argv[0];
    std::string::size_type slash = command.find_last_of('/');
    if(slash != std::string::npos)
	command = command.substr(slash+1);

    // Store all arguments except the first which is the program name
    for(int i=1; i<argc; i++) simple_args.push_back(argv[i]);

    // Check for special arguments, 'list', 'prompt' and 'nodefs'. If found, set respective
    // flags and remove from argument list. Guard against repetition
    std::vector<std::string>::iterator vsi;

    // 'list' means list all input values
    vsi = find_if(simple_args.begin(), simple_args.end(), Check_all_arg("list"));
    if(vsi != simple_args.end()){
	list_values = true;
	simple_args.erase(vsi);
    }
    vsi = find_if(simple_args.begin(), simple_args.end(), Check_all_arg("list"));
    if(vsi != simple_args.end())
	throw Input_Error("Cannot repeat keyword 'list' in argument list");

    // 'prompt' forces prompting of all values 
    vsi = find_if(simple_args.begin(), simple_args.end(), Check_all_arg("prompt"));
    if(vsi != simple_args.end()){
	prompt_all = true;
	simple_args.erase(vsi);
    }
    vsi = find_if(simple_args.begin(), simple_args.end(), Check_all_arg("prompt"));
    if(vsi != simple_args.end())
	throw Input_Error("Cannot repeat keyword 'prompt' in argument list");

    // 'nodefs' means that the default files will not be accessed
    vsi = find_if(simple_args.begin(), simple_args.end(), Check_all_arg("nodefs"));
    if(vsi != simple_args.end()){
	def_files = false;
	simple_args.erase(vsi);
    }
    vsi = find_if(simple_args.begin(), simple_args.end(), Check_all_arg("nodefs"));
    if(vsi != simple_args.end())
	throw Input_Error("Cannot repeat keyword 'nodefs' in argument list");

    // Read the defaults
    if(def_files){
	std::ifstream istr(global_file_.c_str(), std::ios::binary);
	if(istr){
	    try{
		global_defs.read(istr,false);
	    }
	    catch(...){
		throw Input_Error(
		    std::string("Error occurred trying to load global defaults from \"") + global_file_ +
		    std::string("\".\nIt may be corrupted and you might want to delete it."));
	    }
	}
	istr.close();
    
	istr.open(local_file_.c_str(), std::ios::binary);
	if(istr){
	    try{
		local_defs.read(istr,false);
	    }
	    catch(...){
		throw Input_Error(
		    std::string("Error occurred trying to load local defaults from \"") + local_file_ +
		    std::string("\".\nIt may be corrupted and you might want to delete it."));
	    }
	}
	istr.close();
    }
}

/*
 * This constructor works out the standard file names of the local and global default files 
 * from an environment variable and a default if the environment variable does not exist.
 * It expects string arguments, starting with the name of the command.
 * \param args the command line arguments
 * \param directory_environment_var the name of the environment variable specifying the directory where the
 * default files are to be found
 * default_directory the name of the directory to search if the environment variable is not set
 */
Subs::Input::Input(const std::vector<std::string>& args, const std::string& directory_environment_var, const std::string& default_directory)
    : narg(0), set_to_default(false), list_values(false), prompt_all(false), def_files(true) {

    // get the file names
    get_default_files(directory_environment_var.c_str(), default_directory.c_str(), args[0], global_file_, local_file_);

    // Same code as other constructor from now. More efficient to repeat it than to try to use other constructor
    // and copying since that involves unecessary destruction and re-loading.

    // store the command name
    command = args[0];
    std::string::size_type slash = command.find_last_of('/');
    if(slash != std::string::npos)
	command = command.substr(slash+1);

    // Store all arguments except the first which is the program name
    for(size_t i=1; i<args.size(); i++) simple_args.push_back(args[i]);

    // Check for special arguments, 'list', 'prompt' and 'nodefs'. If found, set respective
    // flags and remove from argument list. Guard against repetition
    std::vector<std::string>::iterator vsi;

    // 'list' means list all input values
    vsi = find_if(simple_args.begin(), simple_args.end(), Check_all_arg("list"));
    if(vsi != simple_args.end()){
	list_values = true;
	simple_args.erase(vsi);
    }
    vsi = find_if(simple_args.begin(), simple_args.end(), Check_all_arg("list"));
    if(vsi != simple_args.end())
	throw Input_Error("Cannot repeat keyword 'list' in argument list");

    // 'prompt' forces prompting of all values 
    vsi = find_if(simple_args.begin(), simple_args.end(), Check_all_arg("prompt"));
    if(vsi != simple_args.end()){
	prompt_all = true;
	simple_args.erase(vsi);
    }
    vsi = find_if(simple_args.begin(), simple_args.end(), Check_all_arg("prompt"));
    if(vsi != simple_args.end())
	throw Input_Error("Cannot repeat keyword 'prompt' in argument list");

    // 'nodefs' means that the default files will not be accessed
    vsi = find_if(simple_args.begin(), simple_args.end(), Check_all_arg("nodefs"));
    if(vsi != simple_args.end()){
	def_files = false;
	simple_args.erase(vsi);
    }
    vsi = find_if(simple_args.begin(), simple_args.end(), Check_all_arg("nodefs"));
    if(vsi != simple_args.end())
	throw Input_Error("Cannot repeat keyword 'nodefs' in argument list");

    // Read the defaults
    if(def_files){
	std::ifstream istr(global_file_.c_str(), std::ios::binary);
	if(istr){
	    try{
		global_defs.read(istr,false);
	    }
	    catch(...){
		throw Input_Error(
		    std::string("Error occurred trying to load global defaults from \"") + global_file_ +
		    std::string("\".\nIt may be corrupted and you might want to delete it."));
	    }
	}
	istr.close();
    
	istr.open(local_file_.c_str(), std::ios::binary);
	if(istr){
	    try{
		local_defs.read(istr,false);
	    }
	    catch(...){
		throw Input_Error(
		    std::string("Error occurred trying to load local defaults from \"") + local_file_ +
		    std::string("\".\nIt may be corrupted and you might want to delete it."));
	    }
	}
	istr.close();
    }
}


/**
 * This function must be called once for each variable immediately after the construction 
 * an Input object. It defines the properties of each variable. The variables do not need to
 * be signed in in any order, although standardly I do so in the command order.
 * \param varname      the variable name
 * \param store_where  where the default is stored (local or global)
 * \param prompt       the prompting characteristic of the variable
 */
void Subs::Input::sign_in(const std::string& varname, Var_control1 store_where, Var_control2 prompt){
    Varprop vprop;
    vprop.store_where = store_where;
    if(prompt_all){
	vprop.prompt = PROMPT;
    }else{
	vprop.prompt = prompt;
    }
    std::pair<std::map<std::string,Varprop>::iterator,bool> pins = insert(std::make_pair(toupper(varname),vprop));
    if(!pins.second)
	throw Input_Error("Subs::Input::sign_in(const std::string&, Var_control1, Var_control2): attempt to sign_in " + varname + " twice.");

    // Now look through command arguments for any of the form 'variable=value' associated with this variable
    std::vector<std::string>::iterator vsi;
    vsi = find_if(simple_args.begin(), simple_args.end(), Check_start_arg(varname + "="));
    if(vsi != simple_args.end()){
	std::vector<std::string>::iterator vsi1 = vsi+1, vsi2 = find_if(vsi1, simple_args.end(), Check_start_arg(varname + "="));
	if(vsi2 != simple_args.end())
	    throw Input_Error("bool Subs::Input::sign_in(const std::string&, Var_control1, Var_control2): argument = " + varname + " was specified more than once.");

	// Store in named_args,  remove from simple_args
	named_args.push_back(*vsi);
	simple_args.erase(vsi);
    }
}

/** Destructor writes out the current defaults over the original files.
 * This has the effect of saving the defaults when a program exits correctly.
 */
Subs::Input::~Input(){

    if(def_files){
	std::ofstream ostr(global_file_.c_str(), std::ios::binary);
	if(ostr){
	    try{
		global_defs.write(ostr);
	    }
	    catch(...){
		throw Input_Error(
		    std::string("Error occurred trying to write global defaults to \"") + global_file_ +
		    std::string("\".\nIt may be corrupted and you might want to delete it."));
	    }
	}else{
	    std::cerr << "Failed to open \"" << global_file_ << "\" for storage of global defaults which will"
		      << " therefore not be updated." << std::endl;
	}
	ostr.close();
	
	ostr.open(local_file_.c_str(), std::ios::binary);
	if(ostr){
	    try{
		local_defs.write(ostr);
	    }
	    catch(...){
		throw Input_Error(
		    std::string("Error occurred trying to load local defaults from \"") + local_file_ +
		    std::string("\".\nIt may be corrupted and you might want to delete it."));
	    }
	}else{
	    std::cerr << "Failed to open \"" << local_file_ << "\" for storage of local defaults which will"
		      << " therefore not be updated." << std::endl;
	}
	ostr.close();
    }

}


/**
 * This function saves the current value of the defaults to the standard files. Although
 * this happens anyway when the Input goes out of scope (e.g. when the program exits), there are
 * cases where one often wants to terminate a program with Ctrl-C which prevents the defaults being saved.
 * A call to this function will remedy the problem.
 */
void Subs::Input::save() const {
    if(def_files){
	std::ofstream ostr(global_file_.c_str(), std::ios::binary);
	if(ostr){
	    try{
		global_defs.write(ostr);
	    }
	    catch(...){
		throw Input_Error(
		    std::string("void Subs::Input::save(): error occurred trying to write global defaults to \"") + 
		    global_file_ +
		    std::string("\".\nIt may be corrupted and you might want to delete it."));
	    }
	}else{
	    std::cerr << "void Subs::Input::save(): failed to open \"" << global_file_ 
		      << "\" for storage of global defaults which will therefore not be updated." << std::endl;
	}
	ostr.close();
    
	ostr.open(local_file_.c_str(), std::ios::binary);
	if(ostr){
	    try{
		local_defs.write(ostr);
	    }
	    catch(...){
		throw Input_Error(
		    std::string("void Subs::Input::save(): error occurred trying to load local defaults from \"") + 
		    local_file_ +
		    std::string("\".\nIt may be corrupted and you might want to delete it."));
	    }
	}else{
	    std::cerr << "void Subs::Input::save(): failed to open \"" << local_file_ 
		      << "\" for storage of local defaults which will therefore not be updated." << std::endl;
	}
	ostr.close();
    }
}

/** Gets the value of a string. Allows an intial values, but makes no checks on the nature of the string
 * at all. 
 * \param name name of the variable associated with the string
 * \param value returned string value
 * \param initial_value default value when no default exists
 * \param prompt_string the prompt to use when necessary
 * \exception Subs::Input::Input_Error exceptions can be thrown.
 */
void Subs::Input::get_value(std::string name, std::string& value, const std::string& initial_value, const std::string& prompt_string){

    // Convert variable to upper-case
    name = toupper(name);
  
    // Check to see if the variable has been correctly signed in
    std::map<std::string,Varprop>::const_iterator mp = find(name);
    if(mp == end()) throw Input_Error("Variable = " + name + " has not been signed in");

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
		if(node->has_data()) value = node->value->get_string();

		set_to_default = true;
		
	    }else if(simple_args[narg].find('@') == 0 && is_a_global(simple_args[narg].substr(1))){
		get_from_global(simple_args[narg].substr(1), value);
		
	    }else{
		value = simple_args[narg];
	    }
	    
	}else{
	    
	    // No value specified on the command line, set the default and prompt if required
	    std::string default_value = initial_value;
	    Header::Hnode* node;
	    if(mp->second.store_where == LOCAL)
		node = local_defs.find(name);
	    else
		node = global_defs.find(name);
	    if(node->has_data()) default_value = node->value->get_string();
	    
	    if(mp->second.prompt == PROMPT && !set_to_default){
		std::string entry = "?";
		while(entry == "?"){
		    std::cout << name << " - " << prompt_string << " [" << default_value << "]: " << std::flush;
		    getline(std::cin,entry);
		    if(entry == ""){
			value = default_value;
		    }else if(entry == "\\"){
			value = default_value;
			set_to_default = true;
		    }else if(entry == "?"){
			std::cout << "Parameter = " << name << ", default = " << default_value << std::endl;
		    }else if(entry.find('@') == 0 && is_a_global(entry.substr(1))){
			get_from_global(entry.substr(1), value);
		    }else{
			value = entry;
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
	value = vsi->substr(name.length()+1);
	
    }
    
    // Store new value in defaults
    if(mp->second.store_where == LOCAL){
	local_defs.set(name, new Hstring(value, "Last set by: " + command));
    }else{
	global_defs.set(name, new Hstring(value, "Last set by: " + command));
    }
    
    // Report value if wanted
    if(list_values) std::cout << name << " = " << value << std::endl;
    
}

/** Gets the value of a booleans true/false variable.
 * \param name name of the variable associated with the boolean
 * \param value returned boolean value
 * \param initial_value default value when no default exists
 * \param prompt_string the prompt to use when necessary
 * \exception Subs::Input::Input_Error exceptions can be thrown.
 */
void Subs::Input::get_value(std::string name, bool& value, const bool& initial_value, const std::string& prompt_string){
  
    // Convert variable to upper-case
    name = toupper(name);

    // Check to see if the variable has been correctly signed in
    std::map<std::string,Varprop>::const_iterator mp = find(name);
    if(mp == end()) throw Input_Error("Variable = " + name + " has not been signed in");

    // First scan named_args to see if a variable of this name has been specified. 
    std::vector<std::string>::iterator vsi = find_if(named_args.begin(), named_args.end(), Check_start_arg(name + "="));
    if(vsi == named_args.end()){

	// No named value, see if there is a value by position (only for promted variables)
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
	
	    }else if(Subs::toupper(simple_args[narg]) == "YES" || Subs::toupper(simple_args[narg]) == "Y" ||
		     Subs::toupper(simple_args[narg]) == "TRUE"){
		value = true;
	
	    }else if(Subs::toupper(simple_args[narg]) == "NO" || Subs::toupper(simple_args[narg]) == "N" ||
		     Subs::toupper(simple_args[narg]) == "FALSE"){
		value = false;
	
	    }else{
		throw Input_Error(std::string("Error translating ") + simple_args[narg] + std::string(" for input = ") + name);
	    }
      
	}else{
	
	    // No value specified on the command line, set the default
	
	    bool default_value = initial_value;
	    Header::Hnode *node;
	    if(mp->second.store_where == LOCAL)
		node = local_defs.find(name);
	    else
		node = global_defs.find(name);
	    if(node->has_data()) default_value = node->value->get_bool();
      
	    if(mp->second.prompt == PROMPT && !set_to_default){
		std::string entry = "?";
		while(entry == "?"){
		    if(default_value){
			std::cout << name << " - " << prompt_string << " [yes]: " << std::flush;
		    }else{
			std::cout << name << " - " << prompt_string << " [no]: " << std::flush;
		    }
		    getline(std::cin,entry);
		    if(entry == ""){
			value = default_value;
		    }else if(entry == "\\"){
			value = default_value;
			set_to_default = true;
		    }else if(entry == "?"){
			if(default_value){
			    std::cout << "Parameter = " << name << ", default = yes" << std::endl;
			}else{
			    std::cout << "Parameter = " << name << ", default = no" << std::endl;
			}
		    }else{
			if(toupper(entry) == "YES" || toupper(entry) == "Y" || toupper(entry) == "TRUE"){
			    value = true;
			}else if(toupper(entry) == "NO" || toupper(entry) == "N" || toupper(entry) == "FALSE"){
			    value = false;
			}else{
			    std::cerr << "Must answer one of 'yes', 'y', 'true', 'no', 'n' or 'false'." << std::endl;
			    entry = "?";
			}
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

	std::string entry = Subs::toupper(vsi->substr(name.length()+1));
	if(entry == "YES" || entry == "Y" || entry == "TRUE"){
	    value = true;
	}else if(entry == "NO" || entry == "N" || entry == "FALSE"){
	    value = false;
	}else{
	    throw Input_Error(std::string("Error translating ") + vsi->substr(name.length()+1) + std::string(" for input = ") + name);
	}
    }

    // Store new value in defaults
    if(mp->second.store_where == LOCAL){
	local_defs.set(name, new Hbool(value, "Last set by: " + command));
    }else{
	global_defs.set(name, new Hbool(value, "Last set by: " + command));
    }

    // Report value if wanted
    if(list_values){
	if(value){
	    std::cout << name << " = yes" << std::endl;
	}else{
	    std::cout << name << " = no" << std::endl;
	}
    }

}

/**
 * This gets a character value, checking that it comes from one of a specified set of allowed values.
 * This is useful for specifying one of a range of options. It is case sensitive, so you must specify
 * a character twice as in "aA" to get case-insentive behaviour.
 * \param name variable name for the character value
 * \param value the value returned
 * \param initial_value value in the absence of any default
 * \param allowed_values string of possible characters as in "aAdDfF"
 * \param prompt_string the prompt to use when necessary
 */
void Subs::Input::get_value(std::string name, char& value, const char& initial_value, const std::string& allowed_values, const std::string& prompt_string){

    // Convert variable to upper-case
    name = toupper(name);

    // Check to see if the variable has been correctly signed in
    std::map<std::string,Varprop>::const_iterator mp = find(name);
    if(mp == end()) throw Input_Error("Variable = " + name + " has not been signed in");
   
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
		if(node->has_data()) value = node->value->get_char();

		if(allowed_values.find(value) == std::string::npos) value = initial_value;
		set_to_default = true;

	    }else{

		std::istringstream istr(simple_args[narg]);
		istr >> value;
		if(!istr) throw Input_Error(std::string("Error translating ") + simple_args[narg] + std::string(" for input = ") + name);

	    }

	}else{

	    // No value specified on the command line, set the default
	    char default_value = initial_value;
	    Header::Hnode* node;
	    if(mp->second.store_where == LOCAL)
		node = local_defs.find(name);
	    else
		node = global_defs.find(name);
	    if(node->has_data()) default_value = node->value->get_char();

	    if(allowed_values.find(default_value) == std::string::npos) default_value = initial_value;

	    if(mp->second.prompt == PROMPT && !set_to_default){
		std::string entry = "?";
		while(entry == "?"){
		    std::cout << name << " - " << prompt_string << " [" << default_value << "]: " << std::flush;
		    getline(std::cin,entry);
		    if(entry == ""){
			value = default_value;
		    }else if(entry == "\\"){
			value = default_value;
			set_to_default = true;
		    }else if(entry == "?"){
			std::cout << "Parameter = " << name 
				  << " is a single character variable with the following allowed values: " 
				  << allowed_values << std::endl;
		    }else{
			std::istringstream istr(entry);
			istr >> value;
			if(!istr) throw Input_Error(std::string("Error translating ") + entry +
						    std::string(" for input = ") + name);
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
	std::istringstream istr(vsi->substr(name.length()+1));
	istr >> value;
	if(!istr) throw Input_Error(std::string("Error translating ") + vsi->substr(name.length()+1) +
				    std::string(" for input = ") + name);

    }

    // Check against allowed values
    if(allowed_values.find(value) == std::string::npos)
	throw Input_Error("The only allowed values of the variable " + name + 
			  " are one of the following characters: [" + allowed_values + "]");

    // Store new value in defaults
    if(mp->second.store_where == LOCAL){
	local_defs.set(name, new Hchar(value, "Last set by: " + command));
    }else{
	global_defs.set(name, new Hchar(value, "Last set by: " + command));
    }

    // Report value if wanted
    if(list_values) std::cout << name << " = " << value << std::endl;

}

/** This function supports range-checked user input of a list of doubles with
 * a minimum and maximum that will be tested.
 * \param name the name of the variable (input)
 * \param value its value (returned)
 * \param initial_value its initial value
 * \param low the lowest value it can have
 * \param high the highest value it can have
 * \param nvec number of elements in the vector, 0 if not fixed
 * \param prompt_string the string to use in any prompt
 */
void Subs::Input::get_value(std::string name, std::vector<double>& value, const double& initial_value, 
			    const double& low, const double& high, const std::vector<double>::size_type nvec,
			    const std::string& prompt_string){

    // Convert variable to upper-case
    name = Subs::toupper(name);

    // Check to see if the variable has been correctly signed in
    std::map<std::string,Varprop>::const_iterator mp = find(name);
    if(mp == end()) throw Input_Error("Variable = " + name + " has not been signed in");

    // Scan named_args to see if a variable of this name has been specified
    std::vector<std::string>::iterator vsi = find_if(named_args.begin(), named_args.end(), Check_start_arg(name + "="));
    if(vsi == named_args.end()){
     
	// No named value, see if there is a value by position
	if(mp->second.prompt == PROMPT && narg < simple_args.size()){

	    // A '\' is interpreted meaning take the default for this and all subsequent arguments (unless
	    // it is a named argument
      
	    if(simple_args[narg] == "\\"){
		if(nvec > 0){
		    value.resize(nvec);
		    for(std::vector<double>::size_type i=0; i<value.size(); i++)
			value[i] = initial_value;
		}
		Header::Hnode* node;
		if(mp->second.store_where == LOCAL)
		    node = local_defs.find(name);
		else
		    node = global_defs.find(name);
		if(node->has_data()){
		    node->value->get_value(value);
		}else if(nvec == 0){
		    throw Input_Error(std::string("Error setting input for ") + name + std::string("; need a default for unspecified number of doubles"));
		}

		set_to_default = true;

		if(nvec > 0 && nvec == value.size()){
		    for(std::vector<double>::size_type i=0; i<value.size(); i++)
			value[i] = value[i] < low ? low : (value[i] > high ? high : value[i]);
		}else{
		    throw Input_Error(std::string("Error setting input for ") + name + std::string("; invalid default"));
		}

	    }else{

		// Blank added to avoid apparent bug in istringstream
		std::istringstream ist(simple_args[narg] + " ");
		std::vector<double> temp;
		double d;
		while(!ist.eof()){
		    ist >> d;
		    if(ist.good()){
			temp.push_back(d);
		    }else if(!ist.eof()){
			throw Input_Error(std::string("Error translating ") + simple_args[narg] + std::string(" for input = ") + name);
		    }
		}
		if(nvec > 0 && nvec != temp.size())
		    throw Input_Error("Error setting input " + name + "; expected " + Subs::str(nvec) + " elements.");

		value = temp;
	    }

	}else{

	    // No value specified on the command line, or set to default. Determine the default
	    std::vector<double> default_value;
	    default_value.push_back(initial_value);
	    Header::Hnode* node;
	    if(mp->second.store_where == LOCAL)
		node = local_defs.find(name);
	    else
		node = global_defs.find(name);
	    if(node->has_data()) node->value->get_value(default_value);

	    if(nvec > 0 && nvec != default_value.size()){
		default_value.resize(nvec);
		for(std::vector<double>::size_type i=0; i<default_value.size(); i++)
		    default_value[i] = initial_value;
	    }

	    for(std::vector<double>::size_type i=0; i<default_value.size(); i++)
		default_value[i] = default_value[i] < low ? low : (default_value[i] > high ? high : default_value[i]);

	    if(mp->second.prompt == PROMPT && !set_to_default){
		std::ostringstream ostr;
		for(std::vector<double>::size_type i=0; i<default_value.size()-1; i++)
		    ostr << default_value[i] << " ";
		ostr << default_value[default_value.size()-1];

		std::string entry = "?";
		std::cout << name << " - " << prompt_string << " [" << ostr.str() << "]: " << std::flush;
		while(entry == "?"){
		    getline(std::cin,entry);

		    if(entry == ""){
			value = default_value;
		    }else if(entry == "\\"){
			value = default_value;
			set_to_default = true;
		    }else if(entry == "?"){
			std::cout << "Parameter = " << name << ", range: " << low << " to " << high << std::endl;
		    }else{
			std::istringstream ist(entry + " ");
			std::vector<double> temp;
			double d;
			while(!ist.eof()){
			    ist >> d;
			    if(ist.good()){
				// Range check
				if(d < low || d > high) 
				    throw Input_Error("Value = " + Subs::str(d) + " out of range " + Subs::str(low) + " to " + Subs::str(high) + " of input = " + name);
				temp.push_back(d);
			    }else if(!ist.eof()){
				throw Input_Error("Error translating " + entry + " for input = " + name);
			    }
			}
			if(nvec > 0 && nvec != temp.size())
			    throw Input_Error("Error setting input " + name + "; expected " + Subs::str(nvec) + " elements.");
			value = temp;
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
	std::istringstream ist(entry + " ");
	std::vector<double> temp;
	double d;
	while(!ist.eof()){
	    ist >> d;
	    if(ist.good()){
		// Range check
		if(d < low || d > high) 
		    throw Input_Error("Value = " + Subs::str(d) + " out of range " + Subs::str(low) + " to " + Subs::str(high) + " of input = " + name);
		temp.push_back(d);
	    }else if(!ist.eof()){
		throw Input_Error("Error translating " + entry + " for input = " + name);
	    }
	}
	if(nvec > 0 && nvec != temp.size())
	    throw Input_Error("Error setting input " + name + "; expected " + Subs::str(nvec) + " elements.");
	value = temp;

    }

    // Store new value in defaults
    if(mp->second.store_where == LOCAL){
	local_store(name, value);
    }else{
	global_store(name, value);
    }

    // Report value if wanted
    if(list_values){
	std::ostringstream ostr;
	for(std::vector<double>::size_type i=0; i<value.size()-1; i++)
	    ostr << value[i] << " ";
	ostr << value[value.size()-1]; 
	std::cout << name << " = " << ostr.str() << std::endl;
    }

}


/**
 * Helper routine to create file names for storage and retrieval of command defaults and also checks
 * existence of directory. Needs the name of an environment variable and a default. Essentially
 * it takes the command name and prefixes it with a directory name, either obtained from an environment
 * variable, or failing that, a default name.
 * directory in case the environment variable is not set. 
 * \param ENV_DIR  environment variable which should point to directory where default files are to be stored
 * \param DEF_DIR  alternative directory if ENV_DIR is not set.
 * \param command  command name. To allow for complicated paths, only the part after the last '/' will be used.
 * \param global   name of file for storage of global default values, returned.
 * \param local    name of file for storage of default values local to the program
 */
void Subs::get_default_files(const char* ENV_DIR, const char* DEF_DIR, const std::string& command, std::string& global, std::string& local){

    // Get (and create if ncessary) the directory for default files
    char *cp = getenv(ENV_DIR);
    std::string defaults_dir;
    if(cp == NULL){
	char *home = getenv("HOME");
	if(home == NULL) throw Subs_Error("Can't identify home directory to locate defaults directory");
	defaults_dir = std::string(home) + std::string("/") + std::string(DEF_DIR);

    }else{
	defaults_dir = cp;
    }

    struct stat info;
    if(stat(defaults_dir.c_str(), &info) != 0){
	if(mkdir(defaults_dir.c_str(), 0700) != 0)
	    throw Subs_Error("Failed to create input default directory = " + defaults_dir);
    }
    
    global = defaults_dir + std::string("/GLOBAL.def");
    std::string::size_type n;
    if((n = command.find_last_of('/')) == std::string::npos){
	local  = defaults_dir + std::string("/") + command + std::string(".def");
    }else{
	local  = defaults_dir + std::string("/") + command.substr(n+1) + std::string(".def");
    }
}

/** This function allows you to retrieve a value to the global defaults even if
 * it is not an officially 'signed in' variable. The name is case sensitive. This
 * version is for when you (think that you) know the type of the value to be returned
 * \param name the name of the variable, cased sensitive
 * \return Returns an Hitem pointer.
 * \exception Throws a Subs_Error if the name cannot be found
 */
Subs::Hitem* Subs::Input::get_from_global(const std::string& name) const {
    try{
	return global_defs[name];
    }
    catch(const Subs_Error& err){
	throw Subs_Error("Subs::Hitem* Subs::Input::get_from_global(const std::string&): " + err);
    }
}

/** This function checks on the existence of a variable amongst the global defaults
 * even when it is not a 'signed in' variable
 * \param name the name of the variable, case sensitive
 * \return true if found
 */
bool Subs::Input::is_a_global(const std::string& name) const {
    return global_defs.find(name)->has_data();
}

/** This routine prints out the global defaults to
 * an output stream.
 * \param ostr the output stream the data will be sent to
 */
void Subs::Input::print_global(std::ostream& ostr) const {
    ostr << global_defs;
}

/** This function allows you to copy a value within the global defaults
 * \param name the name of the variable
 * \param copy the name to copy to
 */
void Subs::Input::global_copy(const std::string& name, const std::string& copy){
    try{
	Subs::Hitem *hitem = global_defs[name];
	global_defs.set(copy, hitem->copy());
    }
    catch(const Subs_Error& err){
	throw Subs_Error("void Subs::Input::global_copy(const std::string&, const std::string&): " + err);
    }
}

/** This function allows you to override string defaults loaded from disk files
 * which can be useful if there is a preferable run-time default. Although there
 * is a template this version is neededto avoid having explicitly to put string()
 * around the value argument.
 * \param name the name of the variable (input)
 * \param value the new default
 */
void Subs::Input::set_default(std::string name, const char* value){

    // Convert variable to upper-case  
    name = Subs::toupper(name);

    // Check to see if the variable has been correctly signed in
    std::map<std::string,Varprop>::const_iterator mp = find(name);
    if(mp == end()) throw Input_Error("Variable = " + name + " has not been signed in");
   
    // Set default
    if(mp->second.store_where == LOCAL){
	local_store(name, std::string(value));
    }else{
	global_store(name, std::string(value));
    }
}
