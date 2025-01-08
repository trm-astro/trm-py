#include <ctype.h>
#include <string>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "trm/subs.h"

/** Makes an integer such as '12' into a string
 * of the form '012', with the number of digits
 * fixed by ndig. If ndig is less than the number
 * of digits neede, e.g. ndig=1, con=12, then you
 * get back normal formatting, i.e. '12'
 */

std::string Subs::str(const int& con, int ndig){
    std::ostringstream ss;
    ss.setf(std::ios::right);
    ss << std::setw(ndig) << std::setfill('0') << con;
    return ss.str();
}

/** Makes an integer such as '12' into a string
 * of the form '012', with the number of digits
 * fixed by ndig. If ndig is less than the number
 * of digits neede, e.g. ndig=1, con=12, then you
 * get back normal formatting, i.e. '12'
 */

std::string Subs::str(const long int& con, int ndig){
    std::ostringstream ss;
    ss.setf(std::ios::right);
    ss << std::setw(ndig) << std::setfill('0') << con;
    return ss.str();
}

/** Makes a digit (0-9) into a character
 * \param digit integer from 0 to 9
 * \return equivalent character
 * \exception Thros an exception if the integer is out of range
 */
char Subs::digit_to_char(int digit){

    if(digit < 0 || digit > 9)
        throw Subs_Error("Subs::digit_to_char: integer input = " + Subs::str(digit) + " is out of range 0 - 9");

    switch(digit){
	case 0:
	    return '0';
	case 1:
	    return '1';
	case 2:
	    return '2';
	case 3:
	    return '3';
	case 4:
	    return '4';
	case 5:
	    return '5';
	case 6:
	    return '6';
	case 7:
	    return '7';
	case 8:
	    return '8';
	default:
	    return '9';
    }
}

/* Convert a string to upper case
 * \param str the input string
 * \return the same string converted to upper case
 */
std::string Subs::toupper(const std::string& str){
    std::string temp;
    for(std::string::const_iterator p=str.begin(); p != str.end(); p++)
        temp += char(std::toupper(*p));
    return temp;
}

/* Convert a string to lower case
 * \param str the input string
 * \return the same string converted to lower case
 */
std::string Subs::tolower(const std::string& str){
    std::string temp;
    for(std::string::const_iterator p=str.begin(); p != str.end(); p++)
        temp += char(std::tolower(*p));
    return temp;
}

// writes out a string in a binary format. This requires writing
// first the number of characters and then the characters themselves.

void Subs::write_string(std::ofstream& s, const std::string& str){
    if(!s) return;
    UINT4 n = UINT4(str.length());
    s.write((char*)&n, sizeof(UINT4));
    if(!s) throw Subs_Error("Subs::write(std::ofstream&, const std::string&): error writing number of characters");
    if(n){
        s.write(str.data(), sizeof(char[n]));
        if(!s) throw Subs_Error("Subs::write(std::ofstream&, const std::string&): error writing characters");
    }
}

/** reads a string written in binary format.   
 * \param s the binary input stream
 * \param str the string read in
 */

void Subs::read_string(std::ifstream& s, std::string& str, bool swap_bytes){
    if(!s) return;
    UINT4 n;
    s.read((char*)&n, sizeof(UINT4));
    if(!s) throw Subs_Error("Subs::read(std::ifstream&, std::string&): error reading number of characters");
    if(swap_bytes) n = Subs::byte_swap(n);
    if(n > 0){
        char *p = new char [n];
        if(!p) throw Subs_Error("Subs::read(std::ifstream&, std::string&): error allocating array for characters");
        s.read(p, sizeof(char[n]));
        if(!s) throw Subs_Error("Subs::read(std::ifstream&, std::string&): error reading characters");
        str.assign(p, n);
        delete [] p;
    }
}

// skips a string written in binary format. 

void Subs::skip_string(std::ifstream& s, bool swap_bytes){
    if(!s) return;
    UINT4 n;
    s.read((char*)&n, sizeof(UINT4));
    if(!s) throw Subs_Error("Subs::skip(std::ifstream&): error reading number of characters");
    if(swap_bytes) n = Subs::byte_swap(n);
    if(n > 0){
        s.ignore(sizeof(char[n]));
        if(!s) throw Subs_Error("Subs::skip(std::ifstream&): error skipping characters");
    }
}

/** Reads a line from an input stream, splitting it up into separate lines
 * each of which is split into strings. The lines are split on semi-colons
 * if they are not either escaped or inside quotes. The reading ceases if
 * a newline or the end of the buffer is reached.
 * In the latter case 's' will be left in an error state.
 * The strings are split on blanks unless quoted with double quotes as 
 * in "a string". Escape the " as in \" to have them ignored. An entry such
 * as title="A plot" is preserved as a single string. 
 * \param s the input stream
 * \return a vector of vectors of strings. Leading and trailing blanks are removed unless
 * they occur inside double quotes.
 */

std::vector<std::vector<std::string> > Subs::read_multi_string_line(std::istream& s){

    if(!s) return std::vector<std::vector<std::string> >();
    std::vector<std::vector<std::string> > tvec;
    std::vector<std::string> temp;
    std::string tstr;
    bool started = false, escaped = false, insideQuotes = false;
    char c;
    while(s){
        s.get(c);
        if(!s || c == '\n') break;

        if(!escaped && c == '\\'){
            //	cerr << "not escaped, c = " << c << endl;
	    
            escaped = true;
	    
        }else{
	    
            if(escaped){
                //	std::cerr << "escaped, c = " << c << std::endl;

                escaped = false;
		
                if(started){
                    tstr += c;
                }else{
                    tstr         = c;
                    started      = true;
                }
		
            }else if(!insideQuotes && c == '"'){
                //	std::cerr << "not inside quotes, c = " << c << std::endl;
                insideQuotes = true;
                if(!started){
                    tstr         = "";
                    started      = true;
                }
		
            }else if(insideQuotes){
                //	std::cerr << "inside quotes, c = " << c << std::endl;
                if(c == '"'){
                    insideQuotes = false;
                }else{
                    tstr += c;
                }
		
            }else if(c == ';'){
                //	std::cerr << "semicolon" << std::endl;
                if(started){
                    temp.push_back(tstr);
                    started = false;
                }
                if(temp.size() > 0)
                    tvec.push_back(temp);
                temp.clear();
		
            }else if(started){
                //	std::cerr << "started, c = " << c << std::endl;
                if(c == ' ' || c == '\t'){
                    temp.push_back(tstr);
                    started = false;
                }else{
                    tstr += c;
                }
	
            }else if(!started && (c != ' ' && c != '\t') ){
                // std::cerr << "not started, c = " << c << std::endl;
                started = true;
                tstr = c;
            }
        }
    }
    if(started) 
        temp.push_back(tstr);
    
    if(temp.size() > 0)
        tvec.push_back(temp);
    
    return tvec;
}


/** Reads a line from an input stream. The line is split into strings on blanks 
 * unless quoted with double quotes as in "a string". Escape the " as in \" to have them ignored. An entry such
 * as title="A plot" is preserved as a single string. 
 * \param s the input stream
 * \return a vector of strings. Leading and trailing blanks are removed unless
 * they occur inside double quotes.
 */

std::vector<std::string> Subs::read_line(std::istream& s){

    std::vector<std::string> temp;
    std::string tstr;
    bool started = false, escaped = false, insideQuotes = false;
    char c;
    int nread = 0;
    while(s){
        s.get(c);
        if(!s){
            if(nread)
                break;
            else
                return std::vector<std::string>();
        }
        if(c == '\n') break;
        nread++;

        if(!escaped && c == '\\'){
  
            escaped = true;

        }else{

            if(escaped){

                escaped = false;

                if(started){
                    tstr += c;
                }else{
                    tstr         = c;
                    started      = true;
                }

            }else if(!insideQuotes && c == '"'){

                insideQuotes = true;
                if(!started){
                    tstr         = "";
                    started      = true;
                }
	
            }else if(insideQuotes){

                if(c == '"'){
                    insideQuotes = false;
                }else{
                    tstr += c;
                }
      
            }else if(started){

                if(c == ' ' || c == '\t'){
                    temp.push_back(tstr);
                    started = false;
                }else{
                    tstr += c;
                }
	
            }else if(!started && (c != ' ' && c != '\t') ){

                started = true;
                tstr = c;
            }
        }
    }
    if(started) 
        temp.push_back(tstr);

    return temp;
}

/** Strips trailing white space from a string
 * \param str input string
 * \return same string with no trailing blanks.
 */
std::string Subs::strip_trailing_whitespace(const std::string& str){
    std::string::size_type n = str.length();

    std::string::size_type j = n;
    while(j > 0 && (str[j-1] == ' ' || str[j-1] == '\t' || str[j-1] == '\n'))
        j--;

    if(j == n){
        return str;
    }else{
        std::string temp = str;
        temp.erase(j);
        return temp;
    }
}

/** Strips trailing white space from a string
 * \param str input string
 * \return same string with no trailing blanks.
 */
std::string Subs::strip_leading_whitespace(const std::string& str){
    std::string::size_type n = str.length();

    std::string::size_type j = 0;
    while(j < n && (str[j] == ' ' || str[j] == '\t' || str[j] == '\n'))
        j++;

    if(j == 0){
        return str;
    }else{
        std::string temp = str;
        temp.erase(0,j);
        return temp;
    }
}


/** Returns a double from a string
 * \param entry the string to interpret, which only needs to start with a double
 * \return the double
 * \exception throws a Subs::Subs_Error if the string cannot be read as a double 
 */
double Subs::string_to_double(const std::string& entry){
    double value;
    std::istringstream istr(entry);
    istr >> value;
    if(!istr)
        throw Subs_Error("Subs::string_to_double: could not translate entry = " + entry + " to a double");
    return value;
}

/** Returns an int from a string
 * \param entry the string to interpret, which only needs to start with an int
 * \return the int
 * \exception throws a Subs::Subs_Error if the string cannot be read as an int
 */
int Subs::string_to_int(const std::string& entry){
    int value;
    std::istringstream istr(entry);
    istr >> value;
    if(!istr)
        throw Subs_Error("Subs::string_to_int: could not translate entry = " + entry + " to an int");
    return value;
}

/** Returns a char from a string
 * \param entry the string to interpret, which only needs to start with a char.
 * \return the char
 * \exception throws a Subs::Subs_Error if the string cannot be read as a char
 */
char Subs::string_to_char(const std::string& entry){
    char value;
    std::istringstream istr(entry);
    istr >> value;
    if(!istr)
        throw Subs_Error("Subs::string_to_char: could not translate entry = " + entry + " to a char");
    return value;
}

/** Returns a bool from a string
 * \param entry the string to interpret. Any of '1', 'T', 'true', 'Y' or 'Yes' (case-insensitive)
 * will be translated to 'true'. Similarly, '0', 'F', 'False', 'N' or 'No' will become 'false', anything else will cause an exception.
 * \return the bool
 * \exception throws a Subs::Subs_Error if the string cannot be read as a bool
 */
bool Subs::string_to_bool(const std::string& entry){
    std::string test = Subs::toupper(entry);
    if(test == "T" || test == "TRUE" || test == "Y" || test == "YES" || test == "1"){
        return true;
    }else if(test == "F" || test == "FALSE" || test == "N" || test == "NO" || test == "0"){
        return false;
    }else{
        throw Subs_Error("Subs::string_to_bool: could not translate entry = " + entry + "to a bool");
    }
}





