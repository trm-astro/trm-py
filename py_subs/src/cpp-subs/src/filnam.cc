#include <string>
#include "trm/subs.h"

/** Returns a filename with the supplied extension added if
 * it is not already present. Accounts for any extra
 * space at end of the initial string. Note that 'name' is
 * passed by value because it is modified inside the routine
 * and would therefore require copying if passsed by reference.
 * \param name the filename. Can include the entension already, in which case it is
 * returned unchanged.
 * \param extens the file extension, such as ".dat" or whatever
 * \return Returns the modified name
 */

std::string Subs::filnam(std::string name, const std::string& extens){


  std::string::size_type len = name.length();
  std::string::size_type lsp = name.find_last_not_of(" \t");
  if(lsp < len-1) name.erase(lsp+1);

  std::string::size_type look = name.rfind(extens);
  if(look == std::string::npos || look != name.length() - extens.length()){
    return name + extens;
  }else{
    return name;
  }
}
