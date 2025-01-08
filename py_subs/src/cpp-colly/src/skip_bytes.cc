#include <fstream>
#include "trm/colly.h"

/**
 * skip_bytes skips nbytes of an input stream, ignoring their
 * values. It should be possible to do this with 'ignore' 
 * but it seems to fail on Solaris which seems to regard
 * the byte 'ff' as an end-of-file
 *
 * \param ist input stream
 * \param nbytes number of bytes to be skipped
 */

bool Colly::skip_bytes(std::ifstream& ist, size_t nbytes){
  if(nbytes){
    char c;
    size_t count = 0;
    while(count < nbytes && ist.get(c)) count++;
    if(!ist)
      return false;
    else
      return true;
  }else{
    return true;
  }
}



