#include <fstream>
#include "trm/subs.h"
#include "trm/header.h"
#include "trm/colly.h"

/**
 * skip_molly_data is a routine to skip over molly
 * data assuming that the stream is initially positioned at
 * the start of the data after the header has been read into
 * a Subs::Header object. Thus skip_molly_data can be used 
 * immediately after read_molly_head to jump to the next header.
 *
 * \param ist  input file stream that should have been opened
 * in binary mode and be positioned at the start of the molly data.
 * \param head header read in by read_molly_head
 */

void Colly::skip_molly_data(std::ifstream& ist, Subs::Header& head){
  int  fcode, npix;
  int  nbytes, nstart, nend;

  if(!(ist.read((char *)&nstart,sizeof(int))))
    throw Colly_Error("Colly::skip_molly_data: Error reading number of bytes at start of data");
  fcode = head["Xtra.FCODE"]->get_int();
  npix  = head["Xtra.NPIX"]->get_int();
  switch(fcode){
  case 1: case 4:
    nbytes = npix*sizeof(float);
    break;
  case 2: case 5:
    nbytes = 2*npix*sizeof(float);
    break;
  case 3:
    nbytes = 3*npix*sizeof(float);
    break;
  default:
    std::string message = "Error in Colly::skip_molly_data. Unrecognised molly format code = " + std::to_string(fcode);
    throw Colly_Error(message);
  }

  if(nbytes != nstart){
    std::string message = "Colly::skip_molly_data: estimated number of bytes = ";
    message += nbytes;
    message += " fails to match number expected from start of record = ";
    message += nstart;
    throw Colly_Error(message);
  }
  
  if(!skip_bytes(ist,nbytes))
    throw Colly_Error("Colly::skip_molly_data: Error skipping bytes");

  if(!(ist.read((char *)&nend,sizeof(int))))
    throw Colly_Error("Colly::skip_molly_data: Error reading number of bytes at end of data");

  if(nbytes != nend){
    std::string message = "Colly::skip_molly_data: estimated number of bytes = ";
    message += nbytes;
    message += " fails to match number expected from end of record = ";
    message += nend;
    throw Colly_Error(message);
  }
}





