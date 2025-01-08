#include "trm/subs.h"

/** Determines whether current machine is little endian or not.
 * Needed in order to read various binary files on either type of
 * machine.
 * \return true if the machine is little endian.
 */
bool Subs::is_little_endian() {

  long int one = 1;
  return (*((char *)(&one)) == 1);

}

/** Determines whether current machine is big endian or not.
 * Needed in order to read various binary files on either type of
 * machine.
 * \return true if the machine is big endian.
 */
bool Subs::is_big_endian() {

  long int one = 1;
  return !(*((char *)(&one)) == 1);

}



