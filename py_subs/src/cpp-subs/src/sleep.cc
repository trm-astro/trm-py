#include <ctime>
#include "trm/subs.h"

/** Function to put a program to sleep for a bit
 * \param seconds the time to pause in seconds.
 */

void Subs::sleep(double seconds){
  timespec delay, ret;
  delay.tv_sec  = time_t(int(seconds));
  delay.tv_nsec = long(1.e9*(seconds-int(seconds)));
  nanosleep(&delay,&ret);
}
