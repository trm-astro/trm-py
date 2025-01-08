#include <stdlib.h>
#include "trm/subs.h"
#include "trm/position.h"
#include "trm/colly.h"

/**
 * find_groups recursively allocates groups of spectra
 * by position. The first time its is called, all .group elements 
 * of the vector of spectrum information should be = 0, and the 
 * position should be set equal to that of the first element to ensure a match.
 * It will then call this the first group and start looking for
 * connected members. Each time it finds one of these it calls
 * itself again with the latest member as the start position.
 * By doing so it effectively searches down to the end of every
 * branch to cover the whole group. When it comes back to the
 * top level the group counter should be advanced one. This
 * routine should be called by a loop which goes over every
 * single element of info.
 *
 * The struct Info is defined in molly.h as
 * <pre>
 * struct Info {
 *   Info() : group(0) {}
 *   Subs::Position pos;
 *   Subs::Date date;
 *   string name;
 *   size_t group;
 * };
 *
 * </pre>
 *
 * \param info   vector of basic information as defined above
 * \param pos    initialise to the position of the first spectrum
 * \param ngroup the number of the group
 * \param thresh the cosine threshold that must be exceeded for inclusion in a group
 */

void Colly::find_groups(std::vector<Info>& info, const Subs::Position& pos, 
			const size_t& ngroup, const double& thresh){
  for(size_t i=0; i<info.size(); i++){
    if(!info[i].group && dot(info[i].pos,pos) > thresh){
      info[i].group = ngroup;
      find_groups(info,info[i].pos,ngroup,thresh);
    }
  }
}
