#include <string>
#include "trm/subs.h"

/** Translates a string such as "red" or "black" into an equivalent 
 * enum that can be used to set the colours in the PGPLOT routine
 * cpgsci. Recognised colours are: black, white, red, green, blue,
 * cyan, purple, yellow, orange, light_green dark_green, light_blue,
 * dark_blue, light_grey, pink & none for no plot
 */

Subs::PLOT_COLOUR Subs::what_colour(const std::string& colour){
  if(toupper(colour) == "NONE"){
    return NONE;
  }else if(toupper(colour) == "BLACK"){
    return BLACK;
  }else if(toupper(colour) == "WHITE"){
    return WHITE;
  }else if(toupper(colour) == "RED"){
    return RED;
  }else if(toupper(colour) == "GREEN"){
    return GREEN;
  }else if(toupper(colour) == "BLUE"){
    return BLUE;
  }else if(toupper(colour) == "CYAN"){
    return CYAN;
  }else if(toupper(colour) == "PURPLE"){
    return PURPLE;
  }else if(toupper(colour) == "YELLOW"){
    return YELLOW;
  }else if(toupper(colour) == "ORANGE"){
    return ORANGE;
  }else if(toupper(colour) == "LIGHT_GREEN"){
    return LIGHT_GREEN;
  }else if(toupper(colour) == "DARK_GREEN"){
    return DARK_GREEN;
  }else if(toupper(colour) == "LIGHT_BLUE"){
    return LIGHT_BLUE;
  }else if(toupper(colour) == "DARK_BLUE"){
    return DARK_BLUE;
  }else if(toupper(colour) == "DARK_GREY"){
    return DARK_GREY;
  }else if(toupper(colour) == "LIGHT_GREY"){
    return LIGHT_GREY;
  }else if(toupper(colour) == "PINK"){
    return PINK;
  }else{
    throw Subs_Error("Subs::what_colour(const std::string&): colour = " + colour + " not recognised.");
  }
}

