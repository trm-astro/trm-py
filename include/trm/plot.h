#ifndef TRM_PLOT_H
#define TRM_PLOT_H

#include <string>
#include "plstream.h"
#include "trm/subs.h"

namespace Subs {

  //! A class to represent plots

  /**
   * This class opens a plot device with PGPLOT and remembers
   * the device identifier. It can then be used to reselect the
   * associated plot device and close it later on without the user
   * having to worry when multiple devices are being used. It also
   * closes the device automatically on going out of scope which avoids
   * problems that occur when errors occur.
   */

  class Plot {

  public:

    //! Default constructor. Does not open anything
    Plot() : idev(0) , devname() {}

    //! Standard constructor to open device
    Plot(const std::string& device);

    //! Destructor closes the associated device
    ~Plot();

    //! explicitly close a plot
    void close();

    //! open a plot (first closing it if need be)
    void open(const std::string& device);

    //! Get the pls pointer
    plstream* get_plstream() {return pls;}
    
    //! move focus to a plot
    void focus() const;

    //! Tests whether a plot is open or not.
    bool is_open() const {return (idev > 0);}

    //! Returns device name used for this plot
    std::string get_name() const {return devname;}

    //! Tests whether a plot is in focus or not
    bool is_in_focus() const;

    //! Sets the colors for a plot
    void set_colors(PLINT* r, PLINT* g, PLINT* b, int n);
    void add_colors(PLINT* r, PLINT* g, PLINT* b, int n);

    //! Draws the axes for a plot versus UT 
    void ut_plot(float t1, float t2, float y1, float y2, bool sides = true, bool top = true);

    //! Represents a plot panel. Inside namespace plot
    struct Panel{
      
      //! Default constructor
      Panel(Plot* plotPtr) : 
        plot(plotPtr),
	      x1(0.), x2(1.), y1(0.), y2(1.), just(false), stan(true) {}
      
      //! Constructor specifying the ranges
      Panel(Plot* plotPtr, float x_1, float x_2, float y_1, float y_2, bool scale) :
        plot(plotPtr),
	      x1(x_1), x2(x_2), y1(y_1), y2(y_2), just(scale), stan(true) {}
      
      //! Constructor specifying the ranges and the viewport
      Panel(Plot* plotPtr, float x_1, float x_2, float y_1, float y_2, bool scale, float xv_1, float xv_2, 
	    float yv_1, float yv_2) : 
        plot(plotPtr),
	      x1(x_1), x2(x_2), y1(y_1), y2(y_2), just(scale), stan(false), xv1(xv_1), xv2(xv_2), 
	      yv1(yv_1), yv2(yv_2) {}
    
      //! Sets the ranges
      void set(float x_1, float x_2, float y_1, float y_2, bool scale);
      
      //! Sets the ranges and the viewport
      void set(float x_1, float x_2, float y_1, float y_2, bool scale, float xv_1, 
	       float xv_2, float yv_1, float yv_2);
      
      //! Restore the parameters of this panel.
      void focus() const;
      
      float x1, x2, y1, y2;
      bool just, stan;
      float xv1, xv2, yv1,yv2;
      
      private:
        Plot* plot;
    };

    //! Error class for Plot
    class Plot_Error : public Subs_Error {
    public:
      //! Default constructor
      Plot_Error() : Subs_Error() {}
      
      //! Constructor from a string (e.g. an error message).
      Plot_Error(const std::string& err) : Subs_Error(err) {} 
    };
      
  private:
    plstream* pls;
    int idev;
    std::string devname;
    Plot(const Plot& plot){}; // disable copying because it causes problems when objects go out of scope
    //    Plot& operator=(const Plot& obj){return Plot();} // no assignment
    Plot& operator=(const Plot& obj);
    // Pointer to the start of the RGB colour arrays
    PLINT* rgb[3];
    // Number of colours (array above is 3 x ncol)
    int ncol = 0;
  };    
  

  //! Draws the axes for a plot versus UT 
  void ut_plot(float t1, float t2, float y1, float y2, bool sides = true, bool top = true);
  
  //! Class to handle multiple plots, keeping track of which is in focus
  class Plots {
  public:

    //! Constructor of multiple plots
    Plots(int nplot);

    //! Opens first available plot
    void open(const std::string& device);

    //! Close all plots
    void close();

    //! Close a specific plot
    void close(int nclose);

    //! Focus on whatever plot should be in focus
    void focus();

    //! Focus on a particular plot
    void focus(int i);

    //! Is currently focussed plot open?
    bool is_open() const;

    //! Return the plot number currently in focus (<-1 means none is)
    int get_focus() const {return nfocus_;}

    //! Return the total number of plots
    int get_nplot() const {return nplot_;}
     
    //! Lists the plots and those which are open
    void list(std::ostream& ostr);

    //! Destructor
    ~Plots();
    
  private:
    Plots(const Plots& plot){}; // no copy constructor
    //    Plots& operator=(const Plots& obj){return obj;} // no assignment
    Plots& operator=(const Plots& obj); // no assignment
    int nplot_;
    int nfocus_;
    Plot* plot_;
  };

};

#endif
