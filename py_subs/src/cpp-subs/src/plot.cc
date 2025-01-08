#include "trm/plot.h"

/** Constructor to open a plot device
 * \param device name of plot device to open
 */
Subs::Plot::Plot(const std::string& device) {
	pls = new plstream();
    pls->sdev(device.c_str()); // Set plot device, e.g., "xwin", "png"
	//delay this until colours are set
	//pls->init();               // Initialize the PLplot stream
    devname = device;
}

Subs::Plot::~Plot(){
	if (pls) {
        //pls->end(); // Finalize PLplot stream
        delete pls;
    }
}

/** Closes a plot
 */
void Subs::Plot::close(){
	if (pls) {
        //pls->end(); // Close PLplot
        delete pls;
        pls = nullptr;
    }
}

/** Opens a plot, closing it first if it is already
 * open. The focus of next plotting commands switches to
 * this one
 * \param device plot device to be opened
 */
void Subs::Plot::open(const std::string& device) {
	close();
    pls = new plstream();
    pls->sdev(device.c_str());
    pls->init();
    devname = device;
}

/** Moves the focus of the next plotting commands to the
 * plot 
 */ 
void Subs::Plot::focus() const {
	if (pls) {
        pls->adv(0); // PLplot doesn’t require explicit focus but advance can update
    }
}

/** Tests whether a plot is in focus or not
 */
bool Subs::Plot::is_in_focus() const {
	return pls != nullptr;
}

/** Sets the ranges of the plot panel
 * \param x_1 left-hand X limit
 * \param x_2 right-hand X limit
 * \param y_1 lower Y limit
 * \param y_2 upper Y limit
 * \param scale true for equal physical scales
 */    
void Subs::Plot::Panel::set(float x_1, float x_2, float y_1, float y_2, bool scale) {
	x1   = x_1;
	x2   = x_2;
	y1   = y_1;
	y2   = y_2;
	just = scale;
	stan = true;
}

/** Sets the ranges of the plot panel
 * \param x_1 left-hand X limit
 * \param x_2 right-hand X limit
 * \param y_1 lower Y limit
 * \param y_2 upper Y limit
 * \param scale true for equal physical scales
 * \param xv_1 left-hand X limit of viewport
 * \param xv_2 right-hand X limit of viewport
 * \param yv_1 lower Y limit of viewport
 * \param yv_2 upper Y limit of viewport
 */    
void Subs::Plot::Panel::set(float x_1, float x_2, float y_1, float y_2, bool scale, float xv_1, 
							float xv_2, float yv_1, float yv_2){
	x1   = x_1;
	x2   = x_2;
	y1   = y_1;
	y2   = y_2;
	just = scale;
	stan = false;
	xv1  = xv_1;
	xv2  = xv_2;
	yv1  = yv_1;
	yv2  = yv_2;
}

void Subs::Plot::Panel::focus() const {
	// do we need to set line colour to black here?
	if (plot->pls) {
        if (stan) {
            plot->pls->env(x1, x2, y1, y2, 0, just ? 1 : 0); // 1 for equal scaling
        } else {
            plot->pls->vpor(xv1, xv2, yv1, yv2);
            plot->pls->wind(x1, x2, y1, y2);
        }
    }
}

/** Set the colourmap from input rgb arrays */
void Subs::Plot::set_colors(PLINT* r, PLINT* g, PLINT* b, int n){

	// delete the old colour map (if it exists)
	if(ncol > 0){
		delete[] rgb[0];
	}

	// Update the internal colour map
	PLINT _rgb[n*3]; // Use pointers to access the rgb components
	rgb[0] = _rgb;
	rgb[1] = _rgb + n;
	rgb[2] = _rgb + n*2;
	ncol = n;
	for(int i=0; i<n; i++){
		rgb[0][i] = r[i];
		rgb[1][i] = g[i];
		rgb[2][i] = b[i];
	}

	std::cout << "ncol = " << ncol << std::endl;
	for (int i=0; i<ncol; i++){
		std::cout << "rgb[" << i << "] = " << rgb[0][i] << " " << rgb[1][i] << " " << rgb[2][i] << std::endl;
	}
	
	pls->scmap0n(n);
	pls->scmap0(rgb[0], rgb[1], rgb[2], ncol);
	pls->init();

}

void Subs::Plot::add_colors(PLINT* r, PLINT* g, PLINT* b, int n){
	// Get the current colour map and store
	PLINT new_rgb[ncol*3]; // new array for colours
	// get the old array pointers
	PLINT* _rgb[3]; 
	_rgb[0] = rgb[0];
	_rgb[1] = rgb[1];
	_rgb[2] = rgb[2];
	int nold = ncol; // old number of colours
	// Update the internal colour map number of colours
	ncol = nold + n;
	// Update the internal colour map pointers
	rgb[0] = new_rgb;
	rgb[1] = new_rgb + ncol;
	rgb[2] = new_rgb + ncol*2;

	
	for (int i=0; i<nold; i++){
		new_rgb[i] = _rgb[0][i];
		new_rgb[ncol + i] = _rgb[1][i];
		new_rgb[ncol*2 + i] = _rgb[2][i];
	}
	for (int i=0; i<n; i++){
		new_rgb[nold + i] = r[i];
		new_rgb[ncol + nold + i] = g[i];
		new_rgb[ncol*2 + nold + i] = b[i];
	}
	// Update the colour map
	pls->scmap0n(ncol);
	pls->scmap0(rgb[0], rgb[1], rgb[2], ncol);

	// Free the old array
	delete[] _rgb[0];
}

/** Plots a graph with the UT axis running 21, 22, 23, 0, 1 etc if
 * it has to cross 24. Plot must be open, viewport, colours set etc beforehand.
 * \param t1 start time
 * \param t2 end time
 * \param y1 lower Y limit
 * \param y2 upper Y limit
 * \param sides true to plot axes along the sides
 * \param top true to plot axis along the top
 */
void Subs::Plot::ut_plot(float t1, float t2, float y1, float y2, bool sides, bool top){
	
    float delta = 2.0; // Replace 2.0 as necessary to match the desired tick interval
    if (t2 > 24.0) {
        double xv1, xv2, yv1, yv2;
        pls->gspa(xv1, xv2, yv1, yv2);
        float xvt1 = xv1 + (xv2 - xv1) * (24.0 - t1) / (t2 - t1);

        // Left-hand box
        pls->vpor(xv1, xvt1, yv1, yv2);
        pls->wind(t1, 23.9999, y1, y2);
        pls->box("bcnst", delta, 0, "bnst", delta, 0);

        // Right-hand box
        pls->vpor(xvt1, xv2, yv1, yv2);
        pls->wind(0.0, t2 - 24.0, y1, y2);
        pls->box("bcnst", delta, 0, "cst", delta, 0);

        // Reset viewport and window
        pls->vpor(xv1, xv2, yv1, yv2);
        pls->wind(t1, t2, y1, y2);
    } else {
        pls->env(t1, t2, y1, y2, 0, 0);
        pls->box("bcnst", delta, 0, "bnst", delta, 0);
    }
}

void Subs::ut_plot(float t1, float t2, float y1, float y2, bool sides, bool top){
	std::cout << "Subs::ut_plot deprecated. Please call from Subs::Plot::ut_plot i.e. plot->ut_plot" << std::endl;
}

Subs::Plots::Plots(int nplot) : nplot_(nplot), nfocus_(0) {
	if(nplot < 1)
		throw Subs::Plot::Plot_Error("Subs::Plots(): need at least 1 plot");
	plot_ = new Plot[nplot];
}

/** Opens the plot now in focus, if it is not open already. If
 * the focus is set < 0, then it opens the first one that is not already
 * open.
 * \param device the plot device to open.
 */
void Subs::Plots::open(const std::string& device){
	if(nfocus_ >= 0){
		if(plot_[nfocus_].is_open())
			plot_[nfocus_].focus();
		else
			plot_[nfocus_].open(device);

	}else{
		bool ok = false;
		for(int i=0; i<nplot_; i++){
			if(!plot_[i].is_open()){
				plot_[i].open(device);
				nfocus_ = i;
				ok = true;
				break;
			}
		}
		if(!ok)
			throw Plot::Plot_Error("Subs::Plots::open(const std::string&): no more plots available. You should close some");
	}
}

/** Explicitly close a particular plot.
 * \param nclose the plot to close. -1 to close the current plot.
 */
void Subs::Plots::close(int nclose){
	if(nclose >= 0 && nclose < nplot_){
		plot_[nclose].close();
		nfocus_ = -1;
	}else if(nclose == -1){
		if(nfocus_ >= 0)
			plot_[nfocus_].close();
	}else{
		throw Plot::Plot_Error("Subs::Plots::close(int): plot number out of range -1 to " + Subs::str(nplot_-1));
	}
}

void Subs::Plots::close(){
	for(int i=0; i<nplot_; i++)
		plot_[i].close();
	nfocus_ = -1;
}

bool Subs::Plots::is_open() const {
	if(nfocus_ >= 0 && nfocus_ < nplot_){
		return plot_[nfocus_].is_open();
	}else{
		return false;
	}
}

void Subs::Plots::focus(){
	if(nfocus_ >= 0)
		plot_[nfocus_].focus();
}

void Subs::Plots::focus(int i){
	if(i >= 0 && i < nplot_){
		if(plot_[i].is_open())
			plot_[i].focus();
		nfocus_ = i;
	}else if(i < 0){
		nfocus_ = -1;
	}else{
		throw Plot::Plot_Error("void Subs::Plots::focus(int): attempt to focus on non existent plot number " + Subs::str(i+1));
	}
}
   
/** Lists the plots and those which are open
 * \param ostr the output stream to send the list to
 */
void Subs::Plots::list(std::ostream& ostr){
	for(int i=0; i<nplot_; i++){
		if(plot_[i].is_open()){
			ostr << "Plot " << i+1 << " was opened with device name = " << plot_[i].get_name() << std::endl;
		}else{
			ostr << "Plot " << i+1 << " is not open."  << std::endl;
		}
	}
} 

Subs::Plots::~Plots(){
	for(int i=0; i<nplot_; i++)
		plot_[i].close();
	delete[] plot_;
}


