#include <fstream>
#include <vector>
#include <string>
#include "trm/subs.h"
#include "trm/header.h"
#include "trm/colly.h"

/* write_molly_head writes the headers and arc poly
 * (if set) of a molly spectrum so that molly can still read them.
 * It assumes that the output stream is just positioned at the start of 
 * the header and leaves it just as the data record starts. write_molly_head 
 * assumes that the header is in a form that would be read in by 
 * read_molly_head. Do not rejig the headers beforehand. 
 *
 * \param ost      the output file stream which must have been opened in binary mode.
 * \param head     the header to write out. Although not const, it is not modified.
 * \param original the original order in which header items appeared in molly file. 
 * this can be used to write out in the same order if desired
 * \param old if true then header items will be written in original order
 */


void Colly::write_molly_head(std::ofstream& ost, Subs::Header head, const std::vector<std::string>& original, bool old){
  
    // Start by dumping a number of bytes at start of record as
    // Fortran expects

    int nbytes = 44;
    if(!(ost.write((char *)&nbytes, sizeof(int)))) 
	throw Colly_Error("Colly::write_molly_head: failed to write 4 bytes at start of record 1");

    int fcode;
    head["Xtra.FCODE"]->get_value(fcode);
    if(!(ost.write((char *)&fcode,sizeof(int)))) 
	throw Colly_Error("Colly::write_molly_head: Error writing fcode");

    std::string sunits;
    head["Xtra.UNITS"]->get_value(sunits);
    char hname[16];
    sunits.copy(hname,16);  
    for(std::string::size_type j=sunits.length(); j<16; j++)
	hname[j] = ' ';
    if(!(ost.write(hname,16)))
	throw Colly_Error("Colly::write_molly_head: Error writing units");

    int  npix;
    head["Xtra.NPIX"]->get_value(npix);  
    if(!(ost.write((char *)&npix,sizeof(int))))
	throw Colly_Error("Colly::write_molly_head: Error writing npix");

    int  narc;
    head["Xtra.NARC"]->get_value(narc);  
    if(!(ost.write((char *)&narc,sizeof(int))))
	throw Colly_Error("Colly::write_molly_head: Error writing narc");

    double arc[std::max(1,abs(narc))];
    if(narc != 0){
	Subs::Hitem* hitp = head["Xtra.ARC"];
	for(int k=0; k<abs(narc); k++)
	    arc[k] = hitp->get_dvector()[k];
    }

    head.erase("Xtra");

#ifdef DEBUG
    std::cerr << "FCODE = " << fcode << ", UNITS = " << sunits
	      << ", NPIX = " << npix << ", NARC = " << narc << std::endl;
#endif

    // Count up items

    Subs::Header::iterator hit;
    Subs::Header::Hnode* hnode;
    int nchar = 0, ndoub = 0, nintr = 0, nfloat = 0; 

    if(old){
	for(size_t l=0; l<original.size(); l++){
	    hnode = head.find(original[l]);
	    if(hnode->value->type() == "string"){ 
		nchar++;
	    }else if(hnode->value->type() == "double"){
		ndoub++;
	    }else if(hnode->value->type() == "int"){
		nintr++;
	    }else if(hnode->value->type() == "float"){
		nfloat++;
	    }
	}
    }else{
	for(hit=head.begin(); hit != head.end(); hit++){
	    if(hit->value->type() == "string"){ 
		nchar++;
	    }else if(hit->value->type() == "double"){
		ndoub++;
	    }else if(hit->value->type() == "int"){
		nintr++;
	    }else if(hit->value->type() == "float"){
		nfloat++;
	    }
	}
    }

    // Write out numbers

    if(!(ost.write((char *)&nchar,sizeof(int))))
	throw Colly_Error("Colly::write_molly_head: Error writing nchar");
    if(!(ost.write((char *)&ndoub,sizeof(int))))
	throw Colly_Error("Colly::write_molly_head: Error writing ndoub");
    if(!(ost.write((char *)&nintr,sizeof(int))))
	throw Colly_Error("Colly::write_molly_head: Error writing nintr");
    if(!(ost.write((char *)&nfloat,sizeof(int))))
	throw Colly_Error("Colly::write_molly_head: Error writing nfloat");

    if(!(ost.write((char *)&nbytes, sizeof(int)))) 
	throw Colly_Error("Colly::write_molly_head: failed to write 4 bytes at end of record 1");


#ifdef DEBUG
    std::cerr << "NCHAR = " << nchar << ", NDOUB = " << ndoub
	      << ", NINTR = " << nintr << ", NFLOAT = " << nfloat << std::endl;
#endif

    nbytes = (nchar+ndoub+nintr+nfloat)*16;
    if(!(ost.write((char *)&nbytes, sizeof(int)))) 
	throw Colly_Error("Colly::write_molly_head: failed to write 4 bytes at start of record 2");
  
    // Write header item names. Pad out with blanks.

    if(old){
	for(size_t l=0; l<original.size(); l++){
	    hnode = head.find(original[l]);
	    if(hnode->value->type() == "string"){
		hnode->name.copy(hname,16);
		for(std::string::size_type j=hnode->name.length(); j<16; j++)
		    hname[j] = ' ';
		if(!(ost.write(hname, 16)))
		    throw Colly_Error("Colly::write_molly_head: failed to write header name = " + hnode->name);
	    }
	}
    
	for(size_t l=0; l<original.size(); l++){
	    hnode = head.find(original[l]);
	    if(hnode->value->type() == "double"){
		hnode->name.copy(hname,16);
		for(std::string::size_type j=hnode->name.length(); j<16; j++)
		    hname[j] = ' ';
		if(!(ost.write(hname, 16)))
		    throw Colly_Error("Colly::write_molly_head: failed to write header name = " + hnode->name);
	    }
	}
    
	for(size_t l=0; l<original.size(); l++){
	    hnode = head.find(original[l]);
	    if(hnode->value->type() == "int"){
		hnode->name.copy(hname,16);
		for(std::string::size_type j=hnode->name.length(); j<16; j++)
		    hname[j] = ' ';
		if(!(ost.write(hname, 16)))
		    throw Colly_Error("Colly::write_molly_head: failed to write header name = " + hnode->name);
	    }
	}
    
	for(size_t l=0; l<original.size(); l++){
	    hnode = head.find(original[l]);
	    if(hnode->value->type() == "float"){
		hnode->name.copy(hname,16);
		for(std::string::size_type j=hnode->name.length(); j<16; j++)
		    hname[j] = ' ';
		if(!(ost.write(hname, 16)))
		    throw Colly_Error("Colly::write_molly_head: failed to write header name = " + hnode->name);
	    }
	}

    }else{

	for(hit=head.begin(); hit != head.end(); hit++){
	    if(hit->value->type() == "string"){
		hit->name.copy(hname,16);
		for(std::string::size_type j=hit->name.length(); j<16; j++)
		    hname[j] = ' ';
		if(!(ost.write(hname, 16)))
		    throw Colly_Error("Colly::write_molly_head: failed to write header name = " + hit->name);
	    }
	}
    
	for(hit=head.begin(); hit != head.end(); hit++){
	    if(hit->value->type() == "double"){
		hit->name.copy(hname,16);
		for(std::string::size_type j=hit->name.length(); j<16; j++)
		    hname[j] = ' ';
		if(!(ost.write(hname, 16)))		    
		    throw Colly_Error("Colly::write_molly_head: failed to write header name = " + hit->name);
	    }
	}
    
	for(hit=head.begin(); hit != head.end(); hit++){
	    if(hit->value->type() == "int"){
		hit->name.copy(hname,16);
		for(std::string::size_type j=hit->name.length(); j<16; j++)
		    hname[j] = ' ';
		if(!(ost.write(hname, 16)))
		    throw Colly_Error("Colly::write_molly_head: failed to write header name = " + hit->name);
	    }
	}
    
	for(hit=head.begin(); hit != head.end(); hit++){
	    if(hit->value->type() == "float"){
		hit->name.copy(hname,16);
		for(std::string::size_type j=hit->name.length(); j<16; j++)
		    hname[j] = ' ';
		if(!(ost.write(hname, 16)))
		    throw Colly_Error("Colly::write_molly_head: failed to write header name = " + hit->name);
	    }
	}
    }

    if(!(ost.write((char *)&nbytes, sizeof(int)))) 
	throw Colly_Error("Colly::write_molly_head: failed to write 4 bytes at end of record 2");

#ifdef DEBUG
    std::cerr << "Written header item names" << std::endl;
#endif

    nbytes = nchar*32 + ndoub*sizeof(double) + nintr*sizeof(int) + nfloat*sizeof(float);
    if(!(ost.write((char *)&nbytes, sizeof(int)))) 
	throw Colly_Error("Colly::write_molly_head: failed to write 4 bytes at start of record 3");

    // Now write values

    char citem[32];
    double d;
    int i;
    float f;

    if(old){

	for(size_t l=0; l<original.size(); l++){
	    hnode = head.find(original[l]);
	    if(hnode->value->type() == "string"){
		for(std::string::size_type j=hnode->value->get_string().length(); j<32; j++)
		    citem[j] = ' ';
		hnode->value->get_string().copy(citem,32);
		if(!(ost.write(citem, 32)))
		    throw Colly_Error("Colly::write_molly_head: failed to write header value for item = " + hnode->name);
	    }
	}
    
	for(size_t l=0; l<original.size(); l++){
	    hnode = head.find(original[l]);
	    if(hnode->value->type() == "double"){
		d = hnode->value->get_double();
		if(!(ost.write((char*)&d, sizeof(double))))
		    throw Colly_Error("Colly::write_molly_head: failed to write header value for item = " + hnode->name);
	    }
	}
    
	for(size_t l=0; l<original.size(); l++){
	    hnode = head.find(original[l]);
	    if(hnode->value->type() == "int"){
		i = hnode->value->get_int();
		if(!(ost.write((char*)&i, sizeof(int))))
		    throw Colly_Error("Colly::write_molly_head: failed to write header value for item = " + hnode->name);
	    }
	}
    
	for(size_t l=0; l<original.size(); l++){
	    hnode = head.find(original[l]);
	    if(hnode->value->type() == "float"){
		f = hnode->value->get_float();
		if(!(ost.write((char*)&f, sizeof(float))))
		    throw Colly_Error("Colly::write_molly_head: failed to write header value for item = " + hnode->name);
	    }
	}

    }else{

	for(hit=head.begin(); hit != head.end(); hit++){
	    if(hit->value->type() == "string"){
		for(std::string::size_type j=hit->value->get_string().length(); j<32; j++)
		    citem[j] = ' ';
		hit->value->get_string().copy(citem,32);
		if(!(ost.write(citem, 32)))
		    throw Colly_Error("Colly::write_molly_head: failed to write header value for item = " + hit->name);
	    }
	}
    
	for(hit=head.begin(); hit != head.end(); hit++){
	    if(hit->value->type() == "double"){
		d = hit->value->get_double();
		if(!(ost.write((char*)&d, sizeof(double))))
		    throw Colly_Error("Colly::write_molly_head: failed to write header value for item = " + hit->name);
	    }
	}
    
	for(hit=head.begin(); hit != head.end(); hit++){
	    if(hit->value->type() == "int"){
		i = hit->value->get_int();
		if(!(ost.write((char*)&i, sizeof(int))))
		    throw Colly_Error("Colly::write_molly_head: failed to write header value for item = " + hit->name);
	    }
	}
    
	for(hit=head.begin(); hit != head.end(); hit++){
	    if(hit->value->type() == "float"){
		f = hit->value->get_float();
		if(!(ost.write((char*)&f, sizeof(float))))
		    throw Colly_Error("Colly::write_molly_head: failed to write header value for item = " + hit->name);
	    }
	}
    }

    if(!(ost.write((char *)&nbytes, sizeof(int)))) 
	throw Colly_Error("Colly::write_molly_head: failed to write"
			  " 4 bytes at end of record 3");

#ifdef DEBUG
    std::cerr << "Written header item values" << std::endl;
#endif

    nbytes = abs(narc)*sizeof(double);
    if(!(ost.write((char *)&nbytes, sizeof(int)))) 
	throw Colly_Error("Colly::write_molly_head: failed to write"
			  " 4 bytes at start of record 4");
    if(narc != 0){
	if(!(ost.write((char *)&arc,abs(narc)*sizeof(double))))
	    throw Colly_Error("Colly::write_molly_head: Error writing arc coefficients");

#ifdef DEBUG
	std::cerr << "Written arc coefficients" << std::endl;
#endif

    }

    if(!(ost.write((char *)&nbytes, sizeof(int)))) 
	throw Colly_Error("Colly::write_molly_head: failed to write 4 bytes at end of record 4");
}
