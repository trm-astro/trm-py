#include <iomanip>
#include "trm/header.h"
#include "trm/format.h"

// static members
std::string Subs::Header::start_string = "";

Subs::UINT4 Subs::Header::name_width = 20;

//! Destructor
Subs::Header::~Header() {
    Subs::deleter(head);
}

/** Returns number of items in a Header as would be written to disk
 */
Subs::UINT4 Subs::Header::count() const {
    return Subs::count(head);
}

Subs::UINT4 Subs::count(const Header::Hnode *node){
    Subs::UINT4 n = 0;
    while(node->left != NULL){
	if(node->right != NULL)
	    n += count(node->right) + 1;
	else
	    n++;
	node = node->left;
    }
    return n;
}


void Subs::Header::write(std::ofstream& s) const {

    // Number of header items
    Subs::UINT4 lmap = this->count();
    s.write((char*)&lmap,sizeof(Subs::UINT4));
    
    if(s && lmap) Subs::write(s, head, std::string(""));

}

// Recursive write to a binary file
void Subs::write(std::ofstream& s, const Header::Hnode *node, const std::string& dir) {
    while(node->left != NULL){
	if(dir == "")
	    Subs::write_string(s, node->name);
	else
	    Subs::write_string(s, dir + Header::dir_flag + node->name);
	Subs::Hitem::write_item(s, node->value);
	if(node->right != NULL){
	    if(dir == "")
		write(s, node->right, node->name);
	    else
		write(s, node->right, dir + Header::dir_flag + node->name);
	}
	node = node->left;
    }
}

/**
 * Returns full name of an item with all directories appended.
 */

std::string Subs::Header::Hnode::fullname() const {
    const Subs::Header::Hnode *node = this;
    std::string name = node->name;
    while(node->up != NULL){
	if(node->up->right != NULL && node->up->right == node)
	    name = node->up->name + Subs::Header::dir_flag + name;
	node = node->up;
    }
    return name;
}

void Subs::Header::read(std::ifstream& s, bool swap_bytes) {

    // Read number of header items
    Subs::UINT4 lmap;
    s.read((char*)&lmap,sizeof(Subs::UINT4));
    if(swap_bytes) lmap = Subs::byte_swap(lmap);
    
    if(s && lmap){
	std::string name;
	
	for(Subs::UINT4 i=0; i<lmap; i++){
	    
	    // read the item name
	    Subs::read_string(s,name,swap_bytes);
	    
	    // then its value
	    this->set(name, Subs::Hitem::read_item(s, swap_bytes));
	    
	}
    }
}

// skip header in binary format
void Subs::Header::skip(std::ifstream& s, bool swap_bytes) {
  Subs::UINT4 lmap;
  s.read((char*)&lmap, sizeof(Subs::UINT4));
  if(swap_bytes) lmap = Subs::byte_swap(lmap);

  if(s && lmap){

    for(Subs::UINT4 i=0; i<lmap; i++){

      // skip the item name
      Subs::skip_string(s, swap_bytes);

      // then skip the item
      Subs::Hitem::skip_item(s, swap_bytes);

    }
  }
}

void Subs::Header::set(const std::string& name, Subs::Hitem* hip){
    if(valid_name(name)){
	Subs::set(head, name, hip);
	if(tail->left != NULL) tail = tail->left;
    }else{
	throw Header_Error("Invalid item name in void Subs::Header::set(const std::string&, Subs::Hitem* hip) = [" + name + "]");
    }
}

void Subs::set(Header::Hnode *node, const std::string& name, Hitem *hip) {

    std::string::size_type n = name.find_first_of(Header::dir_flag);
    if(n == std::string::npos){

	while(node->left != NULL){
	    if(node->name == name){
		delete node->value;
		node->value = hip;
		return;
	    }
	    node = node->left;
	}

	// Add the new item to the end.
	node->name     = name;
	node->value    = hip;
	node->left     = new Header::Hnode();
	node->left->up = node;
	if(hip->is_a_dir()){
	    node->right     = new Header::Hnode();
	    node->right->up = node;
	}
	
    }else{
    
	std::string directory = name.substr(0,n);
	bool found = false;
	while(node->left != NULL){
	    if(node->name == directory){
		if(node->right != NULL){
		    found = true;
		    set(node->right, name.substr(n+1), hip);
		}else{
		    throw Subs_Error("Subs::set(const Header::Hnode*, const std::string&, const Hitem*): directory part of name = " + name +
				     " matches an item called " + node->name + " but the latter is not a directory");
		}
	    }
	    node = node->left;
	}

	if(!found)
	    throw Subs_Error("Subs::set(const Header::Hnode*, const std::string&, const Hitem*): did not find directory = " + 
			     directory + " part of name = " + name);
    }
}

void Subs::Header::set_force(const std::string& name, Subs::Hitem* hip){
    if(valid_name(name)){
	Subs::set_force(head, name, hip);
	if(tail->left != NULL) tail = tail->left;
    }else{
	throw Header_Error("Invalid item name in void Subs::Header::set(const std::string&, Subs::Hitem* hip) = [" + name + "]");
    }
}

/** Sets a header item and any directories needed to store it
 */
void Subs::set_force(Header::Hnode *node, const std::string& name, Hitem *hip) {

    std::string::size_type n = name.find_first_of(Header::dir_flag);
    if(n == std::string::npos){

	while(node->left != NULL){
	    if(node->name == name){
		delete node->value;
		node->value = hip;
		return;
	    }
	    node = node->left;
	}

	// Add the new item to the end.
	node->name     = name;
	node->value    = hip;
	node->left     = new Header::Hnode();
	node->left->up = node;
	if(hip->is_a_dir()){
	    node->right     = new Header::Hnode();
	    node->right->up = node;
	}
	
    }else{
    
	std::string directory = name.substr(0,n);
	while(node->left != NULL){
	    if(node->name == directory) break;
	    node = node->left;
	}

	if(node->left == NULL){
	    node->name      = directory;
	    node->value     = new Hdirectory("automatically generated directory");
	    node->left      = new Header::Hnode();
	    node->right     = new Header::Hnode();
	    node->left->up  = node;
	    node->right->up = node;
	}else if(node->right == NULL){
	    throw Subs_Error("Subs::set(const Header::Hnode*, const std::string&, const Hitem*): directory part of name = " + name +
			     " matches an item called " + node->name + " but the latter is not a directory");
	}

	set_force(node->right, name.substr(n+1), hip);

    }
}

/**
 * Checks a header name for validity. This means it must contain no
 * white space, start or end with the directory character flag (e.g. '.')
 * or have two or more of them in a row.
 * \param name the header name to test
 */

bool Subs::Header::valid_name(const std::string &name){ 
    std::string test(2,Header::dir_flag);
    return !(name.find(' ')  != std::string::npos || name.find('\t') != std::string::npos ||
	     name[0] == Header::dir_flag ||          name[name.length()-1] == Header::dir_flag ||
	     name.find(test) != std::string::npos);
}

/* Finds an item matching the string name
 * \param name name to search for. 
*/

Subs::Header::Hnode* Subs::Header::find(const std::string& name) const {

    return Subs::find(name, head);

}


/* Finds a directory called name. 
 * \param dir  directory 
 */
Subs::Header::Hnode* Subs::Header::find_dir(const std::string& dir) const {
    Hnode *node = this->find(dir);
    if(node->right == NULL)
	node = tail;
    return node;
}    


Subs::Header::Hnode* Subs::find(const std::string& name, Header::Hnode* node) {
    
    std::string::size_type n = name.find_first_of(Header::dir_flag);
    if(n == std::string::npos){
	while(node->left != NULL){
	    if(node->name == name) break;
	    node = node->left;
	}
    }else{
	std::string directory = name.substr(0,n);
	while(node->left != NULL){
	    if(node->name == directory){
		if(node->right != NULL)
		    return find(name.substr(n+1), node->right);
		else
		    throw Subs_Error("Header::Hnode* Subs::find(const std::string&, Header::Hnode*): " + directory + 
				     " should have been a directory but was not.");
	    }
	    node = node->left;
	}
    }

    return node;

}

/** This routine searches for a header item with a name that matches the supplied regular expression
 * \param the regular expression (Perl compatible)
 * \param the pointer to the header item if a match is found (the first found)
 * \return 0 if no match is found, 1 is one is found 2 if 2 or more matches exist
 */    
Subs::Header::Hnode* Subs::Header::find(const std::string& regexp, int& nmatch) const {

	int errornumber;
	PCRE2_SIZE erroroffset;
	PCRE2_SPTR pattern = (PCRE2_SPTR)regexp.c_str();

	pcre2_code* re = pcre2_compile(
		pattern,
		PCRE2_ZERO_TERMINATED,
		0,
		&errornumber,
		&erroroffset,
		NULL
	);

	if (re == NULL) {
        PCRE2_UCHAR buffer[256];
        pcre2_get_error_message(errornumber, buffer, sizeof(buffer));
        std::cerr << "PCRE2 compilation failed at offset " << erroroffset
                  << ": " << buffer << std::endl;
        return NULL;
    }
  
    nmatch = 0;
    Hnode *fnode = NULL;
    Subs::find(re, std::string(""), head, nmatch, &fnode);

	pcre2_code_free(re);

    return fnode;
}


/** Recursive regular expression find
 */
void Subs::find(pcre2_code* re, const std::string& dir, Header::Hnode *node, int& nmatch, Header::Hnode **fnode) {
    
    while(node->left != NULL){
        // Prepare the subject string for matching
	    std::string subject = dir + Header::dir_flag + node->name;
        PCRE2_SPTR subject_str = (PCRE2_SPTR)subject.c_str();
		// Create match data block
        pcre2_match_data* match_data = pcre2_match_data_create_from_pattern(re, NULL);
		
		// Perform the match
        int rc = pcre2_match(
            re,
            subject_str,
            subject.size(),
            0,
            0,
            match_data,
            NULL
        );

        if (rc >= 0) {
            nmatch++;
            if (nmatch == 1) *fnode = node;
            if (nmatch > 1) {
                pcre2_match_data_free(match_data);
                break;
            }
        }
		
		// Free match data block after each match attempt
        pcre2_match_data_free(match_data);

		// Recur to the right subtree if it exists
        if (node->right) {
            find(re, dir + Header::dir_flag + node->name, node->right, nmatch, fnode);
            if (nmatch > 1) break;
        }
	
		node = node->left;
    }
}

/** Output operator for Header class.
 * \param s output stream
 * \param obj the Header to list
 */

std::ostream& Subs::operator<<(std::ostream& s, const Subs::Header& obj) {
    obj.print(s);
    return s;
}
    
void Subs::Header::print(std::ostream& s) const {
    Subs::print(s, head, 0, start_string);
}

void Subs::print(std::ostream& s, const Header::Hnode* node, int nlev, const std::string& start) {
    const std::string indent(Header::dir_indent,' ');
    nlev++;

    Format lj;
    while(node->left != NULL){

	if(node->right != NULL){

	    s << start << "\n" << start;
	    for(int i=0; i<nlev-1; i++) s << indent;
	    s << lj(node->name, Header::name_width + Hitem::value_width + 3) << lj(" /directory/",Hitem::type_width) << std::endl;
	    print(s, node->right, nlev, start);

	}else{

	    s << start;
	    for(int i=0; i<nlev-1; i++) s << indent;
	    s << lj(node->name, Header::name_width) << " = " << node->value << std::endl;
	}
	node = node->left;
    }
}


/** Recursive destructor. Build up a chain of
 * pointers to allow reverse traversal from end of a branch
 * back to its head.
 */
void Subs::deleter(Header::Hnode* node) {

    // Wind down to the end of a directory
    while(node->left != NULL){
	if(node->right != NULL){
	    deleter(node->right);
	    node->right = NULL;
	}
	node = node->left;
    }

    // Now wind back up from the bottom
    while(node->up != NULL && node->up->right == NULL){
	delete node->up->value;
	Header::Hnode *old = node;
	node = node->up;
	delete old;
    }

    delete node;
}


//! copy constructor
Subs::Header::Header(const Header& hd) : head(new Hnode()) {
    tail = Subs::copy(hd.head,head);
}
	
//! Assignment
Subs::Header& Subs::Header::operator=(const Header& hd){
    if(this != &hd){
	Subs::deleter(head);
	head = new Hnode();
	tail = Subs::copy(hd.head, head);
    }
    return *this;
}

Subs::Header::Hnode* Subs::copy(const Header::Hnode *node, Header::Hnode *cnode) {

    while(node->left != NULL){
	cnode->name     = node->name;
	cnode->value    = node->value->copy();
	cnode->left     = new Header::Hnode();
	cnode->left->up = cnode;
	if(node->right != NULL){
	    cnode->right     = new Header::Hnode();
	    cnode->right->up = cnode;
	    copy(node->right, cnode->right);
	}
	node  = node->left;
	cnode = cnode->left;
    }
    return cnode;
}

//! Clear a header
void Subs::Header::clear() {
    Subs::deleter(head);
    tail = head = new Hnode();
}


void Subs::Header::erase(const std::string& name) {
    try{
	Hnode *node = find(name);
	this->erase(node);
    }
    catch(const Header_Error& err){
	throw Header_Error("Subs::Header::erase(const std::string&): failed to find item = " + name);
    }
}

void Subs::Header::erase(Hnode* node){

    if(node->has_data()){
	// Delete any sub-directory first
	if(node->right) Subs::deleter(node->right);
	delete node->value;
	if(node->up != NULL){
	    if(node->up->left == node)
		node->up->left  = node->left;
	    else
		node->up->right = node->left;
	}else if(node->left != NULL){
	    head = node->left;
	    head->up = NULL;
	}
	delete node;
    }else{
	throw Header_Error("Subs::Header::erase(Hnode*): node contains no data");
    }

}


void Subs::Header::iterator::operator++(int){

    if(ptr->right != NULL){
	ptr = ptr->right;
    }else if(ptr->left != NULL){
	ptr = ptr->left; 
    }

    ptr = Subs::upstepper(ptr, tail);
}

void Subs::Header::const_iterator::operator++(int){

    if(ptr->right != NULL){
	ptr = ptr->right;
    }else if(ptr->left != NULL){
	ptr = ptr->left; 
    }

    ptr = Subs::upstepper(ptr, tail);
}

Subs::Header::Hnode* Subs::upstepper(Header::Hnode* ptr, Header::Hnode* tail){
    if(ptr->left == NULL && ptr != tail){

	while(ptr->up != NULL && (ptr->up->right == NULL || ptr->up->left == ptr)){
	    ptr = ptr->up;
	}
	if(ptr->up != NULL && ptr->up->right != NULL)
	    ptr = ptr->up->left;

	ptr = upstepper(ptr, tail);
    }
    return ptr;
}

const Subs::Header::Hnode* Subs::upstepper(const Header::Hnode* ptr, const Header::Hnode* tail){
    if(ptr->left == NULL && ptr != tail){

	while(ptr->up != NULL && (ptr->up->right == NULL || ptr->up->left == ptr)){
	    ptr = ptr->up;
	}
	if(ptr->up != NULL && ptr->up->right != NULL)
	    ptr = ptr->up->left;

	ptr = upstepper(ptr, tail);
    }
    return ptr;
}

bool Subs::operator!=(const Header::iterator& it1, const Header::iterator& it2){
    return it1.ptr != it2.ptr;
}

bool Subs::operator!=(const Header::const_iterator& it1, const Header::const_iterator& it2){
    return it1.ptr != it2.ptr;
}

/** This winds down to the bottom of whatever directory node is in
 */
Subs::Header::Hnode* Subs::end_node(Header::Hnode* node){
    while(node->left != NULL){
	node = node->left;
    }
    return node;
}

/** Renames a header item. The old_name must exist and the directory
 * structure for new_name must exist as well, but old_name itself must not
 * exist.
 * \param old_name initial name of item or directory
 * \param new_name new name of item or directory
 */
void Subs::Header::rename(const std::string& old_name, const std::string& new_name){

    if(new_name == old_name) return;

    // Initial checks
    Hnode *onode = this->find(old_name);
    if(!onode->has_data())
	throw Header_Error("Header::rename(const std::string&, const std::string&): failed to find old_name = " + old_name);
    
    Hnode *cnode = this->find(new_name);
    if(cnode->has_data())
	throw Header_Error("Header::rename(const std::string&, const std::string&): new_name = " + new_name + " already present");

    Hnode *nnode = head;
    std::string::size_type n = new_name.find_last_of(Header::dir_flag);
    std::string nname = new_name;
    if(n != std::string::npos){
	nnode = this->find(new_name.substr(0, n));
	if(!nnode->has_data() || nnode->right == NULL)
	    throw Header_Error("Header::rename(const std::string&, const std::string&): no containing directory for new_name = " + new_name);
	nnode = nnode->right;
	nname = new_name.substr(n+1);
    }

    nnode = Subs::end_node(nnode);
 
    // At this point nnode points to the last (dummy) Hnode of the directory that we are going to
    // move the item to. We want to rename the item to be just in front of this.

    // First cut out onode from where it is, accounting for possibility that the node is at the 
    // head of the entire tree
    if(onode->up != NULL){
	if(onode->up->left == onode)
	    onode->up->left  = onode->left;
	else
	    onode->up->right = onode->left;
	onode->left->up = onode->up;
    }else{
	head = onode->left;
	head->up = NULL;
    }

    // Now re-attach to its new place 
    if(nnode->up->left == nnode)
	nnode->up->left = onode;
    else
	nnode->up->right = onode;
    onode->up = nnode->up;
    onode->left = nnode;
    nnode->up = onode;
    onode->name = nname;
}

/** Moves a header item one up in its directory 
 */

void Subs::Header::move_up(const std::string& item){

    // Check that item exists
    Hnode *node = this->find(item);
    if(!node->has_data())
	throw Header_Error("Header::move_up(const std::string&): failed to find item = " + item);
 
    // Check if already at top of tree or directory
    if(node == head || node->up->right == node)
	return;

    // OK, now move. Make sure that whatver was the previous node to node point at whatever
    // node originally pointed to and that node points back to whatever pnode pointed back to 
    // Attach node before the node previous to node to node.
    Hnode *pnode = node->up;
    pnode->left  = node->left;
    node->up     = pnode->up;
    if(pnode->up != NULL){
	if(pnode->up->left == pnode)
	    pnode->up->left  = node;
	else
	    pnode->up->right = node;
    }
    node->left   = pnode;
    pnode->up    = node;

}

/** Moves a header item one up in its directory 
 */

void Subs::Header::move_to_top(const std::string& item){

    // Check that item exists
    Hnode *node = this->find(item);
    if(!node->has_data())
	throw Header_Error("Header::move_up(const std::string&): failed to find item = " + item);
 
    // Check if already at top of tree or directory
    if(node == head || node->up->right == node)
	return;

    // attach node above node to the node following node
    node->left->up = node->up;
    if(node->up != NULL){
	if(node->up->left == node)
	    node->up->left  = node->left;
	else
	    node->up->right = node->left;
    }
    
    // node is now free floating

    // Find the top of the directory
    Hnode* top = node;
    while(top->up != NULL && top->up->right != top){
	top = top->up;
    }
    

    // stick node in there
    if(top->up != NULL)
	top->up->right = node;
    node->left = top;
    node->up   = top->up;
    top->up    = node;
    if(head == top) head = node;

}

