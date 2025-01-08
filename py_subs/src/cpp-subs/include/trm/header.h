#ifndef TRM_HEADER
#define TRM_HEADER

#include "trm/hitem.h"
#define PCRE2_CODE_UNIT_WIDTH 16
#include <pcre2.h>

namespace Subs {

    //! Header class
    /** This stores header items in a tree structure to allow 
     * arbitrary ordering.
     */
    class Header {

    public:

	//! Width to use for header item names
	static Subs::UINT4 name_width;
	
	//! indentation to mark directories
	static const int  dir_indent = 3;   
	
	//! character to indicate a new indentation level
	static const char dir_flag   = '.'; 
	
	//! string to start each line of header output with
	static std::string start_string; 

	//! Hnode structure
	/** Represents one node of a tree-structured header. Each node has a name and associated value
	 * contained in an Hitem together with one or two pointers to other nodes. The 'left-hand' one
	 * stays in the same directory; the right-hand one, if set, dives into sub-directories. Every
	 * directory in a Header ends with a node with the left, value and right pointers set = NULL.
	 */
	struct Hnode {
	    
	    //! Default constructor
	    Hnode() : name(), value(NULL), left(NULL), right(NULL), up(NULL) {}
	    
	    //! Name for the node
	    std::string name;
	    
	    //! Associated value
	    Hitem* value;
	    
	    //! Pointer to left-hand branch in inverted tree
	    Hnode* left;
	    
	    //! Pointer to right-hand branch in inverted tree
	    Hnode* right;

	    //! Pointer back up the tree
	    Hnode* up;
	    
	    //! Checks whether it is a data node
	    bool has_data() const {return left != NULL;}

	    //! Checks whether it is not a data node
	    bool has_no_data() const {return left == NULL;}

	    //! Returns fullname of node
	    std::string fullname() const;
	    
	};

	// Now onto the header functions
	
	//! Default constructor
	Header() : head(new Hnode()), tail(head) {}

	//! Destructor
	~Header();

	//! copy constructor
	Header(const Header& hd);
	
	//! Assignment
	Header& operator=(const Header& hd);

	//! Print to s
	void print(std::ostream& s) const;

	//! write out in binary
	void write(std::ofstream& s) const; 

	//! read from binary
	void read(std::ifstream& s, bool swap_bytes);

	//! skip binary header
	static void skip(std::ifstream& s, bool swap_bytes);

	//! Set a header item
	void set(const std::string& name, Subs::Hitem* hip);

	//! Set a header item, forcing creation of directories if necessary
	void set_force(const std::string& name, Subs::Hitem* hip);

	//! Rename a header item
	void rename(const std::string& old_name, const std::string& new_name); 

	//! Is 'name' a valid header item name
	static bool valid_name(const std::string &name);

	//! Checks for the existence of a directory containing 'name'
//	bool check_for_dir(const std::string& name) const;

	//! Finds Hnode* corresponding to name
	Hnode* find(const std::string& name) const;

	//! Finds Hnode* corresponding to name using regular expression matching
	Hnode* find(const std::string& regexp, int& nmatch) const;

	//! Finds directory dir
	Hnode* find_dir(const std::string& name) const;

	//! Finds value corresponding to name
	Hitem* operator[](const std::string& name) {
	    Hnode *node = find(name);
	    if(node->has_data())
		return node->value;
	    else
		throw Header_Error("Header::operator[](const std::string&): no item called " + name + " could be found");
	}

	//! Clear a header
	void clear();

	//! Erase an item (and any sub-directories it might have)
	void erase(const std::string& name);

	//! Erase an item (and any sub-directories it might have)
	void erase(Hnode* node);

	//! Move an item one up in its directory
	void move_up(const std::string& item);

	//! Move an item to the top of its directory
	void move_to_top(const std::string& item);

	//! Finds value 
	Hitem* operator[](const std::string& name) const {
	    Hnode *node = find(name);
	    if(node->has_data())
		return node->value;
	    else
		throw Header_Error("Header::operator[](const std::string&): no item called " + name + " could be found");
	}

	//! Counts number of header items.
	Subs::UINT4 count() const;

	//! Error class for header class
	class Header_Error : public Subs_Error { 
	public:
	    //! Default constructor
	    Header_Error () : Subs_Error() {}
	    //! Constructor from a string
	    Header_Error (const std::string& str) : Subs_Error(str) {}
	};

	class iterator {
	public:

	    iterator() : ptr(NULL), tail(NULL) {}

	    iterator(Hnode* ptr, Hnode* tail) : ptr(ptr), tail(tail) {}

	    Hnode* operator->(){
		return ptr;
	    }

	    void operator++(int);

	    friend bool operator!=(const iterator& it1, const iterator& it2);

	private:
	    Hnode*  ptr;
	    Hnode*  tail;
	};

	iterator begin() {return iterator(this->head, this->tail);}

	iterator end() {return iterator(this->tail, this->tail);}

	class const_iterator {
	public:

	    const_iterator() : ptr(NULL), tail(NULL) {}

	    const_iterator(Hnode* ptr, Hnode* tail) : ptr(ptr), tail(tail) {}

	    const Hnode* operator->(){
		return ptr;
	    }

	    void operator++(int);

	    friend bool operator!=(const const_iterator& it1, const const_iterator& it2);

	private:
	    const Hnode*  ptr;
	    const Hnode*  tail;
	};

	const_iterator begin() const {return const_iterator(this->head, this->tail);}

	const_iterator end() const {return const_iterator(this->tail, this->tail);}
	    

    private:
	
	// Head of the tree.
	Hnode *head;
	// Tail of the tree
	Hnode *tail;

    };

    //! Counts nodes, starting from node
    Subs::UINT4 count(const Header::Hnode *node);

    //! Sets a node, starting from node, creating it if need be
    void set(Header::Hnode *node, const std::string& name, Hitem* hip);

    //! Sets a node, starting from node, creating it if need be and any directories
    void set_force(Header::Hnode *node, const std::string& name, Hitem* hip);

    //! Prints a node, starting from node
    void print(std::ostream& s, const Header::Hnode* node, int nlev, const std::string& start);

    //! Writes a node in binary format starting from node
    void write(std::ofstream& s, const Header::Hnode *node, const std::string& dir);

    //! Recursive destructor.
    void deleter(Header::Hnode* node);

    //! Returns the end node of a directory
    Header::Hnode* end_node(Header::Hnode* node);

    //! Finds a node
    Header::Hnode* find(const std::string& name, Header::Hnode* node);

    Header::Hnode* find(const std::string& regexp, int& nmatch);

    //! Regular expression find
    void find(pcre2_code* re, const std::string& dir, Header::Hnode *node, int& nmatch, Header::Hnode **fnode);

    //! Copy from node to cnode
    Header::Hnode* copy(const Header::Hnode *node, Header::Hnode *cnode);

    //! ASCII output of a Header
    std::ostream& operator<<(std::ostream& s, const Header& obj);

    Header::Hnode* upstepper(Header::Hnode* ptr, Header::Hnode* tail);

    const Header::Hnode* upstepper(const Header::Hnode* ptr, const Header::Hnode* tail);
};

#endif
