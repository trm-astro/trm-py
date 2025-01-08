#ifndef TRM_CARD
#define TRM_CARD

#include <iostream>
#include <string>
#include "trm/subs.h"

namespace Subs {

  //! Class for representing playing cards

  class Card {

  public:
    
    //! Names of suits
    enum Suit {Clubs, Diamonds, Hearts, Spades};

    //! Names of cards
    enum Number {Ace, Two, Three, Four, Five, Six, Seven, 
		 Eight, Nine, Ten, Jack, Queen, King};

    //! Default constructor
    Card() : number_(Ace), suit_(Clubs) {}

    //! General constructor
    Card(Number number, Suit suit) : number_(number), suit_(suit) {}

    //! Constructor from a number, 0 = Ace of Clubs etc.
    Card(int ncard) {
      ncard %= 52;
      int nsuit = ncard / 13;
      int nnum  = ncard - 13*nsuit;
      suit_   = Suit(nsuit);
      number_ = Number(nnum);
    }

    //! Copy constructor
    Card(const Card& card) : number_(card.number_), suit_(card.suit_) {}

    //! Returns the suit of the card
    Suit suit() const {return suit_;}

    //! Returns the number of the card
    Number number() const {return number_;}

    class Card_Error : public Subs_Error {

    public:
      //! Default constructor
      Card_Error() : Subs::Subs_Error() {};
      //! Constructor from a string
      Card_Error(const std::string& str) : Subs::Subs_Error(str) {}
    };

   private:
    
    Number number_;

    Suit suit_;

  };

};

//! ASCII output

std::ostream& operator<<(std::ostream& ost, const Subs::Card& card);

#endif



















