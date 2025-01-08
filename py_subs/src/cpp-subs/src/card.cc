#include "trm/card.h"

std::ostream& operator<<(std::ostream& ost, const Subs::Card& card){
  
  switch(card.number()){
  case Subs::Card::Ace:
    ost << "Ace of ";
    break;
  case Subs::Card::Two:
    ost << "Two of ";
    break;
  case Subs::Card::Three:
    ost << "Three of ";
    break;
  case Subs::Card::Four:
    ost << "Four of ";
    break;
  case Subs::Card::Five:
    ost << "Five of ";
    break;
  case Subs::Card::Six:
    ost << "Six of ";
    break;
  case Subs::Card::Seven:
    ost << "Seven of ";
    break;
  case Subs::Card::Eight:
    ost << "Eight of ";
    break;
  case Subs::Card::Nine:
    ost << "Nine of ";
    break;
  case Subs::Card::Ten:
    ost << "Ten of ";
    break;
  case Subs::Card::Jack:
    ost << "Jack of ";
    break;
  case Subs::Card::Queen:
    ost << "Queen of ";
    break;
  case Subs::Card::King:
    ost << "King of ";
    break;
  }
  
  switch(card.suit()){
  case Subs::Card::Clubs:
    ost << "Clubs";
    break;
  case Subs::Card::Diamonds:
    ost << "Diamonds";
    break;
  case Subs::Card::Hearts:
    ost << "Hearts";
    break;
  case Subs::Card::Spades:
    ost << "Spades";
    break;
  }

  return ost;
}




















