#include <string>
#include "trm/subs.h"

/** Carries out a byte-swap on an integer as needed when translating between 
 * big- and little-endian machines. 
 *
 * \param i the integer to be byte swapped
 * \return the byte-swapped equivalent to i.
 */
int Subs::byte_swap(int i){
  reverse_bytes((char *)&i, sizeof(int));
  return i;
}

/** Carries out a byte-swap on a long integer as needed when translating between 
 * big- and little-endian machines. 
 *
 * \param i the long integer to be byte swapped
 * \return the byte-swapped equivalent to i.
 */
long int Subs::byte_swap(long int i){
  reverse_bytes((char *)&i, sizeof(long int));
  return i;
}

/** Carries out a byte-swap on a long unsigned integer as needed when translating between 
 * big- and little-endian machines. 
 *
 * \param i the long unsigned integer to be byte swapped
 * \return the byte-swapped equivalent to i.
 */
long unsigned int Subs::byte_swap(long unsigned int i){
  reverse_bytes((char *)&i, sizeof(long unsigned int));
  return i;
}

/** Carries out a byte-swap on an unsigned integer as needed when translating between 
 * big- and little-endian machines. 
 *
 * \param i the unsigned integer to be byte swapped
 * \return the byte-swapped equivalent to i.
 */
unsigned int Subs::byte_swap(unsigned int i){
  reverse_bytes((char *)&i, sizeof(unsigned int));
  return i;
}

/** Carries out a byte-swap on a float as needed when translating between 
 * big- and little-endian machines. 
 *
 * \param f the float to be byte swapped
 * \return the byte-swapped equivalent to i.
 */
float Subs::byte_swap(float f){
  reverse_bytes((char *)&f, sizeof(float));
  return f;
}

/** Carries out a byte-swap on a double as needed when translating between 
 * big- and little-endian machines. 
 *
 * \param f the double to be byte swapped
 * \return the byte-swapped equivalent to i.
 */
double Subs::byte_swap(double d){
  reverse_bytes((char *)&d, sizeof(double));
  return d;
}

/** Reverses order of bytes in the input array
 * \param buffer array of bytes
 * \param nbuff number of bytes. Values supported: 2, 4 and 8
 */
void Subs::reverse_bytes(char* buffer, const int& nbuff) {

  switch(nbuff){
  case 2:
    byte_swap2(buffer);
    break;

  case 4:
    byte_swap4(buffer);
    break;

  case 8:
    byte_swap8(buffer);
    break;
    
  default:
    throw Subs_Error("Subs::reverse_bytes(char*, const int&): unrecognised number of bytes = " + Subs::str(nbuff));
  }
}

// swap bytes in 2 byte buffer

void Subs::byte_swap2(char* buffer) {

  static char c;
  c         = buffer[0];
  buffer[0] = buffer[1];
  buffer[1] = c;
}

// swap bytes in 4 byte buffer

void Subs::byte_swap4(char* buffer) {

  static char c0, c1, c2;

  c0        = buffer[0];
  c1        = buffer[1];
  c2        = buffer[2];
  buffer[0] = buffer[3];
  buffer[1] = c2;
  buffer[2] = c1;
  buffer[3] = c0;

}

// swap bytes in 8 byte buffer

void Subs::byte_swap8(char* buffer) {

  static char c0, c1, c2, c3, c4, c5, c6;
  
  c0  = buffer[0];
  c1  = buffer[1];
  c2  = buffer[2];
  c3  = buffer[3];
  c4  = buffer[4];
  c5  = buffer[5];
  c6  = buffer[6];
  buffer[0] = buffer[7];
  buffer[1] = c6;
  buffer[2] = c5;
  buffer[3] = c4;
  buffer[4] = c3;
  buffer[5] = c2;
  buffer[6] = c1;
  buffer[7] = c0;
}

