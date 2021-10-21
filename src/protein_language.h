/****************************************************************/
/* unifiedcpp - The Bayesian Segmentation program               */
/*                                                              */
/* Please acknowledge the program authors on any publication of */
/* scientific results based in part on use of the program and   */
/* cite the following article in which the program was          */
/* described.                                                   */
/* Liu, J. and C. Lawrence (1999). "Bayesian inference on       */
/* biopolymer models." Bioinformatics 15(1): 38-52.             */
/*                                                              */
/* Copyright (C) 2006   Health Research Inc.                    */
/* HEALTH RESEARCH INCORPORATED (HRI),                          */
/* ONE UNIVERSITY PLACE, RENSSELAER, NY 12144-3455.             */
/* Email:  gibbsamp@wadsworth.org                              */
/*                                                              */
/*                                                              */
/* Copyright (C) 2009   Brown University                        */
/* Brown University                                             */
/* Providence, RI 02912                                         */
/* Email:  gibbs@brown.edu                                     */
/****************************************************************/
/*                                                              */
/* This program is free software; you can redistribute it       */
/* and/or modify it under the terms of the GNU General Public   */
/* License as published by the Free Software Foundation;        */
/* either version 2 of the License, or (at your option)         */
/* any later version.                                           */
/*                                                              */
/* This program is distributed in the hope that it will be      */
/* useful, but WITHOUT ANY WARRANTY; without even the implied   */
/* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR      */
/* PURPOSE. See the GNU General Public License for more         */
/* details.                                                     */
/*                                                              */
/* You should have received a copy of the GNU General Public    */
/* License along with this program; if not, write to the        */
/* Free Software Foundation, Inc., 51 Franklin Street,          */
/* Fifth Floor, Boston, MA  02110-1301, USA.                    */
/****************************************************************/

#ifndef comp_bio_protein_language
#define comp_bio_protein_language
#include <string>
#include <vector>


class protein_language {
public:
  double size;
  protein_language() {
    size = 20.0;
  }
  double get_size() { return size;}
  
  string generateSymbol(int location, int markov) {
    string symbol = "";
    int j, loc_temp;
    if (markov == 0)
      symbol = "0";
    for(j = markov; j > 0; j--)
      {
	loc_temp = location % (int)size;
	switch(loc_temp)
	  {
	  case 0:
	      symbol += "A" ;
	      break;
	    case 1:
	      symbol += "R" ;
	      break;
	    case 2:
	      symbol += "N";
	      break;
	    case 3:
	      symbol += "D";
	      break;
	    case 4:
	      symbol += "C";
	      break;
	    case 5:
	      symbol += "Q";
	      break;
	    case 6:
	      symbol += "E";
	      break;
	    case 7:
	      symbol += "G";
	      break;
	    case 8:
	      symbol += "H";
	      break;
	    case 9:
	      symbol += "I";
	      break;
	    case 10:
	      symbol += "L" ;
	      break;
	    case 11:
	      symbol += "K";
	      break;
	    case 12:
	      symbol += "M";
	      break;
	    case 13:
	      symbol += "F";
	      break;
	    case 14:
	      symbol += "P";
	      break;
	    case 15:
	      symbol += "S";
	      break;
	    case 16:
	      symbol += "T";
	      break;
	    case 17:
	      symbol += "W";
	      break;
	    case 18:
	      symbol += "Y";
	      break;
	    case 19:
	      symbol += "V";
	      break;
	  }
	location /= (int)size;
      }
    return symbol;
  }
  
  int symbol_value(char symbol) {
    int return_val= 0;
    switch(symbol) {
 	  case 'A':
	    case 'a':
	      return_val = 0 ;
	      break;
	    case 'R':
	    case 'r':
	      return_val = 1 ;
	      break;
	    case 'N':
	    case 'n':
	      return_val = 2 ;
	      break;
	    case 'D':
	    case 'd':
	      return_val = 3 ;
	      break;
	    case 'C':
	    case 'c':
	      return_val = 4 ;
	      break;
	    case 'Q':
	    case 'q':
	      return_val = 5 ;
	      break;
	    case 'E':
	    case 'e':
	      return_val = 6 ;
	      break;
	    case 'G':
	    case 'g':
	      return_val = 7 ;
	      break;
	    case 'H':
	    case 'h':
	      return_val = 8 ;
	      break;
	    case 'I':
	    case 'i':
	      return_val = 9 ;
	      break;
	    case 'L':
	    case 'l':
	      return_val = 10 ;
	      break;
	    case 'K':
	    case 'k':
	      return_val = 11 ;
	      break;
	    case 'M':
	    case 'm':
	      return_val = 12 ;
	      break;
	    case 'F':
	    case 'f':
	      return_val = 13 ;
	      break;
	    case 'P':
	    case 'p':
	      return_val = 14 ;
	      break;
	    case 'S':
	    case 's':
	      return_val = 15 ;
	      break;
	    case 'T':
	    case 't':
	      return_val = 16 ;
	      break;
	    case 'W':
	    case 'w':
	      return_val = 17 ;
	      break;
	    case 'Y':
	    case 'y':
	      return_val = 18 ;
	      break;
	    case 'V':
	    case 'v':
	      return_val = 19 ;
	      break;
    }
    return return_val;
  }
  
  vector<string> get_blanks() {
    vector<string> return_vector = vector<string>();
    return_vector.push_back(string("x"));
    return_vector.push_back(string("X"));
    return_vector.push_back(string(" "));
    return_vector.push_back(string("\t"));
    return return_vector;
  }
  int is_language(char c) {
    int return_val = 0;
    switch(c) 
      {
      case 'A':
      case 'a':
      case 'R':
      case 'r':
      case 'N':
      case 'n':
      case 'D':
      case 'd':
      case 'C':
      case 'c':
      case 'Q':
      case 'q':
      case 'E':
      case 'e':
      case 'G':
      case 'g':
      case 'H':
      case 'h':
      case 'I':
      case 'i':
      case 'L': 
      case 'l':
      case 'K':
      case 'k':
      case 'M':
      case 'm':
      case 'F':
      case 'f':
      case 'P':
      case 'p':
      case 'S':
      case 's':
      case 'T':
      case 't':
      case 'W':
      case 'w':
      case 'Y':
      case 'y':
      case 'V':
      case 'v':
	return_val = 1 ;
	break;
      }
    return return_val;
  }
};
#endif
