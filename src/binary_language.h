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

#ifndef comp_bio_binary_language
#define comp_bio_binary_language
#include <string>
#include <vector>


class binary_language {
public:
  double size;
  binary_language() {
    size = 2.0;
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
	    symbol += 'H';
	    break;
	  case 1:
	    symbol += 'T';
	    break;
	  }
	location /= (int)size;
      }
    return symbol;
  }
  
  int symbol_value(char symbol) {
    int return_val= 0;
    switch(symbol) {
    case 'H':
    case 'h':
    case '0':
      return_val = 0 ;
      break;
    case 'T':
    case 't':
    case '1':
      return_val = 1 ;
      break;
    }
    return return_val;
  }
  
  vector<string> get_blanks() {
    vector<string> return_vector = vector<string>();
    return_vector.push_back(string("n"));
    return_vector.push_back(string("N"));
    return_vector.push_back(string("x"));
    return_vector.push_back(string("X"));
    return_vector.push_back(string(" "));
    return_vector.push_back(string("\t"));
    return return_vector;
  }
  int is_language(char c) {
    int return_val = 0;
    switch(c) {
    case 'H':
    case 'h':
    case 'T':
    case 't':
    case '0':
    case '1':
      return_val = 1 ;
      break;
    }
    return return_val;
  }
};
#endif
