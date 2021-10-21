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

#ifndef comp_bio_language
#define comp_bio_language
#include <string>

class language {
public:
  language() {}
  language(string name, int markov,vector<string> options){}
  virtual  void initialize()=0;
  virtual double calculation(int pos, int seq)=0;
  virtual double recursive_calculation(int pos, int seq)=0;
  virtual int get_len(int seq)=0;
  virtual int get_markov(int seq)=0;
  virtual int get_seg_stop(int seq)=0;
  virtual int get_seg_min(int seq)=0;
  virtual int get_seg_max(int seq)=0;
  virtual string get_name(int seq)=0;
  virtual double get_limit(int seq)=0;
  virtual double get_adjustment(int seq)=0;
  virtual void write_details(int random)=0;
  virtual void write_details(int random,int seq)=0;
  virtual void write_seg(vector<vector<double> > &segs, int nbr_segs, int seq)=0;
  virtual void write_res(double evidence, int nbr_segs, int seq)=0;
  virtual void store_prob_details(int start, int stop,int seq)=0;
  virtual void reset_probs(int seq)=0;
  virtual int get_nbr_seq()=0;
};
#endif
