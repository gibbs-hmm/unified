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

#ifndef comp_bio_aligned_dna_fasta
#define comp_bio_aligned_dna_fasta
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sys/ddi.h>
#include "language.h"

class aligned_dna_fasta_sequence {
private:
  double language_size;
  int markov;
  int mark_ct;
  int total_ct;
  int raw;
  vector<double> mark_count;
  vector<double> totals;
  vector<vector<double> > probs;
  vector<int> nbr_cuts;

  //Pseudo counts
  vector<double> pseudo_counts;
  double pseudo_total;

  //Sequence
  int length;
  vector<string> actual;  //Real version of the sequence.
  vector<string> cleaned; //Version of Sequence without N characters
  vector<string> name;

  //Other values
  int seg_stop;
  int seg_min;
  int seg_max;
  int alignment;
  double limit;
  double adjustment; // Log scale

public:
  aligned_dna_fasta_sequence() 
  {
  }
  aligned_dna_fasta_sequence(vector<string> &options) 
  {
    raw = 0;
    mark_ct = 0;
    total_ct = 0;
    adjustment = 0;
    seg_stop = 0;
    seg_min = 0;
    seg_max = 0;
    vector<string>::iterator p;
    for(p = options.begin(); p != options.end(); p++)
      {
	if((*p).compare(0,2,"ss",0,2) == 0) {
	  seg_stop = atoi((*p).data()+2);
	}
	else if((*p).compare(0,2,"sl",0,2) == 0) {
	  seg_min = atoi((*p).data()+2);
	}
	else if((*p).compare(0,2,"sg",0,2) == 0) {
	  seg_max = atoi((*p).data()+2);
	}
	else if((*p).compare(0,2,"ad",0,2) == 0) {
	  adjustment = log(atof((*p).data()+2));
	}
     }
  }
  void initialize(vector<string> seq_name, vector<string> sequence, int markov_val) 
    {
      int i;
      int j;
      alignment = seq_name.size();
      name = seq_name;
      actual = sequence;
      language_size = pow(language.get_size(),(double)alignment);
      markov = markov_val;
      mark_ct = (int)pow(language_size,(double)markov+1);
      total_ct = (int)pow(language_size,(double)markov);
      for(i = 0; i < mark_ct; i++) {
	pseudo_counts.push_back(1.0);
	mark_count.push_back(0.0);
      }
      for(i = 0; i < total_ct;i++)
	totals.push_back(0.0);

      pseudo_total = language_size;
      // Calculate Limit
      limit = pow(10.0,100.0);

      // Calculate Seg_stop
      if(!seg_stop) {
	if(actual[0].length() > 2000)
	  seg_stop = min((int)actual[0].length()/200,199);
	else 
	  {
	    seg_stop = min((int)actual[0].length()/(markov+2),199);
	    seg_stop /=8;
	  }
	seg_stop++;
	if(actual[0].length() < 100)
	  seg_stop += actual[0].length() % seg_stop;
      }
      if(!seg_min) {
	// Calculate seg_min
	seg_min = 1;
      }
      // Clean N values out of string
      cleaned = actual;
      for(i = 0; i < cleaned[0].length();) {
	int remove = 0;
	vector<string> values = language.get_blanks();
	for(int k = 0; k < values.size(); k++) {
	  for(j = 0; j < alignment; j++) {
	    if(values[k] == cleaned[j][i])
	      remove = 1;
	  }
	}
	if(remove) {
	  for(j = 0; j < alignment; j++) {
	    cleaned[j].replace(i,1,"",0,0);
	  }
	}
	else {
	  i++;
	}
	
      }
      if(!seg_max || seg_max > cleaned[0].length()) {
	// Calculate seg_max
	seg_max = cleaned[0].length();
      }
      if (seg_stop < (cleaned[0].length() / seg_max))
	seg_stop = cleaned[0].length() / seg_max;
      
      //Set up probablilty array.
      if(!raw) {
	probs = vector<vector<double> >(cleaned[0].length(), vector<double>((int)mark_ct,0.0));
      }
      else {
	probs = vector<vector<double> >(cleaned[0].length(), vector<double>((int)language_size,0.0));
      }
      nbr_cuts = vector<int>(cleaned[0].length(),0);

      // Calculate Adjustment
      if(!adjustment) {
	adjustment = 0.0;
	double temp_adjustment = 0.0;
	temp_adjustment = calculation(0);
	for(i = 1; i < cleaned[0].length(); i++)
	  {
	    temp_adjustment += recursive_calculation(i);
	  }
	adjustment = (-temp_adjustment -150)/cleaned[0].length();
      }
    }
  double calculation(int pos) 
  {
    int p;
    int location = 0;
    int loc_temp = 1;
    double proba = 0.0;
    double proba_two;
    double proba_three=0.0;


    for(p = 0; p <=markov; p++)
      {
	if((pos-p) >=0) {
	  location = aligned_placement(pos-p);
	  loc_temp *=language_size;
	}
	else
	  return 0.0;
      }
    mark_count[location]++;
    totals[location/language_size]++;
    for(p = 0; p < mark_ct; p++)
      proba += lgamma(mark_count[p] + pseudo_counts[p])-lgamma(pseudo_counts[p]);
    proba_two = lgamma(pseudo_total) * (double)total_ct;
    for(p = 0; p < total_ct; p++)
      proba_three += lgamma(totals[p]+pseudo_total);
    proba_two -= proba_three;
    proba += proba_two;
    proba += adjustment;
    return proba;
  }
  double recursive_calculation(int pos) 
  {
      int p;
    int location = 0;
    int loc_temp = 1;
    double proba = 0.0;

    for(p = 0; p <=markov; p++)
      {
	if((pos-p) >=0) {
	  location = aligned_placement(pos-p);
	  loc_temp *=language_size;
	}
	else
	  return 0.0;
      }
    mark_count[location]++;
    totals[location/language_size]++;
    proba = log((mark_count[location] + pseudo_counts[location] - 1)/(totals[location/language_size]+pseudo_total-1));
    proba+=adjustment;
    return proba;
  }
  int aligned_placement(int p) 
  {
    int i;
    int location = 0;
    int loc_temp = 1;
    for(i = 0; i < alignment; i++) {
      location += language.symbol_value(cleaned[i]][p]) * loc_temp;
      loc_temp *=language.get_size();
    }
    return location;
  }
  bool equality_placement(int c, int a) 
  {
    int i;
    for(i = 0; i < alignment; i++) {
       if(cleaned[i][c] != actual[i][a])
	 return 0;
    }
    return 1;
  }
  int get_markov()
  {
      return markov;
  }
  int get_len() 
  {
    return cleaned[0].length();
  }
  int get_real_len()
  {
    return actual[0].length();
  }
  void set_seg_stop(int stop)
    {
      seg_stop = stop;
    }
  void set_seg_min(int min)
    {
      seg_min = min;
    }
  void set_seg_max(int max)
    {
      seg_max = max;
    }
  void set_name(vector<string> new_name)
  {
    name = new_name;
  }
  void set_limit(double nlimit)
    {
      limit = nlimit;
    }
  void set_adjustment(double adjust)
    {
      adjustment = adjust;
    }
  int get_seg_stop()
    {
      return seg_stop;
    }
  int get_seg_min()
    {
      return seg_min;
    }
  int get_seg_max()
    {
      return seg_max;
    }
  string get_name()
  {
    string temp=name[0];
    for(int i = 1; i < alignment; i++)
      temp += " and " + name[i];
    return temp;
  }
  double get_limit()
    {
      return limit;
    }
  double get_adjustment()
    {
      return adjustment;
    }
  vector<vector<double> > get_probs()
    {
      return probs;
    }
  vector<int> get_nbr_cuts()
    {
      return nbr_cuts;
    }
  void write_details(ofstream &outfile, int random, int seq_nbr)
    {
      int i,j,k,l;
      int n_val;
      outfile.setf(ios::scientific);
      nbr_cuts[cleaned[0].length()-1] = 0;
      random++;
      for(i = 0,k=0; i < markov; i++,k++)
	{
	  if(!equality_placement(i,k)) {
	    for(n_val = 0;!equality_placement(i,k+n_val) && k+n_val < actual[0].length(); n_val++);
	    for(l=0; l<n_val;l++)
	      {
		outfile << "  " << seq_nbr+1 << "\t" << 0;
		for(j = 0; j < mark_ct; j++)
		  outfile << "\t" << 0.0;
		outfile  << "\t" << k+l+1  << "\t" << random  << "\t" << 0.0 / (double) random << endl;
	      }
	    k += n_val;
	  }
	  outfile << "  " << seq_nbr+1 << "\t" <<  nbr_cuts[i];
	  for(j = 0; j < mark_ct; j++)
	    outfile << "\t" << 0.0;
	  outfile  << "\t" << k+1  << "\t" << random  << "\t" << (double)nbr_cuts[i] / (double) random << endl;
	}
      for(i = markov; i < cleaned[0].length(); i++,k++)
	{
	  if(!equality_placement(i,k)) {
	    for(n_val = 0;!equality_placement(i,k+n_val) && k+n_val < actual[0].length(); n_val++);
	    for(l = 0; l < n_val; l++)
	      {
		outfile << "  " << seq_nbr+1 << "\t" << 0;
		for(j = 0; j < mark_ct; j++){
		 if(i == 0)
		   outfile << "\t" << (probs[i][j])/(double)(random);
		 else
		   outfile << "\t" << (probs[i][j] + probs[i-1][j])/(double)(2*random);
		}
		outfile  << "\t" << k+l+1  << "\t" << random  << "\t" << 0.0 / (double) random << endl;
	      }
	    k += n_val;
	  }
	  outfile << "  " << seq_nbr+1 << "\t" <<  nbr_cuts[i];
	  for(j = 0; j < mark_ct; j++)
	    outfile << "\t" << probs[i][j]/(double)random;
	    outfile  << "\t" << k+1  << "\t" << random+1  << "\t" << (double)nbr_cuts[i] / (double) random << endl;
	}
      //Trailing n values.
      if(k != actual[0].length()) {
	for(;k < actual[0].length();k++) {
	  outfile << "  " << seq_nbr+1 << "\t" << 0;
	  for(j = 0; j < mark_ct; j++){
	      outfile << "\t" << (probs[i-1][j])/(double)(random);
	  }
	  outfile  << "\t" << k+l+1  << "\t" << random  << "\t" << 0.0 / (double) random << endl;
	}
      }
      outfile << endl;
    }
  void write_details_raw(ofstream &outfile, int random, int seq_nbr)
    {
      int i,j,k,l;
      int n_val;
      outfile.setf(ios::scientific);
      nbr_cuts[cleaned[0].length()-1] = 0;
      random++;
      for(i = 0,k=0; i < markov; i++,k++)
	{
	  if(!equality_placement(i,k)) {
	    for(n_val = 0;!equality_placement(i,k+n_val) && k+n_val < actual[0].length(); n_val++);
	    for(l=0; l<n_val;l++)
	      {
		outfile << "  " << seq_nbr+1 << "\t" << 0;
		for(j = 0; j < language_size; j++)
		  outfile << "\t" << 0.0;
		outfile  << "\t" << k+l+1  << "\t" << random  << "\t" << 0.0 / (double) random << endl;
	      }
	    k += n_val;
	  }
	  outfile << "  " << seq_nbr+1 << "\t" <<  nbr_cuts[i];
	  for(j = 0; j < language_size; j++)
	    outfile << "\t" << 0.0;
	  outfile  << "\t" << k+1  << "\t" << random  << "\t" << (double)nbr_cuts[i] / (double) random << endl;
	}
      for(i = markov; i < cleaned[0].length(); i++,k++)
	{
	  if(!equality_placement(i,k)) {
	    for(n_val = 0;!equality_placement(i,k+n_val) && k+n_val < actual[0].length(); n_val++);
	    for(l = 0; l < n_val; l++)
	      {
		outfile << "  " << seq_nbr+1 << "\t" << 0;
		for(j = 0; j < language_size; j++){
		 if(i == 0)
		   outfile << "\t" << (probs[i][j])/(double)(random);
		 else
		   outfile << "\t" << (probs[i][j] + probs[i-1][j])/(double)(2*random);
		}
		outfile  << "\t" << k+l+1  << "\t" << random  << "\t" << 0.0 / (double) random << endl;
	      }
	    k += n_val;
	  }
	  outfile << "  " << seq_nbr+1 << "\t" <<  nbr_cuts[i];
	  for(j = 0; j < language_size; j++)
	    outfile << "\t" << probs[i][j]/(double)random;
	    outfile  << "\t" << k+1  << "\t" << random+1  << "\t" << (double)nbr_cuts[i] / (double) random << endl;
	}
      //Trailing n values.
      if(k != actual[0].length()) {
	for(;k < actual[0].length();k++) {
	  outfile << "  " << seq_nbr+1 << "\t" << 0;
	  for(j = 0; j < language_size; j++){
	      outfile << "\t" << (probs[i-1][j])/(double)(random);
	  }
	  outfile  << "\t" << k+l+1  << "\t" << random  << "\t" << 0.0 / (double) random << endl;
	}
      }
      outfile << endl;
    }

 void write_seg(ofstream &outfile, vector<vector<double> > &segs, int nbr_segs, int seq)
  {
    int i;
    for(i = 0; i < seg_stop; i++) {
      outfile << "  " << seq +1 << '\t' << i  << '\t' << nbr_segs  << '\t'; 
      outfile << segs[i][0]   << '\t' <<  segs[i][1]  << '\t' <<  segs[i][2]  << '\t' <<  segs[i][3];
      outfile << endl;
    }
  }
  void write_res(ofstream &outfile, double evidence, int nbr_segs, int seq)
  {
    outfile.setf(ios::scientific);
    outfile << "  " << seq+1 << '\t' << evidence << '\t' << nbr_segs << '\t'; 
    outfile << markov << '\t' << actual[0].length() << endl;
  }
  void store_prob_details(int start, int stop)
  {
    int i;
    int total = 0;
    vector<int> temp_counter;
    nbr_cuts[stop]++;

    temp_counter = vector<int>(mark_ct,0);
    for(i = stop; i > start; i--)
      {
	int location = 0;
	int loc_temp = 1;
	for(int p = 0; p <=markov; p++)
	  {
	    if((i-p) >=0) {
	      location += aligned_placement(i-p) * loc_temp;
	      loc_temp *=language_size;
	    }
	    else 
	      i = start;
	  }
	if(i != start) {
	  temp_counter[location]++;
	  total++;
	}
      }
    int j;
    for(i = stop; i > start; i--)
      {
	for(j = 0; j < mark_ct; j++)
	  {
	    probs[i][j] += ((double)temp_counter[j] + pseudo_counts[j])/ ((double)total + pseudo_total);
	  }
      }

  }
  void store_conditioned_prob_details(int start, int stop, int markov_cond)
  {
    int i,location;
    vector<int> total ;
    vector<int> temp_counter;
    nbr_cuts[stop]++;
    temp_counter = vector<int>( (int)pow(language_size,(double)(markov_cond+1)),0);
    total = vector<int>( (int)pow(language_size,(double)markov_cond),0);

    for(i = stop; i > start; i--)
      {
	location = 0;
	loc_temp = 1;
	for(int p = 0; p <=markov_cond; p++)
	  {
	    if((i-p) >=0) {
	      location += aligned_placement(i-p) * loc_temp;
	      loc_temp *=language_size;
	    }
	    else 
	      i = start;
	  }
	if(i != start) {
	  temp_counter[location]++;
	  total[(int)location/(int)language_size]++;
	}
      }
    int j;
    for(i = stop; i > start; i--)
      {
	location = 0;
	loc_temp = 1;
	for(int p = 1; p <=markov_cond; p++)
	  {
	    if((i-p) >=0) {
	      location += aligned_placement(i-p) * loc_temp;
	      loc_temp *=language_size;
	    }
	    else 
	      i = start;
	  }
	if(i != start) {
	  for(j = 0; j < language_size; j++)
	    {
	      
	      probs[i][j] += ((double)temp_counter[(location * language_size) + j] + pseudo_counts[j])/ ((double)total[location] + pseudo_total);
	    }
	}
      }

  }
  void reset_probs() 
  {
    int pos;
    for(pos = 0; pos < mark_ct; pos++)
      mark_count[pos] = 0.0;
   for(pos = 0; pos < total_ct; pos++)
      totals[pos] = 0.0;
    
  }

};

class aligned_dna_fasta: public language {
private:
  vector<aligned_dna_fasta_sequence> sequences;
  vector<string> options;
  string fileName;
  int markov;
  int nbr_sequences;
  int alignment_value;
public:
  aligned_dna_fasta() {}
  aligned_dna_fasta(string file, int markov_val, vector<string> o) 
  {   
    alignment_value = 2;
    vector<string>::iterator p;
    for(p = options.begin(); p != options.end(); p++)
      {
	if((*p).compare(0,5,"alint",0,5) == 0) {
	  alignment_value = atoi((*p).data()+2);
	}
     }
    options = o;
    markov = markov_val;
    fileName = file;
  }
  void initialize() 
  {
    int seq_counter;
    ifstream infile(fileName.data());
    if(!infile) {
      cout << "File could not be opened." << endl;
      exit(1);
    }
    vector<string>::iterator p;
    for(p = options.begin(); p != options.end(); p++)
      {
	if((*p).compare(0,2,"al",0,2) == 0) {
	  alignment_value = atoi((*p).data()+2);
	}
      }
    char temp[512];
    string temp_name;
    vector<string> names;
    string temp_seq;
    vector<string> seqs;
    string comp = ">";
    int nbr_of_seq = 0;
    while(infile.getline(temp,512))
      {
	if (temp[0] == '>')
	  nbr_of_seq++;
      }
    nbr_sequences = nbr_of_seq/alignment_value;
    nbr_of_seq = 0;
    infile.clear();
    infile.seekg(0);
    infile.getline(temp,512);
    while(temp[0] != '>')
      {
	infile.getline(temp,512);
      }
    temp_name = temp;
    seq_counter = 1;
    while(infile.getline(temp,512))
      {
	if(temp[0] == '>') {
	  if(seq_counter != alignment_value) {
	    names.push_back(temp_name);
	    seqs.push_back(temp_seq);
	    temp_name = temp;
	    seq_counter ++;
	    temp_seq = "";
	  }
	  else {
	    names.push_back(temp_name);
	    seqs.push_back(temp_seq);
	    sequences.push_back(aligned_dna_fasta_sequence(options));
	    sequences[nbr_of_seq].initialize(names,seqs,markov);
	    names.clear();
	    seqs.clear();
	    nbr_of_seq++;
	    temp_name = temp;
	    seq_counter=1;
	    temp_seq = "";
	  }
	}
	else
	  temp_seq = temp_seq + temp;
      }
    names.push_back(temp_name);
    seqs.push_back(temp_seq);
    sequences.push_back(aligned_dna_fasta_sequence(options));
    sequences[nbr_of_seq].initialize(names,seqs,markov);
    names.clear();
    seqs.clear();
  }
  double calculation(int pos, int seq)
  {
    return sequences[seq].calculation(pos);
  }
  double recursive_calculation(int pos, int seq)
  {
    return sequences[seq].recursive_calculation(pos);
  }
  int get_len(int seq)
  {
    return sequences[seq].get_len();
  }
  int get_nbr_seq()
  {
    return nbr_sequences;
  }
  int get_markov(int seq)
  {
    return sequences[seq].get_markov();
  }
  int get_seg_stop(int seq)
  {
    return sequences[seq].get_seg_stop();
  }
  int get_seg_min(int seq)
  {
    return sequences[seq].get_seg_min();
  }
  int get_seg_max(int seq)
  {
    return sequences[seq].get_seg_max();
  }
  string get_name(int seq)
  {
    return sequences[seq].get_name();
  }
  double get_limit(int seq)
  {
    return sequences[seq].get_limit();
  }
  double get_adjustment(int seq)
  {
    return sequences[seq].get_adjustment();
  }
  void set_seg_stop(int seq, int stop)
  {
    sequences[seq].set_seg_stop(stop);
  }
  void set_seg_min(int seq, int min)
  {
    sequences[seq].set_seg_min(min);
  }
  void set_seg_max(int seq, int max)
  {
    sequences[seq].set_seg_max(max);
  }
  void set_name(int seq, vector<string> new_name)
  {
    sequences[seq].set_name(new_name);
  }
  void set_limit(int seq, double limit)
  {
    sequences[seq].set_limit(limit);
  }
  void set_adjustment(int seq, double adjust)
  {
    sequences[seq].set_adjustment(adjust);
  }
  void write_details_raw(int random, int seq, vector<vector<int> > &cuts)
  {
    string temp = fileName + "_info-det";
    ofstream outFile;
    if(seq == 0) {
      outFile.open(temp.data(),ios::out);
    }
    else {
      outFile.open(temp.data(),ios::app);
    }
    if(!outFile)
      {
	cout << "File could not be opened as " << temp << endl;
	exit(1);
      }
    temp = fileName + "_info-name";
    ofstream nameFile;
    if(seq == 0) {
      nameFile.open(temp.data(),ios::out);
    }
    else {
      nameFile.open(temp.data(),ios::app);
    }    
    if(!nameFile)
      {
	cout << "File could not be opened as " << temp << endl;
	exit(1);
      }
    temp = fileName + "_info-inclusive";
    ofstream inclusiveFile;
    if(seq == 0) {
      inclusiveFile.open(temp.data(),ios::out);
    }
    else {
      inclusiveFile.open(temp.data(),ios::app);
    }    
    if(!nameFile)
      {
	cout << "File could not be opened as " << temp << endl;
	exit(1);
      }
    sequences[seq].write_details_raw(outFile, random, seq);
    nameFile << sequences[seq].get_name() << endl;
    for(int i=0; i < cuts[0].size();i++)
      {
	inclusiveFile <<  "  " << seq << "\t" << i << "\t";
	for(int j=cuts.size()-1; j >= 0;j--)
	  {
	    inclusiveFile << cuts[j][i]/(double)random << "\t";
	  }
	inclusiveFile << endl;
      }
  }
  void write_details(int random, int seq)
  {
    string temp = fileName + "_info-det";
    ofstream outFile;
    if(seq == 0) {
      outFile.open(temp.data(),ios::out);
    }
    else {
      outFile.open(temp.data(),ios::app);
    }
    if(!outFile)
      {
	cout << "File could not be opened as " << temp << endl;
	exit(1);
      }
    temp = fileName + "_info-name";
    ofstream nameFile;
    if(seq == 0) {
      nameFile.open(temp.data(),ios::out);
    }
    else {
      nameFile.open(temp.data(),ios::app);
    }    
    if(!nameFile)
      {
	cout << "File could not be opened as " << temp << endl;
	exit(1);
      }
    
    sequences[seq].write_details(outFile, random, seq);
    nameFile << sequences[seq].get_name() << endl;
  }
  void write_details(int random)
  {
    string temp = fileName + "_info-det";
    ofstream outFile(temp.data(),ios::out);
    if(!outFile)
      {
	cout << "File could not be opened as " << temp << endl;
	exit(1);
      }
    temp = fileName + "_info-name";
    ofstream nameFile(temp.data(),ios::out);
    if(!nameFile)
      {
	cout << "File could not be opened as " << temp << endl;
	exit(1);
      }
    
    for(int i = 0; i < nbr_sequences; i++) {
      sequences[i].write_details(outFile, random, i);
      nameFile << sequences[i].get_name() << endl;
    }
  }
  void write_seg(vector<vector<double> > &segs, int nbr_segs, int seq)
  {
    string temp = fileName + "_info-seg";
    ofstream outFile;
    if(seq == 0) {
      outFile.open(temp.data(),ios::out);
    }
    else {
      outFile.open(temp.data(),ios::app);
    }
    if(!outFile)
      {
	cout << "File could not be opened as " << temp << endl;
	exit(1);
      }
      sequences[seq].write_seg(outFile, segs, nbr_segs,seq);
  }
  void write_res(double evidence, int nbr_segs, int seq)
  {
    string temp = fileName + "_info-res";
    ofstream outFile;
    if(seq == 0) {
      outFile.open(temp.data(),ios::out);
    }
    else {
      outFile.open(temp.data(),ios::app);
    }
    if(!outFile)
      {
	cout << "File could not be opened as " << temp << endl;
	exit(1);
      }
      sequences[seq].write_res(outFile, evidence, nbr_segs,seq);
  }
  void store_prob_details(int start, int stop,int seq)
  {
    sequences[seq].store_prob_details(start,stop);
  }
  void store_conditioned_prob_details(int start, int stop,int seq, int markov_cond)
  {
    sequences[seq].store_conditioned_prob_details(start,stop, markov_cond);
  }
  void reset_probs(int seq)
  {
    sequences[seq].reset_probs();
  }
  
};

#endif

