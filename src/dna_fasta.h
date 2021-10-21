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

#ifndef comp_bio_dna_fasta
#define comp_bio_dna_fasta
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sys/ddi.h>
#include "language.h"

class dna_fasta_sequence {
private:
  double language_size;
  int markov;
  int mark_ct;
  int total_ct;
  vector<double> mark_count;
  vector<double> totals;
  vector<vector<double> > probs;
  vector<int> nbr_cuts;

  //Pseudo counts
  vector<double> pseudo_counts;
  double pseudo_total;

  //Sequence
  int length;
  string actual;  //Real version of the sequence.
  string cleaned; //Version of Sequence without N characters
  string name;

  //Other values
  int seg_stop;
  int seg_min;
  int seg_max;
  double limit;
  double adjustment; // Log scale
  int raw;
  int full;
  int segment;

public:
  dna_fasta_sequence() {
    language_size = 4.0;
  }
  dna_fasta_sequence(vector<string> &options) 
  {
    language_size = 4.0;
    mark_ct = 0;
    total_ct = 0;
    adjustment = 0;
    seg_stop = 0;
    seg_min = 0;
    seg_max = 0;
    raw = 0;
    full = 0;
    segment = 0;
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
	else if((*p).compare("full") == 0) {
	  full = 1;
	}
	/*	else if((*p).compare("-i") == 0) {
	  raw = 1;
	  }*/
	
      }
  }
  void initialize(string seq_name, string sequence, int markov_val) 
    {
      int i;
      int j;
      name = seq_name;
      actual = sequence;
      markov = markov_val;
      mark_ct = (int)pow(language_size,(double)markov+1);
      total_ct = (int)pow(language_size,(double)markov);
      pseudo_counts = vector<double>(mark_ct,1.0);
      mark_count = vector<double>(mark_ct,0.0);
      totals = vector<double>(total_ct,0.0);

      pseudo_total = language_size;
      // Calculate Limit
      limit = pow(10.0,50.0);

      // Calculate Seg_stop
      if(!seg_stop) {
	if(actual.length() > 2000)
	  seg_stop = min((int)actual.length()/200,199);
	else 
	  {
	    seg_stop = min((int)actual.length()/(markov+2),199);
	    seg_stop /=8;
	  }
	seg_stop++;
	if(actual.length() < 100)
	  seg_stop += actual.length() % seg_stop;
      }
      if(!seg_min) {
	// Calculate seg_min
	seg_min = 1;
      }
      // Clean X values out of string
      cleaned = actual;
      i = cleaned.find("N");
      while(i < string::npos) {
	cleaned.replace(i,1,"",0,0);
	i = cleaned.find("N",i);
      }
      i = cleaned.find("n");
      while(i < string::npos) {
	cleaned.replace(i,1,"",0,0);
	i = cleaned.find("n",i);
      }

      if(!seg_max || seg_max > cleaned.length()) {
	// Calculate seg_max
	seg_max = cleaned.length();
      }
      if (seg_stop < (cleaned.length() / seg_max))
	seg_stop = cleaned.length() / seg_max;
      
      //Set up probablilty array.
      if(!raw)
	probs = vector<vector<double> >(cleaned.length(), vector<double>((int)mark_ct,0.0));
      else
	probs = vector<vector<double> >(cleaned.length(), vector<double>((int)language_size,0.0));

      nbr_cuts = vector<int>(cleaned.length(),0);

      // Calculate Adjustment
      if(!adjustment) {
	adjustment = 0.0;
	double temp_adjustment = 0.0;
	temp_adjustment = calculation(0);
	for(i = 1; i < cleaned.length(); i++)
	  {
	    temp_adjustment += recursive_calculation(i);
	  }
	adjustment = (-temp_adjustment -150)/cleaned.length();
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
    double temp = 0.0;

    if((pos-markov) >=0) {
      location = generateLocation(cleaned,pos-markov, pos);
      mark_count[location]++;
      totals[location/language_size]++;
      for(p = 0; p < mark_ct; p++)
	proba += lgamma(mark_count[p] + pseudo_counts[p])-lgamma(pseudo_counts[p]);
      proba_two = lgamma(pseudo_total * (double)total_ct);
      for(p = 0; p < total_ct; p++)
	temp += totals[p]+pseudo_total;
      proba_three = lgamma(temp);
      proba = proba + proba_two - proba_three + adjustment;
      return proba;
    }
    else
      return 0.0;
    return 0.0;
  }
  double recursive_calculation(int pos) 
  {
      int p;
    int location = 0;
    int loc_temp = 1;
    double proba = 0.0;
    if((pos-markov) > 0)
      {
	location = generateLocation(cleaned,pos-markov, pos);
	mark_count[location]++;
	totals[location/language_size]++;
	proba = log((double)(mark_count[location] + pseudo_counts[location]-1)/(double)(totals[location/language_size]+pseudo_total-1));
	proba+=adjustment;
	return proba;
      }
    else
      return 0.0;
    return 0.0;
  }
  int generateLocation(string &sequence, int start, int end)
  {
    int i;
    int location = 0;
    for(i = start; i <= end; i++)
      {
	location *= language_size;
	switch(sequence[i]) 
	  {
	  case 'A':
	    case 'a':
	      location += 0 ;
	      break;
	    case 'T':
	    case 't':
	      location += 1 ;
	      break;
	    case 'C':
	    case 'c':
	      location += 2 ;
	      break;
	    case 'G':
	    case 'g':
	      location += 3 ;
	      break;
	  }
      }
    return location;
  }
  int get_markov()
  {
      return markov;
  }
  int get_len() 
  {
    return cleaned.length();
  }
  int get_real_len()
  {
    return actual.length();
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
  void set_name(string new_name)
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
    return name;
  }
  double get_limit()
    {
      return limit;
    }
  double get_adjustment()
    {
      return adjustment;
    }
  int get_full() 
  {
    return full;
  }
  int get_segment() {
    return segment;
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
      nbr_cuts[cleaned.length()-1] = 0;
      
      random++;
      for(i = 0,k=0; i < markov; i++,k++)
	{
	  if(cleaned[i] != actual[k]) {
	    for(n_val = 0;cleaned[i] != actual[k+n_val] && k+n_val < actual.length(); n_val++);
	    for(l=0; l<n_val;l++)
	      {
		outfile << "  " << seq_nbr+1 << "\t" << 0;
		for(j = 0; j < mark_ct; j++)
		  outfile << "\t" << 0.25;
		outfile  << "\t" << k+l+1  << "\t" << random  << "\t" << 0.0 / (double) random << endl;
	      }
	    k += n_val;
	  }
	  outfile << "  " << seq_nbr+1 << "\t" <<  nbr_cuts[i];
	  for(j = 0; j < mark_ct; j++)
	    outfile << "\t" << 0.25;
	  outfile  << "\t" << k+1  << "\t" << random  << "\t" << (double)nbr_cuts[i] / (double) random << endl;
	}
      for(i = markov; i < cleaned.length(); i++,k++)
	{
	  if(cleaned[i] != actual[k]) {
	    for(n_val = 0;cleaned[i] != actual[k+n_val] && k+n_val < actual.length(); n_val++);
	    for(l = 0; l < n_val; l++)
	      {
		outfile << "  " << seq_nbr+1 << "\t" << 0;
		for(j = 0; j < mark_ct; j++) {
		  if(i == 0)
		    outfile << "\t" << (probs[i][j])/(double)(random-1);
		  else
		    outfile << "\t" << (probs[i][j] + probs[i-1][j])/(double)(2*random);
		}
		outfile  << "\t" << k+l+1  << "\t" << random  << "\t" << 0.0 / (double) random << endl;
	      }
	    k += n_val;
	  }
	  outfile << "  " << seq_nbr+1 << "\t" <<  nbr_cuts[i];
	  for(j = 0; j < mark_ct; j++)
	    outfile << "\t" << probs[i][j]/(double)(random-1);
	    outfile  << "\t" << k+1  << "\t" << random+1  << "\t" << (double)nbr_cuts[i] / (double) random << endl;
	}
      outfile << endl;
    }
  void write_details_raw(ofstream &outfile, int random, int seq_nbr)
    {
      int i,j,k,l;
      int n_val;
      outfile.setf(ios::scientific);
      nbr_cuts[cleaned.length()-1] = 0;
      
      random++;
      for(i = 0,k=0; i < markov; i++,k++)
	{
	  if(cleaned[i] != actual[k]) {
	    for(n_val = 0;cleaned[i] != actual[k+n_val] && k+n_val < actual.length(); n_val++);
	    for(l=0; l<n_val;l++)
	      {
		outfile << "  " << seq_nbr+1 << "\t" << 0;
		for(j = 0; j < language_size; j++)
		  outfile << "\t" << 0.25;
		outfile  << "\t" << k+l+1  << "\t" << random  << "\t" << 0.0 / (double) random << endl;
	      }
	    k += n_val;
	  }
	  outfile << "  " << seq_nbr+1 << "\t" <<  nbr_cuts[i];
	  for(j = 0; j < language_size; j++)
	    outfile << "\t" << 0.25;
	  outfile  << "\t" << k+1  << "\t" << random  << "\t" << (double)nbr_cuts[i] / (double) random << endl;
	}
      for(i = markov; i < cleaned.length(); i++,k++)
	{
	  if(cleaned[i] != actual[k]) {
	    for(n_val = 0;cleaned[i] != actual[k+n_val] && k+n_val < actual.length(); n_val++);
	    for(l = 0; l < n_val; l++)
	      {
		outfile << "  " << seq_nbr+1 << "\t" << 0;
		for(j = 0; j < language_size; j++) {
		  if(i == 0)
		    outfile << "\t" << (probs[i][j])/(double)(random-1);
		  else
		    outfile << "\t" << (probs[i][j] + probs[i-1][j])/(double)(2*random);
		}
		outfile  << "\t" << k+l+1  << "\t" << random  << "\t" << 0.0 / (double) random << endl;
	      }
	    k += n_val;
	  }
	  outfile << "  " << seq_nbr+1 << "\t" <<  nbr_cuts[i];
	  for(j = 0; j < language_size; j++)
	    outfile << "\t" << probs[i][j]/(double)(random-1);
	    outfile  << "\t" << k+1  << "\t" << random+1  << "\t" << (double)nbr_cuts[i] / (double) random << endl;
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
    outfile << markov << '\t' << actual.length() << endl;
  }
  void store_prob_details(int start, int stop)
  {
    int i;
    vector<int> total;
    vector<int> temp_counter;
    nbr_cuts[stop]++;
    temp_counter = vector<int>(mark_ct,0);
    total = vector<int>(mark_ct,0);
    for(i = stop; i > start; i--)
      {
	int location;
	if(i-markov >= 0)
	  location = generateLocation(cleaned, i-markov, i);	
	else
	  i = start;	
	if(i != start) {
	  temp_counter[location]++;
	  total[location/language_size]++;
	}
      }
    int j;
    for(i = stop; i > start; i--)
      {
	for(j = 0; j < mark_ct; j++)
	  {
	    probs[i][j] += ((double)temp_counter[j] + pseudo_counts[j])/ ((double)total[language_size] + pseudo_total);
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
	int location = 0;
	if(i-markov_cond >= 0)
	  location = generateLocation(cleaned, i-markov_cond, i);
	else
	  i = start;
	if(i != start) {
	  temp_counter[location]++;
	  total[(int)location/(int)language_size]++;
	}
      }
    int j;
    for(i = stop; i > start; i--)
      {
	if(i-markov_cond >= 0)
	  location =(int)generateLocation(cleaned, i-markov_cond, i) /(int)language_size;	
	else 
	  i = start;
	if(i != start) {
	  for(j = 0; j < language_size; j++)
	    {
	      
	      probs[i][j] += ((double)temp_counter[(location * language_size) + j] + pseudo_counts[j])/ ((double)total[location] + pseudo_total);
	    }
	}
      }

  }
  void write_segment_details(int start, int stop,ofstream &segmentFile) {
    vector<int> temp_counter;
    int i,j,total=0;
    int temp = segment;
    
    temp_counter = vector<int>(mark_ct,0);

    for(i = stop; i >= start; i--)
      {
	int location = 0;
	if(i-markov >= 0)
	  location = generateLocation(cleaned, i-markov, i);	
	else
	  i = start-1;
	if(i != start-1) {
	  temp_counter[location]++;
	  total++;
	}
      }
    if(segment == 0) {
      segmentFile << ">" << name << endl;
segment = 1;
    }
    segmentFile << "Segment " << segment << endl;
      segmentFile << "\tA\tT\tC\tG" << endl;
  for(i =0; i < total_ct; i++)
    {
      segmentFile << generateSymbol(i);
      temp = 0;
      for(j = 0; j < 4 ; j++)
	{
	  temp += temp_counter[(i*4)+j] + 1;
	}
      for(j = 0; j < 4 ; j++)
	{
	  segmentFile << "\t" << (double) (temp_counter[(i*4)+j] +1)/ (double) temp;
	}
      segmentFile << endl;
    }
  segmentFile.close();

  if (start == 0)
    segment = 1;
  else 
    segment++;
  }
string generateSymbol(int location) {
      string symbol = "";
      int j, loc_temp;
      if (markov == 0)
	symbol = "0";
      for(j = markov; j > 0; j--)
	{
	  loc_temp = location % 4;
	  switch(loc_temp)
	    {
	    case 0:
	      symbol += 'A';
	      break;
	    case 1:
	      symbol += 'T';
	      break;
	    case 2:
	      symbol += 'C';
	      break;
	    case 3:
	      symbol += 'G';
	      break;
	    }
	  location /= 4;
	}
      return symbol;
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

class dna_fasta:public language{
private:
  vector<dna_fasta_sequence> sequences;
  vector<string> options;
  string fileName;
  int markov;
  int nbr_sequences;
public:
  dna_fasta() {}
  dna_fasta(string file, int markov_val, vector<string> o) 
  {
    options = o;
    markov = markov_val;
    fileName = file;
  }
  void initialize() 
  {
    ifstream infile(fileName.data());
    if(!infile) {
      cout << "File could not be opened." << endl;
      exit(1);
    }
    char temp[2048];
    string temp_name;
    string temp_seq;
    string comp = ">";
    int nbr_of_seq = 0;
    while(infile.getline(temp,2048))
      {
	if (temp[0] == '>')
	  nbr_of_seq++;
      }
    nbr_sequences = nbr_of_seq;
    sequences = vector<dna_fasta_sequence>(nbr_sequences);
    nbr_of_seq = 0;
    infile.clear();
    infile.seekg(0);
    infile.getline(temp,2048);
    while(temp[0] != '>')
      {
	infile.getline(temp,2048);
      }
    temp_name = temp;
    while(infile.getline(temp,2048))
      {
	if(temp[0] == '>') {
	  if(temp_seq.length() == 0 ) {
	    cout << "Sequence is of length 0. Skipping sequence " << temp_name <<endl; 
	  } else {
	    sequences[nbr_of_seq] = dna_fasta_sequence(options);
	    sequences[nbr_of_seq].initialize(temp_name,temp_seq,markov);
	    nbr_of_seq++;
	  }
	  temp_name = temp;
	  temp_seq = "";
	}
	else
	  temp_seq = temp_seq + temp;
      }
    if(temp_seq.length() == 0 ) {
      cout << "Sequence is of length 0. Skipping sequence " << temp_name <<endl; 
    }
    else {
      sequences[nbr_of_seq] = dna_fasta_sequence(options);
      sequences[nbr_of_seq].initialize(temp_name,temp_seq,markov);
      nbr_of_seq++;
    }
    
    nbr_sequences = nbr_of_seq;
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
  void set_name(int seq, string new_name)
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
  vector<vector<double> > get_probs(int seq)
    {
      return sequences[seq].get_probs();
    }
  vector<int> get_nbr_cuts(int seq)
    {
      return sequences[seq].get_nbr_cuts();
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
    ofstream outFile;
    string temp = fileName + "_info-models";

    sequences[seq].store_prob_details(start,stop);

    if (sequences[seq].get_full() == 1) {
      if(seq == 0 && sequences[seq].get_segment() == 0) {
	outFile.open(temp.data(),ios::out);
      }
      else {
	outFile.open(temp.data(),ios::app);
      }
      sequences[seq].write_segment_details(start,stop,outFile);
    }
  }
  void store_conditioned_prob_details(int start, int stop,int seq, int markov_cond)
  {
    ofstream outFile;
    string temp = fileName + "_info-models";

    sequences[seq].store_conditioned_prob_details(start,stop, markov_cond);

    if (sequences[seq].get_full() == 1) {
      if(seq == 0 && sequences[seq].get_segment() == 0) {
	outFile.open(temp.data(),ios::out);
      }
      else {
	outFile.open(temp.data(),ios::app);
      }
      sequences[seq].write_segment_details(start,stop,outFile);
    }
  }
  void reset_probs(int seq)
  {
    sequences[seq].reset_probs();
  }
  
};

#endif

