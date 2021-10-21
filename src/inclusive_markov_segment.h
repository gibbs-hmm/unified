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

/* inclusive_markov_segment.h
 *
 * last modified: 04-10-06 mjp
 *
 * 03-14-06
 *    - did some minor #ifdef stuff for cygwin compilation - including
 *	protos.h. also see the makefile make.unified.cygwin for a library
 *	module that must be available for linking to happen. the library 
 *	contains long double version of log and exp that this module calls.
 *
 * 03-13-06
 *    - commented out srand48() call in search_details() and placed it in
 *	main.cc 
 *
 * 03-09-06
 *    - finished reconciling. for the most part, i accepted the things in 
 *	ivan's code. we think his code is the most recent. 
 *
 *	i'm going to try running the code and producing an optimal solution.
 *	this is done by specifying the -r value as negative.
 *
 * 03-08-06
 *    - started reconciling this module; any 'debug' code in bill's version
 *	won't be added to this. his debug code was added for vpn bug.
 *
 * 03-06-06
 *    - copied this file from src.ivan/ to src.reconciled/. bill's src has
 *	things in it ivan's doesn't and ivan's src has more things in it that
 *	bill's doesn't. so, i'll start w/ ivan's and add bill changes to this
 *	module. ivan's version might have code in it that computes an optimal
 *	solution.
 *
 * 01-24-06:
 *    - started screwing with file to get it to compile. commented out
 *	the inclusion of sunmath.h; from past experience i've had to use
 *	sunmath.h to get expl(); on linux, expl() is in math.h. i hope 
 *	that's all it was needed for.
 */



#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

// this #define is pre-defined in sun CC compiler; this won't work w/ g++
#ifdef __SUNPRO_CC
#include "sunmath.h"
#endif

// added this from BT version; 
// this is need for compilation under cygwin. it's needed for prototypes for 
// logl, exp. also see the make.unified.cygwin makefile for directory 
// locations and a library file that must be available for cygwin compilation.
#ifdef CYGWIN_COMPILE
#include <protos.h>
#endif

// added this from BT version; not sure it's needed, but since it's 
// conditionally compiled, it can't hurt. it's prototypes for logl, expl etc. 
// for PPC Linux 
#ifdef _PPC64_
extern long double logl ( long double );
extern long double expl ( long double );
#endif


//Function Declarations
template <class T, class D>
void segment(vector<T> &current, int random, D &setup);

template <class T, class D>
void seg_in_two(vector<T> &current, vector<vector<vector<double> > >&proba_k,int seq);

template <class T, class D>
void seg_recu(vector<T> &current, vector<vector<vector<double> > > &proba_k,vector<int> &reajust,int seq);

template <class T, class D>
double find_err(vector<T> &current,vector<vector<double> > &proba_k,int position, int cutpts,int seq,int mark);

template <class T, class D>
int nbr_segments(vector<T> &current, vector<vector<double> > &proba_k,   vector<vector<double> > &segments, vector<int> reajust, double & evidence,int seq, int nbr_mark);

template <class T, class D>
void search_details(vector<T> &current,vector<vector<vector<double> > > &proba_k,  vector<vector<double> > &segments, int random, int seq,vector<vector<int> >  & markov_vals);

template <class T, class D>
void search_cut_details(vector<T> &current,vector<vector<vector<double> > > proba_k, int seg, vector<int > &cuts, int seq, vector<int> &markov_cuts);

template <class T, class D>
void search_maximum_cut_details(vector<T> &current,vector<vector<vector<double> > > proba_k, int seg, vector<int > &cuts, int seq, vector<int> &markov_cuts);


//End Function Declarations
// Set up algorithm to run and report results from each sequence
template <class T, class D>
void segment(vector<T> &current, int random, D &setup) {
  int i;
  int j,k,l;

  // among other things initialize(), reads the .fa file; i believe this
  // occurs for *each' markov level specified,.
  for(i = 0; i < current.size(); i++)
    current[i].initialize();

  vector<vector<vector<D> > > proba_k;
  vector<vector<D> > segments;
  vector<int> reajustments;
  D evidence;
  double adjustment = 0.0;
  int nbr_segs;

  cout << "There are " << current[0].get_nbr_seq() << " sequences in the file" << endl;

  vector<vector<int> > cuts;

  for(i = 0; i < current[0].get_nbr_seq(); i++) {
    adjustment = 0.0;
    for(j = 0; j < current.size(); j++)
      {
	adjustment += current[j].get_adjustment(i);
      }
    for(j = 0; j < current.size(); j++)
      {
	current[j].set_seg_stop(i,current[0].get_seg_stop(i));
	current[j].set_seg_min(i,current[0].get_seg_min(i));
	current[j].set_seg_max(i,current[0].get_seg_max(i));
	current[j].set_limit(i,current[0].get_limit(i));
	current[j].set_adjustment(i,(adjustment/(double)current.size()));
      }
    cout << "Analysing sequence named " << current[0].get_name(i) << endl;
    cout << "Seg Stop: " << current[0].get_seg_stop(i) << endl;
    cout << "Seg Min: " << current[0].get_seg_min(i) << endl;
    cout << "Seg Max: " << current[0].get_seg_max(i) << endl;
    cout << "Adjustment: " << exp(current[0].get_adjustment(i))<< endl;

    proba_k = vector<vector<vector<D> > >(current.size()+1, vector<vector<D> >(current[0].get_seg_stop(i)+1, vector<D>(current[0].get_len(i),0.0)));
    seg_in_two(current, proba_k, i);
    cout << "Full segment " << "=" << proba_k[current.size()][0][current[0].get_len(i)-1] << endl;
    seg_recu(current,proba_k, reajustments,i);
    nbr_segs = nbr_segments(current,proba_k[current.size()],segments,reajustments,evidence, i,current.size());
    // bill's version added some debug code to the following loop; i didn't 
    // add it in. 
    for(l = 0; l < proba_k[j].size(); l++)
      cout << "Segment with  " << l << " cutpoints in it=" << segments[l][3] << endl;
    cout << "Segments = " << nbr_segs << endl;
    cout << "Evidence = " << evidence << endl;

    cuts = vector<vector<int> >(current.size(),vector<int>(current[0].get_len(i),0));    

    search_details(current,proba_k,segments, random,i,cuts);

    // clean up for next run.
    if(random < 0) {
      current[0].write_details_raw((random*-1), i, cuts);
    }
    else {
      current[0].write_details_raw(random, i, cuts);
    }

    proba_k.clear();
    segments.clear();
    reajustments.clear();
  }
}

// Calculate initial probabilities for k = 0 for all markov levels
template <class T, class D>
void seg_in_two(vector<T> &current, vector<vector<vector<D> > > &proba_k,int seq){
  int pos;
  D proba;
  int end;
  int mark;
  for(mark = 0; mark < current.size(); mark++) {
    end = current[mark].get_seg_max(seq);
    current[mark].reset_probs(seq);
    pos = current[0].get_markov(seq);
    proba = exp(current[mark].calculation(pos,seq));
    proba_k[mark][0][pos] = proba;
    pos++;
    for(; pos < end; pos++)
      {
	proba *= exp(current[mark].recursive_calculation(pos,seq));
	proba_k[mark][0][pos] = proba;
	proba_k[current.size()][0][pos] += proba;
      }
  }
}

// Calculate all values for all positions for k > 0 
template <class T, class D>
void seg_recu(vector<T> &current, vector<vector<vector<D> > > &proba_k, vector<int> &reajust,int seq)
{
  int cutpts;
  bool flag = 0;
  int stop = current[0].get_len(seq);
  int mark;
  int end;
  for(cutpts = 1; cutpts <= current[0].get_seg_stop(seq); cutpts++) {
    flag = 0;
    for(mark = 0; mark < current.size(); mark++) {
      if(proba_k[mark][cutpts-1][stop-1] >= pow(current[0].get_limit(seq),1.0/(double)current.size()))
	{
	  flag = 1;
	}
    }
    if(flag == 1) {
      reajust.push_back(cutpts);
    }
    for(mark = 0; mark < current.size(); mark++) {
      int pos;
      if((cutpts*current[0].get_seg_min(seq)) > current[0].get_markov(seq))
	pos = cutpts*current[0].get_seg_min(seq);
      else
	pos = current[0].get_markov(seq);
      if((current[0].get_seg_max(seq) * cutpts) < current[0].get_len(seq))
	end = (current[0].get_seg_max(seq) * cutpts);
      else
	end = current[0].get_len(seq);
      for (; pos < end;pos++)
	{
	  if(flag == 1) {
	    proba_k[mark][cutpts][pos] = find_err(current,proba_k[current.size()],pos,cutpts-1,seq,mark) / (current[0].get_limit(seq));
	    proba_k[current.size()][cutpts][pos] +=  proba_k[mark][cutpts][pos];	    
	  }
	  else {
	    proba_k[mark][cutpts][pos] = find_err(current,proba_k[current.size()],pos,cutpts-1,seq,mark);
	    proba_k[current.size()][cutpts][pos] +=  proba_k[mark][cutpts][pos];
	  }
	}
    }
  }
}

//Calculate probability for a specific position and specific markov level
template <class T, class D>
D find_err(vector<T> &current,vector<vector<D> > &proba_k,int position, int cutpts,int seq, int mark)
{
  current[mark].reset_probs(seq);
  D ret_val = 0;
  D temp;
  int pos = position;
  int end = cutpts;
  if ( end < position - current[mark].get_seg_max(seq))
    end = position - current[mark].get_seg_max(seq);
  temp = exp(current[mark].calculation(pos,seq));
  ret_val = temp * proba_k[cutpts][pos-1];
  for(pos = position -1; (pos - current[mark].get_seg_min(seq)) > end; pos--)
    {
      temp *= exp(current[mark].recursive_calculation(pos,seq));
      ret_val += temp * proba_k[cutpts][pos-1];
    }
  return ret_val;
}

// Find optimal number of segments and provide evidence for the optimal 
// solution in the optimal cutpoint
template <class T, class D>
int nbr_segments(vector<T> &current, vector<vector<D> > &proba_k,   vector<vector<D> > &segments, vector<int> reajust, D & evidence,int seq, int nbr_mark) {
  D total;
  D temp, temp2,pres, val_zero, proba;
  D *table;
  int length, seg, end, ch,kk,pos;
  vector<int>::const_iterator p;
  ch = -2;

  length = current[0].get_len(seq) - current[0].get_markov(seq);
  table = new D[current[0].get_seg_stop(seq)+1];
  pres = 0;
  total = evidence = 0;
  evidence = logl(0.0);
  p = reajust.begin();
  end = current[0].get_len(seq);
  segments = vector<vector<D> >( current[0].get_seg_stop(seq)+1,vector<D>(4,0.0));
  for(seg = 0; seg <= current[0].get_seg_stop(seq); seg++) {
    if(p != reajust.end() && seg == *p)
      {
	pres += logl(current[0].get_limit(seq));
	p++;
      }
    segments[seg][0] = logl(proba_k[seg][end-1]) + pres;
    // Calculate N choose K where N = Length of Sequence - 1
    // N choose K = gamma(N+1) / (gamma(N - k + 1) * gamma(k+1))
    // Since Length of the sequence - 1 = N
    // we do gamma(Length) / (gamma(Length - k) * gamma(k+1))
    temp2 = lgamma(end) - (lgamma(end - seg) + lgamma(seg+1));
    segments[seg][2] = temp2;
    table[seg] = logl(proba_k[seg][end-1])-temp2+logl((double)current[0].get_prior(seq,seg));
    table[seg] += pres;
    if(seg > 0 && evidence < table[seg])
      evidence = table[seg];
    total+=expl(table[seg]);
  }
  if(current[0].get_len(seq) == current[0].get_seg_max(seq))
    val_zero = table[0];
  else
    val_zero = table[(current[0].get_len(seq)/current[0].get_seg_max(seq))+1];
  evidence = expl(val_zero)/(expl(val_zero) + expl(evidence));
  for(seg = 0; seg <= current[0].get_seg_stop(seq); seg++)
    {
      table[seg] = (expl(table[seg]))/total;
      segments[seg][3] = table[seg];
    }
  ch = end/current[0].get_seg_max(seq);
  if(current[0].get_len(seq) % current[0].get_seg_max(seq))
    ch ++;
  else if(current[0].get_len(seq) == current[0].get_seg_max(seq))
    ch = 0;
  else
    ch--;
  for(temp = 0.0, seg = ch; seg <= current[0].get_seg_stop(seq); seg++)
    if(table[seg] >= temp)
      {
	ch = seg + 1;
	temp = table[seg];
      }
  delete []table;
  //Save the seg file.
  return ch;
}

// Do random sample. Take a sampled number of cutpoints, pass to 
// search_cut_details the value of k. Store results from sample. Repeat. 
template <class T, class D>
void search_details(vector<T> &current,vector<vector<vector<D> > > &proba_k,  vector<vector<D> > &segments, int random, int seq,vector<vector<int> > & markov_vals)
{

  vector<D> cut_probs;
  vector<int> cuts;
  vector<int> markov_cuts;
  vector<int> pCuts;
  D cumlative = 0.0;
  double nbr;
  int stop;
  int rep;
  int start = 0;
  int segment;
  int i,j;

  pCuts = vector<int>(current[0].get_seg_stop(seq)+1,0);
  cut_probs = vector<D>(current[0].get_seg_stop(seq)+1);
  for(i = 0; i < (current[0].get_seg_stop(seq)+1);i++)
    {
      cut_probs[i] = segments[i+start][3];
      cumlative += cut_probs[i];
    }
  for(i = 0; i < (current[0].get_seg_stop(seq)+1);i++)
    {
      cut_probs[i] = cut_probs[i] / cumlative;
    }

  // 03-13-06 mjp: moved srand48()  call to main.cc
  // debug:  cout << "\n\nCalling  srand48()" << endl << endl;
  // srand48(time(0));

  if( random < 0) {
    rep = 1;
    cumlative = 0.0;
    nbr = 0.0;
    for(j = 0; j < (current[0].get_seg_stop(seq)+1); j++)
      {
	cumlative = cut_probs[j];
	if(nbr < cumlative ) {
	  rep = j;
	  nbr = cumlative;
	}
      }
    pCuts[rep]++;
    search_maximum_cut_details(current,proba_k,rep,cuts,seq,markov_cuts);
    stop = current[0].get_len(seq)-1;
    vector<int>::const_iterator p,m;
    m = markov_cuts.begin();
    markov_vals[*m][stop]++;
    for(p = cuts.begin(); p != cuts.end();p++,m++)
      {
	current[0].store_conditioned_prob_details(*p, stop,seq,current[(*m)].get_markov(seq));
	for(j = *p; j <stop; j++) {
	  markov_vals[*m][j]++;
	}
	stop = *p;
      }
    current[0].store_conditioned_prob_details(-1, stop,seq,current[(*m)].get_markov(seq));
    for(j = 0; j <stop; j++) {
      markov_vals[*m][j]++;
    }
    cuts.clear();
    markov_cuts.clear();
    random = (random * -1) - 1;
  }

  for(i = 0; i < random; i++) {
    nbr = drand48();
    rep = 1;
    cumlative = 0.0;
    for(j = 0; j < (current[0].get_seg_stop(seq)+1); j++) {
	if(nbr > cumlative && nbr < cumlative + cut_probs[j]) {
	  rep = j;
	  j = current[0].get_seg_stop(seq);
	}
	else
	  cumlative += cut_probs[j];
    }
    pCuts[rep]++;
    search_cut_details(current,proba_k,rep,cuts,seq,markov_cuts);

    stop = current[0].get_len(seq)-1;
    vector<int>::const_iterator p,m;
    m = markov_cuts.begin();
    markov_vals[*m][stop]++;
    for(p = cuts.begin(); p != cuts.end();p++,m++)    {
	current[0].store_conditioned_prob_details(*p, stop,seq,current[(*m)].get_markov(seq));

	for(j = *p; j <stop; j++) {
	  markov_vals[*m][j]++;
	}
	stop = *p;
      }
    current[0].store_conditioned_prob_details(-1, stop,seq,current[(*m)].get_markov(seq));
    for(j = 0; j <stop; j++) {
      markov_vals[*m][j]++;
    }
    cuts.clear();
    markov_cuts.clear();
  }  // for i < random 


  for(j =0; j < pCuts.size();j++) {
    cout << j << " Cuts chosen " << pCuts[j] << " times." << endl;
  }
  
}

//Take details of a given sampled number of cutpoints
template<class T, class D>
void search_cut_details(vector<T> &current,vector<vector<vector<D> > > proba_k, int seg, vector<int > &cuts, int seq, vector<int> &markov_cuts)
{
  int i,j,k,lowest=0;
  D cumlative=0.0;
  D res = 0.0;
  double nbr;
  int extra;
  int stop = current[0].get_len(seq);
  extra = seg - 1;
  for(k = 0; k < current.size(); k++)
    {
      cumlative += proba_k[k][extra+1][stop-1];
    }
  for(k = 0; k < current.size(); k++)
    {
      proba_k[k][extra+1][stop-1] /= cumlative;
    }
  nbr = drand48();
  cumlative = 0.0;
  for(k = 0; k < current.size(); k++)
    {
      if(nbr > cumlative && nbr <= cumlative + proba_k[k][extra+1][stop-1])
	{
	  markov_cuts.push_back(k);
	  k = current.size();
	}
      else
	cumlative += proba_k[k][extra+1][stop-1];
    }
  for(i = extra; i >=0;i--)
    {
      cumlative = 0.0;
      // Choose a markov level
      for(k = 0; k < current.size(); k++) 
	{
	  current[markov_cuts[extra-i]].reset_probs(seq);
	  res = exp(current[markov_cuts[extra-i]].calculation(stop-current[0].get_seg_min(seq), seq));
	  proba_k[k][i][stop-current[0].get_seg_min(seq)] = res * proba_k[k][i][stop-current[0].get_seg_min(seq)];
	  cumlative += proba_k[k][i][stop-current[0].get_seg_min(seq)];
	  if (current[0].get_markov(seq) >(current[0].get_seg_min(seq)*(i+1))) 
	    {
	      lowest = current[0].get_markov(seq);
	    }
	  else
	    {
	      lowest = current[0].get_seg_min(seq)* (i+1);
	    }
	  for(j = (stop-(current[0].get_seg_min(seq)+1)); j >= lowest; j--)
	    {
	      res *= exp(current[markov_cuts[extra-i]].recursive_calculation(j,seq));
	      proba_k[k][i][j] = res * proba_k[k][i][j];
	      cumlative += proba_k[k][i][j];
	    }
	}
      // randomly choose a spot to cut.
      for(k = 0; k < current.size(); k++)
	{
	  if (current[0].get_markov(seq) >(current[0].get_seg_min(seq)*(i+1))) 
	    {
	      j = current[0].get_markov(seq);
	    }
	  else
	    {
	      j = current[0].get_seg_min(seq)* (i+1);
	    }
	  for(; j <= (stop-current[0].get_seg_min(seq)); j++)
	    {
	      proba_k[k][i][j] /= cumlative;
	    }
	}
      nbr = drand48();
      cumlative = 0.0;
      for(k = 0; k < current.size(); k++)
	{
	  if (current[0].get_markov(seq) >(current[0].get_seg_min(seq)*(i+1))) 
	    {
	      j = current[0].get_markov(seq);
	    }
	  else
	    {
	      j = current[0].get_seg_min(seq)* (i+1);
	    }
	  for(; j <= (stop-(current[0].get_seg_min(seq))); j++)
	    {
	      if(nbr > cumlative && nbr <= cumlative + proba_k[k][i][j])
		{
		  cuts.push_back(j);
		  markov_cuts.push_back(k);
		  stop = j;
		  k = current.size();
		}
	      else
		cumlative += proba_k[k][i][j];
	    }
	} 
    }
}

//Calculate the optimal solution from the optimal value of k.
template<class T, class D>
void search_maximum_cut_details(vector<T> &current,vector<vector<vector<D> > > proba_k, int seg, vector<int > &cuts, int seq, vector<int> &markov_cuts)
{
  int i,j,k,lowest=0;
  D cumlative=0.0;
  D res = 0.0;
  double nbr;
  int extra;
  int stop = current[0].get_len(seq);
  int rep,store;
  extra = seg - 1;
  for(k = 0; k < current.size(); k++)
    {
      cumlative += proba_k[k][extra+1][stop-1];
    }
  for(k = 0; k < current.size(); k++)
    {
      proba_k[k][extra+1][stop-1] /= cumlative;
    }
  nbr = -1.0;
  cumlative = 0.0;
  for(k = 0; k < current.size(); k++)
    {
      cumlative = proba_k[k][extra+1][stop-1];
      if(nbr < cumlative )
	{
	  rep = k;
	  nbr = cumlative;
	}
    }
  // Choose markov level of the last segment
  markov_cuts.push_back(rep);
  for(i = extra; i >=0;i--)
    {
      cumlative = 0.0;
      for(k = 0; k < current.size(); k++) 
	{
	  current[markov_cuts[extra-i]].reset_probs(seq);
	  res = exp(current[markov_cuts[extra-i]].calculation(stop-current[0].get_seg_min(seq), seq));
	  proba_k[k][i][stop-current[0].get_seg_min(seq)] = res * proba_k[k][i][stop-current[0].get_seg_min(seq)];
	  cumlative += proba_k[k][i][stop-current[0].get_seg_min(seq)];
	  if (current[0].get_markov(seq) >(current[0].get_seg_min(seq)*(i+1))) 
	    {
	      lowest = current[0].get_markov(seq);
	    }
	  else
	    {
	      lowest = current[0].get_seg_min(seq)*(i+1);
	    }
	  for(j = stop-(1+current[0].get_seg_min(seq)); j >= lowest; j--)
	    {
	      res *= exp(current[markov_cuts[extra-i]].recursive_calculation(j,seq));
	      proba_k[k][i][j] = res * proba_k[k][i][j];
	      cumlative += proba_k[k][i][j];
	    }
	}
      // Choose a maximum markov level
      // choose the maximum spot to cut.
      for(k = 0; k < current.size(); k++)
	{
	  if (current[0].get_markov(seq) >(current[0].get_seg_min(seq)*(i+1))) 
	    {
	      j = current[0].get_markov(seq);
	    }
	  else
	    {
	      j = current[0].get_seg_min(seq)*(i+1);
	    }
	  for(; j < stop-(current[0].get_seg_min(seq)); j++)
	    {
	      proba_k[k][i][j] /= cumlative;
	    }
	}
      nbr = -1.0;
      cumlative = 0.0;
      for(k = 0; k < current.size(); k++)
	{
	  if (current[0].get_markov(seq) >(current[0].get_seg_min(seq)*(i+1))) 
	    {
	      j = current[0].get_markov(seq);
	    }
	  else
	    {
	      j = current[0].get_seg_min(seq)*(i+1);
	    }
	  for(; j < (stop-current[0].get_seg_min(seq)); j++)
	    {
	      cumlative = proba_k[k][i][j];
	      if(nbr < cumlative )
		{
		  store = j;
		  rep = k;
		  nbr = cumlative;
		}
	    }
	} 
      stop = store;
      cuts.push_back(store);
      markov_cuts.push_back(rep);
    } 
}


