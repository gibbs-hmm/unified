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

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "language.h"
//Function Declarations
void seg_in_two(language &current, vector<vector<double> > &proba_k,int seq);

void seg_recu(language &current, vector<vector<double> > &proba_k,vector<int> &reajust,int seq);

double find_err(language &current,vector<vector<double> > &proba_k,int position, int cutpts,int seq);

int nbr_segments(language &current, vector<vector<double> > &proba_k,   vector<vector<double> > &segments, vector<int> reajust, double & evidence,int seq);

void search_details(language &current,vector<vector<double> > &proba_k,  vector<vector<double> > &segments, int random, int seq);

void search_cut_details(language &current,vector<vector<double> > proba_k, int seg, vector<int > &cuts, int seq);

//End Function Declarations

void segment(language &current, int random) {
  int i;
  int j;
  current.initialize();
  vector<vector<double> > proba_k;
  vector<vector<double> > segments;
  vector<int> reajustments;
  double evidence;
  int nbr_segs;
  cout << "There are " << current.get_nbr_seq() << " sequences in the file" << endl;

  for(i = 0; i < current.get_nbr_seq(); i++) {
    cout << "Analysing sequence named " << current.get_name(i) << endl;
    cout << "Seg Stop: " << current.get_seg_stop(i) << endl;
    cout << "Adjustment: " << exp(current.get_adjustment(i))<< endl;
    proba_k = vector<vector<double> >(current.get_seg_stop(i)+1,vector<double>(current.get_len(i),0.0));
 
    seg_in_two(current, proba_k, i);
    cout << "Full segment=" << proba_k[0][current.get_len(i)-1] << endl;
    seg_recu(current,proba_k, reajustments,i);
    nbr_segs = nbr_segments(current,proba_k,segments,reajustments,evidence, i);
    cout << "Segments = " << nbr_segs << endl;
    cout << "Evidence = " << evidence << endl;
    for(j = 0; j < proba_k.size(); j++)
      cout << "Full segment=" << proba_k[j][current.get_len(i)-1] << endl;
    current.write_seg(segments, nbr_segs, i);
    current.write_res(evidence, nbr_segs, i);
    search_details(current,proba_k,segments, random,i);
    current.write_details(random, i);   
    // clean up for next run.
    proba_k.clear();
    segments.clear();
    reajustments.clear();
  }
}


void seg_in_two(language &current, vector<vector<double> > &proba_k,int seq){
  int pos;
  double proba;
  int end;
  end = current.get_len(seq);
  current.reset_probs(seq);
  pos = current.get_markov(seq);
  proba = exp(current.calculation(pos,seq));
  proba_k[0][pos] = proba;
  pos++;
  for(pos; pos < end; pos++)
    {
      proba *= exp(current.recursive_calculation(pos,seq));
      proba_k[0][pos] = proba;
    }
}

void seg_recu(language &current, vector<vector<double> > &proba_k, vector<int> &reajust,int seq)
{
  int cutpts;
  bool flag = 0;
  int stop = current.get_len(seq);
  for(cutpts = 1; cutpts <= current.get_seg_stop(seq); cutpts++) {
    if(proba_k[cutpts-1][stop-1] >= current.get_limit(seq))
      {
	flag = 1;
	reajust.push_back(cutpts);
      }
    int pos;
    if(cutpts > current.get_markov(seq))
      pos = cutpts;
    else
      pos = current.get_markov(seq);
    for (; pos < stop; pos++)
      {
	if(flag)
	  proba_k[cutpts][pos] = find_err(current,proba_k,pos,cutpts-1,seq)/current.get_limit(seq);
	else
	  proba_k[cutpts][pos] = find_err(current,proba_k,pos,cutpts-1,seq);
      }
    flag = 0;
  }
}

double find_err(language &current,vector<vector<double> > &proba_k,int position, int cutpts,int seq)
{
  current.reset_probs(seq);
  double ret_val = 0;
  double temp;
  int pos = position;
  temp = exp(current.calculation(pos,seq));
  ret_val = temp * proba_k[cutpts][pos-1];
  for(pos = position -1; pos > cutpts; pos--)
    {
      temp *= exp(current.recursive_calculation(pos,seq));
      ret_val += temp * proba_k[cutpts][pos-1];
    }
  return ret_val;
}

int nbr_segments(language &current, vector<vector<double> > &proba_k,   vector<vector<double> > &segments, vector<int> reajust, double & evidence,int seq) {
  long double total;
  double temp, temp2,pres;
  double *table;
  int length, seg, end, ch,kk;
  vector<int>::const_iterator p;
  ch = -2;

  length = current.get_len(seq) - current.get_markov(seq);
  table = new double[current.get_seg_stop(seq)+1];
  pres = 0;
  total = evidence = 0;
  evidence = log(0.0);
  p = reajust.begin();
  end = current.get_len(seq);
  segments = vector<vector<double> >( current.get_seg_stop(seq)+1,vector<double>(4,0.0));
  for(seg = 0; seg <= current.get_seg_stop(seq); seg++) {
    for(int i = 0; i < 4; i++)
    if(p != reajust.end() && seg == *p)
      {
	pres += log(current.get_limit(seq));
	p++;
      }
    segments[seg][0] = log(proba_k[seg][end-1]) + pres;
    if(current.get_seg_min(seq) != 1)
      {
	kk = length - (current.get_seg_min(seq) * (seg+1)) + ch;
	if(seg>2)
	  ch++;
	else
	  ch+=2;
      }
    else
      kk = length;
    temp = lgamma(kk) - (lgamma(kk-seg) + lgamma(seg+1));
    temp2 = lgamma(end) - (lgamma(end - seg) + lgamma(seg+1));
    segments[seg][1] = temp;
    segments[seg][2] = temp2;
    table[seg] = log(proba_k[seg][end-1])-temp2;
    table[seg] += pres;
    if(seg > 0 && evidence < table[seg])
      evidence = table[seg];
    total+=exp(table[seg]);
  }
  evidence = exp(table[0])/(exp(table[0]) + exp(evidence));
  for(seg = 0; seg <= current.get_seg_stop(seq); seg++)
    {
      table[seg] = ((long double)exp(table[seg]))/total;
      segments[seg][3] = table[seg];
    }
  ch = end/current.get_seg_max(seq);
  if(current.get_len(seq) % current.get_seg_max(seq))
    ch ++;
  else if(current.get_len(seq) == current.get_seg_max(seq))
    ch = 0;
  else
    ch--;
  for(temp = 0.0, seg = ch; seg <= current.get_seg_stop(seq); seg++)
    if(table[seg] >= temp)
      {
	ch = seg + 1;
	temp = table[seg];
      }
  delete []table;
  //Save the seg file.
  return ch;
}

void search_details(language &current,vector<vector<double> > &proba_k,  vector<vector<double> > &segments, int random, int seq)
{

  vector<double> cut_probs;
  vector<int> cuts;
  double cumlative = 0.0;
  double nbr;
  int stop;
  int rep;
  int start = 0;
  int segment;
  int i,j;
  cut_probs = vector<double>(current.get_seg_stop(seq));
  for(i = 0; i < (current.get_seg_stop(seq));i++)
    {
      cut_probs[i]=segments[i+start][3];
      cumlative += cut_probs[i];
    }
  for(i = 0; i < (current.get_seg_stop(seq));i++)
    {
      cut_probs[i] = cut_probs[i] / cumlative;
    }
  current.store_prob_details(0, current.get_len(seq)-1,seq);
  srand48(time(0));
  for(i = 0; i < random; i++) {
    nbr = drand48();
    rep = 1;
    cumlative = 0.0;
    for(j = 0; j < current.get_seg_stop(seq); j++)
      {
	if(nbr > cumlative && nbr < cumlative + cut_probs[j]) {
	  rep = j;
	  j = current.get_seg_stop(seq);
	}
	cumlative += cut_probs[j];
      }
    search_cut_details(current,proba_k,rep,cuts,seq);
    stop = current.get_len(seq)-1;
    for(vector<int>::const_iterator p = cuts.begin(); p != cuts.end();p++)
      {
	current.store_prob_details(*p, stop,seq);
	stop = *p;
      }
    cuts.clear();
    current.store_prob_details(0, stop,seq);
  }
  
}

void search_cut_details(language &current,vector<vector<double> > proba_k, int seg, vector<int > &cuts, int seq)
{
  int i,j;
  double cumlative=0.0;
  double res = 0.0;
  double nbr;
  int stop = current.get_len(seq);
  if(seg == 1)
    seg--;
  else
    seg -=2;
  for(i = seg; i >=0;i--)
    {
      cumlative = 0.0;
      current.reset_probs(seq);
      res = exp(current.calculation(stop-1, seq));
      proba_k[i][stop-1] *= exp(res);
      cumlative = proba_k[i][stop-1];
      for(j = stop-2; j >= current.get_markov(seq); j--)
	{
	  res += exp(current.recursive_calculation(j,seq));
	  proba_k[i][j] *= exp(res);
	  cumlative += proba_k[i][j];
	}

      // randomly choose a spot to cut.
      for(j = current.get_markov(seq); j < stop; j++)
	{
	  proba_k[i][j] /= cumlative;
	}
      nbr = drand48();
      cumlative = 0.0;
      for(j = current.get_markov(seq); j < stop; j++)
	{
	  if(nbr > cumlative && nbr <= cumlative + proba_k[i][j])
	    {
	      cuts.push_back(j);
	      stop = j;
	    }
	  cumlative += proba_k[i][j];
	}
    }  
}
