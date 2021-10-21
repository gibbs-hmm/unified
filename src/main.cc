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
/* Email:  gibbsamp@wadsworth.org                               */
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

/* main.cc
 *
 * last modified: 04-10-06 mjp
 *
 * 03-14-06
 *    - added compile date string to code. it comes in via the make file.
 *
 * 03-13-06
 *    - moved call to srand48() to this module and incorporated use of -s arg
 *	into the call. also added output of starting seed val to cout.
 *
 * 03-10-06
 *    - creating "reconciled.1"
 *	  - added compiled date output.
 *	  - added --seed | -s <long> arg to specify starting random seed value.
 *
 * 03-06-06
 *    - copied this file from src.ivan/ to src.reconciled. the text at the
 *	explanation of --random was in the ivan version but not the bill 
 *	version.
 *
 * 01-24-06:
 *    - started looking into this program... 
 *
 *
 */

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <cstring>
#include <climits>

using namespace std;

#include "templated_fasta.h"
#include "aligned_template_fasta.h"
#include "dna_language.h"
#include "protein_language.h"
#include "binary_language.h"
#include "inclusive_markov_segment.h"

// this should work under cygwin and linux on ccmb; the #define comes in
// via the make file
#ifdef COMPILE_DATE_STRING
string compileDate = COMPILE_DATE_STRING;
#else
string compileDate = __DATE__;
#endif

// function prototypes:
void fatalError(char *);
static void outputUsage(void);

int main(int argc, char* argv[])
{
  int markov=0; 
  int alphabet=1;
  int inclusive = 0;
  int randomVal = 2000;
  char     *endStr, errBuf[128];
  long int startingSeedArg = 0L, startingSeedVal;
  string fileName;

  string verStr = "reconciled.1";

  vector<string> options;

  {
    char   **p;
    time_t   t = time(NULL);

    cout << "unified " << verStr << endl << endl << "compiled on: " << 
      compileDate << endl << endl << "run: " << ctime(&t) << endl;

    cout << "invoked with:\n  ";

    for ( p = argv; *p != NULL; p++ ) 
      cout << *p << " ";
    cout << endl << endl;

    cout << "For options: unifiedcpp --help" << endl;

    cout << endl << endl;
    cout << flush;   // BT 100807
  }

  // cout << "unified reconciled.1 03/13/06" << endl << endl;

  for(int i = 1; i < argc; i++)
    {
      if(strcmp(argv[i],"--markov") == 0 || strcmp(argv[i],"-m") == 0)
	{
	  markov = atoi(argv[i+1]);
	  i++; // skip option value from decision making
	}
      else if(strcmp(argv[i],"--alphabet") == 0 || strcmp(argv[i],"-a") == 0)
	{
	  alphabet = atoi(argv[i+1]);
	  i++; // skip option value from decision making
	}
      else if(strcmp(argv[i],"--random") == 0 || strcmp(argv[i],"-r") == 0)
	{
	  randomVal = atoi(argv[i+1]);
	  i++; // skip option value from decision making
	}
      else if(strcmp(argv[i],"--seed") == 0 || strcmp(argv[i],"-s") == 0)
	{
	  // starting seed is a long int - allows dec, octal, hex
	  startingSeedArg = strtol(argv[i+1],&endStr,0);
	  // do some error checking
	  if ((startingSeedArg == LONG_MIN) || (startingSeedArg == LONG_MAX)){
	    // do fatal error
	    fatalError("-seed argument out of range");
	  }
	  else if ( strcmp(endStr,"\0") != 0 ) {
	    // do fatal error
	    fatalError("invalid -seed argument");
	  }
	  // else startingSeedArg is set to a valid seed value.
	  i++; // skip option value from decision making
	}
      else if(strcmp(argv[i],"--options") == 0 || strcmp(argv[i],"-o") == 0)
	{
	  options.push_back(string(argv[i+1]));
	  i++; // skip option value from decision making
	}
      else if(strcmp(argv[i],"--inclusive") == 0 || strcmp(argv[i],"-i") == 0)
	{
	  inclusive = 1;
	  options.push_back(string("-i"));
	}
      else if(strcmp(argv[i],"--prior") == 0 || strcmp(argv[i],"-p") == 0)
	{
	  options.push_back(string("p")+string(argv[i+1]));
	}
      else if(strcmp(argv[i],"--help") == 0 || strcmp(argv[i],"-h") == 0)
	{
	  outputUsage();
	  exit(EXIT_FAILURE);
	}
      else if(strncmp(argv[i],"-",1) != 0) {
	// Not an option
	fileName = argv[i];
      }
      else {
	// unrecognized option; fatal error
	sprintf(errBuf,"unrecognized option: %s\n",argv[i]);
	fatalError(errBuf);
      }
    }

  // 03-13-06 mjp
  // moved call to srand48() here (from search_details(), which is in
  // inclusive_markov_segment.h) in order to support -seed arg. also, previous
  // to this move, srand48() was being called before processing every seq in
  // the specified fasta file. that was just plain wrong and (i believe) very
  // undesireable when the user specifies a starting seed. could have put 
  // srand48() call at beginning of segment(), but it would involve passing
  // the seed val as an arg, and there are many calls to segment(). easiest
  // to just do it here.
  if ( (startingSeedVal = startingSeedArg) == 0L )
    startingSeedVal = time(NULL);

  srand48(startingSeedVal);

  // output starting seed value so it can be used again...
  cout << "Starting seed value = " << startingSeedVal << endl;


  if(inclusive) {
    long double temp=0.0;
    if(alphabet == 1){
      cout << "Markov value = " << markov << endl;
      cout << "Language Type is DNA Fasta type. " << endl;

      vector<templated_fasta<dna_language> > current;
      for(int i = markov; i >=0; i--) {
	current.push_back(templated_fasta<dna_language>(fileName,i,options));
       }
      segment(current, randomVal, temp);
    }
    else if(alphabet == 2) {
      cout << "Markov value = " << markov << endl;
      cout << "Language Type is Protein Fasta type. " << endl;
      vector<templated_fasta<protein_language> > current;
      for(int i = markov; i >=0; i--) {
	current.push_back(templated_fasta<protein_language>(fileName,markov,options));
      }
      segment(current, randomVal, temp);
    }
    else if(alphabet == 3) {
      cout << "Markov value = " << markov << endl;
      cout << "Language Type is Binary Fasta type. " << endl;
      vector<templated_fasta<binary_language> > current;
      for(int i = markov; i >=0; i--) {
	current.push_back(templated_fasta<binary_language>(fileName,markov,options));
      }
      segment(current, randomVal, temp);
    }
    else if(alphabet == 4) {
      cout << "Markov value = " << markov << endl;
      cout << "Language Type is Aligned DNA Fasta type. " << endl;
      vector<aligned_templated_fasta<dna_language> > current;
      for(int i = markov; i >=0; i--) {
	current.push_back(aligned_templated_fasta<dna_language>(fileName,markov,options));
      }      
      segment(current, randomVal, temp);
    }
    else if(alphabet == 5) {
      cout << "Markov value = " << markov << endl;
      cout << "Language Type is Aligned Protein Fasta type. " << endl;
      vector<aligned_templated_fasta<protein_language> > current;
      for(int i = markov; i >=0; i--) {
	current.push_back(aligned_templated_fasta<protein_language>(fileName,markov,options));
      }      
      segment(current, randomVal, temp);
    }
    else if(alphabet == 6) {
      cout << "Markov value = " << markov << endl;
      cout << "Language Type is Aligned Binary Fasta type. " << endl;
      vector<aligned_templated_fasta<binary_language> > current;
      for(int i = markov; i >=0; i--) {
	current.push_back(aligned_templated_fasta<binary_language>(fileName,markov,options));
      }
      segment(current, randomVal, temp);
    }
  }
  else{
    double temp = 0.0;
    if(alphabet == 1){
      cout << "Markov value = " << markov << endl;
      cout << "Language Type is DNA Fasta type. " << endl;

      vector<templated_fasta<dna_language> > current;
      current.push_back(templated_fasta<dna_language>(fileName,markov,options));
      segment(current, randomVal, temp);
    }
    else if(alphabet == 2) {
      cout << "Markov value = " << markov << endl;
      cout << "Language Type is Protein Fasta type. " << endl;
      vector<templated_fasta<protein_language> > current;
      current.push_back(templated_fasta<protein_language>(fileName,markov,options));
      segment(current, randomVal, temp);
    }
    else if(alphabet == 3) {
      cout << "Markov value = " << markov << endl;
      cout << "Language Type is Binary Fasta type. " << endl;
      vector<templated_fasta<binary_language> > current;
      current.push_back(templated_fasta<binary_language>(fileName,markov,options));
      segment(current, randomVal, temp);
    }
    else if(alphabet == 4) {
      cout << "Markov value = " << markov << endl;
      cout << "Language Type is Aligned DNA Fasta type. " << endl;
      vector<aligned_templated_fasta<dna_language> > current;
      current.push_back(aligned_templated_fasta<dna_language>(fileName,markov,options));
      segment(current, randomVal, temp);
    }
    else if(alphabet == 5) {
      cout << "Markov value = " << markov << endl;
      cout << "Language Type is Aligned Protein Fasta type. " << endl;
      vector<aligned_templated_fasta<protein_language> > current;
      current.push_back(aligned_templated_fasta<protein_language>(fileName,markov,options));
      segment(current, randomVal, temp);
    }
    else if(alphabet == 6) {
      cout << "Markov value = " << markov << endl;
      cout << "Language Type is Aligned Binary Fasta type. " << endl;
      vector<aligned_templated_fasta<binary_language> > current;
      current.push_back(aligned_templated_fasta<binary_language>(fileName,markov,options));
      segment(current, randomVal, temp);
    }
   }

  cout << endl << "all done. normal exit." << endl;

  return 0;
}


/**************************************************************************
 * fatalError()
 * 
 */

void fatalError(char *msgStr)
{
  cout << "Fatal Error: " << msgStr << endl << endl;

  outputUsage();

  exit(EXIT_FAILURE);

} // fatalError()


/**************************************************************************
 * outputUsage()
 * 
 */

static void outputUsage(void)
{
  cout << "Unified A Bayesian segmentation program" << endl;
  cout << "Usage: unified <options> fasta_filename" << endl;
  cout << endl;

  cout << "Output: The program produces three output files in the same directory as the FASTA file:" << endl;
  cout << "fasta_filename_info-det - a space separated table containing the columns Sequence Sample_Count, A_prob, C_prob, G_prob, T_prob, Position, Samples, Change_Point_probability." << endl;
  cout << "fasta_filename_info-inclusive -  a space separated table containing sequence numbers and positions." << endl;
  cout << "fasta_file_info-name - a list of the FASTA sequence identifiers." << endl;
  cout << endl;
  cout << "fasta_filename_info-det can be used with the -B option of the Gibbs sampler." << endl;
  cout << endl;
  
  cout << "Options: " << endl << "--alphabet num or -a num:Alphabet of sequence. default is 1. 1:DNA, 2:Protein, 3:Binary, 4:Aligned DNA, 5:Aligned Protein, 6:Aligned Binary"<< endl;

  cout << "  --inclusive or -i: Use inclusive markov models. Evaluate all markov levels from 0 to markov, as set by the -m option." << endl;

  cout << "  --markov num or -m num:Set markov level or maximum markov level(for inclusive)" << endl;

  cout << "  --prior a|p or -p a|p: Set prior on k to be either the alternative P(k) uniform for all k >= 0(a) or to the Lawrence and Liu prior where P(k=0) = .5 and P(k) is uniform over for all k > 0 over the remaining .5 probability(p). The default is the alternative prior." << endl;

  cout << "  --random num or -r num: Number of random samples to take. Making the value negative will require the algorithm to include the best solution as the first sampled solution." << endl;

  cout << "  --seed num or -s num: seed for random number generator. specified as long int; dec, octal (0), or hex(0x)  numbers are valid." << endl;

  cout << "  --options <option> or -o <option>: Set a secondary option" << endl;
  cout << "  Secondary Options:" << endl;
  cout << "  ss<num> : Set maximum number of cutpoints to num"<< endl;
  cout << "  sl<num> : Set minimum segment length to num, default 1" << endl;
  cout << "  sg<num> : Set maximum segment length to num, default sequence length." << endl;
  cout << "  pc<num> : Set psuedocount weight value to num. Default 1" << endl;
  cout << "  ad<num> : Set adjustment to num." << endl;
  cout << "  al<num> : Set the alignment level for aligned sequence alphabets." << endl;

} // outputUsage()
