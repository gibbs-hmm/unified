# Unified
Unified - The Baysian Segmentation Program.

This program implements the Bayesain segmentation algorithm described in Liu and Lawrence (1999) https://academic.oup.com/bioinformatics/article/15/1/38/218372.

The main purpose of the proram is to generate a position specific background model for the Gibbs sampler found at https://github.com/gibbs-hmm/Gibbs-Motif-Sampler.

Installation:
A linux binary is provided in the bin directory.<br/>
To compile from source, clone this repo, cd to the src directory and run make.<br/>
A new binary version will be compiled to the ../bin/ direrctory.<br/>

Usage<br/>
Usage: unified < options > fasta_filename<br/>
Output: The program produces three output files in the same directory as the FASTA file:<br/>
  * fasta_filename_info-det - a space separated table containing the columns Sequence Sample_Count, A_prob, C_prob, G_prob, T_prob, Position, Samples, Change_Point_probability.<br/>
  * fasta_filename_info-inclusive -  a space separated table containing sequence numbers and positions.<br/>
  * fasta_file_info-name - a list of the FASTA sequence identifiers.<br/>
<br />
fasta_filename_info-det can be used with the -B option of the Gibbs sampler.<br />
<br />
Options:
<pre><code>
--alphabet num or -a num:Alphabet of sequence. 
	 default is 1. 1:DNA, 2:Protein, 3:Binary, 4:Aligned DNA, 5:Aligned Protein, 6:Aligned Binary
--inclusive or -i: Use inclusive markov models. Evaluate all markov levels from 0 to markov, as set by the -m option.
	--markov num or -m num:Set markov level or maximum markov level(for inclusive)
	--prior a|p or -p a|p: Set prior on k to be either the alternative P(k) uniform for all k >= 0(a) or to the Lawrence and Liu prior where P(k=0) = .5 and P(k) is uniform over for all k > 0 over the remaining .5 probability(p). The default is the alternative prior.
	--random num or -r num: Number of random samples to take. Making the value negative will require the algorithm to include the best solution as the first sampled solution.
	--seed num or -s num: seed for random number generator. specified as long int; dec, octal (0), or hex(0x)  numbers are valid.
	--options < option > or -o < option >: Set a secondary option
</code></pre>
	
<pre><code>
Secondary Options:
	ss< num > : Set maximum number of cutpoints to num
	sl< num > : Set minimum segment length to num, default 1
	sg< num > : Set maximum segment length to num, default sequence length.
	pc< num > : Set psuedocount weight value to num. Default 1
	ad< num > : Set adjustment to num.
	al< num > : Set the alignment level for aligned sequence alphabets.
</code></pre>
<br>
An R function, plot_seq_prob.R, is provided to plot the contents of the fasta_filename_info-det file.
<br \>
