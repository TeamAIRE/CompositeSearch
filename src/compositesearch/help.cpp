/*

        Written by Jananan PATHMANATHAN, 2014-2018
        
        This file is part of CompositeSearch.
        
        CompositeSearch is shared under Creative commons licence: 
        
        Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
        
        See https://creativecommons.org/licenses/by-nc-sa/4.0/

*/

#include <iostream>
#include "help.h"
using namespace std;

void help()
{
	cout << endl;
	cout << "----- DESCRIPTION -----" << endl;
	cout << "CompositeSearch is a memory-efficient and fast software for the detection of " << endl;
	cout << "composite genes and composite gene famalis in very large similarity networks, "<< endl;
	cout << "e.g. with several millions of nodes and hundreds of millions edges." << endl;
	cout << endl;
	cout <<	"Written by Jananan PATHMANATHAN, 2014-2018" << endl;
	cout << "CompositeSearch is shared under Creative commons licence:" << endl;
	cout << "Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)" << endl;
	cout << "See https://creativecommons.org/licenses/by-nc-sa/4.0/" << endl << endl;
	cout << "If you use this program for a scientific paper, please cite the following paper: " << endl;
	cout << "Intorducing CompositeSearch, a generalized network approach for composite gene families " << endl;
	cout << "detection: an application to the analysis of environmental gene remodeling" << endl;
	cout << "JS. Pathmanathan, P. Lopez, F-J. Lapointe and E. Bapteste" << endl;
	cout << "(Article Submitted)" << endl << endl;

	cout << "----- USAGE -----" << endl << endl;

	cout << "compositeSearch -i blastp.out.cleanNetwork -n blastp.out.cleanNetwork.genes -m [fastcomposites or composites] " << endl;
	cout << endl;
	cout << "-- optionnal --" << endl << endl;
	cout << "-e evalue_limit " << endl;
        cout << "-p pident_limit " << endl;
        cout << "-c min_cov " << endl;
        cout << "-l max_overlap " << endl;
        cout << "-t n_cpu " << endl;
        cout << "-g specific_gene_list " << endl;
        cout << "-f gene_families " << endl;
        cout << "-x composite_family_min_size " << endl;
        cout << "-y component_family_min_size " << endl;
        cout << "-v verbose " << endl;
        cout << "-help" << endl << endl;
        
	cout << "--- MANDATORY ARGUMENTS ---" << endl << endl;

        cout << "-i : Input file name. Should be cleaned BLAST output file from cleanBlastp (file.cleanNetwork)." << endl << endl;
	
	cout << "-n : Input file name. Should be gene list of cleaned BLAST from cleanBlastp (file.cleanNetwork.genes)." << endl << endl;

	cout << "-m : Composites detection algorithms." << endl;
	cout << endl;
	cout << "	- fastcomposites : Detection of potential composite genes." << endl;
	cout << "	                   It will only give you a list of potential " << endl;
	cout << "	                   composite genes." << endl;
	cout << endl;
	cout << "	- composites   	 : Detection of potential composite gene families." << endl;
	cout << " 	               	   It will give you a list of potential composite " << endl;
	cout << "	                   and component gene families." << endl << endl;

	cout << "--- OPTIONAL ARGUMENTS ---" << endl << endl;

        cout << "-e : Evalue limit for components detection. [default : 10]" << endl << endl;

        cout << "-p : Minimum percentage of identity used to filter the network. [default : 30]" << endl << endl;

        cout << "-c : Minimum percentage of mutual coverage between two sequences in order to " << endl;
	cout << "     be part of the gene family. This option is used only with the algorithm " << endl;
	cout << "     composites. [default : 80]" << endl << endl;

        cout << "-l : Maximum overlap authorized between two components aligned on a potential " << endl;
	cout << "     composite gene. [default : 20]" << endl << endl;

        cout << "-t : Number of threads (CPUs) to use" << endl << endl;

	cout << "-g : File with a list of specific genes to analyze (one gene per line)." << endl << endl;

	cout << "-f : File with a list of gene and it family ID (if you do not want to use CompositeSearch family detection algorithm)." << endl;
	cout << "     Genes should be sorted by their family ID." << endl << endl;

	cout << "-x : Minimum size of composite gene famimlies [default: 1]" << endl << endl;

	cout << "-x : Minimum size of component gene famimlies [default: 1]" << endl << endl;

	cout << "-v : verbose [0 or 1, default : 1]" << endl << endl;
}
