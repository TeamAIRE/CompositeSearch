/*

	Written by Jananan PATHMANATHAN, 2014-2018
	
	This file is part of CompositeSearch.
	
	CompositeSearch is shared under Creative commons licence: 
	
	Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
	
	See https://creativecommons.org/licenses/by-nc-sa/4.0/

*/

#include "parseHeader.h"
#include "functions.h"
#include "loadNetwork.h"
#include "dfscc.h"
#include "loadFamilies.h"
#include "composites.h"
#include "compositeFamilies.h"
#include "fastComposites.h"
#include "mergeFiles.h"
#include "getTime.h"
#include "help.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <list>
#include <iterator>
#include <unistd.h>
#include <stdlib.h>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <map>
#include <thread>
#include <cmath>
#include <ctime>
#include <set>
#include <vector>

using namespace std;

int main(int argc, char * argv[])
{
	// ALL TIME (start)
	time_t Astart = time(NULL);
	// *************** Default Values ***************
	// Blastp file header
	// Clean Blastp file header
	string header = "qseqid sseqid qstart qend qlen sstart send slen pident evalue";
	// Evalue threshold
	double evalueLimit = 10;
	// Percentage of identity threshold
	float pidentLimit = 30;
	// Minimum mutual coverage between two sequences
	unsigned short int minCov = 80;
	// Maximum overlap between two components
	unsigned short int MAX_OVERLAP = 20;
	// Number of CPU(s)
	unsigned int nCpu = 1;
	// composite family minimum size 
	unsigned int compositeFamMinSize = 1;
	// component family minimu size 
	unsigned int componentFamMinSize = 1;
	// Verbose mode
	unsigned short int verbose = 1;
	// Check cmd line 
	if(argc > 1 && string(argv[1]) == "-help")
	{
		help();
		exit(1);
	}
	if( argc == 1 || argc < 7 || argc > 21 )
	{
		cout << endl;
		cout << "Usage : " << argv[0] << " -i blastp.out.cleanNetwork -n blastp.out.cleanNetwork.genes -m [fastcomposites or composites] " << endl << endl;
		cout << "-- optional -- " << endl << endl;
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
		cout << "-- Default values --" << endl << endl;
		cout << "evalue_limit               : " << evalueLimit << endl;
		cout << "pident_limit               : " << pidentLimit << endl;
		cout << "min_cov                    : " << minCov << endl;
		cout << "max_overlap                : " << MAX_OVERLAP << endl;
		cout << "n_cpu                      : " << nCpu << endl;
		cout << "composite_family_min_size  : " << compositeFamMinSize << endl;
		cout << "component_family_min_size  : " << componentFamMinSize << endl;
		cout << "verbose                    : " << verbose << endl << endl;
		exit(1);
	}
	// Time
	time_t start, end;
	int c;
	// Blastp.out file
	string fileIn;
	// Blastp genes list
	string blastpGenesFile;
	// FusedTriplets output file
	string fileOut;
	// File containing specific genes to analyse
	string geneList;
	// File containing pre-computed gene families
	string geneFamiliesFile;
	// triplets, fusedtriplets or fusedgenes
	string found;
	// select algorithm
	unsigned short int methodSelected = 0;
	// blastp genes list given
	unsigned short int blastpGenes = 0;
	// Get arguments
	while( (c = getopt (argc, argv, "i:n:m:e:p:c:l:t:g:f:x:y:v:")) != -1)
	{
		switch(c)
		{
			case 'i':
				fileIn = optarg;
				break;
			case 'n':
				blastpGenesFile = optarg;
				blastpGenes = 1;
				break;
			case 'm':
				found = optarg;
				methodSelected = 1;
				if( found != "fastcomposites" && found != "composites" )
				{
					cout << "Compute option [-c] " << found << " unknown." << endl;
					cout << "Please select one of those : fastcomposites or composites" << endl;
					exit(1);
				}
				break;
			case 'e':
				evalueLimit = atof(optarg);
				break;
			case 'p':
				pidentLimit = atof(optarg);
				break;
			case 'c':
				minCov = atoi(optarg);
				break;
			case 'l':
				MAX_OVERLAP = atoi(optarg);
				break;
			case 't':
				nCpu = atoi(optarg);
				break;
			case 'g':
				geneList = optarg;
				break;
			case 'f':
				geneFamiliesFile = optarg;
				break;
			case 'x':
				compositeFamMinSize = atoi(optarg);
				break;
			case 'y':
				componentFamMinSize = atoi(optarg);
				break;
			case 'v':
				verbose = atoi(optarg);
				break;
			default:
				cout << endl;
				cout << "Usage : " << argv[0] << " -i blastp.out.cleanNetwork -n blastp.out.cleanNetwork.genes -m [fastcomposites or composites] " << endl << endl;
				cout << "-- optionnal -- " << endl << endl;
				cout << "-e evalue_limit " << endl;
				cout << "-p pident_limit " << endl;
				cout << "-c min_cov " << endl;
				cout << "-l max_overlap " << endl;
				cout << "-t n_cpu " << endl;
				cout << "-t n_cpu " << endl;
				cout << "-g specific_gene_list " << endl;
				cout << "-f gene_families " << endl;
				cout << "-x composite_family_min_size " << endl;
				cout << "-y component_family_min_size " << endl;
				cout << "-v verbose " << endl;
				cout << "-help" << endl << endl;
				cout << "-- Default values --" << endl << endl;
				cout << "evalue_limit               : " << evalueLimit << endl;
				cout << "pident_limit               : " << pidentLimit << endl;
				cout << "min_cov                    : " << minCov << endl;
				cout << "max_overlap                : " << MAX_OVERLAP << endl;
				cout << "n_cpu                      : " << nCpu << endl;
				cout << "composite_family_min_size  : " << compositeFamMinSize << endl;
				cout << "component_family_min_size  : " << componentFamMinSize << endl;
				cout << "verbose                    : " << verbose << endl << endl;
				exit(1);
				break;
		}
	}
	if( methodSelected == 0 && blastpGenes == 0 )
	{
		cout << "Please select one of those method : fastcomposites or composites." << endl;
		cout << "The file containing the list of blastp genes is missing : -n blastp_genes.list" << endl;
		exit(1);
	}
	else if(methodSelected == 0 && blastpGenes == 1 )
	{
		cout << "Please select one of those method : fastcomposites or composites." << endl;
		exit(1);
	}
	else if( methodSelected == 1 && blastpGenes == 0 )
	{
		cout << "The file containing the list of blastp genes is missing : -n blastp_genes.list" << endl;
		exit(1);
	}
	// input file basename
	string basename = fileIn.substr(fileIn.find_last_of("\\/")+1);
	replace(basename.begin(), basename.end(), '.', '_');
	string folder = basename + "_" + found;
	// Output directory
	string folderTime = asctime(localtime(&Astart));
	replace(folderTime.begin(), folderTime.end(), ' ', '_');
	replace(folderTime.begin(), folderTime.end(), ':', '_');
	folderTime.pop_back();
	string outputDir = folder + "_" + folderTime;
	system(("mkdir " + outputDir).c_str());
        // Ouput file basename
        fileOut = outputDir + "/" + basename;
	// Print running command line
	if( verbose == 1 )
	{
		if(geneFamiliesFile == "")
		{
			cout << "---------- PARAMETERS ----------" << endl;
			cout << "Input         : " << fileIn << endl;
			cout << "Ouput         : " << outputDir << endl;
			cout << "Compute       : " << found << endl;
			cout << "E-value       : " << evalueLimit << endl;
			cout << "Pident        : " << pidentLimit << endl;
			cout << "MinCov        : " << minCov << endl;
			cout << "Overlap       : " << MAX_OVERLAP << endl;
			cout << "Nb CPUs       : " << nCpu << endl;
			cout << "Composite FMS : " << compositeFamMinSize << endl;
			cout << "Component FMS : " << componentFamMinSize << endl;
			cout << "Verbose       : " << verbose << endl;
			cout << "--------------------------------" << endl << endl;
		}
		else
		{
			cout << "---------- PARAMETERS ----------" << endl;
                        cout << "Input         : " << fileIn << endl;
                        cout << "Ouput         : " << outputDir << endl;
                        cout << "Compute       : " << found << endl;
                        cout << "E-value       : " << evalueLimit << endl;
                        cout << "Pident        : " << pidentLimit << endl;
                        cout << "GeneFamilies  : " << geneFamiliesFile << endl;
                        cout << "Overlap       : " << MAX_OVERLAP << endl;
                        cout << "Nb CPUs       : " << nCpu << endl;
                        cout << "Composite FMS : " << compositeFamMinSize << endl;
                        cout << "Component FMS : " << componentFamMinSize << endl;
                        cout << "Verbose       : " << verbose << endl;
                        cout << "--------------------------------" << endl << endl;
		}
	}
	// Create all maps and their reference
	// edges[query,subject]:qstart, qenid, sstart, send, evalue and pident
	map<pair<unsigned long long int, unsigned long long int>, edgeValues> edges;
	map<pair<unsigned long long int, unsigned long long int>, edgeValues>& refEdges = edges;
	// all_neighbors[node]: list of all neighbors
	map<unsigned long long int, geneInfo > genes;
	map<unsigned long long int, geneInfo >& refGenes = genes;
	// node's degree
	vector<pair<unsigned long long int,unsigned long long int>> nodesDegree;
	vector<pair<unsigned long long int,unsigned long long int>>& refND = nodesDegree;
	// Check header
	map<string, unsigned short int> positionsList;
	map<string, unsigned short int>& refPositionsList = positionsList;
	// family dectection isVisted
	map<unsigned long long int, bool> isVisited;
	map<unsigned long long int, bool>& refIV = isVisited;
	// for loading network
	map<unsigned long long int, bool> geneIsVisited;
	map<unsigned long long int, bool>& refGIV = geneIsVisited;
	// gene families info
	map<unsigned long long int, familyInfo> geneFamilies;
	map<unsigned long long int, familyInfo>& refGeneFamilies = geneFamilies;
	
	// Step 1 : get header position
	positionsList = getPositions(header);
	
	// Step 2 : Initialize isVisited and geneIsVisited
	start = time(NULL);
	if( verbose == 1 )
	{
		cout << "----- GENE LIST INITIALIZATION -----" << endl << endl;
		cout << "START - " << asctime(localtime(&start)) << endl;
	}
	initializeGeneList(blastpGenesFile,refIV,refGIV);
	end = time(NULL);
	if( verbose == 1 )
	{
		cout << "END   - " << asctime(localtime(&end)) << endl;
	}
	
	// Step 3 : Load network
	start = time(NULL);
	if( verbose == 1 )
	{
		cout << "----- NETWORK LOADING -----" << endl << endl;
		cout << "START - " << asctime(localtime(&start)) << endl;
	}
	runLoadNetwork(fileIn, geneList, pidentLimit, evalueLimit, refPositionsList, refEdges, refGenes, refND, refGIV, nCpu, folderTime);
	end = time(NULL);
	if( verbose == 1 )
	{
		cout << "END   - " << asctime(localtime(&end)) << endl;
	}
	duration( start, end, verbose );
	cout << endl;
	if( verbose == 1 )
	{
		cout << "---------- NETWORK INFO ----------" << endl << endl;
		cout << "Number of nodes :\t" << genes.size() << endl << endl;
		cout << "Number of edges :\t" << edges.size() << endl << endl;
		cout << "----------------------------------" << endl << endl;
	}
	// Step 4 : compute composite genes [fastcomposites] or composite gene families [composites]
	if( found == "composites" )
	{
		start = time(NULL);
		if( geneFamiliesFile == "")
		{
			if( verbose == 1 )
			{
				cout << "----- GENE FAMILIES DETECTION -----" << endl << endl;
				cout << "START - " << asctime(localtime(&start)) << endl;
			}
			computeFamilies(refGenes, refEdges, refGeneFamilies, refIV, minCov, outputDir);
			end = time(NULL);
		}
		else
		{
			if( verbose == 1 )
                        {
                                cout << "----- LOADING GENE FAMILIES -----" << endl << endl;
                                cout << "START - " << asctime(localtime(&start)) << endl;
                        }
                        loadFamilies(geneFamiliesFile, refGenes, refGeneFamilies, refIV, outputDir);
                        end = time(NULL);
		}
		if( verbose == 1 )
		{
			cout << "END   - " << asctime(localtime(&end)) << endl;
		}
		duration( start, end, verbose );
		cout << endl;
		start = time(NULL);
		if( verbose == 1 )
		{
			cout << "----- COMPOSITES DETECTION -----" << endl << endl;
			cout << "START - " << asctime(localtime(&start)) << endl;
		}
		// Detect composites
		runComposites(refGenes, refEdges, refGeneFamilies, refND, MAX_OVERLAP, minCov, compositeFamMinSize, componentFamMinSize, nCpu, folderTime);
		// Compute composite Families
		runCompositeFamilies(refGenes, refGeneFamilies, nCpu, folderTime);
		end = time(NULL);
		// Merge all output files in one
		mergeFiles(nCpu, "composites", fileOut, folderTime);
		if( verbose == 1 )
		{
			cout << "END   - " << asctime(localtime(&end)) << endl;;
		}
		duration( start, end, verbose);
		cout << endl;
	}
	else if ( found == "fastcomposites" )
	{
		start = time(NULL);
		if( verbose == 1 )
		{
			cout << "----- FAST COMPOSITES DETECTION -----" << endl << endl;
			cout << "START - " << asctime(localtime(&start));
		}
		// Run FusedGenes in MultiThread mode
		runFastComposites( refGenes, refEdges, refND, MAX_OVERLAP, nCpu, folderTime);
		end = time(NULL);
		// Merge all output files in one
		mergeFiles( nCpu, "fastcomposites", fileOut, folderTime);
		if( verbose == 1 )
		{
			cout << "END   - " << asctime(localtime(&end));
		}
		duration( start, end, verbose);
		cout << endl;
	}
	// ALL TIME (end)
	time_t Aend = time(NULL);
	// REAL TIME
	if( verbose == 1 )
	{
		cout << "----- TOTAL COMPUTING TIME -----" << endl;
	}
	duration( Astart, Aend, verbose);
	cout << endl;
}
//END
