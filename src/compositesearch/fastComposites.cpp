/*

        Written by Jananan PATHMANATHAN, 2014-2018
        
        This file is part of CompositeSearch.
        
        CompositeSearch is shared under Creative commons licence: 
        
        Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
        
        See https://creativecommons.org/licenses/by-nc-sa/4.0/

*/

#include "functions.h"
#include <map>
#include <string>
#include <iostream>
#include <cstring>
#include <list>
#include <set>
#include <stdlib.h>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <thread>
#include <vector>

using namespace std;

void fastComposites( map<unsigned int, list<unsigned long long int> >& refSubNodes, map<unsigned long long int, geneInfo>& genes, map<pair<unsigned long long int, unsigned long long int>, edgeValues>& edges, unsigned short int MAX_OVERLAP, unsigned int l, string timeInfo )
{
	// temporary output file
	string fileOut = to_string(l) + ".fastcomposites_" + timeInfo;
	ofstream output(fileOut.c_str());
	// subnodes list
	// progress
	string progressInfo = to_string(l) + ".progress_" + timeInfo;
	unsigned long long int totalNodes = refSubNodes[l].size();
	unsigned long long int processedNodes = 0;
	while(!refSubNodes[l].empty())
	{
		// Composite candidate genes
		unsigned long long int composite = refSubNodes[l].front();
		refSubNodes[l].pop_front();
		// List of potential components
		vector<unsigned long long int> components = genes[composite].getNeighbors();
		// FusedGene -> 0 (not detected) or 1 (detected)
		unsigned int FusedGene = 0;
		// Start checking compistes
		while( components.size() > 1 && FusedGene != 1)
		{
			// Get neighbor (last)
			unsigned long long int component_1 = components.back();
			// Get the list of all neighbors of component_1
			vector<unsigned long long int> all_c1_neighbors = genes[component_1].getNeighbors();
			vector<unsigned long long int>& refAc1n = all_c1_neighbors;
			// Delete component_1 form the neighbors list
			components.pop_back();
			// Check if the selected gene is a composite
			vector<unsigned long long int>::iterator j;
			for( j = components.begin(); j != components.end(); ++j)
			{
				unsigned long long int component_2 = *j;
				// Check if component_1 and component_2 have an homology
				if(!FoundIn(component_2, refAc1n))
				{
					//cout << component_1 << "\t" << composite << "\t" << component_2 << endl;
					edgeValues edge_1 = getEdgeValues(edges, composite, component_1);
					edgeValues edge_2 = getEdgeValues(edges, composite, component_2);
					if( checkTriplet(edge_1, edge_2, MAX_OVERLAP) )
					{
						//cout << component_1 << "\t" << composite << "\t" << component_2 << endl;
						//Print fusedgene
						output << composite << endl;
						// Fusedgene found
						FusedGene = 1;
					}
				}
				// If fusedgene found (1) exit the for and while loop
				if( FusedGene == 1 ) break;
			}
		}
		processedNodes++;
		unsigned int cpuP = ((float)processedNodes/(float)totalNodes)*100.0;
		ofstream progress(progressInfo.c_str());
		progress << "FastComposites - CPU\t" << l << "\t" << cpuP << "%" << endl;
	}
}

void runFastComposites( map<unsigned long long int, geneInfo>& genes, map<pair<unsigned long long int, unsigned long long int>, edgeValues>& edges, vector<pair<unsigned long long int,unsigned long long int>>& nodesDegree, unsigned short int MAX_OVERLAP, unsigned int nCpu, string timeInfo)
{
	//Map used to split nodes list into sublist
	map<unsigned int, list<unsigned long long int> > subList;
	map<unsigned int, list<unsigned long long int> >& refSubList = subList;
	unsigned int k = 0;
	//Split nodes to distribute
	unsigned int n;
	for(n = 0; n < nodesDegree.size(); n++)
	{
		if( k == nCpu )
		{
			k = 0;
		}
		subList[k].push_back(nodesDegree[n].second);
		k++;
	}
        //cout << subList.size() << endl;
	// MultiThreading
	vector<thread> threads;
	for(unsigned int i = 0; i < nCpu; ++i)
	{
		threads.push_back(thread(fastComposites, std::ref(refSubList), std::ref(genes), std::ref(edges), MAX_OVERLAP, i, timeInfo));
	}
	for(auto& thread : threads)
	{
		thread.join();
	}
	system(("rm *.progress_" + timeInfo).c_str());
}
//end
