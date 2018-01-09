/*

        Written by Jananan PATHMANATHAN, 2014-2018
        
        This file is part of CompositeSearch.
        
        CompositeSearch is shared under Creative commons licence: 
        
        Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
        
        See https://creativecommons.org/licenses/by-nc-sa/4.0/

*/

#include "functions.h"
#include "dfscc.h"
#include "getTime.h"
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
#include <iomanip>

using namespace std;

void compositeFamilies( map<unsigned int, list<unsigned long long int> >& familiesToCheck, map<unsigned long long int, geneInfo>& genes, map<unsigned long long int, familyInfo>& geneFamilies, unsigned int l, string timeInfo)
{
	string fileOut = to_string(l) + ".compositefamilies_" + timeInfo;
	ofstream output(fileOut.c_str());
	string fileInfo = to_string(l) + ".compositefamiliesinfo_" + timeInfo;
	ofstream outputInfo(fileInfo.c_str());
	// progress
	string progressInfo = to_string(l) + ".progress_" + timeInfo;
	unsigned long long int totalFamilies = familiesToCheck[l].size();
	unsigned long long int processedFamilies = 0;
	while(!familiesToCheck[l].empty())
	{
		// Family id
		unsigned long long int familyID = familiesToCheck[l].back();
		// Delete CCid from vector
		familiesToCheck[l].pop_back();
		if(geneFamilies[familyID].isComposite() == true)
		{
			// Nb composite
			unsigned long long int NbComposites = 0;
			// nOverlap
			float sumNoOverlap = 0.0;
			// Nb Domain
			unsigned int sumNbDomains = 0;
			vector<unsigned long long int> genesID = geneFamilies[familyID].getGenesID();
			// Component families ID list
			map<unsigned long long int, unsigned long long int> componentFamilies;
			while(!genesID.empty())
			{
				unsigned long long int gene = genesID.back();
				genesID.pop_back();
				if(genes[gene].isComposite() == true)
				{
					NbComposites++;
					vector<unsigned long long int> geneComponentFamiliesID = genes[gene].getComponentFamiliesID();
					while(!geneComponentFamiliesID.empty())
					{
						unsigned long long int componentFamilyID = geneComponentFamiliesID.back();
						geneComponentFamiliesID.pop_back();
						componentFamilies[componentFamilyID]++;
					}
					sumNbDomains+=genes[gene].getNbDomains();
					sumNoOverlap+=genes[gene].getNoOverlap();
				}
			}
			// composite
			output << ">CF" << familyID << endl;
			map<unsigned long long int, unsigned long long int>::iterator it;
			unsigned long long int singleNodeFamilies = 0;
			for(it = componentFamilies.begin(); it != componentFamilies.end(); ++it)
			{
				//output << it->first << "\t" << setprecision (2) << fixed << geneFamilies[it->first].getConnectivity() << "\t" << geneFamilies[it->first].getSize() << "\t" << it->second << "\t" << setprecision (2) << fixed << (float)it->second/(float)geneFamilies[familyID].getSize() << endl;
				output << "F" << it->first << "\t" << setprecision (2) << fixed << ((float)it->second/(float)geneFamilies[familyID].getSize())*100.00 << endl;
				if(geneFamilies[it->first].getSize() == 1) singleNodeFamilies++;
			}
			output << endl;
			// info
			outputInfo << "F" << familyID << "\t" << geneFamilies[familyID].getSize() << "\t" << NbComposites << "\t" << setprecision (2) << fixed << (float)NbComposites/(float)geneFamilies[familyID].getSize() << "\t" << componentFamilies.size() << "\t" << singleNodeFamilies << "\t" <<  setprecision (2) << fixed << geneFamilies[familyID].getConnectivity() << "\t" << round((float)sumNbDomains/(float)NbComposites)  << "\t" << setprecision (2) << fixed << sumNoOverlap/(float)NbComposites << endl;
		}
		processedFamilies++;
		unsigned int cpuP = ((float)processedFamilies/(float)totalFamilies)*100.0;
		ofstream progress(progressInfo.c_str());
		progress << "CompositeFamilies - CPU\t" << l << "\t" << cpuP << "%" << endl;
	}
}

void runCompositeFamilies( map<unsigned long long int, geneInfo>& genes, map<unsigned long long int, familyInfo>& geneFamilies, unsigned int nCpu, string timeInfo)
{
	//Map used to split nodes list into sublist
	map<unsigned int, list<unsigned long long int> > subList;
	map<unsigned int, list<unsigned long long int> >& refSubList = subList;
	unsigned int k = 0;
	//Split nodes to distribute
	map<unsigned long long int, familyInfo>::iterator n;
	for( n = geneFamilies.begin(); n != geneFamilies.end(); n++)
	{
		if( k == nCpu )
		{
			k = 0;
		}
		subList[k].push_back(n->first);
		k++;
	}
	// MultiThreading
	vector<thread> threads;
	for(unsigned int i = 0; i < nCpu; ++i)
	{
		threads.push_back(thread(compositeFamilies, std::ref(refSubList), std::ref(genes), std::ref(geneFamilies), i, timeInfo));
	}
	for(auto& thread : threads)
	{
		thread.join();
	}
	system(("rm *.progress_" + timeInfo).c_str());
}
//end
