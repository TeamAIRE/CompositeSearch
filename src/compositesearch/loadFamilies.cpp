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
#include <vector>
#include <iomanip>
#include <iterator>


using namespace std;

// Composites Families detection algo

void loadFamilies( string geneFamiliesFile, map<unsigned long long int, geneInfo>& genes, map<unsigned long long int, familyInfo>& geneFamilies, map<unsigned long long int, bool>& isVisited, string outputDir)
{
	// family original ID
	map<unsigned long long int,string> famOriginalID;
	// family info file
	string fileOutCCinfo = outputDir + "/family.info";
	ofstream outputCCinfo(fileOutCCinfo.c_str());
	outputCCinfo << "#familyID\toriginalID\tnodes\tedges\tconnectivity" << endl;
	// nodes file
	string fileOutCCnodes = outputDir + "/family.nodes";
	ofstream outputCCnodes(fileOutCCnodes.c_str());
	// edges file
	string fileOutCCedges = outputDir + "/family.edges";
	ofstream outputCCedges(fileOutCCedges.c_str());
	// Read gene families file
        ifstream geneFam(geneFamiliesFile.c_str());
	string gF;
	string CID;
	// Number of connected components and new family ID
	unsigned long long int nCC = 0;
	while(getline(geneFam,gF))
	{
		istringstream gInfo(gF);
		vector<string> gFid;
		copy(istream_iterator<string>(gInfo),istream_iterator<string>(),back_inserter<vector<string> >(gFid));
		unsigned long long int gID;
		string ID;
		gID = atoll(gFid[0].c_str());
		ID = gFid[1].c_str();
		//cout << gID << "\t" << ID << endl;
		if(nCC == 0)
		{
			nCC++;
			geneFamilies[nCC]=familyInfo();
			CID = ID;
			famOriginalID[nCC]=ID;
			//cout << "First CC\t" << nCC << endl;
		}
		else if(CID != ID)
		{
			nCC++;
			geneFamilies[nCC]=familyInfo();
			CID = ID;
			famOriginalID[nCC]=ID;
			//cout << "Next CC\t" << nCC << endl;
		}
		geneFamilies[nCC].insertGeneID(gID);
                genes[gID].setFamilyID(nCC);
		//cout << gID << "\t" << nCC << endl;
	}
	//
	nCC = 0;
        // Compute connected components
	map<unsigned long long int, geneInfo>::iterator it;
	for( it = genes.begin(); it != genes.end(); it++)
	{
		// Check if node has been visited, if true go to next one
		if(isVisited[it->first] == true ) continue;
		// Start
		// cc family number
		nCC++;
		// Nb edges for connectivity
		unsigned long long int NbEdges = 0;
		unsigned long long int NbNodes = 0;
		// Temporary list containing nodes belonging to CC
		set<unsigned long long int> tmpNodes;
		set<unsigned long long int>::iterator val;
		// Insert CC starting node
		tmpNodes.insert(it->first);
		// family ID
		unsigned long long int fID = genes[it->first].getFamilyID();
		// output nodes
                outputCCnodes << ">F" << fID << endl;
                // outptu edges
                outputCCedges << ">F" << fID << endl;
		// Get CC until tmpNodes is not empty
		while(!tmpNodes.empty())
		{
			// Get node
			val = tmpNodes.begin();
			unsigned long long int node = *val;
			// Delete node from tmpNodes
			tmpNodes.erase(val);
			// If not visited check node
			if( isVisited[node] == false )
			{
				// Change node status to visited
				isVisited[node] = true;
				NbNodes++;
				outputCCnodes << node << endl;
				// Current node's neighnors
				vector<unsigned long long int > currentNodeNeighbors = genes[node].getNeighbors();
				//for(i = 0; i < currentNodeNeighbors.size(); i++ )
				while(!currentNodeNeighbors.empty())
				{
					// Get neighbor
					unsigned long long int neighbor = currentNodeNeighbors.back();
					currentNodeNeighbors.pop_back();
					// Check neighbor's status
					if(isVisited[neighbor] == false)
					{
						if(genes[node].getFamilyID() == genes[neighbor].getFamilyID())
						{
							// Add neighbor to tmpNodes
							tmpNodes.insert(neighbor);
							// Count edges
							NbEdges++;
							if(node < neighbor)
							{
								outputCCedges << node << "\t" << neighbor << endl;
							}
							else
							{
								outputCCedges << neighbor << "\t" << node << endl;
							}
						}
					}
				}
			}
		}
		// Check CC integrity
		if( NbNodes < 3 )
		{
			//info
			geneFamilies[fID].setConnectivity(1.0);
			geneFamilies[fID].setComposite(false);
			outputCCinfo << fID << "\t" << famOriginalID[fID] << "\t" << NbNodes << "\t" << NbEdges << "\t1.00" << endl;
			geneFamilies[fID].sortGenes();
		}
		else
		{
			float connectivity = (2.0*(float)NbEdges)/((float)NbNodes*((float)NbNodes-1.0));
                        outputCCinfo << fID << "\t" << famOriginalID[fID] << "\t" << NbNodes << "\t" << NbEdges << "\t" << setprecision (2) << fixed << connectivity << endl;
                        //
                       	geneFamilies[fID].setConnectivity(connectivity);
			geneFamilies[fID].setComposite(false);
			geneFamilies[fID].sortGenes();
		}
	}
	// Clear isVisted map container
	isVisited.clear();
}

//end
