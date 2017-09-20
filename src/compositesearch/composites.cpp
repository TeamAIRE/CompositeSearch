/*

        Written by Jananan PATHMANATHAN, 2014-2016
        
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

void composites( map<unsigned int, list<unsigned long long int> >& refSubNodes, map<unsigned long long int, geneInfo >& genes, map<pair<unsigned long long int, unsigned long long int>, edgeValues>& edges, map<unsigned long long int, familyInfo >& geneFamilies, unsigned short int MAX_OVERLAP, unsigned short int minCov, unsigned int compositeFamMinSize, unsigned int componentFamMinSize, unsigned int l)
{
	// temporary output file
	string fileOut = to_string(l) + ".composites";
	ofstream output(fileOut.c_str());
	string fileInfo = to_string(l) + ".compositesinfo";
	ofstream info(fileInfo.c_str());
	// progress
	string progressInfo = to_string(l) + ".progress";
	unsigned long long int totalNodes = refSubNodes[l].size();
	unsigned long long int processedNodes = 0;
	// Check if genes are composites
	while(!refSubNodes[l].empty())
	{
		// Composite candidate genes
		unsigned long long int composite = refSubNodes[l].front();
		refSubNodes[l].pop_front();
		// Composite candidate genes neighbors
		vector<unsigned long long int> compositeNeighbors = genes[composite].getNeighbors();
		vector<unsigned long long int> components;
		vector<unsigned long long int>& refCP = components;
		if(geneFamilies[genes[composite].getFamilyID()].getSize() >= compositeFamMinSize)
		{
			while(!compositeNeighbors.empty())
			{
				// Get compositeNeighbor's last gene
				unsigned long long int neighbor = compositeNeighbors.back();
				// Remove the gene from compositeNieghbor
				compositeNeighbors.pop_back();
				// Get composite and neighbor hit evalue
				if(genes[composite].getFamilyID() != genes[neighbor].getFamilyID() && geneFamilies[genes[neighbor].getFamilyID()].getSize() >= componentFamMinSize)
				{
					components.push_back(neighbor);
				}
			}
		}
		compositeNeighbors.clear();
		// components
		if(components.size() > 1)
		{
			// Cluster components into families
			map<unsigned long long int, vector<unsigned long long int> > componentFamilies;
			map<unsigned long long int, vector<unsigned long long int> >& refCF = componentFamilies;
			map<unsigned long long int, vector<unsigned long long int> >::iterator cf;
			computeComponentFamilies(refCP, refCF, genes);
			//
			map<unsigned long long int, pair<unsigned int,unsigned int>> familiesZone;
			// area - familyID
			vector<pair<unsigned long long int,unsigned long long int>> orderZones;
			// is good Domain
			map<unsigned long long int,bool> isGoodDomain;
			// mean pident by family
			map<unsigned long long int,float> familyPident;
			// sequence percentage by family
			map<unsigned long long int,float> familySeqPercent;
			// Decteted components family zone
			for(cf = componentFamilies.begin(); cf != componentFamilies.end(); cf++)
			{
				//vistedFamily[cf->first]=false;
				isGoodDomain[cf->first]=false;
				//cout << "Family " << cf->first << endl;
				vector<unsigned int> startVal;
				vector<unsigned int>& refSV = startVal;
				vector<unsigned int> endVal;
				vector<unsigned int>& refEV = endVal;
				vector<unsigned long long int>::iterator c;
				for(c = componentFamilies[cf->first].begin(); c != componentFamilies[cf->first].end(); c++)
				{
					edgeValues edge = getEdgeValues(edges, composite, *c);
					startVal.push_back(edge.getQstart());
					endVal.push_back(edge.getQend());
					//cout << composite << "\t" << *c << "\t" << edge.getPident() << endl;
				}
				//
				unsigned int startMedian = computeMedian(refSV);
				unsigned int endMedian = computeMedian(refEV);
				familiesZone[cf->first]=make_pair(startMedian,endMedian);
				//
				//+ to_string(endMedian);
				string zone = to_string(startMedian);
				unsigned long long int zoneInt = atoll(zone.c_str()); 
				orderZones.push_back(make_pair(zoneInt,cf->first));
			}
			sort(orderZones.begin(),orderZones.end());
			//
			unsigned int z;
			// find good component families
			for(z = 0; z < orderZones.size(); z++)
			{
				unsigned int domainCount = 0;
				unsigned int inter1Start;
				unsigned int inter1End;
				unsigned int inter2Start;
				unsigned int inter2End;
				inter1Start = familiesZone[orderZones.at(z).second].first;
				inter1End = familiesZone[orderZones.at(z).second].second;
				unsigned int u;
				unsigned int next = z + 1;
				for(u = next; u < orderZones.size(); u++)
				{
					inter2Start = familiesZone[orderZones.at(u).second].first;
					inter2End = familiesZone[orderZones.at(u).second].second;
					if(checkOverlap(MAX_OVERLAP,inter1Start,inter1End,inter2Start,inter2End))
					{
						//cout << inter1Start << " " << inter1End << " " << inter2Start << " " << inter2End << endl;
						isGoodDomain[orderZones.at(z).second]=true;
						isGoodDomain[orderZones.at(u).second]=true;
						inter1Start = inter2Start;
						inter1End = inter2End;
					}
				}	
			}
			// check overlaps between component families
			map<unsigned long long int, unsigned long long int> getClusterID;
			map<unsigned long long int,vector<unsigned long long int>> clusters;
			vector<unsigned long long int> goodDomains;
			for(z = 0; z < orderZones.size(); z++)
			{
				if(isGoodDomain[orderZones.at(z).second] == true)
				{
					getClusterID[orderZones.at(z).second]=0;
					goodDomains.push_back(orderZones.at(z).second);
				}
			}
			unsigned long long int clusterID = 0;
			for(z = 0; z < goodDomains.size(); z++)
			{
				unsigned long long int id = goodDomains.at(z);
				if(getClusterID[id] == 0)
				{
					clusterID++;
					getClusterID[id]=clusterID;
					clusters[clusterID].push_back(id);
				}
				unsigned int inter1Start;
				unsigned int inter1End;
				unsigned int inter2Start;
				unsigned int inter2End;
				inter1Start = familiesZone[id].first;
				inter1End = familiesZone[id].second;
				unsigned int u;
				for(u = 0; u < goodDomains.size(); u++)
				{
					unsigned long long int nextId = goodDomains.at(u);
					inter2Start = familiesZone[nextId].first;
					inter2End = familiesZone[nextId].second;
					if(getClusterID[id] == getClusterID[nextId] && getClusterID[id] != 0)
					{
						continue;
					}
					else if(nextId == id)
					{
						continue;
					}
					else if(inter2Start >= inter1End)
					{
						break;
					}
					else if(inter2End <= inter1Start)
					{
						continue;
					}
					if(mergeRegion(inter1Start,inter1End,inter2Start,inter2End,minCov,MAX_OVERLAP))
					{
						if(getClusterID[nextId] != 0 && getClusterID[id] == 0)
						{
							unsigned long long int NC = getClusterID[nextId];
							clusters[NC].push_back(id);
							getClusterID[id]=NC;
						}
						else if(getClusterID[id] != 0 && getClusterID[nextId] == 0)
						{
							unsigned long long int NC = getClusterID[id];
							clusters[NC].push_back(nextId);
							getClusterID[nextId]=NC;
						}
					}
				}
			}
			//
			if(clusters.size() > 1)
			{
				unsigned int NbDomains = 0;
				unsigned long long int overlap = 0;
				unsigned long long int familyPairs = 0;
				map<unsigned long long int,unsigned long long int> regionStart;
				map<unsigned long long int,unsigned long long int> regionEnd;
				map<unsigned long long int,unsigned long long int> regionID;
				map<unsigned long long int,vector<unsigned long long int>>::iterator cid;
				for(cid = clusters.begin(); cid != clusters.end(); cid++)
				{
					unsigned long long int rStart = 0;
					unsigned long long int rEnd = 0;
					NbDomains++;
					vector<unsigned long long int>::iterator g;
					for(g = clusters[cid->first].begin(); g != clusters[cid->first].end(); g++)
					{
						unsigned long long int tmpStart;
			                        unsigned long long int tmpEnd;
						tmpStart = familiesZone[*g].first;
		                                tmpEnd = familiesZone[*g].second;
						if(rStart > tmpStart || rStart == 0)
						{
							rStart = tmpStart;
						}
						if(rEnd < tmpEnd || rEnd == 0)
						{
							rEnd = tmpEnd;
						}
					}
					regionStart[NbDomains]=rStart;
					regionEnd[NbDomains]=rEnd;
					regionID[NbDomains]=cid->first;
				}
				//
				genes[composite].setComposite(true);
		                genes[composite].setNbDomains(NbDomains);
		                geneFamilies[genes[composite].getFamilyID()].setComposite(true);
				unsigned int nbD;
				unsigned int singleNodeFamilies = 0;
                                unsigned int NbFamilies = 0;
				output << ">C" << composite << endl;
				for(nbD = 1; nbD < NbDomains+1; nbD++)
				{
					if(nbD == NbDomains)
					{
						output << "[Domain " << nbD << "]\t" << regionStart[nbD] << "\t" << regionEnd[nbD] << endl;
					}
					else
					{
						unsigned long long int End1 = regionEnd[nbD];
						unsigned long long int Start2 = regionStart[nbD+1];
						if(End1 <= Start2)
						{
							output << "[Domain " << nbD << "]\t" << regionStart[nbD] << "\t" << regionEnd[nbD] << endl;
						}
						else
						{
							overlap++;
							unsigned long long int diff = End1 - Start2 + 1;
							unsigned long long int mean = floor((float)diff/2.0);
							regionEnd[nbD] = Start2 + mean;
							regionStart[nbD+1] = Start2 + mean;
							output << "[Domain " << nbD << "]\t" << regionStart[nbD] << "\t" << regionEnd[nbD] << endl;
						}
					}
					vector<unsigned long long int>::iterator g;
					for(g = clusters[regionID[nbD]].begin(); g != clusters[regionID[nbD]].end(); g++)
					{
						NbFamilies++;
						if(geneFamilies[*g].getSize() == 1) singleNodeFamilies++;
						genes[composite].insertComponentFamilyID(*g);
						vector<unsigned long long int>::iterator c;
						for(c = componentFamilies[*g].begin(); c != componentFamilies[*g].end(); c++)	
						{
							edgeValues edge = getEdgeValues(edges, composite, *c);
							output << "F" << *g << "\t" << *c << "\t" << edge.getQstart() << "\t" << edge.getQend() << "\t" << edge.getPident() << endl;
						}
					}
				}
				familyPairs = NbDomains - 1;
				float noOverlap = 1.0-((float)overlap/(float)familyPairs);
				genes[composite].setNoOverlap(noOverlap);
				//
				info << "C" << composite << "\t" << genes[composite].getLength() << "\tF" << genes[composite].getFamilyID() << "\t" << NbFamilies << "\t" << singleNodeFamilies << "\t" << genes[composite].getNbDomains() << "\t" << setprecision (2) << fixed << genes[composite].getNoOverlap() << endl;
		      	 	orderZones.clear();
                      	  	familiesZone.clear();
                        	componentFamilies.clear();
                        	isGoodDomain.clear();
                        	regionStart.clear();
                       		regionEnd.clear();
                        	regionID.clear();
                        	getClusterID.clear();
                        	clusters.clear();
			}
		}
		processedNodes++;
		unsigned int cpuP = ((float)processedNodes/(float)totalNodes)*100.0;
		ofstream progress(progressInfo.c_str());
		progress << "Composites - CPU\t" << l << "\t" << cpuP << "%" << endl;
	}
}

void runComposites( map<unsigned long long int, geneInfo >& genes, map<pair<unsigned long long int, unsigned long long int>, edgeValues>& edges, map<unsigned long long int, familyInfo>& geneFamilies, vector<pair<unsigned long long int,unsigned long long int>>& nodesDegree, unsigned short int MAX_OVERLAP, unsigned short int minCov, unsigned int compositeFamMinSize, unsigned int componentFamMinSize, unsigned int nCpu)
{
	//Map used to split nodes list into sublist
	map<unsigned int, list<unsigned long long int> > subList;
	map<unsigned int, list<unsigned long long int> >& refSubList = subList;
	unsigned int k = 0;
	//Split nodes to distribute
	unsigned int n;
	for( n = 0; n < nodesDegree.size(); n++)
	{
		if( k == nCpu )
		{
			k = 0;
		}
		subList[k].push_back(nodesDegree[n].second);
		k++;
	}
	// MultiThreading
	vector<thread> threads;
	for(unsigned int i = 0; i < nCpu; ++i)
	{
		threads.push_back(thread(composites, std::ref(refSubList), std::ref(genes), std::ref(edges), std::ref(geneFamilies), MAX_OVERLAP, minCov, compositeFamMinSize, componentFamMinSize, i));
	}
	for(auto& thread : threads){
		thread.join();
	}
	// delete progress file
	system("rm *.progress");
}
//end
