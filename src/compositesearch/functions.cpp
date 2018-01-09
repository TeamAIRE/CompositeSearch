/*

        Written by Jananan PATHMANATHAN, 2014-2018
        
        This file is part of CompositeSearch.
        
        CompositeSearch is shared under Creative commons licence: 
        
        Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
        
        See https://creativecommons.org/licenses/by-nc-sa/4.0/

*/

#include "functions.h"
#include <iostream>
#include <map>
#include <list>
#include <vector>
#include <math.h>
#include <algorithm>
#include <cmath>

using namespace std;

// Classe geneInfo
geneInfo::geneInfo() : m_length(0), m_composite(false), m_NbDomains(0), m_noOverlap(0), m_familyID(0)
{

}

geneInfo::geneInfo(unsigned int length, bool composite, unsigned int NbDomains, float noOverlap, unsigned long long int familyID) : m_length(length), m_composite(composite), m_NbDomains(NbDomains), m_noOverlap(noOverlap), m_familyID(familyID)
{

}

// Set or insert values
void geneInfo::setLength(unsigned int length)
{
	m_length = length;
}

void geneInfo::setComposite(bool composite)
{
	m_composite = composite;
}

void geneInfo::setNbDomains(unsigned int NbDomains)
{
	m_NbDomains = NbDomains;
}

void geneInfo::setNoOverlap( float noOverlap)
{
	m_noOverlap = noOverlap;
}

void geneInfo::setFamilyID(unsigned long long int familyID)
{
	m_familyID = familyID;
}

void geneInfo::insertComponentFamilyID(unsigned long long int ID)
{
	m_componentFamiliesID.push_back(ID);
}

void geneInfo::insertNeighbor(unsigned long long int neighbor)
{
	m_neighbors.push_back(neighbor);
}

void geneInfo::sortNeighbors()
{
	sort(m_neighbors.begin(),m_neighbors.end());
}

// get values

unsigned int geneInfo::getLength()
{
	return m_length;
}

bool geneInfo::isComposite()
{
	return m_composite;
}

unsigned int geneInfo::getNbDomains()
{
	return m_NbDomains;
}

float geneInfo::getNoOverlap()
{
	return m_noOverlap;
}

unsigned long long int geneInfo::getFamilyID()
{
	return m_familyID;
}

vector<unsigned long long int> geneInfo::getComponentFamiliesID()
{
	return m_componentFamiliesID;
}

vector<unsigned long long int> geneInfo::getNeighbors()
{
	return m_neighbors;
}

unsigned long long int geneInfo::getNeighborsSize()
{
	return m_neighbors.size();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Class edgeValues
edgeValues::edgeValues() : m_qstart(0), m_qend(0), m_sstart(0), m_send(0), m_pident(0), m_evalue(0)
{

}

edgeValues::edgeValues(unsigned short int qstart, unsigned short int qend, unsigned short int sstart, unsigned short int send, float pident, double evalue) : m_qstart(qstart), m_qend(qend), m_sstart(sstart), m_send(send), m_pident(pident), m_evalue(evalue)
{

}

void edgeValues::setQstart(unsigned short int qstart)
{
	m_qstart = qstart;
}

void edgeValues::setQend(unsigned short int qend)
{
	m_qend = qend;
}

void edgeValues::setSstart(unsigned short int sstart)
{
	m_sstart = sstart;
}

void edgeValues::setSend(unsigned short int send)
{
	m_send = send;
}

void edgeValues::setPident(float pident)
{
        m_pident = pident;
}

void edgeValues::setEvalue(double evalue)
{
	m_evalue = evalue;
}

unsigned short int edgeValues::getQstart()
{
	return m_qstart;
}

unsigned short int edgeValues::getQend()
{
	return m_qend;
}

unsigned short int edgeValues::getSstart()
{
	return m_sstart;
}

unsigned short int edgeValues::getSend()
{
	return m_send;
}

float edgeValues::getPident()
{
        return m_pident;
}

double edgeValues::getEvalue()
{
	return m_evalue;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
familyInfo::familyInfo() : m_connectivity(0), m_composite(false)
{

}

familyInfo::familyInfo(float connectivity, bool composite) : m_connectivity(connectivity), m_composite(composite)
{

}

// set and insert info
void familyInfo::setConnectivity(float connectivity)
{
	m_connectivity = connectivity;
}

void familyInfo::setComposite(bool composite)
{
	m_composite = composite;
}

void familyInfo::insertGeneID(unsigned long long int ID)
{
	m_genes.push_back(ID);
}

void familyInfo::sortGenes()
{
	sort(m_genes.begin(),m_genes.end());
}


// get info
float familyInfo::getConnectivity()
{
	return m_connectivity;
}

bool familyInfo::isComposite()
{
	return m_composite;
}

vector<unsigned long long int> familyInfo::getGenesID()
{
        return m_genes;
}

unsigned long long int familyInfo::getSize()
{
	return m_genes.size();
}


// Get edge values for a given hit
edgeValues getEdgeValues(map<pair<unsigned long long int,unsigned long long int>, edgeValues> &edges, unsigned long long int a, unsigned long long int b)
{
	if( a > b )
	{
		edgeValues toReverse = edges[make_pair(b,a)];
		return reverseEdgeValues(toReverse);
	}
	else
	{
		return edges[make_pair(a,b)];
	}
}

// Function reversing edge values
edgeValues reverseEdgeValues(edgeValues edge)
{
        unsigned short int qstart, qend, qlen, sstart, send, slen;
        // Initial values
        qstart = edge.getSstart();
        qend = edge.getSend();
        sstart = edge.getQstart();
        send = edge.getQend();
        // Set reverse values
        edge.setQstart(qstart);
        edge.setQend(qend);
        edge.setSstart(sstart);
        edge.setSend(send);
        // Return
        return edge;
}


// Function checking the overlap between two components
bool checkTriplet(edgeValues edge1, edgeValues edge2, unsigned int MAX_OVERLAP)//, unsigned int& noOverLap)
{
	if( edge1.getQend() < edge2.getQstart() + MAX_OVERLAP )
	{
		//cout << edge1.getQend() << " e1e < " << edge2.getQstart() << " e2s + " << MAX_OVERLAP << endl;
		return true;
	}
	else if( edge2.getQend() < edge1.getQstart() + MAX_OVERLAP )
	{
		//cout << edge2.getQend() << " e1e < " << edge1.getQstart() << " e2s + " << MAX_OVERLAP << endl;
		return true;
	}
	else
	{
		//cout << edge1.getQend() << " e1e " << edge2.getQstart() << " e2s  and mxo " << MAX_OVERLAP << endl;
		return false;
	}
}

// Dichotomic search to check if a value is present in a list
bool FoundIn(unsigned long long int value, vector<unsigned long long int>& List)
{
	bool found;
	unsigned long long int start;
	unsigned long long int end;
	unsigned long long int middle;
	// initialization
	found = false;
	start = 0;
	end = List.size();
	// start dichotomic search
	while(!found && ((end - start) > 1))
	{
		middle = (start + end)/2;
		found = (List.at(middle) == value);
		if(List.at(middle) > value)
		{
			end = middle;
		}
		else
		{
			start = middle;
		}
	}
	if(List.at(start) == value)
	{
		 return true;
	}
	else
	{
		return false;
	}
}

bool checkCoverage(unsigned long long int node1, unsigned long long int node2, unsigned short int minCov, map<unsigned long long int, geneInfo>& genes, map<pair<unsigned long long int, unsigned long long int>, edgeValues>& edges)
{
        // Edge value
        edgeValues edge = getEdgeValues(edges, node1, node2);
        // Query coverage
        unsigned short int qcov = floor((((float)edge.getQend()-(float)edge.getQstart()+1.0)*100.0)/(float)genes[node1].getLength());
        // Subject coverage
        unsigned short int scov = floor((((float)edge.getSend()-(float)edge.getSstart()+1.0)*100.0)/(float)genes[node2].getLength());
        // Return true if both cov are higher or equal to the minimum coverage value 
        // fixed by the user or the default one (80%)
        if(qcov >= minCov && scov >= minCov)
        {
                return true;
        }
        else
        {
                return false;
        }
}

unsigned int computeMedian(vector<unsigned int>& val)
{
	// med
	unsigned int median;
	// Sort
	sort(val.begin(),val.end());
	// Size
	unsigned int N = val.size();
	// Median
	if( N%2 == 0 )
	{
		unsigned int n = N/2;
		median = floor(((float)val[n-1] + (float)val[n])/2.0);
	}
	else
	{
		median = val[floor((float)N/2.0)];
	}
	return median;
}

bool checkOverlap(unsigned int MAX_OVERLAP, unsigned int inter1Start, unsigned int inter1End, unsigned int inter2Start, unsigned int inter2End)
{
	if( inter1End < inter2Start + MAX_OVERLAP )
	{
		return true;
	}
	else if( inter2End < inter1Start + MAX_OVERLAP )
	{
		return true;
	}
	else
	{
		return false;
	}
}


bool mergeRegion(unsigned int inter1Start,unsigned int inter1End, unsigned int inter2Start, unsigned int inter2End, unsigned short int minCov, unsigned int MAX_OVERLAP)
{
	unsigned int START;
	unsigned int END;
	if(inter1Start < inter2Start)
	{
		START = inter2Start;
	}
	else
	{
		START = inter1Start;
	}
	if(inter1End < inter2End)
	{
		END = inter1End;
	}
	else
	{
		END = inter2End;
	}
	//
	unsigned int length1 = inter1End - inter1Start + 1;
	unsigned int length2 = inter2End - inter2Start + 1;
	// Overlap
	unsigned int overlapLength = END - START + 1;
	unsigned short int cscov = floor(((float)overlapLength / (float)length1)*100.0);
	unsigned short int nscov = floor(((float)overlapLength / (float)length2)*100.0);
	if(cscov >= minCov || nscov >= minCov || overlapLength > MAX_OVERLAP)
	{
		return true;
	}
	else
	{
		return false;
	}
}

//END
